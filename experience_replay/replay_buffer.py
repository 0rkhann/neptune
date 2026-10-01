"""
Some code is based on the implementation from https://github.com/MolecularAI/Reinvent.

Supports:
- Saturn (SMILES): Bemis-Murcko scaffolds
- Neptune (HELM): Order-independent fingerprints (scaffold_type="subclass"/"monomer")
- Atom-level (SMILES): ECFP4 sidechain fingerprints (scaffold_type="sidechain")
"""

from typing import Tuple, List
import logging
import numpy as np
import pandas as pd
from copy import deepcopy
from rdkit import Chem, DataStructs
from utils.chemistry_utils import (
    randomize_smiles_batch,
)

from experience_replay.dataclass import ExperienceReplayParameters
from utils.diversity_utils import (
    to_smiles,
    get_diversity_fingerprint,
    HELM_AVAILABLE,
    SIDECHAIN_FP_AVAILABLE,
    get_sidechain_fingerprint,
    get_subclass_fingerprint,
    get_monomer_fingerprint,
    subclass_fingerprint_from_monomers,
    monomer_fingerprint_from_monomers,
)

# Oracle is called if seeding molecules into the Replay Buffer at the start of the generative experiment
from oracles.oracle import Oracle


class ReplayBuffer:
    """
    Replay buffer class which stores the top N highest reward sequences generated so far.

    Supports:
    - Saturn (SMILES): Stores SMILES, uses Bemis-Murcko scaffolds for purging
    - Neptune (HELM): Stores HELM sequences, uses order-independent fingerprints for purging
    - Atom-level (SMILES): ECFP4 sidechain fingerprints for Tanimoto-based purging
    """

    def __init__(self, parameters: ExperienceReplayParameters):
        self.parameters = parameters
        self.memory_size = parameters.memory_size
        self.sample_size = parameters.sample_size
        # Stores the top N highest reward sequences generated so far
        self.memory = pd.DataFrame(columns=["smiles", "reward"])

        # Scaffold type: "subclass", "monomer", "sidechain", or None
        self.scaffold_type = parameters.scaffold_type

        # SynthesizabilityChecker and threshold for sidechain mode (injected by RL agent via set_synth_checker)
        self._synth_checker = None
        self._tanimoto_threshold = 0.65  # Default; overridden by set_synth_checker

        if self.scaffold_type == "sidechain":
            if not SIDECHAIN_FP_AVAILABLE:
                raise ImportError(
                    "scaffold_type='sidechain' requires GenAI4Peptidomimetic_native. "
                    "Install with: pip install -e <GenAI4Peptidomimetic_native_dir>"
                )
            self._fingerprint_fn = None  # Not used; sidechain mode uses Tanimoto comparison
            self._monomer_fingerprint_fn = None
            logging.info(
                "[ReplayBuffer] Using XOR-ECFP sidechain fingerprints for selective memory purge"
            )
        elif self.scaffold_type in ("subclass", "monomer"):
            if not HELM_AVAILABLE:
                raise ImportError(
                    "HELM fingerprint support requires subclass_fingerprint. "
                    "Please ensure these modules are available."
                )
            if self.scaffold_type == "subclass":
                self._fingerprint_fn = get_subclass_fingerprint
                self._monomer_fingerprint_fn = subclass_fingerprint_from_monomers
            else:
                self._fingerprint_fn = get_monomer_fingerprint
                self._monomer_fingerprint_fn = monomer_fingerprint_from_monomers
            logging.info(
                f"[ReplayBuffer] Using {self.scaffold_type} fingerprints"
            )
        else:
            self._fingerprint_fn = None
            self._monomer_fingerprint_fn = None
            logging.info(
                "[ReplayBuffer] Using Bemis-Murcko scaffolds for SMILES (Saturn)"
            )

    def set_synth_checker(self, synth_checker, tanimoto_threshold: float = 0.65) -> None:
        """Inject SynthesizabilityChecker (shared with DiversityFilter)."""
        self._synth_checker = synth_checker
        self._tanimoto_threshold = tanimoto_threshold
        logging.info(
            f"[ReplayBuffer] SynthesizabilityChecker injected "
            f"(tanimoto_threshold={tanimoto_threshold})"
        )

    def add(self, smiles: np.ndarray[str], rewards: np.ndarray[float]) -> None:
        df = pd.DataFrame({"smiles": smiles, "reward": rewards})
        self.memory = pd.concat([self.memory, df]) if len(self.memory) > 0 else df
        # Keep only the top N (by reward)
        self.purge_memory()

    def sample_memory(self) -> Tuple[np.ndarray[str], np.ndarray[float]]:
        sample_size = min(len(self.memory), self.sample_size)
        if sample_size > 0:
            sampled = self.memory.sample(sample_size)
            smiles = sampled["smiles"].values
            rewards = sampled["reward"].values
            return np.array(smiles), np.array(rewards)
        else:
            return [], []

    def augmented_memory_replay(self, prior) -> Tuple[List[str], np.array]:
        """
        Augmented Memory's key operation for sample efficiency:
        Randomizes all SMILES in the memory and returns the randomized SMILES and their corresponding rewards.
        """
        if len(self.memory) != 0:
            smiles = self.memory["smiles"].values
            # Randomize the smiles
            randomized_smiles = randomize_smiles_batch(smiles, prior)
            rewards = self.memory["reward"].values
            return randomized_smiles, rewards
        else:
            return [], []

    def purge_memory(self):
        """
        Removes duplicate sequences in the memory and keeps the top N by reward.
        Keep only non-zero rewards sequences.
        """
        unique_df = self.memory.drop_duplicates(subset=["smiles"])
        sorted_df = unique_df.sort_values("reward", ascending=False)
        self.memory = sorted_df.head(self.memory_size)
        self.memory = self.memory.loc[self.memory["reward"] != 0.0]

    def selective_memory_purge(
        self, sequences: np.ndarray[str], rewards: np.ndarray[float]
    ) -> None:
        """
        Augmented Memory's key operation to prevent mode collapse and promote diversity:
        Purges the memory of sequences that have penalized rewards (0.0) *before* executing Augmented Memory updates.
        Intuitively, this operation prevents penalized sequences from directing the Agent's chemical space navigation.

        - Neptune (HELM): uses order-independent fingerprints (matching diversity filter)
        - Saturn (SMILES): uses Bemis-Murcko scaffolds
        - Sidechain mode: uses Tanimoto similarity between XOR-ECFP fingerprints

        # NOTE: Consider a MPO objective task using a product aggregator. If one of the OracleComponent's reward is 0,
        #       then the aggregated reward may be 0. But other OracleComponents may have a non-zero reward. We do not
        #       want to purge the memory of these scaffolds. This is already handled because 0 reward sequences are not
        #       added to the memory in the first place. Selective Memory Purge *only* removes scaffolds that are
        #       penalized by the Diversity Filter.
        """
        if self.scaffold_type == "sidechain":
            self._sidechain_memory_purge(sequences, rewards)
            return

        zero_reward_indices = np.where(rewards == 0.0)[0]
        if len(zero_reward_indices) > 0:
            sequences_to_purge = sequences[zero_reward_indices]
            scaffolds_to_purge = [
                get_diversity_fingerprint(
                    s, self._fingerprint_fn, self._monomer_fingerprint_fn,
                    self._synth_checker, "[ReplayBuffer]"
                )
                for s in sequences_to_purge
            ]
            purged_memory = deepcopy(self.memory)
            purged_memory["scaffolds"] = purged_memory["smiles"].apply(
                lambda s: get_diversity_fingerprint(
                    s, self._fingerprint_fn, self._monomer_fingerprint_fn,
                    self._synth_checker, "[ReplayBuffer]"
                )
            )
            purged_memory = purged_memory.loc[
                ~purged_memory["scaffolds"].isin(scaffolds_to_purge)
            ]
            purged_memory.drop("scaffolds", axis=1, inplace=True)
            self.memory = purged_memory
        else:
            # If no scaffolds are penalized, do nothing
            return

    def _sidechain_memory_purge(
        self, sequences: np.ndarray[str], rewards: np.ndarray[float]
    ) -> None:
        """
        Tanimoto-based selective memory purge for SidechainFingerprint mode.

        Computes ECFP4 fingerprints for penalized sequences (reward=0.0),
        then removes buffer entries whose FP has Tanimoto > threshold with any
        penalized FP.
        """
        zero_idx = np.where(rewards == 0.0)[0]
        if len(zero_idx) == 0:
            return

        # Compute sidechain FPs for penalized sequences
        penalized_fps = []
        for i in zero_idx:
            smiles = to_smiles(sequences[i])
            fp = get_sidechain_fingerprint(smiles, self._synth_checker) if smiles else None
            if fp is not None:
                penalized_fps.append(fp)

        if not penalized_fps:
            return

        # Filter memory: remove entries similar to any penalized FP
        keep_mask = []
        for _, row in self.memory.iterrows():
            smiles = to_smiles(row["smiles"])
            mem_fp = get_sidechain_fingerprint(smiles, self._synth_checker) if smiles else None
            if mem_fp is None:
                keep_mask.append(True)  # Keep if FP can't be computed
                continue
            similar = any(
                DataStructs.TanimotoSimilarity(mem_fp, pfp) > self._tanimoto_threshold
                for pfp in penalized_fps
            )
            keep_mask.append(not similar)

        self.memory = self.memory.loc[keep_mask]

    def prepopulate_buffer(self, oracle: Oracle) -> Oracle:
        """
        Seeds the Replay Buffer with a set of sequences (SMILES or HELM).
        Useful if there are known high-reward molecules to pre-populate the Replay Buffer with.

        Oracle is returned here because seeding updates the Oracle's cache with the seeded sequences.

        NOTE: With more sequences to seed with, the generative experiment will become more like
              transfer learning rather than reinforcement learning (at the start). Continuing
              the run will more and more leverage reinforcement learning to find other diverse
              solutions. Therefore, while seeding will quick-start the Agent's learning, there
              are implications on the diversity of the solutions found.
        """
        if len(self.parameters.sequences) > 0:
            from src.genai_utils.helm import is_helm_notation, convert_helm_to_smiles

            # Convert sequences (SMILES or HELM) to canonical SMILES for oracle evaluation
            canonical_smiles = []
            for seq in self.parameters.sequences:
                if is_helm_notation(seq):
                    # Convert HELM to SMILES for oracle evaluation
                    smiles = convert_helm_to_smiles(seq)
                    if smiles:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            canonical_smiles.append(Chem.MolToSmiles(mol))
                else:
                    # It's already SMILES
                    mol = Chem.MolFromSmiles(seq)
                    if mol:
                        canonical_smiles.append(Chem.MolToSmiles(mol))

            assert len(self.parameters.sequences) == len(
                canonical_smiles
            ), f"Warmup sequences are not all valid: {len(self.parameters.sequences)} sequences vs {len(canonical_smiles)} valid SMILES"

            # Create RDKit molecules from canonical SMILES for oracle evaluation
            mols = [Chem.MolFromSmiles(s) for s in canonical_smiles]

            oracle_components_df = pd.DataFrame()
            rewards = np.empty((len(oracle.oracle), len(mols)))
            for idx, oracle_component in enumerate(oracle.oracle):
                raw_property_values, component_rewards = (
                    oracle_component.calculate_reward(mols, oracle_calls=0)
                )
                oracle_components_df[f"{oracle_component.name}_raw_values"] = (
                    raw_property_values
                )
                oracle_components_df[f"{oracle_component.name}_reward"] = (
                    component_rewards
                )
                rewards[idx] = component_rewards

            aggregated_rewards = oracle.aggregator(rewards, oracle.oracle_weights)

            # Add the original sequences to the Replay Buffer (HELM for Neptune, SMILES for Saturn)
            # The model needs to be able to tokenize these sequences
            self.add(smiles=self.parameters.sequences, rewards=aggregated_rewards)

            # Update the Oracle Cache with the canonical SMILES (to avoid re-evaluating if rediscovered)
            # NOTE: We do NOT add warmup molecules to oracle_history.csv or increment oracle.calls
            # because they are pre-selected seeds, not generated by the model
            oracle.update_oracle_cache(canonical_smiles, aggregated_rewards)

        return oracle
