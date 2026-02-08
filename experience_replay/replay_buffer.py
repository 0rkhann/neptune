"""
Some code is based on the implementation from https://github.com/MolecularAI/Reinvent.

Supports:
- Saturn (SMILES): Bemis-Murcko scaffolds
- Neptune (HELM): Bigram scaffolds (when use_bigram_scaffold=True)
"""

from typing import Tuple, List
import logging
import numpy as np
import pandas as pd
from copy import deepcopy
from rdkit import Chem
from utils.chemistry_utils import (
    randomize_smiles_batch,
    get_bemis_murcko_scaffold,
)

from experience_replay.dataclass import ExperienceReplayParameters

# Oracle is called if seeding molecules into the Replay Buffer at the start of the generative experiment
from oracles.oracle import Oracle

# Import bigram scaffold utilities (for Neptune/HELM only)
try:
    from utils.bigram_diversity_utils import get_bigram_scaffold
    from utils.helm import is_helm_notation

    BIGRAM_AVAILABLE = True
except ImportError:
    BIGRAM_AVAILABLE = False
    get_bigram_scaffold = None
    is_helm_notation = None


class ReplayBuffer:
    """
    Replay buffer class which stores the top N highest reward sequences generated so far.

    Supports:
    - Saturn (SMILES): Stores SMILES, uses Bemis-Murcko scaffolds for purging
    - Neptune (HELM): Stores HELM sequences, uses bigram scaffolds for purging
    """

    def __init__(self, parameters: ExperienceReplayParameters):
        self.parameters = parameters
        self.memory_size = parameters.memory_size
        self.sample_size = parameters.sample_size
        # Stores the top N highest reward sequences generated so far
        self.memory = pd.DataFrame(columns=["smiles", "reward"])

        # Bigram scaffold support (Neptune/HELM only)
        self.use_bigram_scaffold = parameters.use_bigram_scaffold
        if self.use_bigram_scaffold:
            if not BIGRAM_AVAILABLE:
                raise ImportError(
                    "Bigram scaffold support requires bigram_diversity_utils. "
                    "Please ensure these modules are available."
                )
            logging.info(
                "[ReplayBuffer] Using bigram scaffolds for HELM sequences (Neptune)"
            )
        else:
            logging.info(
                "[ReplayBuffer] Using Bemis-Murcko scaffolds for SMILES (Saturn)"
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

    def _get_scaffold(self, sequence: str) -> str:
        """
        Get scaffold for a sequence.

        - Neptune (HELM with use_bigram_scaffold=True): bigram scaffolds
        - Saturn (SMILES): Bemis-Murcko scaffolds

        NOTE: Bigram scaffolds are ONLY supported for HELM (Neptune).
        """
        if self.use_bigram_scaffold:
            # Bigram scaffolds are ONLY for HELM (Neptune)
            if not is_helm_notation(sequence):
                raise ValueError(
                    f"use_bigram_scaffold=True requires HELM sequences (Neptune), "
                    f"but received SMILES: {sequence[:60]}... "
                    f"Use use_bigram_scaffold=False for Saturn (SMILES)."
                )

            scaffold = get_bigram_scaffold(sequence)
            if not scaffold:
                logging.error(
                    f"[ReplayBuffer] Failed to extract bigram scaffold from HELM: {sequence[:60]}..."
                )
                raise ValueError(
                    f"Failed to extract bigram scaffold from HELM: {sequence[:60]}..."
                )
            return scaffold
        else:
            # Bemis-Murcko scaffolds for SMILES (Saturn)
            return get_bemis_murcko_scaffold(sequence)

    def selective_memory_purge(
        self, sequences: np.ndarray[str], rewards: np.ndarray[float]
    ) -> None:
        """
        Augmented Memory's key operation to prevent mode collapse and promote diversity:
        Purges the memory of sequences that have penalized rewards (0.0) *before* executing Augmented Memory updates.
        Intuitively, this operation prevents penalized sequences from directing the Agent's chemical space navigation.

        - Neptune (HELM): uses bigram scaffolds (matching bigram diversity filter)
        - Saturn (SMILES): uses Bemis-Murcko scaffolds

        # NOTE: Consider a MPO objective task using a product aggregator. If one of the OracleComponent's reward is 0,
        #       then the aggregated reward may be 0. But other OracleComponents may have a non-zero reward. We do not
        #       want to purge the memory of these scaffolds. This is already handled because 0 reward sequences are not
        #       added to the memory in the first place. Selective Memory Purge *only* removes scaffolds that are
        #       penalized by the Diversity Filter.
        """
        zero_reward_indices = np.where(rewards == 0.0)[0]
        if len(zero_reward_indices) > 0:
            sequences_to_purge = sequences[zero_reward_indices]
            scaffolds_to_purge = [self._get_scaffold(s) for s in sequences_to_purge]
            purged_memory = deepcopy(self.memory)
            purged_memory["scaffolds"] = purged_memory["smiles"].apply(
                self._get_scaffold
            )
            purged_memory = purged_memory.loc[
                ~purged_memory["scaffolds"].isin(scaffolds_to_purge)
            ]
            purged_memory.drop("scaffolds", axis=1, inplace=True)
            self.memory = purged_memory
        else:
            # If no scaffolds are penalized, do nothing
            return

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
