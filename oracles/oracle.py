from typing import List, Tuple, Optional
import os
import json
import pandas as pd
import numpy as np
import logging
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import Mol
from utils.chemistry_utils import canonicalize_smiles_batch, get_bemis_murcko_scaffold

from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters, OracleConfiguration
from oracles.oracle_utils import extract_oracle_metadata
from oracles.reward_aggregator.reward_aggregator import RewardAggregator
from diversity_filter.diversity_filter import DiversityFilter

from oracles.utils import construct_oracle_component

# HELM to SMILES conversion for Neptune
try:
    from utils.helm import (
        convert_helm_to_smiles,
        is_helm_notation,
        load_valid_monomers,
        check_helm_validness,
    )

    HELM_CONVERSION_AVAILABLE = True
except ImportError:
    HELM_CONVERSION_AVAILABLE = False
    convert_helm_to_smiles = None
    is_helm_notation = None
    load_valid_monomers = None
    check_helm_validness = None


class Oracle:
    """
    The Oracle function to optimize for in the generative experiment.
    Can be composed of multiple OracleComponents, each of which handles a specific property.
    Aggregating the rewards from each OracleComponent returns a scalar reward (return in RL terminology) for Agent update.
    """

    def __init__(self, oracle_configuration: OracleConfiguration):
        self.oracle_configuration = oracle_configuration

        # Construct the oracle function which can be composed of >1 individual oracles (multi-parameter optimization)
        self.oracle = self.construct_oracle(oracle_configuration.components)
        # Preliminary oracles can be executed as a first pass to filter out poor candidates
        self.preliminary_oracles = [
            oracle for oracle in self.oracle if oracle.preliminary_check
        ]
        self.oracle_weights = [oracle.weight for oracle in self.oracle]
        self.aggregator = RewardAggregator(oracle_configuration.aggregator)

        # Track oracle budget
        self.budget = oracle_configuration.budget
        self.allow_oracle_repeats = oracle_configuration.allow_oracle_repeats
        self.calls = 0

        # Cache dictionary to store the results of previous oracle calls
        self.cache = dict()

        # Oracle history to assess sample efficiency via Generative Yield and Oracle Burden metrics
        self.oracle_history = pd.DataFrame(
            {
                "oracle_calls": [],
                "scaffold": [],
                "smiles": [],
                "reward": [],
                "penalized_reward": [],
            }
        )
        # Add oracle components' raw value and reward to the oracle history DataFrame
        for oracle in self.oracle:
            self.oracle_history[f"{oracle.name}_raw_values"] = []
            self.oracle_history[f"{oracle.name}_reward"] = []
            # Composite oracles (e.g. GEAM) expose sub-component columns
            if hasattr(oracle, "get_component_breakdown"):
                for key in oracle.get_component_breakdown():
                    self.oracle_history[f"{oracle.name}_{key}"] = []

        # Track how many times the same SMILES is sampled
        self.repeated_sampled_smiles = {}
        self.repeated_hallucinated_smiles = {}

        # HELM validity checking: Load valid monomers if available
        # Extract building blocks path from constructed oracle components (e.g., from TANGO's enforced_structures)
        self.valid_monomers = None
        if HELM_CONVERSION_AVAILABLE:
            enforced_structures_path = None
            for oracle_component in self.oracle:
                if (
                    oracle_component.specific_parameters
                    and "enforced_structures" in oracle_component.specific_parameters
                ):
                    enforced_structures_path = oracle_component.specific_parameters[
                        "enforced_structures"
                    ]
                    break

            if enforced_structures_path is not None:
                try:
                    from pathlib import Path

                    bb_path = Path(enforced_structures_path)
                    if bb_path.exists():
                        self.valid_monomers = load_valid_monomers(bb_path)
                        logging.info(
                            f"[Oracle] Loaded {len(self.valid_monomers)} valid monomers for HELM validation from {enforced_structures_path}"
                        )
                except Exception as e:
                    logging.warning(
                        f"[Oracle] Could not load valid monomers for HELM validation: {e}"
                    )

    def __call__(
        self,
        sequences: np.ndarray[str],  # Can be SMILES or HELM sequences
        diversity_filter: DiversityFilter,
        is_hallucinated_batch: bool = False,
        step: Optional[int] = None,  # RL step for logging
    ) -> Tuple[np.ndarray[str], np.ndarray[float]]:
        """
        The Oracle is called at every generation epoch and performs the following:
            1. Calls each oracle component in the Oracle
            2. Aggregates the oracle feedback into a single scalar reward
            3. Penalizes the reward based on the Diversity Filter
            4. Updates the Diversity Filter
            5. Updates the Oracle History which tracks oracle calls, rewards, and penalized rewards
            6. Updates the Oracle Cache to store the results of previous oracle calls

        Returns the original sequences (SMILES or HELM) and the penalized rewards.
        """
        # 1. Convert HELM to SMILES if needed, validate, deduplicate
        smiles, is_helm_input, smiles_to_helm_map = self._convert_helm_sequences(
            sequences
        )

        # 2. Check oracle cache for repeat SMILES
        repeat_smiles, cached_rewards, new_smiles = self.rewards_from_oracle_cache(
            smiles, is_hallucinated_batch
        )
        if is_helm_input:
            for s in repeat_smiles:
                if (
                    s in self.cache
                    and isinstance(self.cache[s], dict)
                    and "helm" in self.cache[s]
                ):
                    smiles_to_helm_map[s] = self.cache[s]["helm"]

        # 3. Evaluate new molecules through oracle components
        oracle_components_df = pd.DataFrame()
        if len(new_smiles) > 0:
            new_mols = np.vectorize(Chem.MolFromSmiles)(new_smiles)
            new_smiles, new_mols = self.execute_preliminary_check(
                new_smiles, new_mols
            )

            if len(new_smiles) > 0:
                new_smiles, new_mols, oracle_components_df, rewards = (
                    self._evaluate_oracle_components(new_mols, new_smiles)
                )
                aggregated_rewards = self._aggregate_and_enrich(
                    new_smiles, oracle_components_df, rewards, is_helm_input
                )
            else:
                aggregated_rewards = np.array([0.0])
        else:
            aggregated_rewards = np.array([0.0])

        # 4. Finalize: penalize, update history/cache, return
        return self._finalize(
            repeat_smiles,
            cached_rewards,
            new_smiles,
            aggregated_rewards,
            oracle_components_df,
            diversity_filter,
            is_helm_input,
            smiles_to_helm_map,
            step,
        )

    # ── Extracted private methods ────────────────────────────────────────────

    def _convert_helm_sequences(
        self, sequences: np.ndarray[str]
    ) -> Tuple[np.ndarray[str], bool, dict]:
        """
        Convert HELM sequences to SMILES, validate, and deduplicate.

        Returns:
            smiles: Validated and deduplicated SMILES array
            is_helm_input: Whether the input was HELM notation
            smiles_to_helm_map: Reverse mapping from SMILES to HELM
        """
        helm_to_smiles_map = {}
        is_helm_input = False
        smiles = sequences

        if len(sequences) > 0:
            logging.info(
                f"[Oracle] Received {len(sequences)} sequences. First sequence: {sequences[0][:80]}..."
            )
            logging.info(
                f"[Oracle] HELM_CONVERSION_AVAILABLE: {HELM_CONVERSION_AVAILABLE}"
            )
            if HELM_CONVERSION_AVAILABLE:
                is_helm = is_helm_notation(sequences[0])
                logging.info(
                    f"[Oracle] is_helm_notation(sequences[0]): {is_helm}"
                )

        if HELM_CONVERSION_AVAILABLE and len(sequences) > 0:
            if is_helm_notation(sequences[0]):
                is_helm_input = True
                converted_smiles = []
                logging.info(
                    f"[Oracle] Processing {len(sequences)} HELM sequences"
                )
                logging.info(
                    f"[Oracle] valid_monomers loaded: {self.valid_monomers is not None}"
                )
                if self.valid_monomers is not None:
                    logging.info(
                        f"[Oracle] Number of valid monomers: {len(self.valid_monomers)}"
                    )
                else:
                    logging.warning(
                        "[Oracle] valid_monomers is None! HELM validation will be skipped!"
                    )

                # For HELM: Validate using monomer library check
                invalid_count = 0
                conversion_fail_count = 0
                conversion_empty_count = 0

                # Show first 3 HELM sequences for debugging
                for i, helm_seq in enumerate(sequences[:3]):
                    logging.info(f"[Oracle] Sample HELM #{i+1}: {helm_seq}")

                for helm_seq in sequences:
                    # Check if HELM is valid (all monomers in library)
                    if self.valid_monomers is not None:
                        is_valid, valid_count, total_count = check_helm_validness(
                            helm_seq, self.valid_monomers
                        )
                        if not is_valid:
                            invalid_count += 1
                            if invalid_count <= 5:
                                logging.warning(
                                    f"[Oracle] Invalid HELM (unknown monomers {valid_count}/{total_count}): {helm_seq[:80]}..."
                                )
                            continue

                    # Convert valid HELM to SMILES
                    try:
                        smi = convert_helm_to_smiles(helm_seq)
                        if smi:
                            converted_smiles.append(smi)
                            helm_to_smiles_map[helm_seq] = smi
                            if (
                                conversion_fail_count
                                + conversion_empty_count
                                + len(converted_smiles)
                                <= 3
                            ):
                                logging.info(
                                    f"[Oracle] Converted HELM to SMILES: {helm_seq[:60]}... -> {smi[:60]}..."
                                )
                        else:
                            conversion_empty_count += 1
                            if conversion_empty_count <= 5:
                                logging.warning(
                                    f"[Oracle] HELM->SMILES returned None: {helm_seq[:80]}..."
                                )
                    except Exception as e:
                        conversion_fail_count += 1
                        if conversion_fail_count <= 3:
                            logging.warning(
                                f"[Oracle] HELM->SMILES exception for {helm_seq[:60]}...: {type(e).__name__}: {e}"
                            )
                        continue

                logging.info("[Oracle] HELM Conversion Summary:")
                logging.info(
                    f"[Oracle]   - Input HELM sequences: {len(sequences)}"
                )
                logging.info(
                    f"[Oracle]   - Successfully converted: {len(converted_smiles)}"
                )
                logging.info(
                    f"[Oracle]   - Invalid (unknown monomers): {invalid_count}"
                )
                logging.info(
                    f"[Oracle]   - Conversion returned None: {conversion_empty_count}"
                )
                logging.info(
                    f"[Oracle]   - Conversion exceptions: {conversion_fail_count}"
                )
                smiles = np.array(converted_smiles)

        # Validate SMILES with RDKit (HELM already validated via monomer check)
        if not is_helm_input:
            valid_mask = np.array(
                [Chem.MolFromSmiles(s) is not None for s in smiles]
            )
            smiles = smiles[valid_mask]

        # De-duplicate
        smiles = self.de_duplicate_smiles(smiles)

        # Build reverse mapping: SMILES → HELM
        smiles_to_helm_map = (
            {v: k for k, v in helm_to_smiles_map.items()} if is_helm_input else {}
        )

        return smiles, is_helm_input, smiles_to_helm_map

    def _evaluate_oracle_components(
        self, new_mols: np.ndarray[Mol], new_smiles: np.ndarray[str]
    ) -> Tuple[np.ndarray[str], np.ndarray[Mol], pd.DataFrame, np.ndarray]:
        """
        Call each oracle component and filter out molecules with failed calculations.

        Returns:
            new_smiles: Filtered SMILES (NaN rows removed)
            new_mols: Filtered Mol objects
            oracle_components_df: DataFrame with raw values and rewards per component
            rewards: 2D array (n_oracles × n_molecules) of component rewards
        """
        oracle_components_df = pd.DataFrame()
        rewards = np.empty((len(self.oracle), len(new_mols)))
        for idx, oracle in enumerate(self.oracle):
            raw_property_values, component_rewards = oracle.calculate_reward(
                new_mols, self.calls
            )
            oracle_components_df[f"{oracle.name}_raw_values"] = raw_property_values
            oracle_components_df[f"{oracle.name}_reward"] = component_rewards
            rewards[idx] = component_rewards
            # Composite oracles expose sub-component columns
            if hasattr(oracle, "get_component_breakdown"):
                for key, values in oracle.get_component_breakdown().items():
                    oracle_components_df[f"{oracle.name}_{key}"] = values

        # Filter out molecules with NaN values (failed oracle calculations)
        raw_value_cols = [
            col
            for col in oracle_components_df.columns
            if col.endswith("_raw_values")
        ]
        if raw_value_cols:
            valid_mask = (
                ~oracle_components_df[raw_value_cols].isna().any(axis=1).values
            )
            if not np.all(valid_mask):
                num_filtered = np.sum(~valid_mask)
                logging.warning(
                    f"[Oracle] Filtering out {num_filtered} molecules with failed oracle calculations (NaN values)"
                )
                new_smiles = new_smiles[valid_mask]
                new_mols = new_mols[valid_mask]
                oracle_components_df = oracle_components_df[
                    valid_mask
                ].reset_index(drop=True)
                rewards = rewards[:, valid_mask]

        return new_smiles, new_mols, oracle_components_df, rewards

    def _aggregate_and_enrich(
        self,
        new_smiles: np.ndarray[str],
        oracle_components_df: pd.DataFrame,
        rewards: np.ndarray,
        is_helm_input: bool,
    ) -> np.ndarray[float]:
        """
        Aggregate component rewards and enrich oracle_components_df with metadata.

        Returns:
            aggregated_rewards: Single reward per molecule
        """
        if len(new_smiles) == 0:
            return np.array([])

        aggregated_rewards = self.aggregator(rewards, self.oracle_weights)

        # Extract and load metadata from oracles (ETL pattern)
        metadata = extract_oracle_metadata(self.oracle)
        expected_len = len(new_smiles)

        if (
            metadata["bb_sequences"]
            and len(metadata["bb_sequences"]) == expected_len
        ):
            oracle_components_df["bb_sequences"] = metadata["bb_sequences"]

        if (
            metadata["synth_factors"]
            and len(metadata["synth_factors"]) == expected_len
            and "synthesizability_factor_raw_values"
            not in oracle_components_df.columns
        ):
            oracle_components_df["synthesizability_factor_raw_values"] = (
                metadata["synth_factors"]
            )

        # For HELM inputs: valid HELM sequences are synthesizable by construction
        if is_helm_input:
            oracle_components_df["synthesizability_factor_raw_values"] = 1.0

        return aggregated_rewards

    def _finalize(
        self,
        repeat_smiles: np.ndarray[str],
        cached_rewards: np.ndarray[float],
        new_smiles: np.ndarray[str],
        aggregated_rewards: np.ndarray[float],
        oracle_components_df: pd.DataFrame,
        diversity_filter: DiversityFilter,
        is_helm_input: bool,
        smiles_to_helm_map: dict,
        step: Optional[int],
    ) -> Tuple[np.ndarray[str], np.ndarray[float]]:
        """
        Finalize oracle evaluation: penalize, update history/cache, and return.

        Returns:
            return_sequences: Original sequences (HELM if input was HELM)
            penalized_all_rewards: Diversity-penalized rewards
        """
        # Increment oracle calls
        self.calls += len(new_smiles)

        # Concatenate repeated and new SMILES
        all_smiles = np.concatenate([repeat_smiles, new_smiles])
        all_rewards = np.concatenate([cached_rewards, aggregated_rewards])

        # Penalize rewards and update the Diversity Filter in a single pass
        n_repeat = len(repeat_smiles)
        if is_helm_input:
            all_helm = np.array(
                [smiles_to_helm_map.get(s, s) for s in all_smiles]
            )
            penalized_all_rewards = diversity_filter.penalize_and_update(
                all_helm, all_rewards
            )
        else:
            penalized_all_rewards = diversity_filter.penalize_and_update(
                all_smiles, all_rewards
            )
        penalized_new_rewards = penalized_all_rewards[n_repeat:]

        # Update Oracle History
        if len(new_smiles) > 0:
            if not self.allow_oracle_repeats:
                oracle_history_smiles = new_smiles
                oracle_history_rewards = aggregated_rewards
                oracle_history_penalized_rewards = penalized_new_rewards
            elif self.allow_oracle_repeats:
                oracle_history_smiles = all_smiles
                oracle_history_rewards = all_rewards
                oracle_history_penalized_rewards = penalized_all_rewards

            if is_helm_input:
                original_sequences = np.array(
                    [smiles_to_helm_map.get(s, s) for s in oracle_history_smiles]
                )
            else:
                original_sequences = oracle_history_smiles

            self.update_oracle_history(
                smiles=oracle_history_smiles,
                original_sequences=original_sequences,
                scaffolds=np.vectorize(get_bemis_murcko_scaffold)(
                    oracle_history_smiles
                ),
                rewards=oracle_history_rewards,
                penalized_rewards=oracle_history_penalized_rewards,
                oracle_components_df=oracle_components_df,
                step=step,
            )

        # Update Oracle Cache
        if is_helm_input:
            all_helm = np.array(
                [smiles_to_helm_map.get(s, s) for s in all_smiles]
            )
            self.update_oracle_cache(
                all_smiles, penalized_all_rewards, helm=all_helm
            )
        else:
            self.update_oracle_cache(all_smiles, penalized_all_rewards)

        # Return original sequences (HELM if input was HELM)
        if is_helm_input:
            return_sequences = np.array(
                [smiles_to_helm_map.get(s, s) for s in all_smiles]
            )
            return return_sequences, penalized_all_rewards
        else:
            return all_smiles, penalized_all_rewards

    # ── Oracle construction and caching ──────────────────────────────────────

    def construct_oracle(
        self, oracle_components: List[OracleComponentParameters]
    ) -> List[OracleComponent]:
        """
        Construct the oracle function which can be composed of multiple individual oracle components.

        Components that need a SynthesizabilityChecker (PeptideLength, SynthesizabilityFactor)
        share a single instance: the first component creates it, subsequent ones receive it via DI.
        The shared instance is also exposed as self.synth_checker for the DiversityFilter.
        """
        oracle = []
        shared_synth_checker = None
        for component in oracle_components:
            oracle_component = construct_oracle_component(
                OracleComponentParameters(**component),
                synth_checker=shared_synth_checker,
            )
            oracle.append(oracle_component)
            # Capture the checker from the first component that creates one
            if shared_synth_checker is None and getattr(oracle_component, 'synth_checker', None) is not None:
                shared_synth_checker = oracle_component.synth_checker

        self.synth_checker = shared_synth_checker
        return oracle

    def rewards_from_oracle_cache(
        self, smiles: np.ndarray[str], is_hallucinated_batch: bool
    ) -> Tuple[np.ndarray[str], np.ndarray[float], np.ndarray[str]]:
        """
        Checks if there are any Cached rewards in a sampled batch of SMILES.
        Also updates trackers for repeated SMILES. If Oracle repeats are permitted, directly return.
        """
        if not self.allow_oracle_repeats:
            # Canonicalize the SMILES before checking Cache
            canonical_smiles = canonicalize_smiles_batch(smiles)
            repeat_indices = []
            cached_rewards = []
            for idx, s in enumerate(canonical_smiles):
                if s in self.cache:
                    repeat_indices.append(idx)
                    # Cache format: {'helm': helm_seq, 'rewards': [r1, r2, ...]} or just [r1, r2, ...]
                    cache_entry = self.cache[s]
                    if isinstance(cache_entry, dict):
                        # Dict format with helm
                        cached_rewards.append(np.mean(cache_entry["rewards"]))
                    else:
                        # Simple list format (for SMILES-only runs)
                        cached_rewards.append(np.mean(cache_entry))

            if len(repeat_indices) != 0:
                # Track the repeated SMILES and their rewards
                repeated_smiles = smiles[repeat_indices]
                for idx, s in enumerate(repeated_smiles):
                    if is_hallucinated_batch:
                        if s not in self.repeated_hallucinated_smiles:
                            self.repeated_hallucinated_smiles[s] = (
                                1,
                                cached_rewards[idx],
                            )
                        else:
                            self.repeated_hallucinated_smiles[s] = (
                                self.repeated_hallucinated_smiles[s][0] + 1,
                                cached_rewards[idx],
                            )
                    else:
                        if s not in self.repeated_sampled_smiles:
                            self.repeated_sampled_smiles[s] = (1, cached_rewards[idx])
                        else:
                            self.repeated_sampled_smiles[s] = (
                                self.repeated_sampled_smiles[s][0] + 1,
                                cached_rewards[idx],
                            )

                return (
                    repeated_smiles,
                    np.array(cached_rewards),
                    np.delete(smiles, repeat_indices),
                )
            else:
                return np.array([]), np.array([]), smiles

        else:
            return np.array([]), np.array([]), smiles

    def update_oracle_cache(
        self,
        smiles: np.ndarray[str],
        rewards: np.ndarray[float],
        helm: np.ndarray[str] = None,
    ) -> None:
        """
        Updates the Oracle Cache to store the results of previous oracle calls.

        Cache format:
        - HELM runs: {canonical_smiles: {'helm': helm_seq, 'rewards': [r1, r2, ...]}}
        - SMILES runs: {canonical_smiles: [r1, r2, ...]}

        Args:
            smiles: SMILES strings
            rewards: Corresponding rewards
            helm: HELM sequences (optional, for HELM-based runs)
        """
        # Canonicalize the SMILES before adding to Cache
        canonical_smiles = canonicalize_smiles_batch(smiles)
        for idx, (s, r) in enumerate(zip(canonical_smiles, rewards)):
            helm_seq = helm[idx] if helm is not None and idx < len(helm) else None

            # If the same SMILES is sampled, all rewards are tracked for two reasons:
            #   1. Potential stochasticity in the oracle feedback
            #   2. Penalized rewards (by the Diversity Filter) should be reflected so the Agent is steered away from these scaffolds
            if s not in self.cache:
                # New entry
                if helm_seq is not None:
                    self.cache[s] = {"helm": helm_seq, "rewards": [r]}
                else:
                    self.cache[s] = [r]
            elif r == 0.0:
                # Zero reward: reset cache entry
                if helm_seq is not None:
                    self.cache[s] = {"helm": helm_seq, "rewards": [r]}
                else:
                    self.cache[s] = [r]
            else:
                # Update existing entry
                if isinstance(self.cache[s], dict):
                    # Dict format: append reward
                    self.cache[s]["rewards"].append(r)
                else:
                    # List format: append reward
                    self.cache[s].append(r)

    def execute_preliminary_check(
        self, smiles: np.ndarray[str], mols: np.ndarray[Mol]
    ) -> Tuple[np.ndarray[str], np.ndarray[Mol]]:
        """
        Executes a preliminary check (if applicable). Each oracle component has a preliminary_check flag that can be set to True.
        Components set to True will be executed first to check that the molecule satisfies that component based on a reward threshold.
        If the molecule does not satisfy the threshold, it is removed from the sampled batch.
        """
        # FIXME: Set a threshold for each component. If not using Step transformation, rewards are not necessarily 0
        THRESHOLD = 0.05

        if len(self.preliminary_oracles) > 0:
            filtered_indices = []
            # TODO: Vectorize the batch of Mols
            for idx, mol in enumerate(mols):
                for oracle in self.preliminary_oracles:
                    _, reward = oracle.calculate_reward(np.array([mol]), self.calls)
                    if reward < THRESHOLD:
                        # If the reward is below the threshold for at least one oracle component, add the index to the filtered indices
                        filtered_indices.append(idx)
                        break

            return np.delete(smiles, filtered_indices), np.delete(
                mols, filtered_indices
            )

        else:
            return smiles, mols

    def update_oracle_history(
        self,
        scaffolds: np.ndarray[str],
        smiles: np.ndarray[str],
        original_sequences: np.ndarray[str],  # HELM if input was HELM, else SMILES
        rewards: np.ndarray[float],
        penalized_rewards: np.ndarray[float],
        oracle_components_df: pd.DataFrame,
        step: Optional[int] = None,
    ) -> None:
        """
        This method performs the following on every generation epoch:
        1. Increments the number of oracle calls so far
        2. Updates the Oracle History that tracks the generative sampling as a function of oracle calls

        # NOTE: If self.allow_oracle_repeats = True, the Oracle History tracks every single SMILES generated and not just the unique set.
        #       This can be useful to interrogate the stochasticity of the oracle to inform downstream molecule prioritization.
        """
        # Extract commonly tracked metrics from oracle_components_df
        # Note: All oracle component raw values are preserved via concatenation below

        # Synthesizability factor
        synth_factor = oracle_components_df.get(
            "synthesizability_factor_raw_values", pd.Series([None] * len(smiles))
        )

        # Calculate is_synthesizable (s >= 0.999 means fully synthesizable)
        if "synthesizability_factor_raw_values" in oracle_components_df.columns:
            is_synth = synth_factor >= 0.999
        else:
            is_synth = pd.Series([False] * len(smiles))

        # Determine which sequence to store:
        # - bb_sequences: comma-separated building blocks from SMILES decomposition (Saturn/SMILES models)
        # - original_sequences: HELM notation for Neptune models, or SMILES if no decomposition available
        if "bb_sequences" in oracle_components_df.columns:
            sequences = oracle_components_df["bb_sequences"]
        else:
            sequences = original_sequences

        # Track generated SMILES + reward as a function of oracle calls
        df = pd.DataFrame(
            {
                "step": np.full_like(
                    smiles, step if step is not None else -1, dtype=int
                ),
                "oracle_calls": np.full_like(smiles, self.calls),
                "smiles": smiles,
                "is_valid": np.full(
                    len(smiles), True
                ),  # All entries here passed RDKit validation
                "synthesizability_factor": synth_factor,
                "is_synthesizable": is_synth,
                "sequence": sequences,  # HELM or comma-separated BB list
                "reward": rewards,
                "scaffold": scaffolds,
                "penalized_reward": penalized_rewards,
            }
        )

        # Concatenate with oracle component columns (which includes ALL raw values and rewards)
        # This preserves all objective metrics: logp, qed, sa_score, docking_score, etc.
        df = pd.concat([df, oracle_components_df], axis=1)

        self.oracle_history = (
            pd.concat([self.oracle_history, df]) if len(self.oracle_history) > 0 else df
        )

    @staticmethod
    def de_duplicate_smiles(smiles: np.ndarray[str]) -> np.ndarray[str]:
        """
        De-duplicate a batch of SMILES.
        """
        smiles_copy = deepcopy(smiles)
        smiles_copy = canonicalize_smiles_batch(smiles_copy)
        _, unique_indices = np.unique(smiles_copy, return_index=True)
        unique_indices.sort()

        return smiles[unique_indices]

    def budget_exceeded(self) -> bool:
        """Check if the oracle budget has been exceeded."""
        return self.calls >= self.budget

    def write_out_oracle_history(self, path: str) -> None:
        """Write out the oracle history as a CSV."""
        import csv

        self.oracle_history.to_csv(
            os.path.join(path, "oracle_history.csv"),
            index=False,
            quoting=csv.QUOTE_NONNUMERIC,  # Quote all non-numeric fields (prevents comma issues in SMILES)
        )

    def write_out_repeat_history(self, path: str) -> None:
        """Write out the repeated SMILES histories as JSON."""
        # FIXME: Reproduce json dump error
        try:
            with open(
                os.path.join(path, "repeated_sampled_smiles_history.json"), "w"
            ) as f:
                json.dump(self.repeated_sampled_smiles, f, indent=2)
            with open(
                os.path.join(path, "repeated_hallucinated_smiles_history.json"), "w"
            ) as f:
                json.dump(self.repeated_hallucinated_smiles, f, indent=2)
        except Exception:
            print("Failed to write out repeat histories.")
