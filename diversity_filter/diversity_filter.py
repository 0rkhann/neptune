"""
Some code is based on the implementation from https://github.com/MolecularAI/Reinvent

Supports all combinations of tokenization × diversity filter:
- IdenticalMurckoScaffold: Bemis-Murcko scaffolds (SMILES input)
- SubclassFingerprint: groups by subclass (l_n_38, l_n_22 → l_n)
- MonomerFingerprint: keeps full monomer identity (l_n_38 ≠ l_n_22)
- SidechainFingerprint: ECFP4 whole-molecule fingerprint with Tanimoto clustering
- MAP4CFingerprint: MAP4C (MinHashed Atom-Pair with Chirality) with Jaccard clustering

Cross-tokenization is handled transparently:
- SMILES + SubclassFingerprint/MonomerFingerprint: decomposes SMILES → monomer IDs via SynthesizabilityChecker
- HELM + SidechainFingerprint/MAP4CFingerprint: converts HELM → SMILES via convert_helm_to_smiles
"""

import numpy as np
import logging
from diversity_filter.dataclass import DiversityFilterParameters

from utils.diversity_utils import (
    to_smiles,
    get_diversity_fingerprint,
    SIDECHAIN_FP_AVAILABLE,
    SynthesizabilityChecker,
    get_sidechain_fingerprint,
    get_subclass_fingerprint,
    get_monomer_fingerprint,
    subclass_fingerprint_from_monomers,
    monomer_fingerprint_from_monomers,
    MAP4C_AVAILABLE,
    get_map4c_fingerprint,
    map4c_jaccard_similarity,
)


class ScaffoldMemory:
    """
    Similarity-based cluster memory for fingerprint diversity filtering.

    Groups fingerprints into clusters using a configurable similarity function —
    the first fingerprint in each cluster becomes the centroid (single-linkage clustering).

    A cluster is "saturated" when its count exceeds bucket_size, triggering reward
    penalization (reward = 0.0) for subsequent molecules with similar patterns.

    Supports:
    - SidechainFingerprint (ECFP4): ExplicitBitVect + DataStructs.TanimotoSimilarity
    - MAP4CFingerprint (MAP4C): np.ndarray (uint64 MinHash) + Jaccard similarity
    """

    def __init__(self, bucket_size: int = 10, tanimoto_threshold: float = 0.65,
                 similarity_fn=None):
        self.bucket_size = bucket_size
        self.tanimoto_threshold = tanimoto_threshold
        self._similarity_fn = similarity_fn
        # Each cluster: (centroid_fp, count)
        self._clusters = []

    def _find_cluster(self, fp) -> int:
        """Find the index of the matching cluster (similarity > threshold), or -1."""
        if self._similarity_fn is None:
            from rdkit import DataStructs
            self._similarity_fn = DataStructs.TanimotoSimilarity

        for i, (centroid, _count) in enumerate(self._clusters):
            if self._similarity_fn(fp, centroid) > self.tanimoto_threshold:
                return i
        return -1

    def is_saturated(self, fp) -> bool:
        """Check if the cluster matching this fingerprint has exceeded bucket_size."""
        idx = self._find_cluster(fp)
        if idx < 0:
            return False
        return self._clusters[idx][1] > self.bucket_size

    def update(self, fp) -> None:
        """Add fingerprint to its matching cluster, or create a new cluster."""
        idx = self._find_cluster(fp)
        if idx >= 0:
            centroid, count = self._clusters[idx]
            self._clusters[idx] = (centroid, count + 1)
        else:
            self._clusters.append((fp, 1))

    @property
    def num_clusters(self) -> int:
        """Number of distinct sidechain pattern clusters."""
        return len(self._clusters)


class DiversityFilter:
    """
    Implements Diversity Filter as described in the paper:
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00473-0

    Supports:
    - Saturn (SMILES): Uses Bemis-Murcko scaffolds (IdenticalMurckoScaffold)
    - Neptune (HELM): Order-independent fingerprints:
      - SubclassFingerprint: groups monomers by subclass identity
      - MonomerFingerprint: keeps exact monomer identity
    - Atom-level (SMILES): ECFP4 whole-molecule fingerprint (SidechainFingerprint)
      - Standard Morgan radius=2 fingerprint with chirality
      - Uses Tanimoto similarity clustering (ScaffoldMemory) instead of exact-match buckets
    - All tokenizations: MAP4C MinHash fingerprint (MAP4CFingerprint)
      - Chirality-aware atom-pair substructure fingerprint
      - Uses Jaccard similarity clustering (ScaffoldMemory) instead of exact-match buckets
    """

    def __init__(self, parameters: DiversityFilterParameters, synth_checker=None):
        self.parameters = parameters
        self.name = parameters.name
        self.bucket_size = parameters.bucket_size

        # Fingerprint-based clustering modes
        self.use_sidechain_fingerprint = False
        self._synth_checker = None
        self._scaffold_memory = None

        if self.name == "SidechainFingerprint":
            if not SIDECHAIN_FP_AVAILABLE:
                raise ImportError(
                    "SidechainFingerprint requires GenAI4Peptidomimetic_native. "
                    "Install with: pip install -e <GenAI4Peptidomimetic_native_dir>"
                )
            self._scaffold_memory = ScaffoldMemory(
                bucket_size=self.bucket_size,
                tanimoto_threshold=parameters.tanimoto_threshold,
            )
            self.use_sidechain_fingerprint = True
            self._fp_mode = "sidechain"
            self.bucket_history = None  # Not used in sidechain mode
            logging.info(
                f"[DiversityFilter] Using ECFP4 fingerprints "
                f"(tanimoto_threshold={parameters.tanimoto_threshold})"
            )
            self._fingerprint_fn = None
            self._monomer_fingerprint_fn = None
            # Store synth_checker if provided (for configure_replay_buffer)
            if synth_checker is not None:
                self._synth_checker = synth_checker
            return

        if self.name == "MAP4CFingerprint":
            if not MAP4C_AVAILABLE:
                raise ImportError(
                    "MAP4CFingerprint requires mapchiral. "
                    "Install from: ../mapchiral/"
                )
            self._scaffold_memory = ScaffoldMemory(
                bucket_size=self.bucket_size,
                tanimoto_threshold=parameters.tanimoto_threshold,
                similarity_fn=map4c_jaccard_similarity,
            )
            self.use_sidechain_fingerprint = True  # reuse same code path
            self._fp_mode = "map4c"
            self.bucket_history = None
            logging.info(
                f"[DiversityFilter] Using MAP4C fingerprints "
                f"(jaccard_threshold={parameters.tanimoto_threshold})"
            )
            self._fingerprint_fn = None
            self._monomer_fingerprint_fn = None
            if synth_checker is not None:
                self._synth_checker = synth_checker
            return

        # Track the number of times a given fingerprint has been generated
        self.bucket_history = dict()

        if self.name == "SubclassFingerprint":
            self._fingerprint_fn = get_subclass_fingerprint
            self._monomer_fingerprint_fn = subclass_fingerprint_from_monomers
            logging.info(
                "[DiversityFilter] Using subclass fingerprints"
            )
        elif self.name == "MonomerFingerprint":
            self._fingerprint_fn = get_monomer_fingerprint
            self._monomer_fingerprint_fn = monomer_fingerprint_from_monomers
            logging.info(
                "[DiversityFilter] Using monomer fingerprints"
            )
        else:
            self._fingerprint_fn = None
            self._monomer_fingerprint_fn = None
            logging.info("[DiversityFilter] Using Bemis-Murcko scaffolds for SMILES (Saturn)")

        # Store synth_checker for cross-tokenization (SMILES → monomer decomposition)
        if synth_checker is not None:
            self._synth_checker = synth_checker
        elif self._fingerprint_fn is not None and parameters.bb_csv_path and SIDECHAIN_FP_AVAILABLE:
            self._synth_checker = SynthesizabilityChecker(parameters.bb_csv_path)
            logging.info(
                "[DiversityFilter] SynthesizabilityChecker created for cross-tokenization support"
            )

    def _compute_similarity_fp(self, smiles: str):
        """Compute fingerprint for similarity-based filtering (SidechainFP or MAP4C)."""
        if getattr(self, '_fp_mode', None) == "map4c":
            return get_map4c_fingerprint(smiles, self._synth_checker)
        return get_sidechain_fingerprint(smiles, self._synth_checker)

    def configure_replay_buffer(self, replay_buffer) -> None:
        """Inject shared SynthesizabilityChecker and similarity config into replay buffer."""
        if self._synth_checker is not None:
            replay_buffer.set_synth_checker(
                self._synth_checker,
                tanimoto_threshold=self.parameters.tanimoto_threshold,
            )

    def update(
        self, sequences: np.ndarray[str]  # SMILES for Saturn, HELM for Neptune
    ) -> None:
        """
        Update the bucket history based on the sampled (or hallucinated) batch.
        """
        if len(sequences) == 0:
            return

        if self.use_sidechain_fingerprint:
            for seq in sequences:
                smiles = to_smiles(seq)
                fp = self._compute_similarity_fp(smiles) if smiles else None
                if fp is not None:
                    self._scaffold_memory.update(fp)
            if self._scaffold_memory.num_clusters % 50 == 0:
                logging.debug(
                    f"[DiversityFilter] {self._scaffold_memory.num_clusters} sidechain clusters tracked"
                )
            return

        scaffolds = [
            get_diversity_fingerprint(
                seq, self._fingerprint_fn, self._monomer_fingerprint_fn,
                self._synth_checker, "[DiversityFilter]"
            )
            for seq in sequences
        ]

        for scaf in scaffolds:
            if scaf and scaf.strip():
                if scaf in self.bucket_history:
                    self.bucket_history[scaf] += 1
                else:
                    self.bucket_history[scaf] = 1

        # Log diversity filter statistics periodically
        if len(self.bucket_history) % 50 == 0:
            logging.debug(
                f"[DiversityFilter] {len(self.bucket_history)} unique fingerprints tracked"
            )

    def penalize_reward(
        self,
        sequences: np.ndarray[str],  # SMILES for Saturn, HELM for Neptune
        rewards: np.ndarray[float],
    ) -> np.ndarray[float]:
        """
        Penalize sampled (or hallucinated) sequences based on the bucket history.
        Does NOT update the bucket history — call update() separately if needed.
        """
        if len(sequences) == 0:
            return np.array([])

        if self.use_sidechain_fingerprint:
            penalized_rewards = []
            num_penalized = 0
            for idx, seq in enumerate(sequences):
                smiles = to_smiles(seq)
                fp = self._compute_similarity_fp(smiles) if smiles else None
                if fp is not None and self._scaffold_memory.is_saturated(fp):
                    penalized_rewards.append(0.0)
                    num_penalized += 1
                else:
                    penalized_rewards.append(rewards[idx])
            if num_penalized > 0:
                logging.debug(
                    f"[DiversityFilter] Penalized {num_penalized}/{len(sequences)} sequences (sidechain saturation)"
                )
            return np.array(penalized_rewards)

        scaffolds = [
            get_diversity_fingerprint(
                seq, self._fingerprint_fn, self._monomer_fingerprint_fn,
                self._synth_checker, "[DiversityFilter]"
            )
            for seq in sequences
        ]

        penalized_rewards = []
        num_penalized = 0
        for idx, scaf in enumerate(scaffolds):
            if (
                scaf in self.bucket_history
                and self.bucket_history[scaf] > self.bucket_size
            ):
                # Fingerprint seen too many times: penalize
                penalized_rewards.append(0.0)
                num_penalized += 1
            else:
                # Fingerprint not seen or within bucket size: keep original reward
                penalized_rewards.append(rewards[idx])

        if num_penalized > 0:
            logging.debug(
                f"[DiversityFilter] Penalized {num_penalized}/{len(sequences)} sequences"
            )

        return np.array(penalized_rewards)

    def penalize_and_update(
        self,
        sequences: np.ndarray[str],
        rewards: np.ndarray[float],
    ) -> np.ndarray[float]:
        """
        Compute fingerprints once, penalize all against pre-update state, then update.

        Semantically equivalent to penalize_reward() followed by update(),
        but computes each fingerprint/scaffold only once.
        """
        if len(sequences) == 0:
            return np.array([])

        if self.use_sidechain_fingerprint:
            # Phase 1: compute all fingerprints
            fingerprints = []
            for seq in sequences:
                smiles = to_smiles(seq)
                fp = self._compute_similarity_fp(smiles) if smiles else None
                fingerprints.append(fp)

            # Phase 2: penalize (all checked against current state, before any updates)
            penalized_rewards = []
            num_penalized = 0
            for idx, fp in enumerate(fingerprints):
                if fp is not None and self._scaffold_memory.is_saturated(fp):
                    penalized_rewards.append(0.0)
                    num_penalized += 1
                else:
                    penalized_rewards.append(rewards[idx])

            # Phase 3: update memory with pre-computed fingerprints
            for fp in fingerprints:
                if fp is not None:
                    self._scaffold_memory.update(fp)

            if num_penalized > 0:
                logging.debug(
                    f"[DiversityFilter] Penalized {num_penalized}/{len(sequences)} sequences (sidechain saturation)"
                )
            if self._scaffold_memory.num_clusters % 50 == 0:
                logging.debug(
                    f"[DiversityFilter] {self._scaffold_memory.num_clusters} sidechain clusters tracked"
                )
            return np.array(penalized_rewards)

        # Subclass/Monomer/Murcko: compute scaffolds once
        scaffolds = [
            get_diversity_fingerprint(
                seq, self._fingerprint_fn, self._monomer_fingerprint_fn,
                self._synth_checker, "[DiversityFilter]"
            )
            for seq in sequences
        ]

        # Penalize (all checked against current bucket_history)
        penalized_rewards = []
        num_penalized = 0
        for idx, scaf in enumerate(scaffolds):
            if (
                scaf in self.bucket_history
                and self.bucket_history[scaf] > self.bucket_size
            ):
                penalized_rewards.append(0.0)
                num_penalized += 1
            else:
                penalized_rewards.append(rewards[idx])

        # Update bucket_history with pre-computed scaffolds
        for scaf in scaffolds:
            if scaf and scaf.strip():
                if scaf in self.bucket_history:
                    self.bucket_history[scaf] += 1
                else:
                    self.bucket_history[scaf] = 1

        if num_penalized > 0:
            logging.debug(
                f"[DiversityFilter] Penalized {num_penalized}/{len(sequences)} sequences"
            )
        if self.bucket_history and len(self.bucket_history) % 50 == 0:
            logging.debug(
                f"[DiversityFilter] {len(self.bucket_history)} unique fingerprints tracked"
            )

        return np.array(penalized_rewards)
