"""
Some code is based on the implementation from https://github.com/MolecularAI/Reinvent

Supports:
- Saturn (SMILES): Bemis-Murcko scaffolds (IdenticalMurckoScaffold)
- Neptune (HELM): Bigram scaffolds (IdenticalBigramScaffold)
"""

import numpy as np
import logging
from utils import chemistry_utils
from diversity_filter.dataclass import DiversityFilterParameters

# Import HELM utilities for bigram scaffold support (Neptune only)
try:
    from utils.bigram_diversity_utils import get_bigram_scaffold
    from utils.helm import is_helm_notation

    HELM_DIVERSITY_AVAILABLE = True
except ImportError:
    HELM_DIVERSITY_AVAILABLE = False
    get_bigram_scaffold = None
    is_helm_notation = None


class DiversityFilter:
    """
    Implements Diversity Filter as described in the paper:
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00473-0

    Supports:
    - Saturn (SMILES): Uses Bemis-Murcko scaffolds (IdenticalMurckoScaffold)
    - Neptune (HELM): Uses bigram scaffolds (IdenticalBigramScaffold)
    
    NOTE: Bigram scaffolds are ONLY supported for HELM sequences (Neptune).
          Saturn uses Bemis-Murcko scaffolds exclusively.
    """

    def __init__(self, parameters: DiversityFilterParameters):
        self.parameters = parameters
        self.name = parameters.name
        # Track the number of times a given scaffold has been generated
        self.bucket_history = dict()
        self.bucket_size = parameters.bucket_size

        # Determine scaffold type based on name
        self.use_bigram_scaffold = (
            HELM_DIVERSITY_AVAILABLE and self.name == "IdenticalBigramScaffold"
        )

        if self.use_bigram_scaffold:
            logging.info(
                "[DiversityFilter] Using bigram scaffolds for HELM sequences (Neptune)"
            )
        else:
            logging.info("[DiversityFilter] Using Bemis-Murcko scaffolds for SMILES (Saturn)")

    def _get_scaffold(self, sequence: str) -> str:
        """
        Get scaffold for a sequence.
        
        - HELM (Neptune with IdenticalBigramScaffold): bigram scaffolds
        - SMILES (Saturn): Bemis-Murcko scaffolds
        """
        if self.use_bigram_scaffold:
            # Bigram scaffolds are ONLY for HELM (Neptune)
            if not is_helm_notation(sequence):
                raise ValueError(
                    f"IdenticalBigramScaffold diversity filter requires HELM sequences (Neptune), "
                    f"but received SMILES: {sequence[:60]}... "
                    f"Use IdenticalMurckoScaffold for Saturn (SMILES)."
                )
            scaffold = get_bigram_scaffold(sequence)
            if not scaffold:
                logging.error(
                    f"[DiversityFilter] Failed to extract bigram scaffold from HELM: {sequence[:60]}..."
                )
                raise ValueError(f"Failed to extract bigram scaffold from HELM: {sequence[:60]}...")
            return scaffold
        else:
            # Bemis-Murcko scaffolds for SMILES (Saturn)
            return chemistry_utils.get_bemis_murcko_scaffold(sequence)

    def update(
        self, sequences: np.ndarray[str]  # SMILES for Saturn, HELM for Neptune
    ) -> None:
        """
        Update the bucket history based on the sampled (or hallucinated) batch.
        """
        if len(sequences) == 0:
            return

        scaffolds = [self._get_scaffold(seq) for seq in sequences]

        for scaf in scaffolds:
            if scaf and scaf.strip():
                if scaf in self.bucket_history:
                    self.bucket_history[scaf] += 1
                else:
                    self.bucket_history[scaf] = 1

        # Log diversity filter statistics periodically
        if len(self.bucket_history) % 50 == 0:
            logging.debug(
                f"[DiversityFilter] {len(self.bucket_history)} unique scaffolds tracked"
            )

    def penalize_reward(
        self,
        sequences: np.ndarray[str],  # SMILES for Saturn, HELM for Neptune
        rewards: np.ndarray[float],
    ) -> np.ndarray[float]:
        """
        Penalize sampled (or hallucinated) sequences based on the bucket history.
        """
        if len(sequences) == 0:
            return np.array([])

        scaffolds = [self._get_scaffold(seq) for seq in sequences]

        penalized_rewards = []
        num_penalized = 0
        for idx, scaf in enumerate(scaffolds):
            if (
                scaf in self.bucket_history
                and self.bucket_history[scaf] > self.bucket_size
            ):
                # Scaffold seen too many times: penalize
                penalized_rewards.append(0.0)
                num_penalized += 1
            else:
                # Scaffold not seen or within bucket size: keep original reward
                penalized_rewards.append(rewards[idx])

        if num_penalized > 0:
            logging.debug(
                f"[DiversityFilter] Penalized {num_penalized}/{len(sequences)} sequences"
            )

        return np.array(penalized_rewards)
