"""
Diversity utilities for Saturn.
Implements bigram scaffold extraction for diversity filtering with HELM sequences and SMILES.

IMPORTANT: Requires GenAI4Peptidomimetic_native to be installed as a package.
Install with: pip install -e /path/to/GenAI4Peptidomimetic_native
"""

from typing import Optional
import logging

# Import HELM utilities from GenAI4Peptidomimetic_native package
try:
    from src.genai_utils.helm import decode_helm_string, is_helm_notation
    from src.genai_utils.synthesizability import SynthesizabilityChecker

    HELM_DIVERSITY_AVAILABLE = True
except ImportError as e:
    HELM_DIVERSITY_AVAILABLE = False
    raise ImportError(
        "GenAI4Peptidomimetic_native package is required for HELM diversity filtering. "
        f"Original error: {e}\n"
        "Please install it with: pip install -e /path/to/GenAI4Peptidomimetic_native"
    )

# Global synthesizability checker instance (lazy initialization)
_synth_checker: Optional[SynthesizabilityChecker] = None


def _get_synth_checker() -> SynthesizabilityChecker:
    """Get or create synthesizability checker instance."""
    global _synth_checker
    if _synth_checker is None:
        _synth_checker = SynthesizabilityChecker()
    return _synth_checker


def get_bigram_scaffold(sequence: str) -> str:
    """
    Extract bigram fingerprint from HELM sequence or SMILES for diversity filter.

    Supports both HELM and SMILES:
    - HELM: Directly extracts monomers from HELM notation
    - SMILES: Decomposes SMILES into building blocks first, then extracts bigrams

    For 2-4 mer catalysts, bigrams capture:
    - Sequence order (Ala-Gly ≠ Gly-Ala)
    - Catalytic dyad/triad geometry
    - Local pairwise interactions
    - N/C-cap effects (important for catalysis!)

    Args:
        sequence: HELM notation string or SMILES string

    Returns:
        String representation of sorted bigrams, or empty string if extraction fails

    Example (HELM):
        Input:  "PEPTIDE1{[na5].[l_a_1].[l_g_1].[l_s_1].[ca2]}$$$$"
        Output: "[('l_a_1', 'l_g_1'), ('l_g_1', 'l_s_1'), ('l_s_1', 'ca2'), ('na5', 'l_a_1')]"

    Example (SMILES):
        Input:  "CC(=O)N[C@@H](C)C(=O)N"
        Output: "[('l_a_1', 'cp1'), ('na1', 'l_a_1')]"  # After decomposition

    Note: Includes N/C-caps in bigrams because they affect:
        - Catalytic mechanism (N-cap can be part of active site)
        - Electronic properties
        - Solubility and stability
        - Steric effects near catalytic residues
    """
    monomers = []

    # Check if it's HELM notation
    if is_helm_notation(sequence):
        # Extract monomers from HELM
        try:
            decoded = decode_helm_string(sequence, return_short_names=True)
            monomers = [name for name, block_type in decoded["monomers"]]
        except Exception as e:
            logging.warning(f"Failed to decode HELM string: {sequence[:60]}... -> {e}")
            return ""
    else:
        # It's SMILES - decompose into building blocks first
        try:
            checker = _get_synth_checker()
            sequence_list, valid_count, total_count = checker.decompose_partial(
                sequence
            )

            if not sequence_list:
                # Decomposition failed, return empty string
                return ""

            # Filter out "UNKNOWN" markers - they can't form meaningful bigrams
            # But keep them in sequence for position information
            monomers = sequence_list

            # If all are UNKNOWN, can't form bigrams
            if all(m == "UNKNOWN" for m in monomers):
                return ""

        except Exception as e:
            logging.warning(f"Failed to decompose SMILES: {sequence[:60]}... -> {e}")
            return ""

    if not monomers:
        return ""

    if len(monomers) < 2:
        # Single monomer, no bigrams possible
        return str(monomers)

    # Extract bigrams (consecutive pairs) - INCLUDING caps
    # Note: UNKNOWN markers are included in bigrams to maintain sequence position
    bigrams = []
    for i in range(len(monomers) - 1):
        bigram = (monomers[i], monomers[i + 1])
        bigrams.append(bigram)

    # Sort for consistency and convert to string
    return str(sorted(bigrams))
