"""
Shared diversity fingerprint utilities for DiversityFilter and ReplayBuffer.

Core logic lives in GenAI4Peptidomimetic_native (src.genai_utils) and is
re-exported here, following the same pattern as utils/subclass_fingerprint.py
and utils/helm/__init__.py.

When GenAI4Peptidomimetic_native is not installed, only Bemis-Murcko scaffold
mode is available (SMILES-only Saturn runs).
"""

from utils import chemistry_utils

# ── Re-export from GenAI4Peptidomimetic_native ──────────────────────────────

try:
    from src.genai_utils.helm import (
        is_helm_notation,
        to_smiles,
        get_subclass_fingerprint,
        get_monomer_fingerprint,
        subclass_fingerprint_from_monomers,
        monomer_fingerprint_from_monomers,
    )
    from src.genai_utils.diversity import get_diversity_fingerprint

    HELM_AVAILABLE = True
except ImportError as _helm_import_error:
    import warnings

    HELM_AVAILABLE = False
    is_helm_notation = None
    to_smiles = None
    get_subclass_fingerprint = None
    get_monomer_fingerprint = None
    subclass_fingerprint_from_monomers = None
    monomer_fingerprint_from_monomers = None

    # UserWarning, not ImportWarning: Python ignores ImportWarning by default, so
    # the degraded mode would otherwise be invisible.
    warnings.warn(
        "GenAI4Peptidomimetic_native is not importable "
        f"({_helm_import_error}). Diversity filtering falls back to "
        "Bemis-Murcko scaffolds, which is valid only for SMILES-only runs over "
        "ring-containing molecules. Install it with 'pip install -e .' from a "
        "checkout of https://github.com/schwallergroup/GenAI4Peptidomimetic",
        UserWarning,
        stacklevel=2,
    )

    # Fallback: Murcko-only mode for SMILES-only Saturn runs
    def get_diversity_fingerprint(sequence, fingerprint_fn, monomer_fingerprint_fn,
                                  synth_checker, log_prefix=""):
        scaffold = chemistry_utils.get_bemis_murcko_scaffold(sequence)
        if not scaffold:
            # A linear peptide has no ring system, so its Murcko scaffold is the
            # empty string. Every such molecule would collide into one bucket and
            # silently disable the diversity filter. Fail instead.
            raise RuntimeError(
                f"{log_prefix} Bemis-Murcko fallback produced an empty scaffold for "
                f"{sequence!r}. This molecule has no ring system, so the fallback "
                "cannot distinguish it from any other ring-free molecule. "
                "GenAI4Peptidomimetic_native is required for peptide diversity "
                "filtering; install it with 'pip install -e .' from a checkout of "
                "https://github.com/schwallergroup/GenAI4Peptidomimetic"
            )
        return scaffold

try:
    from src.genai_utils.sidechain_fingerprint import get_sidechain_fingerprint
    from src.genai_utils.synthesizability import SynthesizabilityChecker

    SIDECHAIN_FP_AVAILABLE = True
except ImportError:
    SIDECHAIN_FP_AVAILABLE = False
    get_sidechain_fingerprint = None
    SynthesizabilityChecker = None

try:
    from src.genai_utils.map4c_fingerprint import (
        get_map4c_fingerprint,
        map4c_jaccard_similarity,
        MAP4C_AVAILABLE,
    )
except ImportError:
    MAP4C_AVAILABLE = False
    get_map4c_fingerprint = None
    map4c_jaccard_similarity = None

