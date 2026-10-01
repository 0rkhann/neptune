"""
Order-independent diversity fingerprints for Neptune (HELM sequences).

Two granularity levels, both order-independent:

1. **Subclass fingerprint** (get_subclass_fingerprint): Strips numeric IDs so
   monomers sharing subclass codes are grouped (l_n_38 and l_n_22 → l_n).
   Most aggressive diversity — forces truly different subclass compositions.

2. **Monomer fingerprint** (get_monomer_fingerprint): Keeps the full monomer name
   including numeric ID (l_n_38 ≠ l_n_22). Less aggressive — different specific
   building blocks within the same subclass still count as diverse.

Core logic (monomer_to_subclass, get_subclass_fingerprint, get_monomer_fingerprint)
lives in GenAI4Peptidomimetic_native (src.genai_utils.helm) and is re-exported here.

IMPORTANT: Requires GenAI4Peptidomimetic_native to be installed as a package.
Install with: pip install -e . from a checkout of that repository.
"""

# Import and re-export from GenAI4Peptidomimetic_native
try:
    from src.genai_utils.helm import (
        monomer_to_subclass,
        get_subclass_fingerprint,
        get_monomer_fingerprint,
        subclass_fingerprint_from_monomers,
        monomer_fingerprint_from_monomers,
    )

    HELM_DIVERSITY_AVAILABLE = True
except ImportError as e:
    HELM_DIVERSITY_AVAILABLE = False
    raise ImportError(
        "GenAI4Peptidomimetic_native package is required for HELM diversity filtering. "
        f"Original error: {e}\n"
        "Install it with 'pip install -e .' from a checkout of "
        "https://github.com/schwallergroup/GenAI4Peptidomimetic"
    )