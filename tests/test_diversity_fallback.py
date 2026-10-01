"""Guard the diversity-fingerprint path against silently degrading to a no-op.

Background: when GenAI4Peptidomimetic_native is not importable,
utils/diversity_utils.py falls back to Bemis-Murcko scaffolds. A linear peptide
has no ring system, so its Murcko scaffold is the empty string, which collapses
every linear peptide into one diversity bucket and disables the filter without
raising anything. These tests pin both halves of that: the empty-scaffold fact
that makes the fallback dangerous, and the fact that we are not silently running
in fallback mode.
"""

import pytest

from utils import chemistry_utils
from utils import diversity_utils

# Three distinct linear tripeptides. None contains a ring.
LINEAR_PEPTIDES = [
    "CC(N)C(=O)NCC(=O)NC(CO)C(=O)O",                    # Ala-Gly-Ser
    "CC(C)C(N)C(=O)NC(CC(C)C)C(=O)NC(C(C)O)C(=O)O",     # Val-Leu-Thr
    "NCC(=O)NCC(=O)NCC(=O)O",                           # Gly-Gly-Gly
]

CYCLIC_PEPTIDE = "O=C1NC(Cc2ccccc2)C(=O)NCC(=O)NC1"


def test_murcko_scaffold_is_empty_for_linear_peptides():
    """The premise: Murcko cannot distinguish ring-free peptides."""
    scaffolds = [chemistry_utils.get_bemis_murcko_scaffold(s) for s in LINEAR_PEPTIDES]
    assert scaffolds == ["", "", ""], (
        "Expected every linear peptide to yield an empty Murcko scaffold. "
        f"Got {scaffolds!r}"
    )


def test_murcko_scaffold_is_non_empty_for_cyclic_peptide():
    """Contrast case: the fallback is only degenerate for ring-free molecules."""
    assert chemistry_utils.get_bemis_murcko_scaffold(CYCLIC_PEPTIDE)


def test_helm_backend_is_available():
    """Fail loudly if the run would silently use the Murcko fallback.

    HELM_AVAILABLE is False only when GenAI4Peptidomimetic_native is missing.
    Any peptide result produced in that state has no working diversity filter.
    """
    assert diversity_utils.HELM_AVAILABLE, (
        "GenAI4Peptidomimetic_native is not importable, so diversity filtering "
        "would fall back to Bemis-Murcko scaffolds. Install it with "
        "'pip install -e .' from a checkout of that repository."
    )


@pytest.mark.skipif(
    diversity_utils.HELM_AVAILABLE,
    reason="Fallback is not installed; nothing to check.",
)
def test_fallback_rejects_ring_free_molecules():
    """When the fallback *is* active, it must raise rather than return ''."""
    with pytest.raises(RuntimeError, match="empty scaffold"):
        diversity_utils.get_diversity_fingerprint(
            LINEAR_PEPTIDES[0], None, None, None
        )
