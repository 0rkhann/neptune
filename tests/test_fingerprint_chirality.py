"""Guard every Morgan fingerprint call against silently discarding stereochemistry.

RDKit's `useChirality` defaults to False, so a fingerprint call that omits it
treats D and L forms of a residue as the same molecule. For peptidomimetics that
erases the distinction the whole project rests on, and nothing errors when it
happens.
"""

import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs import TanimotoSimilarity

from utils import chemistry_utils

REPO_ROOT = Path(__file__).resolve().parent.parent

# L- and D-alanine: identical constitution, opposite configuration.
L_ALA = "C[C@@H](N)C(=O)O"
D_ALA = "C[C@H](N)C(=O)O"


def test_enantiomers_differ_when_chirality_is_on():
    """The premise: chirality-aware fingerprints separate D from L."""
    fps = [
        GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(s), radius=3, nBits=2048, useChirality=True
        )
        for s in (L_ALA, D_ALA)
    ]
    assert TanimotoSimilarity(*fps) < 1.0


def test_enantiomers_collide_when_chirality_is_off():
    """The failure mode, pinned so the reason for the flag stays visible."""
    fps = [
        GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), radius=3, nBits=2048)
        for s in (L_ALA, D_ALA)
    ]
    assert TanimotoSimilarity(*fps) == 1.0


def test_shared_helper_keeps_chirality():
    assert (
        TanimotoSimilarity(
            chemistry_utils.construct_morgan_fingerprint(L_ALA),
            chemistry_utils.construct_morgan_fingerprint(D_ALA),
        )
        < 1.0
    )


def test_no_morgan_call_omits_usechirality():
    """Every GetMorganFingerprintAsBitVect call in the repo must set it.

    Catches the defect returning anywhere, not just where it was found.
    """
    offenders = []
    for path in REPO_ROOT.rglob("*.py"):
        if "__pycache__" in path.parts or path.name == Path(__file__).name:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for call in re.finditer(r"GetMorganFingerprintAsBitVect\(([^)]*)\)", text):
            if "useChirality" not in call.group(1):
                line = text[: call.start()].count("\n") + 1
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{line}")
    assert not offenders, (
        "Morgan fingerprint calls without useChirality (D and L collapse): "
        + ", ".join(offenders)
    )
