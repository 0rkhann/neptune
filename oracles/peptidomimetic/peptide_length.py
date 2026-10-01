"""
Peptide Length Oracle — binary constraint on backbone building block count.

Returns 1 if the number of backbone building blocks (amino acids + superbases,
excluding N-caps and C-caps) is ≤ a user-specified maximum, 0 otherwise.

Uses the SynthesizabilityChecker to decompose the molecule into its
constituent building blocks and then counts only the non-cap residues.

Config example
--------------
{
    "name": "peptide_length",
    "weight": 1.0,
    "reward_shaping_function_parameters": {
        "transformation_function": "no_transformation",
        "parameters": {}
    },
    "specific_parameters": {
        "max_length": 3,
        "enforced_structures": "/path/to/bb.csv"
    }
}
"""

import numpy as np
from rdkit import Chem

from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters

try:
    from src.genai_utils.synthesizability import SynthesizabilityChecker
    from src.genai_utils.building_blocks import default_building_blocks_path
except ImportError as e:
    raise ImportError(
        "Failed to import GenAI4Peptidomimetic_native. Install it with "
        "'pip install -e .' from a checkout of "
        "https://github.com/schwallergroup/GenAI4Peptidomimetic"
    ) from e


class PeptideLength(OracleComponent):
    """
    Binary oracle that enforces a maximum number of backbone residues.

    Counts only amino acid and superbase building blocks (not caps).
    Returns 1.0 if count ≤ max_length, 0.0 otherwise.
    """

    def __init__(self, parameters: OracleComponentParameters, synth_checker=None):
        super().__init__(parameters)

        self.max_length = int(parameters.specific_parameters.get("max_length", 3))

        enforced_structures = parameters.specific_parameters.get("enforced_structures")
        if enforced_structures is None:
            enforced_structures = str(default_building_blocks_path())

        print("[PeptideLength] Initializing with:")
        print(f"  - max_length: {self.max_length}")
        print(f"  - Building blocks: {enforced_structures}")

        if synth_checker is not None:
            self.synth_checker = synth_checker
        else:
            self.synth_checker = SynthesizabilityChecker(enforced_structures)

        # Build sets of cap names for fast lookup
        self._ncap_names = set(self.synth_checker.ncaps.keys())
        self._ccap_names = set(self.synth_checker.ccaps.keys())

    def __call__(self, mols: np.ndarray) -> np.ndarray:
        scores = np.zeros(len(mols), dtype=np.float64)
        for idx, mol in enumerate(mols):
            if mol is None:
                continue
            scores[idx] = self._score_molecule(mol)
        return scores

    def _score_molecule(self, mol: Chem.Mol) -> float:
        """Return 1.0 if backbone residue count ≤ max_length, else 0.0."""
        smiles = Chem.MolToSmiles(mol)
        sequence, valid_count, total_count = (
            self.synth_checker.decompose_partial(smiles)
        )

        if total_count == 0:
            return 0.0

        # Count only backbone residues (AA + superbase), not caps
        backbone_count = sum(
            1 for name in sequence
            if name not in self._ncap_names and name not in self._ccap_names
            and name != "UNKNOWN"
        )

        if backbone_count == 0:
            return 0.0

        return 1.0 if backbone_count <= self.max_length else 0.0
