"""
Synthesizability Factor Oracle (s-factor).

Returns the fraction of building blocks that are in a provided library:
    s = (# of building blocks in library) / (total # of building blocks)

IMPORTANT: This oracle requires GenAI4Peptidomimetic_native to be installed as a package.
Install with: pip install -e . from a checkout of that repository.
"""

import numpy as np
from rdkit import Chem

from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters

# Import synthesizability checker from GenAI4Peptidomimetic_native package
try:
    from src.genai_utils.synthesizability import SynthesizabilityChecker
    from src.genai_utils.building_blocks import default_building_blocks_path
except ImportError as e:
    raise ImportError(
        "Failed to import GenAI4Peptidomimetic_native. Install it with "
        "'pip install -e .' from a checkout of "
        "https://github.com/schwallergroup/GenAI4Peptidomimetic"
    ) from e


class SynthesizabilityFactor(OracleComponent):
    """
    Synthesizability Factor Oracle (s-factor).

    Returns the fraction of building blocks that can be found in a provided library:
        s = (# of building blocks in library) / (total # of building blocks)

    This is designed to be composed with other oracles using Saturn's aggregators.

    Example config (combining with slogp using product aggregator):
    {
        "aggregator": "product",
        "components": [
            {
                "name": "slogp",
                "weight": 1.0,
                "reward_shaping_function_parameters": {
                    "transformation_function": "sigmoid",
                    "parameters": {"low": -2.0, "high": 18.0, "k": 0.8}
                }
            },
            {
                "name": "synthesizability_factor",
                "weight": 1.0,
                "specific_parameters": {
                    "enforced_structures": "/path/to/bb.csv",
                    "power": 1.0
                }
            }
        ]
    }

    Config parameters:
        enforced_structures: Path to building blocks CSV (default: uses GenAI4Peptidomimetic_native)
        power: Exponent for the s-factor (default: 1.0) -> returns s^power
    """

    def __init__(self, parameters: OracleComponentParameters, synth_checker=None):
        """
        Initialize synthesizability factor oracle.

        Args:
            parameters: Oracle component parameters (from Saturn's config)
            synth_checker: Optional shared SynthesizabilityChecker instance
        """
        super().__init__(parameters)

        # Get building blocks library path from config or use default
        enforced_structures = parameters.specific_parameters.get("enforced_structures")
        if enforced_structures is None:
            enforced_structures = str(default_building_blocks_path())

        # Power parameter (exponent for s-factor)
        self.power = parameters.specific_parameters.get("power", 1.0)

        print("[SynthesizabilityFactor] Initializing with:")
        print(f"  - Building blocks library: {enforced_structures}")
        print(f"  - Power (exponent): {self.power}")

        # Use shared checker if provided, otherwise create a new one
        if synth_checker is not None:
            self.synth_checker = synth_checker
        else:
            self.synth_checker = SynthesizabilityChecker(enforced_structures)

        # Store sequences for tracking (updated in __call__)
        self.last_sequences = []

    def __call__(self, mols: np.ndarray) -> np.ndarray:
        """
        Calculate synthesizability factor for each molecule.

        Args:
            mols: Array of RDKit Mol objects

        Returns:
            Array of s-factors in [0, 1], where s = (valid_bb_count / total_bb_count)^power
        """
        scores = np.zeros(len(mols), dtype=np.float32)
        self.last_sequences = []  # Reset sequences for this batch

        for idx, mol in enumerate(mols):
            if mol is None:
                scores[idx] = 0.0
                self.last_sequences.append("")
                continue

            try:
                # Get SMILES and decompose to building blocks
                smiles = Chem.MolToSmiles(mol)
                sequence, valid_count, total_count = (
                    self.synth_checker.decompose_partial(smiles)
                )

                # Store sequence as comma-separated string
                if sequence:
                    self.last_sequences.append(",".join(sequence))
                else:
                    self.last_sequences.append("")

                # Calculate s-factor
                if total_count > 0:
                    s = valid_count / total_count
                    scores[idx] = s**self.power
                else:
                    scores[idx] = 0.0

            except Exception as e:
                print(f"[SynthesizabilityFactor] Error scoring molecule: {e}")
                scores[idx] = 0.0
                self.last_sequences.append("")

        return scores
