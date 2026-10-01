"""
Topological Potential Field (TPF) Oracle — Level 2 Cooperative Reward.

Calculates a Continuous Cooperative Reward (CCR) based on the topological
(graph-distance) arrangement of H-Bond Donors (HBDs) relative to the
superbase center in a peptidomimetic catalyst.

The reward favours placing multiple HBDs at a "Goldilocks" distance of
~10 bonds from the superbase imine nitrogen, using a Gaussian potential:

    C_i = exp( -(d_i - μ)² / (2σ²) )
    P   = Σ C_i

Molecules with ≠ 1 superbase centre receive a reward of 0.

Config example
--------------
{
    "name": "topological_potential_field",
    "weight": 1.0,
    "reward_shaping_function_parameters": {
        "transformation_function": "no_transformation",
        "parameters": {}
    },
    "specific_parameters": {
        "mu": 10,
        "sigma": 2
    }
}
"""

import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops

from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters

# ---------------------------------------------------------------------------
# SMARTS definitions
# ---------------------------------------------------------------------------

# Superbase imine nitrogen: P=N (phosphazene) or C=N in guanidine context.
# TMG_Dap superbases contain C(=N)N(C)C — the imine N is =N.
# Phosphazene superbases contain P=N.
# We match the nitrogen atom that is double-bonded to C or P.
SUPERBASE_SMARTS = [
    "[#7;X2]=[#15]",  # P=N  (phosphazene)
    "[#7;X2]=[#6](-[#7])-[#7]",  # C(=N)(N)N  guanidine imine N (TMG)
]

# ---------------------------------------------------------------------------
# HBD SMARTS — comprehensive set covering all building-block chemistries.
#
# Sources of HBDs in peptidomimetics:
#   • Backbone amide N-H  (from peptide/amide coupling)
#   • Urea N-H            (from isocyanate N-caps)
#   • Thiourea N-H        (from isothiocyanate N-caps)
#   • Sulfonamide N-H     (from sulfonylating N-caps)
#   • Free amine N-H      (uncapped or side-chain)
#   • Hydroxyl O-H        (serine/threonine side-chains, hydroxamic acid C-cap)
#   • Aromatic N-H        (indole, pyrrole — aromatic backbone AAs)
#   • Hydrazide N-H       (hydrazine C-cap)
#   • Carbamate N-H       (carbamylating N-caps)
# ---------------------------------------------------------------------------
HBD_SMARTS = [
    "[N;!a;!H0]",  # any non-aromatic nitrogen with ≥1 H (amide, urea, thiourea,
    #   sulfonamide, amine, hydrazide, carbamate NH)
    "[n;H1]",  # aromatic NH (indole, pyrrole backbone)
    "[O;H1]",  # hydroxyl (Ser/Thr side-chain, hydroxamic acid)
]


class TopologicalPotentialField(OracleComponent):
    """
    Topological Potential Field oracle for peptidomimetic catalyst optimisation.

    Rewards molecules that place H-Bond Donors at an optimal graph-distance
    from a single superbase centre.  Molecules with ≠1 superbase receive 0.
    """

    def __init__(self, parameters: OracleComponentParameters):
        super().__init__(parameters)

        # Gaussian parameters
        self.mu = float(parameters.specific_parameters.get("mu", 10))
        self.sigma = float(parameters.specific_parameters.get("sigma", 2))

        print("[TopologicalPotentialField] Initializing with:")
        print(f"  - mu (optimal distance): {self.mu}")
        print(f"  - sigma (width):         {self.sigma}")

        # Pre-compile SMARTS
        self._superbase_queries = [Chem.MolFromSmarts(s) for s in SUPERBASE_SMARTS]
        self._hbd_queries = [Chem.MolFromSmarts(s) for s in HBD_SMARTS]

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def __call__(self, mols: np.ndarray) -> np.ndarray:
        scores = np.zeros(len(mols), dtype=np.float64)
        for idx, mol in enumerate(mols):
            if mol is None:
                continue
            scores[idx] = self._score_molecule(mol)
        return scores

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------
    def _get_superbase_atoms(self, mol: Chem.Mol) -> list[int]:
        """Return atom indices of superbase imine nitrogens."""
        hits: set[int] = set()
        for query in self._superbase_queries:
            for match in mol.GetSubstructMatches(query):
                hits.add(match[0])  # the nitrogen is first atom in each SMARTS
        return list(hits)

    def _get_hbd_atoms(self, mol: Chem.Mol) -> list[int]:
        """Return unique atom indices of all HBD atoms."""
        hits: set[int] = set()
        for query in self._hbd_queries:
            for match in mol.GetSubstructMatches(query):
                hits.add(match[0])
        return list(hits)

    def _score_molecule(self, mol: Chem.Mol) -> float:
        """
        P = Σ exp( -(d_i - μ)² / (2σ²) )

        Returns 0 if the molecule does not contain exactly 1 superbase centre.
        """
        # --- superbase check ---
        sb_atoms = self._get_superbase_atoms(mol)
        if len(sb_atoms) != 1:
            return 0.0

        sb_idx = sb_atoms[0]

        # --- HBD identification ---
        hbd_atoms = self._get_hbd_atoms(mol)
        if not hbd_atoms:
            return 0.0

        # Remove superbase N itself from HBD list (it is not a donor)
        hbd_atoms = [a for a in hbd_atoms if a != sb_idx]
        if not hbd_atoms:
            return 0.0

        # --- topological distances ---
        dist_matrix = rdmolops.GetDistanceMatrix(mol)
        two_sigma_sq = 2.0 * self.sigma * self.sigma

        potential = 0.0
        for hbd_idx in hbd_atoms:
            d_i = dist_matrix[sb_idx][hbd_idx]
            potential += math.exp(-((d_i - self.mu) ** 2) / two_sigma_sq)

        return potential
