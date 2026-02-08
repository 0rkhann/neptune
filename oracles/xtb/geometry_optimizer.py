from typing import Tuple, List
import os
import tempfile
import shutil
import subprocess
from rdkit.Chem import Mol
from rdkit import Chem
from rdkit.Chem import AllChem


class GeometryOptimizer:
    def __init__(self):
        # TODO: can include xTB specific parameters here
        pass

    def optimize_geometry(self, mol: Mol) -> Tuple[str, str, List[str]]:
        """
        Optimize geometry using RDKit ETKDGv3 + MMFF94, then single-point GFN2-xTB.

        Steps:
        1. Add hydrogens to molecule
        2. Generate 3D conformer using ETKDGv3
        3. Optimize with MMFF94 for 50 iterations
        4. Run single-point GFN2-xTB calculation (fast)

        Returns:
            1. the path to the temp directory containing the geometry and xTB output
            2. the path to the geometry file
            3. the xTB output
        """
        # Make temp folder to generate and store geometries
        temp_dir = tempfile.mkdtemp()

        # Add hydrogens
        mol_with_h = Chem.AddHs(mol)

        # Generate 3D conformer using ETKDGv3
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # For reproducibility
        embed_result = AllChem.EmbedMolecule(mol_with_h, params)

        if embed_result == -1:
            # Embedding failed, try without ETKDG
            embed_result = AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
            if embed_result == -1:
                raise ValueError("Failed to generate 3D conformer")

        # Optimize with MMFF94 for 50 iterations
        AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=50)

        # Write to XYZ file for xTB
        xyz_path = os.path.join(temp_dir, "rdkit_conformer.xyz")
        AllChem.MolToXYZFile(mol_with_h, xyz_path)

        # Run single-point GFN2-xTB calculation (--sp flag for single-point, no optimization)
        xtb_output = subprocess.run(
            [
                "xtb",
                xyz_path,
                "--gfn",
                "2",
                "--sp",  # Single-point calculation (fast!)
                "--namespace",
                f"{temp_dir}/temp",
            ],
            capture_output=True,
        )

        return temp_dir, xyz_path, xtb_output.stdout.decode().splitlines()

    @staticmethod
    def clean_up_temp_dir(path: str):
        shutil.rmtree(path)
