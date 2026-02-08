import numpy as np
from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters
from rdkit import Chem
from rdkit.Chem import Mol

from utils.chemistry_utils import construct_morgan_fingerprints_batch
from oracles.synthesizability.utils.utils import (
    extract_functional_groups,
    get_node_reward,
)
from oracles.synthesizability.s_factor import SynthesizabilityIntegration


class Tango(OracleComponent):
    """
    TANGO (TANimoto Group Overlap) Oracle — Per-Fragment Scoring.

    Implements TANGO as described in:
    "It Takes Two to Tango: Directly Optimizing for Constrained
    Synthesizability in Generative Molecular Design" (Guo & Schwaller, 2024).

    The molecule is first decomposed into building block fragments using
    retrosynthetic analysis. TANGO is then computed for each fragment
    against the enforced building block library. The final score is the
    mean of per-fragment TANGO scores.

    For fully synthesizable molecules (all fragments match library BBs),
    TANGO → 1.0. For partially synthesizable molecules, TANGO reflects
    the average fragment-level similarity to the library.

    """

    def __init__(self, parameters: OracleComponentParameters):
        super().__init__(parameters)

        # --- Reward type ---
        self.reward_type = parameters.specific_parameters.get("reward_type", None)
        assert self.reward_type in [
            "tango_fg",
            "tango_fms",
            "tango_all",
        ], "Reward Type must be one of tango_fg, tango_fms, or tango_all"

        # --- TANGO weights ---
        self.tango_weights = self.parameters.specific_parameters.get(
            "tango_weights", None
        )
        assert self.tango_weights is not None, "Please provide TANGO weights."

        # --- Enforced structures (BB library for similarity computation) ---
        self.enforced_structures = parameters.specific_parameters.get(
            "enforced_structures"
        )
        assert (
            self.enforced_structures is not None
        ), "Enforced Structures must be provided"
        assert type(self.enforced_structures) in [
            list,
            str,
        ], "Enforced Structures must be a list of SMILES or a path to a file"

        if isinstance(self.enforced_structures, list):
            self.enforced_structures_smiles = self.enforced_structures
        else:
            file_path = self.enforced_structures
            if file_path.endswith(".csv"):
                import pandas as pd

                df = pd.read_csv(file_path)
                if "smiles" in df.columns:
                    self.enforced_structures_smiles = df["smiles"].dropna().tolist()
                elif "SMILES" in df.columns:
                    self.enforced_structures_smiles = df["SMILES"].dropna().tolist()
                else:
                    raise ValueError(
                        f"CSV file {file_path} must contain a 'smiles' or 'SMILES' column. "
                        f"Found columns: {list(df.columns)}"
                    )
            else:
                self.enforced_structures_smiles = [
                    s.strip() for s in open(file_path, "r").readlines() if s.strip()
                ]

        # --- Retrosynthetic decomposition (required for per-fragment TANGO) ---
        # TANGO must score individual fragments, not the full assembled molecule.
        # The SynthesizabilityChecker decomposes the molecule into BB fragments.
        synth_config = parameters.specific_parameters.get("synthesizability_config", {})
        synth_config["enforced_structures"] = self.enforced_structures
        self.synth_integration = SynthesizabilityIntegration(synth_config)

        # --- Extend enforced structures with caps and superbases ---
        # The CSV only contains amino acids. Caps and superbases are loaded
        # from token_definitions.py via BuildingBlockLibrary. We need them
        # in the enforced set so get_node_reward finds exact matches for ALL
        # library building blocks (not just AAs).
        bb_lib = self.synth_integration.checker.bb_library
        caps_and_sb_smiles = []
        for info in bb_lib.ncaps.values():
            smi = info.get("smiles")
            if smi and smi not in self.enforced_structures_smiles:
                caps_and_sb_smiles.append(smi)
        for info in bb_lib.ccaps.values():
            smi = info.get("smiles")
            if smi and smi not in self.enforced_structures_smiles:
                caps_and_sb_smiles.append(smi)
        for info in bb_lib.superbases.values():
            smi = info.get("smiles")
            if smi and smi not in self.enforced_structures_smiles:
                caps_and_sb_smiles.append(smi)

        if caps_and_sb_smiles:
            self.enforced_structures_smiles.extend(caps_and_sb_smiles)
            print(
                f"[Tango] Extended enforced structures with {len(caps_and_sb_smiles)} "
                f"caps/superbases (total: {len(self.enforced_structures_smiles)})"
            )

        self.enforced_structures_fps = construct_morgan_fingerprints_batch(
            self.enforced_structures_smiles
        )
        self.enforced_structures_functional_groups = extract_functional_groups(
            self.enforced_structures_smiles
        )
        print("[Tango] Per-fragment TANGO scoring via retrosynthetic decomposition")

    def __call__(self, mols: np.ndarray[Mol]) -> np.ndarray[float]:
        """
        Compute per-fragment TANGO reward for each molecule.

        For each molecule:
        1. Decompose into building block fragments via retrosynthesis
        2. Compute TANGO for each fragment against the BB library
        3. Return the mean of per-fragment TANGO scores

        Fully synthesizable molecules → each fragment matches a BB → TANGO = 1.0
        Non-synthesizable molecules → UNKNOWN fragments score low → TANGO < 1.0
        """
        rewards = []

        for query_mol in mols:
            query_smiles = Chem.MolToSmiles(query_mol, canonical=True)

            # Decompose molecule into fragments with their SMILES
            sequence, fragment_smiles, valid_count, total_count = (
                self.synth_integration.checker.decompose_partial(
                    query_smiles, unknown_marker="UNKNOWN", return_details=True
                )
            )

            # Per-fragment TANGO scoring
            # Each fragment is scored against the full BB library
            # (AAs + caps + superbases). Matched BBs → 1.0, UNKNOWN → similarity.
            node_rewards = []
            for i, (bb_name, frag_smi) in enumerate(zip(sequence, fragment_smiles)):
                if Chem.MolFromSmiles(frag_smi) is None:
                    print(
                        f"[Tango] WARNING: Chem.MolFromSmiles failed for fragment '{frag_smi}' "
                        f"(bb={bb_name}) from molecule '{query_smiles}'",
                        flush=True,
                    )
                    node_rewards.append(0.0)
                    continue
                nr = get_node_reward(
                    reward_type=self.reward_type,
                    query_smiles=frag_smi,
                    enforce_blocks_fps=self.enforced_structures_fps,
                    enforced_blocks_functional_groups=self.enforced_structures_functional_groups,
                    tango_weights=self.tango_weights,
                )
                node_rewards.append(nr)

            tango = float(np.mean(node_rewards)) if node_rewards else 0.0
            rewards.append(tango)

        return np.array(rewards, dtype=np.float32)
