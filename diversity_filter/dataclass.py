from dataclasses import dataclass
from typing import Optional


@dataclass
class DiversityFilterParameters:
    # Based on REINVENT 3.2: https://github.com/MolecularAI/Reinvent
    name: str = "IdenticalMurckoScaffold"  # "IdenticalMurckoScaffold" for Saturn, "SubclassFingerprint"/"MonomerFingerprint" for Neptune, "SidechainFingerprint" for ECFP4
    bucket_size: int = 10
    bb_csv_path: Optional[str] = None  # Required for SidechainFingerprint (building blocks CSV)
    tanimoto_threshold: float = 0.65  # Similarity threshold for SidechainFingerprint (Tanimoto)