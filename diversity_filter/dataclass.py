from dataclasses import dataclass
from typing import Optional


@dataclass
class DiversityFilterParameters:
    # Based on REINVENT 3.2: https://github.com/MolecularAI/Reinvent
    name: str = "IdenticalMurckoScaffold"  # "IdenticalMurckoScaffold" for Saturn, "IdenticalBigramScaffold" for Neptune
    bucket_size: int = 10
    bb_csv_path: Optional[str] = None  # Not used for diversity filtering (kept for backward compatibility)