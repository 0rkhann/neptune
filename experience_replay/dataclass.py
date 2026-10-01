from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class ExperienceReplayParameters:
    memory_size: int = 100
    sample_size: int = 10
    sequences: List[str] = field(
        default_factory=list
    )  # Can be SMILES (Saturn) or HELM (Neptune)
    scaffold_type: Optional[str] = None  # "subclass", "monomer", "sidechain", or None (Bemis-Murcko for SMILES)
