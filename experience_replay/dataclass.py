from dataclasses import dataclass, field
from typing import List


@dataclass
class ExperienceReplayParameters:
    memory_size: int = 100
    sample_size: int = 10
    sequences: List[str] = field(
        default_factory=list
    )  # Can be SMILES (Saturn) or HELM (Neptune)
    use_bigram_scaffold: bool = (
        False  # Use bigram scaffolds for selective memory purge (Neptune/HELM only)
    )
