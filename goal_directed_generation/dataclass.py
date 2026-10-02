from dataclasses import dataclass
from typing import Optional
from experience_replay.dataclass import ExperienceReplayParameters
from hallucinated_memory.dataclass import HallucinatedMemoryParameters
from beam_enumeration.dataclass import BeamEnumerationParameters
from diversity_filter.dataclass import DiversityFilterParameters

@dataclass
class ReinforcementLearningParameters:
    prior: str
    agent: str
    batch_size: int
    learning_rate: float = 0.0001
    sigma: float = 128.0
    augmented_memory: bool = True
    augmentation_rounds: int = 10
    selective_memory_purge: bool = True
    # Overrides the value stored in the checkpoint, which counts TOKENS. The same
    # token budget means different molecule sizes per tokenizer: at a 1190-monomer
    # vocabulary 128 tokens is ~128 residues, at a 34-character vocabulary it is
    # ~8. Set this to compare tokenizers at a matched molecule length instead of a
    # matched token count. None keeps the checkpoint value.
    max_sequence_length: Optional[int] = None

@dataclass
class GoalDirectedGenerationConfiguration:
    seed: int
    model_architecture: str
    reinforcement_learning: ReinforcementLearningParameters
    experience_replay: ExperienceReplayParameters
    diversity_filter: DiversityFilterParameters
    hallucinated_memory: HallucinatedMemoryParameters
    beam_enumeration: BeamEnumerationParameters
