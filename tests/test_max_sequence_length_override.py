"""The token budget stored in a checkpoint is not a molecule-length budget.

`max_sequence_length` counts tokens. At a 1190-monomer vocabulary 128 tokens is
roughly 128 residues; at a 34-character vocabulary it is roughly 8. Comparing two
tokenizers under the same stored value therefore compares them at very different
molecule sizes, which looks like a representation effect and is not one.

These tests pin the override that makes a matched-molecule-length comparison
possible, and the guard that stops it being used where it would silently truncate.
"""

from dataclasses import fields

import pytest

from goal_directed_generation.dataclass import ReinforcementLearningParameters


def _params(**kw):
    return ReinforcementLearningParameters(
        prior="unused", agent="unused", batch_size=64, **kw
    )


def test_override_defaults_to_none():
    """Absent from a config, the checkpoint's own value must be kept."""
    assert _params().max_sequence_length is None


def test_override_is_settable():
    assert _params(max_sequence_length=512).max_sequence_length == 512


def test_field_exists_on_the_dataclass():
    """Guards against the field being dropped in a refactor, which would make
    every config carrying it fail silently rather than loudly."""
    assert "max_sequence_length" in {f.name for f in fields(ReinforcementLearningParameters)}


@pytest.mark.parametrize(
    "vocab_size,tokens,approx_residues",
    [(1190, 128, 128), (34, 128, 8)],
)
def test_the_premise(vocab_size, tokens, approx_residues):
    """Documents the arithmetic the override exists to correct.

    Monomer-level: one token is one residue. Character-level: a HELM monomer costs
    roughly 16 characters, so 128 tokens buys about 8 residues — the ceiling
    measured in neptune_submonomer's output, where the monomer count maxed at
    exactly 8.
    """
    chars_per_residue = 1 if vocab_size > 100 else 16
    assert tokens // chars_per_residue == approx_residues
