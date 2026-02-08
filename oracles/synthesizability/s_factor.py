"""
Shared utility for integrating synthesizability checking into oracles.

This module provides a unified interface for synthesizability checking
that can be used by any oracle (TANGO, QED, etc.) without code duplication.

Requires GenAI4Peptidomimetic_native to be installed:
    git clone https://github.com/yourorg/GenAI4Peptidomimetic_native.git
    cd GenAI4Peptidomimetic_native
    pip install -e .
"""

from pathlib import Path
from typing import Protocol

try:
    from src.genai_utils.synthesizability import SynthesizabilityChecker
except ImportError as e:
    raise ImportError(
        "GenAI4Peptidomimetic_native is required but not found. "
        "Please install it:\n"
        "  git clone <repo_url>/GenAI4Peptidomimetic_native.git\n"
        "  cd GenAI4Peptidomimetic_native\n"
        "  pip install -e .\n"
    ) from e


class RetrosynthesisChecker(Protocol):
    """
    Protocol for retrosynthetic checkers that decompose molecules to building blocks.

    Any checker that implements decompose_partial() can be used:
    - GenAI4Peptidomimetic SynthesizabilityChecker
    - Syntheseus-based checkers (future)
    - Custom checkers
    """

    def decompose_partial(
        self, smiles: str, unknown_marker: str = "UNKNOWN"
    ) -> tuple[list[str], int, int]:
        """
        Decompose SMILES into building blocks.

        Returns:
            (sequence, valid_count, total_count)
        """
        ...


class SynthesizabilityIntegration:
    """
    Manages synthesizability checking for oracles.

    This is a shared utility to avoid duplicating synthesizability
    initialization logic across multiple oracle classes.

    Usage:
        config = {
            "enforced_structures": "/path/to/building_blocks.csv",
            "max_sequence_length": 5,
            "max_search_attempts": 10000,
            "include_caps": True,
            "n_workers": 4,
            "timeout_per_molecule": 180
        }
        synth = SynthesizabilityIntegration(config)
        s_factor = synth.get_synthesizability_factor(smiles, power=1.0)
    """

    def __init__(self, config: dict):
        """
        Initialize synthesizability checker from config.

        Args:
            config: Dictionary with synthesizability configuration.
                Required:
                    - enforced_structures: Path to building blocks CSV file (absolute path)

                Optional (GenAI SynthesizabilityChecker parameters):
                    - max_sequence_length: Maximum number of building blocks per molecule (default: 5)
                    - max_search_attempts: Maximum decomposition attempts (default: 10000)
                    - include_caps: Whether to include N-caps and C-caps (default: True)
                    - n_workers: Number of parallel workers (default: CPU count)
                    - timeout_per_molecule: Timeout in seconds per molecule (default: 180)

        Raises:
            ValueError: If required config parameters are missing
            FileNotFoundError: If building blocks CSV not found
        """
        self.checker: RetrosynthesisChecker = self._initialize_checker(config)

    def _initialize_checker(self, config: dict) -> RetrosynthesisChecker:
        """
        Initialize retrosynthetic checker from GenAI4Peptidomimetic_native package.

        Raises:
            ValueError: If required config is missing
            FileNotFoundError: If paths not found
        """
        # Get building blocks library path
        enforced_structures = config.get("enforced_structures")
        if not enforced_structures:
            raise ValueError(
                "Missing required config: 'enforced_structures'. "
                "Provide path to building blocks CSV file."
            )

        # Validate path exists
        bb_csv = Path(enforced_structures)
        if not bb_csv.exists():
            raise FileNotFoundError(
                f"Building blocks CSV not found: {bb_csv}\n"
                f"Please provide correct path to building blocks CSV file."
            )

        # Initialize checker with all optional parameters
        checker = SynthesizabilityChecker(
            bb_csv_path=str(bb_csv),
            max_sequence_length=config.get("max_sequence_length", 5),
            max_search_attempts=config.get("max_search_attempts", 10000),
            include_caps=config.get("include_caps", True),
            n_workers=config.get("n_workers", None),  # None = auto-detect CPU count
            timeout_per_molecule=config.get("timeout_per_molecule", 180),
        )

        print(
            f"[SynthesizabilityIntegration] Initialized SynthesizabilityChecker "
            f"with {len(checker.building_blocks)} building blocks "
            f"(max_seq={config.get('max_sequence_length', 5)}, "
            f"timeout={config.get('timeout_per_molecule', 180)}s)"
        )

        return checker

    def get_synthesizability_factor(
        self, smiles: str, power: float = 1.0
    ) -> tuple[float, str]:
        """
        Get synthesizability factor and building block sequence for a SMILES string.

        The synthesizability factor 's' is calculated as:
            s = (# of building blocks in library) / (# total building blocks)

        For example:
            - Molecule with 5 BBs, 4 in library: s = 4/5 = 0.8
            - Molecule with 5 BBs, all in library: s = 5/5 = 1.0
            - Molecule with 5 BBs, none in library: s = 0/5 = 0.0

        Args:
            smiles: SMILES string to check
            power: Exponent to apply to synthesizability factor (s^power)

        Returns:
            Tuple of (s^power, sequence) where:
            - s^power is synthesizability factor ∈ [0, 1]
            - sequence is comma-separated building block names (e.g., "na1,l_a_1,cp1")

        Power control:
            - power = 1.0: Linear penalty (s)
            - power = 2.0: s² (steeper penalty for partial matches)
            - power = 3.0: s³ (even steeper)
            - power = 0.5: s^0.5 (softer penalty)

        Raises:
            Exception: If decomposition fails
        """
        try:
            # Use decompose_partial to get sequence, valid_count and total_count
            sequence, valid_count, total_count = self.checker.decompose_partial(
                smiles, unknown_marker="UNKNOWN"
            )

            # Convert sequence list to comma-separated string
            sequence_str = ",".join(sequence) if sequence else ""

            # Calculate synthesizability factor
            if total_count == 0:
                s = 0.0  # Invalid molecule
            else:
                s = valid_count / total_count

            # Apply power: s^power
            s_scaled = s**power
            return s_scaled, sequence_str
        except Exception as e:
            # Re-raise with context
            raise RuntimeError(
                f"Synthesizability check failed for SMILES: {smiles[:50]}..."
            ) from e
