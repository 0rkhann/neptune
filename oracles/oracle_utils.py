"""
Utility functions for extracting and processing oracle metadata.

This module provides helper functions for extracting various types of metadata
from oracle components using the Extract-Transform-Load (ETL) pattern.
"""

from typing import Dict, List, Any, Optional


def extract_oracle_metadata(oracles: List[Any]) -> Dict[str, Optional[Any]]:
    """
    Extract all relevant metadata from oracle components.

    Uses Extract-Transform-Load pattern to encapsulate the logic for
    finding and extracting metadata from various oracle types.

    Priority order:
    1. Direct SynthesizabilityFactor oracle (highest priority)
    2. TANGO oracle with synthesizability enabled
    3. Other oracles with metadata attributes

    Args:
        oracles: List of oracle component objects

    Returns:
        dict: Contains extracted metadata with keys:
            - 'bb_sequences': Building block sequences (comma-separated)
            - 'synth_factors': Synthesizability factor values (s^power)
            - 'tango_similarities': Raw TANGO similarity scores
            - 'property_values': Other raw property values
    """
    metadata = {
        "bb_sequences": None,
        "synth_factors": None,
        "tango_similarities": None,
        "property_values": {},
    }

    for oracle in oracles:
        # Priority 1: Direct Synthesizability Oracle
        if oracle.name == "synthesizability_factor":
            metadata["bb_sequences"] = getattr(oracle, "last_sequences", None)
            # SynthesizabilityFactor already populates synthesizability_factor_raw_values
            # via the oracle component mechanism
            break  # Exit early if primary source found

        # Priority 2: TANGO Oracle with synthesizability enabled
        if oracle.name == "tango" and getattr(oracle, "use_synthesizability", False):
            metadata["bb_sequences"] = getattr(oracle, "last_sequences", None)
            metadata["synth_factors"] = getattr(oracle, "last_synth_factors", None)
            # TANGO similarities are captured via oracle component raw_values
            break

    # Note: Other oracle component values (slogp_raw_values, qed_raw_values, etc.)
    # are already captured by the oracle component mechanism in Oracle.__call__()
    # This function only extracts metadata that requires special handling

    return metadata


def validate_metadata_length(
    metadata: Dict[str, Any], expected_length: int, strict: bool = False
) -> Dict[str, bool]:
    """
    Validate that metadata arrays have the expected length.

    Args:
        metadata: Dictionary of metadata arrays
        expected_length: Expected number of elements
        strict: If True, raise ValueError on mismatch. If False, return validation status.

    Returns:
        dict: Validation status for each metadata key

    Raises:
        ValueError: If strict=True and any metadata has wrong length
    """
    validation = {}

    for key, value in metadata.items():
        if value is None:
            validation[key] = True  # None is valid (optional data)
        elif isinstance(value, (list, tuple)):
            is_valid = len(value) == expected_length
            validation[key] = is_valid

            if strict and not is_valid:
                raise ValueError(
                    f"Metadata '{key}' has length {len(value)}, "
                    f"expected {expected_length}"
                )
        else:
            validation[key] = True  # Non-list values don't need length validation

    return validation


def merge_metadata_into_dataframe(
    df_dict: Dict[str, Any],
    metadata: Dict[str, Any],
    expected_length: int,
    skip_existing: bool = True,
) -> None:
    """
    Merge extracted metadata into a dataframe dictionary (in-place).

    Args:
        df_dict: Dictionary that will be used to create a DataFrame
        metadata: Extracted metadata from extract_oracle_metadata()
        expected_length: Expected length for validation
        skip_existing: If True, skip keys that already exist in df_dict
    """
    # Map metadata keys to dataframe column names
    column_mapping = {
        "bb_sequences": "bb_sequences",
        "synth_factors": "synthesizability_factor_raw_values",
    }

    for meta_key, column_name in column_mapping.items():
        value = metadata.get(meta_key)

        # Skip if None, wrong length, or already exists
        if value is None:
            continue
        if len(value) != expected_length:
            continue
        if skip_existing and column_name in df_dict:
            continue

        df_dict[column_name] = value
