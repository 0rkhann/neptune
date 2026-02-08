"""
HELM utilities package for Saturn.

Provides HELM tokenization and imports conversion utilities from GenAI4Peptidomimetic_native.

Requirements:
    GenAI4Peptidomimetic_native must be installed:
        git clone <repo_url>/GenAI4Peptidomimetic_native.git
        cd GenAI4Peptidomimetic_native
        pip install -e .
"""

# Import local tokenizers
from .monomer_tokenizer import HELMTokenizer
from .submonomer_tokenizer import SubMonomerTokenizer, HelmDictionary


# Import HELM utilities from GenAI4Peptidomimetic_native package
try:
    from src.genai_utils.helm import (
        create_helm_string,
        decode_helm_string,
        is_helm_notation,
        deduplicate_helm_strings,
        convert_helm_to_smiles,
        load_valid_monomers,
        check_helm_validness,
        HELMConversionError,
    )

    HELM_CONVERSION_AVAILABLE = True
except ImportError as e:
    # Provide fallback with helpful error message
    import warnings

    warnings.warn(
        f"Could not import HELM utilities from GenAI4Peptidomimetic_native: {e}. "
        "Only tokenizers will be available. "
        "To enable HELM conversion, install GenAI4Peptidomimetic_native:\n"
        "  git clone <repo_url>/GenAI4Peptidomimetic_native.git\n"
        "  cd GenAI4Peptidomimetic_native\n"
        "  pip install -e .",
        ImportWarning,
    )
    create_helm_string = None
    decode_helm_string = None
    is_helm_notation = None
    deduplicate_helm_strings = None
    convert_helm_to_smiles = None
    load_valid_monomers = None
    check_helm_validness = None

    class HELMConversionError(Exception):
        """HELM conversion error (fallback)."""

        pass


__all__ = [
    # Tokenizers (local)
    "HELMTokenizer",
    "SubMonomerTokenizer",
    "HelmDictionary",
    # HELM utilities (from GenAI4Peptidomimetic_native)
    "create_helm_string",
    "decode_helm_string",
    "is_helm_notation",
    "deduplicate_helm_strings",
    "convert_helm_to_smiles",
    "load_valid_monomers",
    "check_helm_validness",
    "HELMConversionError",
]
