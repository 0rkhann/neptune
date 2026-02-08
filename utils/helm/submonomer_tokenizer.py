"""
Sub-monomer tokenizer for Saturn using HelmDictionary (character-level).

This tokenizer tokenizes HELM sequences at a character level using HelmDictionary,
matching HELM-GPT's original tokenization approach. This provides fine-grained
character-level tokenization compared to monomer-level tokenization.

Adapted from helm-gpt/utils/dataset.py (HelmDictionary)
"""

from typing import List


class HelmDictionary:
    """
    Character-level dictionary for HELM sequences.

    This is the same tokenizer used by HELM-GPT. It tokenizes at the character level
    with a fixed vocabulary of ~76 characters, including special replacements for
    common HELM patterns.
    """

    PAD, BEGIN, END = " ", "@", "\n"

    def __init__(self) -> None:
        self.char_idx = {
            self.PAD: 0,
            self.BEGIN: 1,
            self.END: 2,
            "A": 3,
            "R": 4,
            "N": 5,
            "D": 6,
            "C": 7,
            "E": 8,
            "Q": 9,
            "G": 10,
            "H": 11,
            "I": 12,
            "L": 13,
            "K": 14,
            "M": 15,
            "F": 16,
            "P": 17,
            "S": 18,
            "T": 19,
            "W": 20,
            "Y": 21,
            "V": 22,
            "X": 23,  # Natural amino acids
            "$": 24,
            "(": 25,
            ")": 26,
            ",": 27,
            "-": 28,
            ".": 29,
            ":": 30,
            "[": 31,
            "]": 32,
            "{": 33,
            "|": 34,
            "}": 35,  # Common symbols in HELM
            "0": 36,
            "1": 37,
            "2": 38,
            "3": 39,
            "4": 40,
            "5": 41,
            "6": 42,
            "7": 43,
            "8": 44,
            "9": 45,
            ">": 46,
            "B": 47,
            "O": 48,
            "_": 49,
            "a": 50,
            "b": 51,
            "c": 52,
            "d": 53,
            "e": 54,
            "f": 55,
            "g": 56,
            "h": 57,
            "i": 58,
            "l": 59,
            "m": 60,
            "n": 61,
            "o": 62,
            "p": 63,
            "r": 64,
            "s": 65,
            "t": 66,
            "u": 67,
            "v": 68,
            "x": 69,
            "y": 70,
            "z": 71,  # characters
            "/": 72,
            "*": 73,
            "\t": 74,
            "&": 75,  # special tokens
        }
        self.idx_char = {v: k for k, v in self.char_idx.items()}

        # Special replacements for common HELM patterns
        self.encode_dict = {"PEPTIDE": "/", "me": "*", "am": "\t", "ac": "&"}
        self.decode_dict = {v: k for k, v in self.encode_dict.items()}

    def get_vocab_size(self) -> int:
        """Get the vocabulary size."""
        return len(self.idx_char)

    def encode(self, helm: str) -> str:
        """
        Replace multi-char tokens with single tokens in HELM string.

        Args:
            helm: HELM notation string

        Returns:
            HELM string with special patterns replaced
        """
        temp_helm = helm
        for symbol, token in self.encode_dict.items():
            temp_helm = temp_helm.replace(symbol, token)
        return temp_helm

    def decode(self, helm: str) -> str:
        """
        Replace special tokens with their multi-character equivalents.

        Args:
            helm: Encoded HELM string

        Returns:
            Decoded HELM string
        """
        temp_helm = helm
        for symbol, token in self.decode_dict.items():
            temp_helm = temp_helm.replace(symbol, token)
        return temp_helm

    def string_to_indices(
        self, helm: str, with_begin_and_end: bool = True
    ) -> List[int]:
        """
        Convert HELM string to list of character indices.

        Args:
            helm: HELM notation string
            with_begin_and_end: Whether to add BEGIN and END tokens

        Returns:
            List of character indices
        """
        # Apply special replacements
        encoded_helm = self.encode(helm)

        # Convert to indices
        indices = []
        if with_begin_and_end:
            indices.append(self.char_idx[self.BEGIN])

        for char in encoded_helm:
            if char in self.char_idx:
                indices.append(self.char_idx[char])
            else:
                # Unknown character - use PAD as fallback
                indices.append(self.char_idx[self.PAD])

        if with_begin_and_end:
            indices.append(self.char_idx[self.END])

        return indices

    def indices_to_string(
        self, indices: List[int], remove_begin_end: bool = True
    ) -> str:
        """
        Convert list of character indices back to HELM string.

        Args:
            indices: List of character indices
            remove_begin_end: Whether to remove BEGIN and END tokens

        Returns:
            HELM notation string
        """
        chars = []
        for idx in indices:
            if idx in self.idx_char:
                char = self.idx_char[idx]
                if remove_begin_end and char in [self.BEGIN, self.END]:
                    continue
                chars.append(char)

        helm_string = "".join(chars)
        # Decode special replacements
        return self.decode(helm_string)

    @property
    def begin_idx(self) -> int:
        """Get BEGIN token index."""
        return self.char_idx[self.BEGIN]

    @property
    def end_idx(self) -> int:
        """Get END token index."""
        return self.char_idx[self.END]

    @property
    def pad_idx(self) -> int:
        """Get PAD token index."""
        return self.char_idx[self.PAD]


class SubMonomerTokenizer:
    """
    Sub-monomer tokenizer wrapper for Saturn using HelmDictionary.

    This tokenizer uses character-level tokenization (HelmDictionary) to tokenize
    HELM sequences, matching HELM-GPT's original approach. It provides a
    compatible interface with Saturn's tokenizer (tokenize/untokenize methods).

    Example:
        tokenizer = SubMonomerTokenizer()
        tokens = tokenizer.tokenize("PEPTIDE1{[m1].[m2]}$$$$")
        # Returns: ['@', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '1', '{', '[', 'm', '1', ']', ...]
    """

    def __init__(self):
        """
        Initialize sub-monomer tokenizer with HelmDictionary.
        """
        self.helm_dict = HelmDictionary()

    def tokenize(self, helm_string: str, with_begin_and_end: bool = True) -> List[str]:
        """
        Tokenize HELM string at character level using HelmDictionary.

        Args:
            helm_string: HELM notation string
            with_begin_and_end: Whether to add BEGIN ('@') and END ('\\n') tokens

        Returns:
            List of character tokens
        """
        # Apply special replacements (PEPTIDE -> /, etc.)
        encoded_helm = self.helm_dict.encode(helm_string)

        # Convert to character tokens
        tokens = []
        if with_begin_and_end:
            tokens.append(self.helm_dict.BEGIN)

        for char in encoded_helm:
            tokens.append(char)

        if with_begin_and_end:
            tokens.append(self.helm_dict.END)

        return tokens

    def untokenize(self, tokens: List[str]) -> str:
        """
        Reconstruct HELM string from character tokens.

        Args:
            tokens: List of character tokens (may include '@' and '\\n')

        Returns:
            HELM notation string
        """
        # Remove begin/end tokens
        filtered_tokens = [
            t for t in tokens if t not in [self.helm_dict.BEGIN, self.helm_dict.END]
        ]

        # Join tokens back into string
        helm_string = "".join(filtered_tokens)

        # Decode special replacements (/, *, \t, & -> PEPTIDE, me, am, ac)
        return self.helm_dict.decode(helm_string)

    def get_vocab_size(self) -> int:
        """Get the vocabulary size."""
        return self.helm_dict.get_vocab_size()

    def get_begin_token(self) -> str:
        """Get the BEGIN token (for compatibility with generator)."""
        return self.helm_dict.BEGIN  # '@'

    def get_end_token(self) -> str:
        """Get the END token (for compatibility with generator)."""
        return self.helm_dict.END  # '\n'

    def __repr__(self):
        return f"SubMonomerTokenizer(vocab_size={self.get_vocab_size()})"
