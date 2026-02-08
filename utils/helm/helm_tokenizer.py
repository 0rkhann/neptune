"""
Enhanced HELM Tokenizer for Full HELM Notation Support

Handles complete HELM syntax including:
- Linear peptides: PEPTIDE1{[m1].[m2].[m3]}$$$$
- Cyclic peptides: PEPTIDE1{[m1].[m2].[m3]}$$PEPTIDE1,PEPTIDE1,3:R3-1:R1$$$$
- Branched structures: PEPTIDE1{[m1].[m2]}|PEPTIDE2{[m3]}$PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$$
- Multiple chains with connections

This tokenizer maintains monomer-level tokenization while preserving connection information.
"""
from typing import List
import re


class HELMTokenizer:
    """
    HELM-aware monomer tokenizer that preserves full HELM syntax.
    
    Tokenization strategy:
    1. Extract and tokenize chain sequences (monomer-level)
    2. Preserve connection annotations as structured tokens
    3. Maintain HELM delimiters ($, |, {, }, etc.)
    
    Examples:
        Linear: PEPTIDE1{[m1].[m2]}$$$$
        -> ['PEPTIDE1', '{', '[m1]', '.', '[m2]', '}', '$$$$']
        
        Cyclic: PEPTIDE1{[m1].[m2].[m3]}$$PEPTIDE1,PEPTIDE1,3:R3-1:R1$$$$
        -> ['PEPTIDE1', '{', '[m1]', '.', '[m2]', '.', '[m3]', '}', 
            '$$', 'PEPTIDE1,PEPTIDE1,3:R3-1:R1', '$$$$']
    """
    
    def __init__(self):
        # HELM component patterns
        self.patterns = {
            'chain_name': re.compile(r'(PEPTIDE\d+|RNA\d+|CHEM\d+)'),
            'monomer': re.compile(r'\[[^\]]+\]'),
            'connection': re.compile(r'\$\$([^$]+)\$\$\$\$'),  # Between $$ and $$$$
            'simple_connection': re.compile(r'([A-Z]+\d+,[A-Z]+\d+,\d+:[A-Z]\d+-\d+:[A-Z]\d+)'),
        }
    
    def tokenize(self, helm_string: str, with_begin_and_end: bool = True) -> List[str]:
        """
        Tokenize HELM string preserving monomer-level granularity and connections.
        
        Args:
            helm_string: Full HELM notation string
            with_begin_and_end: Whether to add ^ and $ tokens
            
        Returns:
            List of tokens
        """
        tokens = []
        
        # Check if this is a simple linear peptide (no connections)
        if helm_string.endswith('$$$$') and '$$$$' in helm_string:
            # Check if there's connection info (something between $$ and $$$$)
            parts = helm_string.split('$$$$')
            if len(parts) == 2 and parts[0].endswith('$$'):
                # Has connection annotation
                tokens = self._tokenize_with_connections(helm_string)
            elif '$$' not in helm_string[:-4]:  # Only trailing $$$$
                # Simple linear peptide
                tokens = self._tokenize_simple_linear(helm_string)
            else:
                # Complex structure
                tokens = self._tokenize_with_connections(helm_string)
        else:
            # Fallback for any other format
            tokens = self._tokenize_with_connections(helm_string)
        
        if with_begin_and_end:
            tokens = ['^'] + tokens + ['$']
        
        return tokens
    
    def _tokenize_simple_linear(self, helm_string: str) -> List[str]:
        """
        Tokenize simple linear peptide: PEPTIDE1{[m1].[m2].[m3]}$$$$
        """
        tokens = []
        
        # Extract chain name (e.g., PEPTIDE1)
        chain_match = self.patterns['chain_name'].match(helm_string)
        if chain_match:
            tokens.append(chain_match.group(1))
            remaining = helm_string[chain_match.end():]
        else:
            remaining = helm_string
        
        # Process remaining: {[m1].[m2].[m3]}$$$$
        i = 0
        current_token = ""
        inside_bracket = False
        
        while i < len(remaining):
            char = remaining[i]
            
            if char == '[':
                # Start of monomer
                if current_token:
                    tokens.append(current_token)
                    current_token = ""
                current_token = '['
                inside_bracket = True
            elif char == ']':
                # End of monomer
                current_token += ']'
                tokens.append(current_token)
                current_token = ""
                inside_bracket = False
            elif char == '.':
                # Dot separator
                if current_token:
                    tokens.append(current_token)
                    current_token = ""
                tokens.append('.')
            elif char in ['{', '}']:
                # Braces
                if current_token:
                    tokens.append(current_token)
                    current_token = ""
                tokens.append(char)
            elif char == '$':
                # Dollar signs - collect all consecutive $
                if current_token:
                    tokens.append(current_token)
                    current_token = ""
                dollar_run = '$'
                j = i + 1
                while j < len(remaining) and remaining[j] == '$':
                    dollar_run += '$'
                    j += 1
                tokens.append(dollar_run)
                i = j - 1
            else:
                # Regular character
                if inside_bracket or not current_token or current_token[0] != '[':
                    current_token += char
            
            i += 1
        
        if current_token:
            tokens.append(current_token)
        
        return tokens
    
    def _tokenize_with_connections(self, helm_string: str) -> List[str]:
        """
        Tokenize HELM with connection annotations (cyclic, branched).
        Example: PEPTIDE1{[m1].[m2].[m3]}$$PEPTIDE1,PEPTIDE1,3:R3-1:R1$$$$
        """
        tokens = []
        
        # Split into main structure and connections
        # Pattern: {...}$$<connections>$$$$
        if '$$' in helm_string:
            # Find the connection part
            match = self.patterns['connection'].search(helm_string)
            if match:
                # Has connection annotation
                structure_part = helm_string[:match.start() + 2]  # Include first $$
                connection_part = match.group(1)  # Connection string
                end_part = '$$$$'
                
                # Tokenize structure part (without connection)
                structure_tokens = self._tokenize_simple_linear(structure_part[:-2] + '$$$$')
                # Remove the trailing $$$$ from structure tokens
                if structure_tokens and structure_tokens[-1] == '$$$$':
                    structure_tokens = structure_tokens[:-1]
                
                tokens.extend(structure_tokens)
                tokens.append('$$')
                
                # Add connection as a structured token
                tokens.append(connection_part)
                tokens.append('$$$$')
                
                return tokens
        
        # No special connection handling needed
        return self._tokenize_simple_linear(helm_string)
    
    def untokenize(self, tokens: List[str]) -> str:
        """
        Reconstruct HELM string from tokens.
        
        Args:
            tokens: List of tokens
            
        Returns:
            HELM string
        """
        helm_string = ""
        for token in tokens:
            if token == "$":
                break
            if token != "^":
                helm_string += token
        return helm_string
    
    def __repr__(self):
        return "HELMTokenizer(monomer_level=True, full_helm_support=True)"


# For backward compatibility
MonomerTokenizer = HELMTokenizer
