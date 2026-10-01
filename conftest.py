"""Put the repository root on sys.path so tests can import the package modules.

Without this, tests only run from their own directory (see the sys.path.append
in tests/oracles/syntheseus/test_syntheseus.py).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
