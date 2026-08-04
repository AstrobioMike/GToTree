"""
Shared filesystem paths for the test suite

`DATA_DIR` is anchored to this file's directory

Importable at module level (unlike a fixture), so it works for module-level constants:

    from gtotree.tests.paths import DATA_DIR
    MOCK_PFAM_HMM = DATA_DIR / "mock-pfams.hmm"
"""

from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"

MOCK_PFAM_HMM = DATA_DIR / "mock-pfams.hmm"
MOCK_PFAM_INFO = DATA_DIR / "mock-pfamA.txt"
