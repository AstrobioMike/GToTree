import pytest  # type: ignore

from gtotree.tests.paths import DATA_DIR, MOCK_PFAM_HMM, MOCK_PFAM_INFO


@pytest.fixture(scope="session")
def data_dir():
    """
    The test data directory
    """
    return DATA_DIR


@pytest.fixture(scope="session")
def mock_pfam_hmm():
    """
    The small mock Pfam HMM used across the gen-scg-hmms tests
    """
    return MOCK_PFAM_HMM


@pytest.fixture(scope="session")
def mock_pfam_info():
    """
    The small mock Pfam metadata table used across the gen-scg-hmms tests
    """
    return MOCK_PFAM_INFO
