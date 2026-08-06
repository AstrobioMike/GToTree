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


@pytest.fixture(scope="session")
def nt_fixture_dir(tmp_path_factory):
    from gtotree.tests import nt_fixtures

    directory = tmp_path_factory.mktemp("nt-fixtures")
    nt_fixtures.write_genomes(directory)
    nt_fixtures.write_hmm(directory)
    return directory


@pytest.fixture(scope="session")
def nt_genomes(nt_fixture_dir):
    """Paths to the generated nucleotide genome FASTAs, in order."""
    return sorted(nt_fixture_dir.glob("nt-mock-*.fna"))


@pytest.fixture(scope="session")
def nt_hmm(nt_fixture_dir):
    """Path to the generated SCG profile matching `nt_genomes`."""
    return nt_fixture_dir / "nt-mock.hmm"
