import pytest  # type: ignore
from gtotree.tests.paths import DATA_DIR, MOCK_PFAM_HMM, MOCK_PFAM_INFO
from gtotree.utils.misc.context import log_file_var

RUN_LOG_NAME = "gtotree-runlog.txt"


@pytest.fixture(autouse=True)
def isolated_run_log(tmp_path_factory):
    """
    Point the run log somewhere disposable for every test.

    `log_file_var` defaults to the bare relative name "gtotree-runlog.txt", and
    `capture_stdout_to_log` opens it with a plain `open(path, 'a')`. Only
    preflight_checks rewrites that to an absolute path inside the output directory, so
    any test that reaches a decorated reporter without going through preflight appends
    to the directory pytest was invoked from -- the repo root. Autouse rather than
    opt-in because the write is a side effect of *reporting*, which means the tests
    that trip it are rarely the ones that look like they touch logging.
    """
    token = log_file_var.set(str(tmp_path_factory.mktemp("run-log") / RUN_LOG_NAME))
    try:
        yield
    finally:
        log_file_var.reset(token)


@pytest.fixture(scope="session", autouse=True)
def no_run_log_left_in_rootdir(pytestconfig):
    """
    Backstop for the above: fail loudly if a run log still lands in the repo root.

    The autouse redirect covers the main thread, but `log_file_var` is a ContextVar and
    worker threads start with a fresh context, so a reporter called from inside a pool
    would fall back to the relative default and slip past it.
    """
    log_path = pytestconfig.rootpath / RUN_LOG_NAME
    preexisting = log_path.exists()

    yield

    if log_path.exists() and not preexisting:
        log_path.unlink()
        pytest.fail(
            f"the test run wrote {RUN_LOG_NAME} into {pytestconfig.rootpath} -- "
            "something reached a @capture_stdout_to_log reporter while log_file_var "
            "still held its relative default",
            pytrace=False,
        )


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
