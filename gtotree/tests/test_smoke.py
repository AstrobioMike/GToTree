import shutil
import pytest  # type: ignore
from gtotree.tests.smoke import run_smoke_test


REQUIRED_BINARIES = ["muscle", "hmmsearch", "trimal", "FastTree"]


def _missing_binaries():
    return [b for b in REQUIRED_BINARIES if shutil.which(b) is None]


def test_smoke_end_to_end(tmp_path, monkeypatch):
    missing = _missing_binaries()
    if missing:
        pytest.skip(f"external binaries not on PATH: {', '.join(missing)}")

    monkeypatch.chdir(tmp_path)
    assert run_smoke_test([]) == 0


def test_smoke_run_leaves_no_files_behind(tmp_path, monkeypatch):

    missing = _missing_binaries()
    if missing:
        pytest.skip(f"external binaries not on PATH: {', '.join(missing)}")

    monkeypatch.chdir(tmp_path)
    assert run_smoke_test([]) == 0
    leftovers = sorted(p.name for p in tmp_path.iterdir())
    assert leftovers == [], f"smoke test left files behind: {leftovers}"
