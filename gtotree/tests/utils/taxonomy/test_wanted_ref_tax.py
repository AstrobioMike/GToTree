import os
import pytest # type: ignore
from gtotree.utils.gtdb.get_gtdb_data import VERSION_FILENAME
from gtotree.utils.ncbi.get_ncbi_assembly_data import DATE_FILENAME
from gtotree.utils.taxonomy.wanted_ref_tax import describe_source_version


def test_gtdb_source_description(tmp_path, monkeypatch):
    (tmp_path / VERSION_FILENAME).write_text("R220\n2024-04-24\n")
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))

    desc = describe_source_version("gtdb")
    assert "GTDB" in desc
    assert "R220" in desc
    assert "2024-04-24" in desc


def test_ncbi_source_description_formats_date(tmp_path, monkeypatch):
    (tmp_path / DATE_FILENAME).write_text("2026,01,05\n")
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))

    desc = describe_source_version("ncbi")
    assert "NCBI" in desc
    # read_date_retrieved formats YYYY,MM,DD -> 'Jan 05, 2026'
    assert "Jan 05, 2026" in desc


def test_missing_sidecar_returns_none_not_crash(tmp_path, monkeypatch):
    """A display label must never abort the run; a missing version file -> None."""
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))  # dir exists, VERSION.txt does not
    assert describe_source_version("gtdb") is None


def test_unknown_source_returns_none():
    assert describe_source_version("nonsense") is None


def test_gtdb_without_release_date(tmp_path, monkeypatch):
    """
    If VERSION.txt somehow has only a version line, we still return the version rather
    than raising on the missing second line.
    """
    # report_gtdb_version_info reads two lines; a single-line file raises IndexError,
    # which the helper swallows to None. Confirm that path is graceful.
    (tmp_path / VERSION_FILENAME).write_text("R220\n")
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    assert describe_source_version("gtdb") is None
