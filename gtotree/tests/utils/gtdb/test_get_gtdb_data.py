"""
Tests for the GTDB-specific half of the reference-table setup.

As with the NCBI twin, the download/verify/cleanup machinery is shared and tested in
gtotree/tests/utils/misc/test_hosted_parquet_asset.py against both assets. This module
covers what's GTDB's alone: the asset spec, the wrappers, the fact that
get_gtdb_data() hands back the directory (its NCBI counterpart returns None), and
report_gtdb_version_info().
"""

import pytest  # type: ignore
from unittest.mock import patch

from gtotree.utils.gtdb.get_gtdb_data import (
    GTDB_ASSET,
    GTDB_DATA_URL,
    GTDB_VERSION_URL,
    PARQUET_FILENAME,
    VERSION_FILENAME,
    check_gtdb_location_var_is_set,
    check_if_gtdb_data_present,
    get_gtdb_data,
    gtdb_data_table_path,
    report_gtdb_version_info,
)

MODPATH = "gtotree.utils.gtdb.get_gtdb_data"

_VERSION_BODY = "r220\n2024-04-24\n"


# --- the asset spec --------------------------------------------------------

class TestAssetSpec:

    def test_it_points_at_the_gtdb_variable_and_files(self):
        assert GTDB_ASSET.env_var == "GTDB_DIR"
        assert GTDB_ASSET.parquet_filename == PARQUET_FILENAME == "gtdb-data.parquet"
        assert GTDB_ASSET.sidecar_filename == VERSION_FILENAME == "VERSION.txt"

    def test_the_module_urls_match_the_asset(self):
        assert GTDB_DATA_URL == GTDB_ASSET.data_url
        assert GTDB_VERSION_URL == GTDB_ASSET.sidecar_url

    def test_the_asset_is_pulled_from_bits_rolling_release(self):
        assert GTDB_DATA_URL.startswith("https://github.com/AstrobioMike/bit/releases")
        assert GTDB_DATA_URL.endswith(PARQUET_FILENAME)


# --- the wrappers delegate -------------------------------------------------

class TestWrappersDelegateToTheAsset:

    def test_location_var_returns_the_path(self, monkeypatch, tmp_path):
        monkeypatch.setenv("GTDB_DIR", str(tmp_path))
        assert check_gtdb_location_var_is_set() == str(tmp_path)

    def test_location_var_exits_if_missing(self, monkeypatch):
        monkeypatch.delenv("GTDB_DIR", raising=False)
        with pytest.raises(SystemExit) as excinfo:
            check_gtdb_location_var_is_set()
        assert excinfo.value.code == 1

    def test_table_path_derives_from_the_filename_constant(self, monkeypatch,
                                                           tmp_path):
        monkeypatch.setenv("GTDB_DIR", str(tmp_path))
        assert gtdb_data_table_path() == str(tmp_path / PARQUET_FILENAME)
        assert gtdb_data_table_path("/somewhere") == f"/somewhere/{PARQUET_FILENAME}"

    def test_presence_check_delegates(self, tmp_path):
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / VERSION_FILENAME).write_text(_VERSION_BODY)
        assert check_if_gtdb_data_present(str(tmp_path)) is True


# --- routing ---------------------------------------------------------------

class TestRouting:
    """
    Unlike get_ncbi_assembly_data(), this one returns the directory -- callers in
    get_accessions_from_gtdb and dl_ncbi_assemblies use the return value.
    """

    def test_a_present_asset_is_not_re_downloaded_and_the_dir_comes_back(
            self, monkeypatch, tmp_path):
        monkeypatch.setenv("GTDB_DIR", str(tmp_path))
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / VERSION_FILENAME).write_text(_VERSION_BODY)

        with patch(f"{MODPATH}.get_slim_gtdb_tab") as mock_dl:
            assert get_gtdb_data(force_update=False) == str(tmp_path)

        mock_dl.assert_not_called()

    def test_an_absent_asset_is_downloaded_and_the_dir_still_comes_back(
            self, monkeypatch, tmp_path):
        monkeypatch.setenv("GTDB_DIR", str(tmp_path))

        with patch(f"{MODPATH}.get_slim_gtdb_tab") as mock_dl:
            assert get_gtdb_data() == str(tmp_path)

        mock_dl.assert_called_once()

    def test_force_update_downloads_even_if_present(self, monkeypatch, tmp_path):
        monkeypatch.setenv("GTDB_DIR", str(tmp_path))
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / VERSION_FILENAME).write_text(_VERSION_BODY)

        with patch(f"{MODPATH}.get_slim_gtdb_tab") as mock_dl:
            get_gtdb_data(force_update=True)

        mock_dl.assert_called_once()


# --- report_gtdb_version_info ----------------------------------------------

class TestReportGtdbVersionInfo:

    def test_the_two_lines_come_back_as_a_pair(self, tmp_path):
        (tmp_path / VERSION_FILENAME).write_text(_VERSION_BODY)
        assert report_gtdb_version_info(str(tmp_path)) == ("r220", "2024-04-24")

    def test_a_released_prefix_is_stripped(self, tmp_path):
        # bit publishes the date line as "Released <date>"
        (tmp_path / VERSION_FILENAME).write_text("R220\nReleased 2024-04-24\n")
        assert report_gtdb_version_info(str(tmp_path)) == ("R220", "2024-04-24")

    def test_blank_lines_are_skipped(self, tmp_path):
        (tmp_path / VERSION_FILENAME).write_text("\nR220\n\n2024-04-24\n\n")
        assert report_gtdb_version_info(str(tmp_path)) == ("R220", "2024-04-24")

    def test_a_one_line_file_raises_rather_than_returning_junk(self, tmp_path):
        # the validator rejects this shape on download, so reaching it means the file
        # was hand-placed; better to fail loudly than to report a bogus release date
        (tmp_path / VERSION_FILENAME).write_text("R220\n")
        with pytest.raises(IndexError):
            report_gtdb_version_info(str(tmp_path))
