"""
Tests for the NCBI-specific half of the assembly-info setup.

The download, verification, cleanup and present/absent routing all live in
gtotree/utils/misc/hosted_parquet_asset.py now and are tested there, against both this
asset and the GTDB one. What's left here is what's genuinely NCBI's: that the asset is
configured with the right variable, filenames and URLs, that the wrappers the rest of
the codebase imports still delegate to it, and read_date_retrieved(), which has no GTDB
counterpart.
"""

import pytest  # type: ignore
from unittest.mock import patch

from gtotree.utils.ncbi.get_ncbi_assembly_data import (
    DATE_FILENAME,
    NCBI_ASSET,
    NCBI_DATA_URL,
    NCBI_DATE_URL,
    PARQUET_FILENAME,
    check_if_data_present,
    check_ncbi_assembly_info_location_var_is_set,
    get_ncbi_assembly_data,
    get_ncbi_assembly_summary_tab,
    ncbi_data_table_path,
    read_date_retrieved,
)

MODPATH = "gtotree.utils.ncbi.get_ncbi_assembly_data"

# the asset's date stamp shape: a single 'YYYY,MM,DD' line
_DATE_BODY = "2026,01,05\n"


# --- the asset spec --------------------------------------------------------

class TestAssetSpec:

    def test_it_points_at_the_ncbi_variable_and_files(self):
        assert NCBI_ASSET.env_var == "NCBI_ASSEMBLY_DATA_DIR"
        assert NCBI_ASSET.parquet_filename == PARQUET_FILENAME == "ncbi-data.parquet"
        assert NCBI_ASSET.sidecar_filename == DATE_FILENAME == "date-retrieved.txt"

    def test_the_module_urls_match_the_asset(self):
        # these two names are imported by tests elsewhere, so they have to stay in step
        assert NCBI_DATA_URL == NCBI_ASSET.data_url
        assert NCBI_DATE_URL == NCBI_ASSET.sidecar_url

    def test_the_asset_is_pulled_from_bits_rolling_release(self):
        assert NCBI_DATA_URL.startswith("https://github.com/AstrobioMike/bit/releases")
        assert NCBI_DATA_URL.endswith(PARQUET_FILENAME)


# --- the wrappers delegate -------------------------------------------------

class TestWrappersDelegateToTheAsset:

    def test_location_var_returns_the_path(self, monkeypatch, tmp_path):
        monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
        assert check_ncbi_assembly_info_location_var_is_set() == str(tmp_path)

    def test_location_var_exits_if_missing(self, monkeypatch):
        monkeypatch.delenv("NCBI_ASSEMBLY_DATA_DIR", raising=False)
        with pytest.raises(SystemExit):
            check_ncbi_assembly_info_location_var_is_set()

    def test_table_path_derives_from_the_filename_constant(self, monkeypatch,
                                                           tmp_path):
        monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
        assert ncbi_data_table_path() == str(tmp_path / PARQUET_FILENAME)
        assert ncbi_data_table_path("/somewhere") == f"/somewhere/{PARQUET_FILENAME}"

    def test_summary_tab_resolves_the_env_var_on_demand(self, monkeypatch, tmp_path):
        # not at import time, so importing this module doesn't require the variable
        monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
        assert get_ncbi_assembly_summary_tab() == str(tmp_path / PARQUET_FILENAME)

    def test_presence_check_delegates(self, tmp_path):
        # the exhaustive present/absent/cleanup cases live in the shared module's tests
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
        assert check_if_data_present(str(tmp_path)) is True
        assert check_if_data_present(str(tmp_path)) is NCBI_ASSET.is_present(
            str(tmp_path))


# --- routing ---------------------------------------------------------------

class TestRouting:
    """
    get_ncbi_assembly_data() returns None rather than the directory (unlike its GTDB
    counterpart), so its routing is its own and worth pinning here.
    """

    def test_a_present_asset_is_not_re_downloaded(self, monkeypatch, tmp_path):
        monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            assert get_ncbi_assembly_data(force_update=False) is None

        mock_dl.assert_not_called()

    def test_an_absent_asset_is_downloaded(self, monkeypatch, tmp_path):
        monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            get_ncbi_assembly_data()

        mock_dl.assert_called_once()

    def test_force_update_downloads_even_if_present(self, monkeypatch, tmp_path):
        monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            get_ncbi_assembly_data(force_update=True)

        mock_dl.assert_called_once()


# --- read_date_retrieved ---------------------------------------------------

class TestReadDateRetrieved:

    def test_a_stamp_is_formatted_for_humans(self, tmp_path):
        (tmp_path / DATE_FILENAME).write_text("2026,01,05\n")
        assert read_date_retrieved(str(tmp_path)) == "Jan 05, 2026"

    def test_an_unparseable_stamp_comes_back_raw(self, tmp_path):
        # better a slightly odd date in the output than a crash on a cosmetic line
        (tmp_path / DATE_FILENAME).write_text("weird-stamp\n")
        assert read_date_retrieved(str(tmp_path)) == "weird-stamp"

    def test_an_out_of_range_date_comes_back_raw(self, tmp_path):
        (tmp_path / DATE_FILENAME).write_text("2026,13,45\n")
        assert read_date_retrieved(str(tmp_path)) == "2026,13,45"

    def test_trailing_content_after_the_stamp_is_ignored(self, tmp_path):
        (tmp_path / DATE_FILENAME).write_text("2026,01,05\nsomething else\n")
        assert read_date_retrieved(str(tmp_path)) == "Jan 05, 2026"
