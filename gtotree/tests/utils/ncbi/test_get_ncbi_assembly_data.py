"""
Tests for the Parquet-based NCBI assembly-info setup (get_ncbi_assembly_data).

These replace the old tarball/rebuild tests: GToTree now consumes bit's hosted
ncbi-data.parquet directly, so there is no NCBI-rebuild path to exercise -- only
download, Parquet-footer verification, atomic date-file write, and present/absent
routing.
"""

import socket
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore
from unittest.mock import patch

from gtotree.utils.ncbi.get_ncbi_assembly_data import (
    PARQUET_FILENAME,
    DATE_FILENAME,
    NCBI_DATA_URL,
    NCBI_DATE_URL,
    check_ncbi_assembly_info_location_var_is_set,
    ncbi_data_table_path,
    get_ncbi_assembly_summary_tab,
    check_if_data_present,
    get_slim_ncbi_assembly_data,
    get_ncbi_assembly_data,
    read_date_retrieved,
)

MODPATH = "gtotree.utils.ncbi.get_ncbi_assembly_data"


# --- helpers ---------------------------------------------------------------

def _valid_parquet(path):
    pq.write_table(pa.table({"assembly_accession": pa.array(["GCA_1"])}), str(path))


# the asset's date stamp shape: a single 'YYYY,MM,DD' line
_DATE_BODY = "2026,01,05\n"


def _fake_downloader(date_body=_DATE_BODY, break_parquet=False):
    """
    Stand-in for download_with_tqdm(url, label, filename, ...): serves the parquet URL
    to a real (small) Parquet file and the date URL to `date_body`.
    """
    def _dl(url, label, filename=None, **kw):
        if url == NCBI_DATA_URL:
            if break_parquet:
                with open(filename, "w") as fh:
                    fh.write("not a parquet")
            else:
                _valid_parquet(filename)
        elif url == NCBI_DATE_URL:
            with open(filename, "w") as fh:
                fh.write(date_body)
        else:
            raise AssertionError(f"unexpected download URL: {url}")
    return _dl


# --- location var ----------------------------------------------------------

def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    assert check_ncbi_assembly_info_location_var_is_set() == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("NCBI_ASSEMBLY_DATA_DIR", raising=False)
    with pytest.raises(SystemExit):
        check_ncbi_assembly_info_location_var_is_set()


def test_table_path_derives_from_filename_constant(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    assert ncbi_data_table_path() == str(tmp_path / PARQUET_FILENAME)
    assert ncbi_data_table_path("/somewhere") == f"/somewhere/{PARQUET_FILENAME}"


def test_summary_tab_accessor_is_lazy(monkeypatch, tmp_path):
    # get_ncbi_assembly_summary_tab() resolves via the env var on demand
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    assert get_ncbi_assembly_summary_tab() == str(tmp_path / PARQUET_FILENAME)


# --- check_if_data_present -------------------------------------------------

def test_present_when_both_files_nonempty(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
    assert check_if_data_present(str(tmp_path)) is True


def test_absent_when_table_missing(tmp_path):
    (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
    assert check_if_data_present(str(tmp_path)) is False


def test_absent_when_date_missing(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    assert check_if_data_present(str(tmp_path)) is False


def test_half_present_pair_is_cleaned_up(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / PARQUET_FILENAME).exists()


def test_empty_files_count_as_absent_and_are_removed(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("")
    (tmp_path / DATE_FILENAME).write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / DATE_FILENAME).exists()


# --- get_slim_ncbi_assembly_data (download path) ---------------------------

def test_download_writes_both_files(tmp_path):
    with patch(f"{MODPATH}.download_with_tqdm", side_effect=_fake_downloader()):
        get_slim_ncbi_assembly_data(str(tmp_path))
    assert (tmp_path / PARQUET_FILENAME).exists()
    assert (tmp_path / DATE_FILENAME).read_text().startswith("2026,01,05")
    assert pq.ParquetFile(str(tmp_path / PARQUET_FILENAME)).metadata.num_rows == 1


def test_download_verifies_parquet_and_bails_on_corruption(tmp_path):
    with patch(f"{MODPATH}.download_with_tqdm",
               side_effect=_fake_downloader(break_parquet=True)), \
         patch(f"{MODPATH}.report_early_exit", side_effect=SystemExit(1)):
        with pytest.raises(SystemExit):
            get_slim_ncbi_assembly_data(str(tmp_path))
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / DATE_FILENAME).exists()


def test_download_network_error_is_translated(tmp_path):
    def _boom(*a, **k):
        raise socket.timeout("no route")
    with patch(f"{MODPATH}.download_with_tqdm", side_effect=_boom), \
         patch(f"{MODPATH}.report_early_exit", side_effect=SystemExit(1)):
        with pytest.raises(SystemExit):
            get_slim_ncbi_assembly_data(str(tmp_path))
    assert not (tmp_path / PARQUET_FILENAME).exists()


def test_bad_date_file_is_rejected(tmp_path):
    # a date file that isn't a 'YYYY,MM,DD' stamp fails validation -> cleanup + exit
    with patch(f"{MODPATH}.download_with_tqdm",
               side_effect=_fake_downloader(date_body="not-a-date\n")), \
         patch(f"{MODPATH}.report_early_exit", side_effect=SystemExit(1)):
        with pytest.raises(SystemExit):
            get_slim_ncbi_assembly_data(str(tmp_path))
    assert not (tmp_path / DATE_FILENAME).exists()


# --- get_ncbi_assembly_data (routing) --------------------------------------

def test_routing_skips_download_when_present(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
    with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
        get_ncbi_assembly_data(force_update=False)
    mock_dl.assert_not_called()


def test_routing_downloads_when_absent(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
        get_ncbi_assembly_data()
    mock_dl.assert_called_once()


def test_routing_force_update_downloads_even_if_present(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_ASSEMBLY_DATA_DIR", str(tmp_path))
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
    with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
        get_ncbi_assembly_data(force_update=True)
    mock_dl.assert_called_once()


# --- read_date_retrieved ---------------------------------------------------

def test_read_date_retrieved_formats_stamp(tmp_path):
    (tmp_path / DATE_FILENAME).write_text("2026,01,05\n")
    assert read_date_retrieved(str(tmp_path)) == "Jan 05, 2026"


def test_read_date_retrieved_returns_raw_if_unparseable(tmp_path):
    (tmp_path / DATE_FILENAME).write_text("weird-stamp\n")
    assert read_date_retrieved(str(tmp_path)) == "weird-stamp"
