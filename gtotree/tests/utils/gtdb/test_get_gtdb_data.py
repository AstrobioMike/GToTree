import io
import tarfile
from pathlib import Path
from unittest.mock import patch
import pytest # type: ignore

import gtotree.utils.gtdb.get_gtdb_data as mod
from gtotree.utils.gtdb.get_gtdb_data import (
    check_gtdb_location_var_is_set,
    check_if_gtdb_data_present,
    get_slim_gtdb_tab,
    get_gtdb_data,
    GTDB_TABLE_FILENAME,
    GTDB_VERSION_FILENAME,
)


MODPATH = "gtotree.utils.gtdb.get_gtdb_data"


def _build_tarball(path, *, include_table=True, include_version=True,
                   extra_member=None):
    table = "accession\tdomain\nGB_GCA_1\tBacteria\n"
    version = "r220\n"
    with tarfile.open(path, "w:gz") as tar:
        def _add(name, content):
            b = content.encode(); ti = tarfile.TarInfo(name); ti.size = len(b)
            tar.addfile(ti, io.BytesIO(b))
        if include_table:
            _add(GTDB_TABLE_FILENAME, table)
        if include_version:
            _add(GTDB_VERSION_FILENAME, version)
        if extra_member:
            _add(extra_member, "junk\n")


def _mock_tarball_download(src):
    def _dl(url, label, dest):
        import shutil
        shutil.copy(src, dest)
        return dest
    return _dl


################################################################################
# env var + presence checks
################################################################################

def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("GTDB_dir", str(tmp_path))
    assert check_gtdb_location_var_is_set() == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("GTDB_dir", raising=False)
    with pytest.raises(SystemExit):
        check_gtdb_location_var_is_set()


def test_check_present_both_files_nonempty(tmp_path):
    (tmp_path / GTDB_TABLE_FILENAME).write_text("data")
    (tmp_path / GTDB_VERSION_FILENAME).write_text("r220")
    assert check_if_gtdb_data_present(str(tmp_path)) is True


def test_check_present_table_missing(tmp_path):
    (tmp_path / GTDB_VERSION_FILENAME).write_text("r220")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / GTDB_VERSION_FILENAME).exists()


def test_check_present_version_missing(tmp_path):
    (tmp_path / GTDB_TABLE_FILENAME).write_text("data")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / GTDB_TABLE_FILENAME).exists()


def test_check_present_empty_files_removed(tmp_path):
    (tmp_path / GTDB_TABLE_FILENAME).write_text("")
    (tmp_path / GTDB_VERSION_FILENAME).write_text("")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / GTDB_TABLE_FILENAME).exists()
    assert not (tmp_path / GTDB_VERSION_FILENAME).exists()


################################################################################
# fast path (get_slim_gtdb_tab)
################################################################################

def test_slim_path_extracts_both_files(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_tarball(src, extra_member="README-junk.txt")
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_tarball_download(src)):
        get_slim_gtdb_tab(str(tmp_path))
    assert (tmp_path / GTDB_TABLE_FILENAME).exists()
    assert (tmp_path / GTDB_VERSION_FILENAME).read_text().startswith("r220")
    assert not (tmp_path / "README-junk.txt").exists()
    assert not (tmp_path / "GTDB-slim.tar.gz").exists()


def test_slim_path_no_url_falls_back(tmp_path):
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", ""), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    mock_gen.assert_called_once_with(str(tmp_path))


def test_slim_path_download_error_falls_back(tmp_path):
    def boom(url, label, dest):
        raise ConnectionError("server down")
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=boom), \
         patch(f"{MODPATH}.time.sleep"), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    mock_gen.assert_called_once_with(str(tmp_path))
    assert not (tmp_path / "GTDB-slim.tar.gz").exists()


def test_slim_path_missing_file_in_archive_falls_back(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_tarball(src, include_version=False)
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_tarball_download(src)), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    mock_gen.assert_called_once_with(str(tmp_path))


def test_slim_path_corrupt_tarball_falls_back(tmp_path):
    def write_garbage(url, label, dest):
        Path(dest).write_bytes(b"not a gzip stream")
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=write_garbage), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    mock_gen.assert_called_once_with(str(tmp_path))
    assert not (tmp_path / "GTDB-slim.tar.gz").exists()


def test_slim_path_retries_then_succeeds(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_tarball(src)
    import shutil
    import socket
    seen = {"n": 0}
    def flaky(url, label, dest):
        if seen["n"] < 2:
            seen["n"] += 1
            raise socket.timeout("stall")
        shutil.copy(src, dest)
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=flaky), \
         patch(f"{MODPATH}.time.sleep"), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    mock_gen.assert_not_called()
    assert (tmp_path / GTDB_TABLE_FILENAME).exists()


def test_slim_path_404_not_retried(tmp_path):
    import urllib.error
    calls = {"n": 0}
    def not_found(url, label, dest):
        calls["n"] += 1
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=not_found), \
         patch(f"{MODPATH}.time.sleep"), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    assert calls["n"] == 1
    mock_gen.assert_called_once_with(str(tmp_path))


################################################################################
# routing (get_gtdb_data)
################################################################################

def test_get_data_already_present_skips(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_dir", str(tmp_path))
    (tmp_path / GTDB_TABLE_FILENAME).write_text("data")
    (tmp_path / GTDB_VERSION_FILENAME).write_text("r220")
    with patch(f"{MODPATH}.get_slim_gtdb_tab") as mock_slim:
        result = get_gtdb_data(force_update=False)
    mock_slim.assert_not_called()
    assert result == str(tmp_path)


def test_get_data_default_uses_slim_path(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_dir", str(tmp_path))
    with patch(f"{MODPATH}.get_slim_gtdb_tab") as mock_slim:
        get_gtdb_data()
    mock_slim.assert_called_once_with(str(tmp_path))


def test_get_data_force_update_refetches(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_dir", str(tmp_path))
    (tmp_path / GTDB_TABLE_FILENAME).write_text("data")
    (tmp_path / GTDB_VERSION_FILENAME).write_text("r220")
    with patch(f"{MODPATH}.get_slim_gtdb_tab") as mock_slim:
        get_gtdb_data(force_update=True)
    mock_slim.assert_called_once_with(str(tmp_path))
