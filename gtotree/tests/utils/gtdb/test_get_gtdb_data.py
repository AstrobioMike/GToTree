import io
import tarfile
from pathlib import Path
import gzip
from unittest.mock import patch
import pytest # type: ignore
import pandas as pd # type: ignore
import gtotree.utils.gtdb.get_gtdb_data as mod
from gtotree.utils.gtdb.get_gtdb_data import (
    check_gtdb_location_var_is_set,
    check_if_gtdb_data_present,
    get_slim_gtdb_tab,
    get_gtdb_data,
    gen_gtdb_tab,
    GTDB_KEPT_COLUMNS,
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
    def _dl(url, label, dest, *args, **kwargs):
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
    def boom(url, label, dest, *args, **kwargs):
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
    def write_garbage(url, label, dest, *args, **kwargs):
        Path(dest).write_bytes(b"not a gzip stream")
    with patch.object(mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=write_garbage), \
         patch(f"{MODPATH}.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path))
    mock_gen.assert_called_once_with(str(tmp_path))
    assert not (tmp_path / "GTDB-slim.tar.gz").exists()


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


# a minimal upstream GTDB metadata schema: the taxonomy column plus a handful of
# GTDB_KEPT_COLUMNS members (and one column that should be dropped by slimming).
# The first column is renamed to "accession" by gen_gtdb_tab.
BASE_COLS = [
    "accession", "gtdb_taxonomy", "ncbi_genbank_assembly_accession",
    "ncbi_taxid", "gtdb_representative", "checkm2_completeness",
    "genome_size", "some_unused_column",
]

ARC_TAX = "d__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;f__Haloferacaceae;g__Haloferax;s__Haloferax volcanii"
BAC_TAX = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"


def _row(acc, tax):
    return [acc, tax, "GCA_x", "562", "t", "99.5", "4600000", "junk"]


def _write_gz_tsv(path, rows, cols=BASE_COLS):
    lines = ["\t".join(cols)]
    lines += ["\t".join(map(str, r)) for r in rows]
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def _mock_download(arc_rows, bac_rows, arc_cols=BASE_COLS, bac_cols=BASE_COLS):
    """
    Return a download_with_tqdm side effect that writes gzipped arc/bac fixtures
    to whichever dest path is requested (keyed on the ar53/bac120 link).
    """
    def _dl(link, label, dest):
        if "ar53" in link:
            _write_gz_tsv(dest, arc_rows, arc_cols)
        elif "bac120" in link:
            _write_gz_tsv(dest, bac_rows, bac_cols)
        else:
            raise AssertionError(f"unexpected download link: {link}")
        return dest
    return _dl


def _run_gen(tmp_path, arc_rows, bac_rows, **kw):
    """Run gen_gtdb_tab with mocked downloads + VERSION fetch; return the dir."""
    with patch(f"{MODPATH}.download_with_tqdm",
               side_effect=_mock_download(arc_rows, bac_rows, **kw)), \
         patch("urllib.request.urlretrieve") as mock_ver:
        mock_ver.side_effect = lambda url, path: Path(path).write_text("r220\n")
        gen_gtdb_tab(str(tmp_path))
    return tmp_path


def _read_out(tmp_path):
    p = tmp_path / GTDB_TABLE_FILENAME
    return pd.read_csv(p, sep="\t", dtype=str)


################################################################################
# taxonomy splitting
################################################################################

def test_splits_taxonomy_into_seven_ranks(tmp_path):
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX)])
    df = _read_out(tmp_path)

    arc = df[df["accession"] == "GB_GCA_001"].iloc[0]
    assert arc["domain"] == "Archaea"
    assert arc["phylum"] == "Halobacteriota"
    assert arc["class"] == "Halobacteria"
    assert arc["order"] == "Halobacteriales"
    assert arc["family"] == "Haloferacaceae"
    assert arc["genus"] == "Haloferax"
    assert arc["species"] == "Haloferax volcanii"


def test_strips_rank_prefixes(tmp_path):
    """Each rank value must have its 'x__' prefix removed (the [3:] slice)."""
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX)])
    df = _read_out(tmp_path)
    for col in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        for val in df[col]:
            assert not val.startswith("d__")
            assert "__" not in val   # no rank prefix survived


################################################################################
# combining arc + bac
################################################################################

def test_combines_archaea_and_bacteria(tmp_path):
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX), _row("RS_GCF_002", BAC_TAX)])
    df = _read_out(tmp_path)
    accs = set(df["accession"])
    assert accs == {"GB_GCA_001", "RS_GCF_001", "RS_GCF_002"}


################################################################################
# slimming to GTDB_KEPT_COLUMNS
################################################################################

def test_output_columns_are_subset_of_kept_columns(tmp_path):
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX)])
    df = _read_out(tmp_path)
    # every output column is in the kept set (order preserved as in KEPT_COLUMNS)
    assert list(df.columns) == [c for c in GTDB_KEPT_COLUMNS if c in df.columns]
    # the unused upstream column was dropped
    assert "some_unused_column" not in df.columns
    # the raw taxonomy column was dropped in favor of the 7 split ranks
    assert "gtdb_taxonomy" not in df.columns


def test_missing_kept_column_is_skipped_not_errored(tmp_path):
    """A GTDB release lacking one of GTDB_KEPT_COLUMNS should not error; that
    column is simply absent from the output (intersect-with-present behavior)."""
    cols = [c for c in BASE_COLS if c != "checkm2_completeness"]
    def row(acc, tax):
        return [acc, tax, "GCA_x", "562", "t", "4600000", "junk"]
    _run_gen(tmp_path,
             arc_rows=[row("GB_GCA_001", ARC_TAX)],
             bac_rows=[row("RS_GCF_001", BAC_TAX)],
             arc_cols=cols, bac_cols=cols)
    df = _read_out(tmp_path)
    assert "checkm2_completeness" not in df.columns
    assert "genome_size" in df.columns          # other kept columns still present


################################################################################
# version-info file
################################################################################

def test_writes_version_info_file(tmp_path):
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX)])
    assert (tmp_path / GTDB_VERSION_FILENAME).read_text().strip() == "r220"


def test_removes_downloaded_gz_intermediates(tmp_path):
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX)])
    assert not (tmp_path / "ar53_metadata.tsv.gz").exists()
    assert not (tmp_path / "bac120_metadata.tsv.gz").exists()


################################################################################
# malformed lineage guard
################################################################################

def test_malformed_lineage_aborts(tmp_path):
    """A lineage without exactly 7 ranks should abort (sys.exit) rather than
    silently produce a misaligned table."""
    bad_tax = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria"  # only 3
    with pytest.raises(SystemExit):
        _run_gen(tmp_path,
                 arc_rows=[_row("GB_GCA_001", ARC_TAX)],
                 bac_rows=[_row("RS_GCF_bad", bad_tax)])


def test_valid_lineage_does_not_abort(tmp_path):
    """Control for the guard test: well-formed lineages complete normally."""
    _run_gen(tmp_path,
             arc_rows=[_row("GB_GCA_001", ARC_TAX)],
             bac_rows=[_row("RS_GCF_001", BAC_TAX)])
    assert (tmp_path / GTDB_TABLE_FILENAME).exists()
