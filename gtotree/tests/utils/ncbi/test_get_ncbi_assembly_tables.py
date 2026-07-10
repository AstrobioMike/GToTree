import io
import tarfile
from datetime import date
from pathlib import Path
from unittest.mock import patch
import pytest # type: ignore

import gtotree.utils.ncbi.get_ncbi_assembly_tables as mod
from gtotree.utils.ncbi.get_ncbi_assembly_tables import (
    check_ncbi_assembly_info_location_var_is_set,
    check_if_data_present,
    download_ncbi_assembly_summary_data,
    get_slim_ncbi_assembly_data,
    get_ncbi_assembly_data,
    get_ncbi_assembly_summary_tab,
)


NCBI_COLUMNS = [
    "assembly_accession", "bioproject", "biosample", "wgs_master",
    "refseq_category", "taxid", "species_taxid", "organism_name",
    "infraspecific_name", "isolate", "version_status", "assembly_level",
    "release_type", "genome_rep", "seq_rel_date", "asm_name", "asm_submitter",
    "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq",
    "relation_to_type_material", "asm_not_live_date", "assembly_type", "group",
    "genome_size", "genome_size_ungapped", "gc_percent", "replicon_count",
    "scaffold_count", "contig_count", "annotation_provider", "annotation_name",
    "annotation_date", "total_gene_count", "protein_coding_gene_count",
    "non_coding_gene_count", "pubmed_id",
]
SLIM_COLUMNS = [
    "assembly_accession", "taxid", "organism_name", "infraspecific_name",
    "version_status", "assembly_level", "asm_name", "ftp_path",
]

MODPATH = "gtotree.utils.ncbi.get_ncbi_assembly_tables"


def _ncbi_row(acc):
    f = [f"c{i}" for i in range(38)]
    f[0] = acc; f[5] = "562"; f[7] = "Escherichia coli"; f[15] = "ASM"
    f[19] = "https://ftp.ncbi.nlm.nih.gov/genomes/all/X/" + acc
    return "\t".join(f)


def _write_ncbi_summary(path, accs):
    Path(path).write_text(
        "##provenance\n#" + "\t".join(NCBI_COLUMNS) + "\n" +
        "\n".join(_ncbi_row(a) for a in accs) + "\n")


def _mock_ncbi_download(link, label, dest):
    """download_with_tqdm side effect writing fake NCBI-format summaries."""
    if "genbank" in link:
        _write_ncbi_summary(dest, ["GCA_000001.1", "GCA_000002.1"])
    else:
        _write_ncbi_summary(dest, ["GCF_000001.1"])


def _build_tarball(path, *, include_table=True, include_date=True,
                   extra_member=None):
    """Build a slim assembly-info tarball like the hosted asset."""
    table = "\t".join(SLIM_COLUMNS) + "\nGCA_000001.1\t562\tE. coli\tNA\tlatest\tComplete\tASM\tna\n"
    datestr = "2026,06,25\n"
    with tarfile.open(path, "w:gz") as tar:
        def _add(name, content):
            b = content.encode(); ti = tarfile.TarInfo(name); ti.size = len(b)
            tar.addfile(ti, io.BytesIO(b))
        if include_table:
            _add("ncbi-assembly-info.tsv", table)
        if include_date:
            _add("date-retrieved.txt", datestr)
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
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    assert check_ncbi_assembly_info_location_var_is_set() == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("NCBI_assembly_data_dir", raising=False)
    with pytest.raises(SystemExit):
        check_ncbi_assembly_info_location_var_is_set()


def test_summary_tab_accessor_is_lazy(monkeypatch, tmp_path):
    """get_ncbi_assembly_summary_tab resolves on call, not at import; with the
    env var set it returns <dir>/ncbi-assembly-info.tsv."""
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    assert get_ncbi_assembly_summary_tab() == str(tmp_path / "ncbi-assembly-info.tsv")


def test_check_if_data_present_both_files_nonempty(tmp_path):
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    assert check_if_data_present(str(tmp_path)) is True


def test_check_if_data_present_table_missing(tmp_path):
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "date-retrieved.txt").exists()


def test_check_if_data_present_date_file_missing(tmp_path):
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "ncbi-assembly-info.tsv").exists()


def test_check_if_data_present_both_missing(tmp_path):
    assert check_if_data_present(str(tmp_path)) is False


def test_check_if_data_present_empty_files_removed(tmp_path):
    (tmp_path / "ncbi-assembly-info.tsv").write_text("")
    (tmp_path / "date-retrieved.txt").write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "ncbi-assembly-info.tsv").exists()
    assert not (tmp_path / "date-retrieved.txt").exists()


################################################################################
# rebuild path (download_ncbi_assembly_summary_data) -> slim format
################################################################################

def test_rebuild_produces_slim_combined_table(tmp_path):
    with patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_ncbi_download):
        download_ncbi_assembly_summary_data(str(tmp_path))
    out = (tmp_path / "ncbi-assembly-info.tsv").read_text().splitlines()
    assert out[0].split("\t") == SLIM_COLUMNS
    accs = [l.split("\t")[0] for l in out[1:]]
    assert "GCA_000001.1" in accs and "GCA_000002.1" in accs   # genbank
    assert "GCF_000001.1" in accs                               # refseq


def test_rebuild_removes_temp_files(tmp_path):
    with patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_ncbi_download):
        download_ncbi_assembly_summary_data(str(tmp_path))
    assert not (tmp_path / "assembly_summary_genbank.txt").exists()
    assert not (tmp_path / "assembly_summary_refseq.txt").exists()


def test_rebuild_writes_date_file(tmp_path):
    fixed = date(2024, 6, 15)
    with patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_ncbi_download), \
         patch(f"{MODPATH}.date") as mock_date:
        mock_date.today.return_value = fixed
        download_ncbi_assembly_summary_data(str(tmp_path))
    assert (tmp_path / "date-retrieved.txt").read_text().strip() == "2024,06,15"


def test_rebuild_download_failure_exits(tmp_path):
    with patch(f"{MODPATH}.download_with_tqdm", side_effect=Exception("network error")), \
         patch(f"{MODPATH}.report_early_exit") as mock_exit:
        mock_exit.side_effect = SystemExit(1)
        with pytest.raises(SystemExit):
            download_ncbi_assembly_summary_data(str(tmp_path))
    mock_exit.assert_called_once()


################################################################################
# fast path (get_slim_ncbi_assembly_data)
################################################################################

def test_slim_path_extracts_both_files(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_tarball(src, extra_member="README-junk.txt")
    with patch.object(mod, "NCBI_ASSEMBLY_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_tarball_download(src)):
        get_slim_ncbi_assembly_data(str(tmp_path))
    assert (tmp_path / "ncbi-assembly-info.tsv").exists()
    assert (tmp_path / "date-retrieved.txt").read_text().startswith("2026,06,25")
    assert not (tmp_path / "README-junk.txt").exists()          # only wanted files
    assert not (tmp_path / "ncbi-assembly-info.tar.gz").exists()  # tarball cleaned up


def test_slim_path_no_url_falls_back_to_rebuild(tmp_path):
    with patch.object(mod, "NCBI_ASSEMBLY_TARBALL_URL", ""), \
         patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb:
        get_slim_ncbi_assembly_data(str(tmp_path))
    mock_rb.assert_called_once_with(str(tmp_path))


def test_slim_path_download_error_falls_back(tmp_path):
    def boom(url, label, dest, *args, **kwargs):
        raise ConnectionError("server down")
    with patch.object(mod, "NCBI_ASSEMBLY_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=boom), \
         patch(f"{MODPATH}.time.sleep"), \
         patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb:
        get_slim_ncbi_assembly_data(str(tmp_path))
    mock_rb.assert_called_once_with(str(tmp_path))
    assert not (tmp_path / "ncbi-assembly-info.tar.gz").exists()


def test_slim_path_missing_file_in_archive_falls_back(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_tarball(src, include_date=False)   # table only, no date file
    with patch.object(mod, "NCBI_ASSEMBLY_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=_mock_tarball_download(src)), \
         patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb:
        get_slim_ncbi_assembly_data(str(tmp_path))
    mock_rb.assert_called_once_with(str(tmp_path))


def test_slim_path_corrupt_tarball_falls_back(tmp_path):
    def write_garbage(url, label, dest, *args, **kwargs):
        Path(dest).write_bytes(b"not a gzip stream")
    with patch.object(mod, "NCBI_ASSEMBLY_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch(f"{MODPATH}.download_with_tqdm", side_effect=write_garbage), \
         patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb:
        get_slim_ncbi_assembly_data(str(tmp_path))
    mock_rb.assert_called_once_with(str(tmp_path))
    assert not (tmp_path / "ncbi-assembly-info.tar.gz").exists()


################################################################################
# routing (get_ncbi_assembly_data)
################################################################################

def test_get_data_already_present_skips_everything(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    with patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb, \
         patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_slim:
        get_ncbi_assembly_data(force_update=False)
    mock_rb.assert_not_called()
    mock_slim.assert_not_called()


def test_get_data_default_uses_slim_path(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_slim, \
         patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb:
        get_ncbi_assembly_data()
    mock_slim.assert_called_once_with(str(tmp_path))
    mock_rb.assert_not_called()


def test_get_data_force_update_forces_slim_fetch(tmp_path, monkeypatch):
    """-f re-fetches the hosted asset even when local data already exists; it
    does not rebuild from NCBI directly (that's only the fallback)."""
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_slim, \
         patch(f"{MODPATH}.download_ncbi_assembly_summary_data") as mock_rb:
        get_ncbi_assembly_data(force_update=True)
    mock_slim.assert_called_once_with(str(tmp_path))
    mock_rb.assert_not_called()
