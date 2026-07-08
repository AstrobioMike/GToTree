from pathlib import Path
import pytest # type: ignore

from gtotree.tests.utils.utils import FakeRunData
from gtotree.utils.ncbi.parse_assembly_summary_file import (
    parse_assembly_summary,
    build_base_link,
    resolve_base_link,
    sanitize_assembly_name,
)


SLIM_COLUMNS = [
    "assembly_accession", "taxid", "organism_name", "infraspecific_name",
    "version_status", "assembly_level", "asm_name", "ftp_path",
]

# NCBI full-file positions used to build a headerless legacy fixture row
LEGACY_N = 38


def _slim_row(acc, taxid="562", org="Escherichia coli", infra="", ver="latest",
              level="Complete Genome", asm="ASM584v2",
              ftp="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"):
    return "\t".join([acc, taxid, org, infra, ver, level, asm, ftp])


def _write_slim(path, rows):
    Path(path).write_text("\t".join(SLIM_COLUMNS) + "\n" + "\n".join(rows) + "\n")


def _legacy_row(acc, taxid="562", org="Escherichia coli", infra="", ver="latest",
                level="Complete Genome", asm="ASM584v2",
                ftp="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"):
    """Full 38-col NCBI row with values in the legacy positions, no header."""
    f = [f"c{i}" for i in range(LEGACY_N)]
    f[0] = acc; f[5] = taxid; f[7] = org; f[8] = infra
    f[10] = ver; f[11] = level; f[15] = asm; f[19] = ftp
    return "\t".join(f)


def _read_subtable(rd):
    lines = Path(rd.ncbi_sub_table_path).read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    return header, rows


################################################################################
# sanitize / build / resolve base link (pure functions)
################################################################################

def test_sanitize_assembly_name_replaces_specials():
    assert sanitize_assembly_name("ASM 9 (weird/name),#x") == "ASM_9_weird_name_x"


def test_build_base_link_nine_digit_accession():
    url, base = build_base_link("GCF_000005845.2", "ASM584v2")
    assert url == ("https://ftp.ncbi.nlm.nih.gov/genomes/all/"
                   "GCF/000/005/845/GCF_000005845.2_ASM584v2/")
    assert base == "GCF_000005845.2_ASM584v2"


def test_resolve_base_link_prefers_real_ftp_normalizes_https():
    out = resolve_base_link(
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2",
        "GCF_000005845.2", "ASM584v2")
    assert out == ("https://ftp.ncbi.nlm.nih.gov/genomes/all/"
                   "GCF/000/005/845/GCF_000005845.2_ASM584v2")
    assert not out.endswith("/")


@pytest.mark.parametrize("missing", ["na", "NA", "Na", ""])
def test_resolve_base_link_rebuilds_when_ftp_missing(missing):
    out = resolve_base_link(missing, "GCF_000005845.2", "ASM584v2")
    assert out == ("https://ftp.ncbi.nlm.nih.gov/genomes/all/"
                   "GCF/000/005/845/GCF_000005845.2_ASM584v2")


def test_resolve_base_link_na_when_unrebuildable():
    assert resolve_base_link("na", "na", "na") == "na"


################################################################################
# parse: header-aware (slim) path
################################################################################

def test_parse_slim_header_extracts_correct_columns(tmp_path):
    summary = tmp_path / "slim.tsv"
    _write_slim(summary, [_slim_row("GCF_000005845.2", taxid="562",
                                    org="Escherichia coli", infra="strain=K-12")])
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    header, rows = _read_subtable(rd)
    assert header == ["input_accession", "found_accession", "assembly_name",
                      "taxid", "organism_name", "infraspecific_name",
                      "version_status", "assembly_level", "http_base_link"]
    r = rows[0]
    assert r["found_accession"] == "GCF_000005845.2"
    assert r["assembly_name"] == "ASM584v2"
    assert r["taxid"] == "562"
    assert r["organism_name"] == "Escherichia coli"
    assert r["infraspecific_name"] == "strain=K-12"
    assert r["http_base_link"].endswith("GCF_000005845.2_ASM584v2")
    assert not r["http_base_link"].endswith("/")   # canonical: no trailing slash


def test_parse_slim_and_legacy_agree(tmp_path):
    """The header-aware and headerless-legacy paths must produce identical
    output for the same underlying row (this is the whole point of the dual
    path -- the slim table and the raw NCBI file are interchangeable)."""
    slim = tmp_path / "slim.tsv"
    legacy = tmp_path / "legacy.txt"
    _write_slim(slim, [_slim_row("GCF_000005845.2")])
    Path(legacy).write_text(_legacy_row("GCF_000005845.2") + "\n")

    rd_slim = FakeRunData(["GCF_000005845.2"], tmp_path / "a")
    (tmp_path / "a").mkdir()
    rd_legacy = FakeRunData(["GCF_000005845.2"], tmp_path / "b")
    (tmp_path / "b").mkdir()

    parse_assembly_summary(str(slim), rd_slim)
    parse_assembly_summary(str(legacy), rd_legacy)

    _, rows_slim = _read_subtable(rd_slim)
    _, rows_legacy = _read_subtable(rd_legacy)
    assert rows_slim == rows_legacy


def test_parse_rebuilds_link_when_ftp_na(tmp_path):
    summary = tmp_path / "slim.tsv"
    _write_slim(summary, [_slim_row("GCF_000005845.2", asm="ASM584v2", ftp="na")])
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    _, rows = _read_subtable(rd)
    # rebuilt from acc + assembly name, 9-digit split, no trailing slash
    assert rows[0]["http_base_link"] == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCF/000/005/845/GCF_000005845.2_ASM584v2")


################################################################################
# parse: found / not-found bookkeeping
################################################################################

def test_parse_marks_found_and_not_found(tmp_path):
    summary = tmp_path / "slim.tsv"
    _write_slim(summary, [_slim_row("GCF_000005845.2")])
    rd = FakeRunData(["GCF_000005845.2", "GCF_999999999.1"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    by_id = {gd.id: gd for gd in rd.ncbi_accs}
    assert by_id["GCF_000005845.2"].acc_was_found is True
    assert by_id["GCF_999999999.1"].acc_was_found is False
    assert by_id["GCF_999999999.1"].removed is True

    not_found = Path(tmp_path / "ncbi-accessions-not-found.txt").read_text().split()
    assert not_found == ["GCF_999999999.1"]


def test_parse_matches_on_root_accession_ignoring_version(tmp_path):
    """A wanted accession should match the summary row even if the version
    suffix differs (matching is on the pre-dot root)."""
    summary = tmp_path / "slim.tsv"
    _write_slim(summary, [_slim_row("GCF_000005845.3")])   # .3 in table
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)         # .2 wanted

    parse_assembly_summary(str(summary), rd)

    _, rows = _read_subtable(rd)
    assert rows[0]["input_accession"] == "GCF_000005845.2"   # original preserved
    assert rows[0]["found_accession"] == "GCF_000005845.3"   # what NCBI had


def test_parse_idempotent_if_subtable_already_set(tmp_path):
    summary = tmp_path / "slim.tsv"
    _write_slim(summary, [_slim_row("GCF_000005845.2")])
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)
    rd.ncbi_sub_table_path = "already/done.tsv"

    out = parse_assembly_summary(str(summary), rd)
    assert out.ncbi_sub_table_path == "already/done.tsv"   # unchanged, no work
