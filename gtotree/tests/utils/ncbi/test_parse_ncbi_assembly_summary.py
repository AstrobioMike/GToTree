from pathlib import Path
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore

from gtotree.tests.utils.utils import FakeRunData
from gtotree.utils.ncbi.parse_ncbi_assembly_summary import (
    parse_assembly_summary,
    build_base_link,
    resolve_base_link,
    sanitize_assembly_name,
)


# columns the parser reads from the hosted NCBI Parquet asset
_PARQUET_COLUMNS = [
    "assembly_accession", "asm_name", "taxid", "organism_name",
    "infraspecific_name", "version_status", "assembly_level", "ftp_path",
]


def _row(acc, taxid="562", org="Escherichia coli", infra="", ver="latest",
         level="Complete Genome", asm="ASM584v2",
         ftp="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"):
    return {
        "assembly_accession": acc, "asm_name": asm, "taxid": taxid,
        "organism_name": org, "infraspecific_name": infra,
        "version_status": ver, "assembly_level": level, "ftp_path": ftp,
    }


def _write_parquet(path, rows):
    """Write a schema-faithful slim NCBI Parquet asset with just the read columns."""
    cols = {c: [r.get(c, "") for r in rows] for c in _PARQUET_COLUMNS}
    pq.write_table(
        pa.table({c: pa.array(cols[c], type=pa.string()) for c in _PARQUET_COLUMNS}),
        str(path))



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
# parse: reading the Parquet asset
################################################################################

def test_parse_extracts_correct_columns(tmp_path):
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2", taxid="562",
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


def test_parse_rebuilds_link_when_ftp_na(tmp_path):
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2", asm="ASM584v2", ftp="na")])
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
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2")])
    rd = FakeRunData(["GCF_000005845.2", "GCF_999999999.1"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    by_id = {gd.id: gd for gd in rd.ncbi_accs}
    assert by_id["GCF_000005845.2"].acc_was_found is True
    assert by_id["GCF_999999999.1"].acc_was_found is False
    assert by_id["GCF_999999999.1"].removed is True

    rows = Path(tmp_path / "removed-genomes.tsv").read_text().splitlines()
    assert [r.split("\t")[0] for r in rows[1:]] == ["GCF_999999999.1"]


def test_parse_matches_on_root_accession_ignoring_version(tmp_path):
    """A wanted accession should match the summary row even if the version
    suffix differs (matching is on the pre-dot root)."""
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.3")])   # .3 in table
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)         # .2 wanted

    parse_assembly_summary(str(summary), rd)

    _, rows = _read_subtable(rd)
    assert rows[0]["input_accession"] == "GCF_000005845.2"   # original preserved
    assert rows[0]["found_accession"] == "GCF_000005845.3"   # what NCBI had


def test_parse_idempotent_if_subtable_already_set(tmp_path):
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2")])
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)
    rd.ncbi_sub_table_path = "already/done.tsv"

    out = parse_assembly_summary(str(summary), rd)
    assert out.ncbi_sub_table_path == "already/done.tsv"   # unchanged, no work


################################################################################
# taxid propagation onto GenomeData
################################################################################

def test_parse_records_the_taxid_on_the_genome_data(tmp_path):
    """
    The sub-table lives in the working directory, which a finished run deletes, so
    the summary table reads the taxid off the GenomeData instead -- where it also
    survives into run-data.json and is still there for `gtt update-headers`.
    """
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2", taxid="511145")])
    rd = FakeRunData(["GCF_000005845.2"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    by_id = {gd.id: gd for gd in rd.ncbi_accs}
    assert by_id["GCF_000005845.2"].taxid == "511145"


def test_the_taxid_on_genome_data_matches_the_sub_table(tmp_path):
    """
    The two are written from the same row, and the summary table's taxid column
    switched from one source to the other -- so they must not drift.
    """
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2", taxid="511145"),
                             _row("GCA_000009065.1", taxid="83333")])
    rd = FakeRunData(["GCF_000005845.2", "GCA_000009065.1"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    _header, rows = _read_subtable(rd)
    from_sub_table = {r["input_accession"]: r["taxid"] for r in rows}
    from_genome_data = {gd.id: gd.taxid for gd in rd.ncbi_accs}

    assert from_genome_data == from_sub_table


def test_a_not_found_accession_has_no_taxid(tmp_path):
    summary = tmp_path / "ncbi-data.parquet"
    _write_parquet(summary, [_row("GCF_000005845.2")])
    rd = FakeRunData(["GCF_000005845.2", "GCF_999999999.1"], tmp_path)

    parse_assembly_summary(str(summary), rd)

    by_id = {gd.id: gd for gd in rd.ncbi_accs}
    assert by_id["GCF_999999999.1"].taxid is None
