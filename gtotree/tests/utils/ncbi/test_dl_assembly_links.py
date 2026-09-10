"""
Link construction for `gtt dl-ncbi-assemblies`.

The subcommand builds a per-format file URL from a directory URL plus the directory
basename. `resolve_base_link` (ftp_path present) and `build_base_link` (fallback) used
to disagree on the trailing slash, so the ftp_path path -- which is nearly every row --
produced a separator-less URL and NCBI answered 404 for every genome. These pin the
joined result rather than the trailing-slash convention itself, so either convention
stays fine as long as the final URL is right.
"""

from pathlib import Path

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.ncbi.dl_assembly_links import (parse_ncbi_assembly_summary,
                                                  _resolve_links,
                                                  FORMAT_EXTENSIONS)
from gtotree.utils.ncbi.dl_ncbi_assemblies import RunData
from gtotree.utils.taxonomy.tax_ranks import RANKS, accession_core
from gtotree.utils.taxonomy.lineage_lookup import NO_LINEAGE


_PARQUET_COLUMNS = [
    "assembly_accession", "asm_name", "taxid", "organism_name",
    "infraspecific_name", "version_status", "assembly_level", "ftp_path",
]

_FTP_DIR = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
            "GCF_000005845.2_ASM584v2")


def _row(acc="GCF_000005845.2", asm="ASM584v2", ftp=_FTP_DIR):
    return {
        "assembly_accession": acc, "asm_name": asm, "taxid": "562",
        "organism_name": "Escherichia coli", "infraspecific_name": "",
        "version_status": "latest", "assembly_level": "Complete Genome",
        "ftp_path": ftp,
    }


def _write_parquet(tmp_path, rows):
    table = pa.table({col: [str(r[col]) for r in rows] for col in _PARQUET_COLUMNS})
    path = tmp_path / "ncbi-assembly-summary.parquet"
    pq.write_table(table, path)
    return path


def _run(tmp_path, rows, accs, wanted_format="fasta"):
    parquet = _write_parquet(tmp_path, rows)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    run_data = RunData(
        wanted_format=wanted_format,
        output_dir=str(out_dir),
        wanted_accs=list(accs),
        num_wanted=len(accs),
        ncbi_sub_table_path=out_dir / "downloaded-assemblies-info.tsv",
        not_found_path=out_dir / "ncbi-accessions-not-found.txt",
        not_downloaded_path=out_dir / "ncbi-accessions-not-downloaded.tsv",
    )
    parse_ncbi_assembly_summary(parquet, run_data)
    lines = Path(run_data.ncbi_sub_table_path).read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def test_target_link_has_separator_between_dir_and_filename(tmp_path):
    """The ftp_path branch: the bug that made every download 404."""
    (record,) = _run(tmp_path, [_row()], ["GCF_000005845.2"])

    assert record["target_link"] == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
        "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    )
    assert "ASM584v2GCF" not in record["target_link"]


def test_target_link_built_when_ftp_path_missing(tmp_path):
    """The fallback branch must produce the identical URL."""
    (record,) = _run(tmp_path, [_row(ftp="na")], ["GCF_000005845.2"])

    assert record["target_link"] == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
        "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    )


def test_ftp_path_with_trailing_slash_does_not_double_up(tmp_path):
    (record,) = _run(tmp_path, [_row(ftp=_FTP_DIR + "/")], ["GCF_000005845.2"])

    assert "//GCF_000005845.2_ASM584v2_genomic" not in record["target_link"]
    assert record["target_link"].endswith(
        "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz")


@pytest.mark.parametrize("wanted_format", sorted(FORMAT_EXTENSIONS))
def test_every_format_gets_a_well_formed_link(tmp_path, wanted_format):
    ncbi_ext, local_ext = FORMAT_EXTENSIONS[wanted_format]
    (record,) = _run(tmp_path, [_row()], ["GCF_000005845.2"],
                     wanted_format=wanted_format)

    assert record["target_link"] == f"{_FTP_DIR}/GCF_000005845.2_ASM584v2{ncbi_ext}"
    assert record["local_destination"].endswith(f"GCF_000005845.2{local_ext}")


def test_both_branches_agree_on_trailing_slash_convention():
    """
    Whatever the convention, the two branches must share it -- disagreeing is what
    let the separator go missing on only one of them.
    """
    from_ftp, base_from_ftp = _resolve_links(
        "GCF_000005845.2", "ASM584v2", _FTP_DIR)
    built, base_built = _resolve_links("GCF_000005845.2", "ASM584v2", "")

    assert from_ftp.endswith("/") == built.endswith("/")
    assert from_ftp == built
    assert base_from_ftp == base_built == "GCF_000005845.2_ASM584v2"


def test_unresolvable_row_yields_na_link(tmp_path):
    (record,) = _run(tmp_path, [_row(asm="", ftp="na")], ["GCF_000005845.2"])

    assert record["target_link"] == "NA"


# --- lineage columns in the info table ------------------------------------
# --add-ncbi-tax / --add-gtdb-tax are independent of each other and of --source, and
# both default off. NCBI lineage rides along on the asset scan; GTDB lineage comes
# from a map the caller builds and hangs off run_data.

_LINEAGE = ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
            "Enterobacteriaceae", "Escherichia", "Escherichia coli")


def _lineage_summary(tmp_path, with_ranks=True):
    """A one-row Parquet fixture, optionally carrying the asset's lineage columns."""
    row = _row()
    cols = {c: [str(row[c])] for c in _PARQUET_COLUMNS}
    if with_ranks:
        for rank, value in zip(RANKS, _LINEAGE):
            cols[rank] = [value]
    path = tmp_path / "ncbi-summary-lineage.parquet"
    pq.write_table(pa.table(cols), path)
    return path


def _lineage_run_data(tmp_path):
    out_dir = tmp_path / "lineage-out"
    out_dir.mkdir()
    return RunData(
        wanted_format="fasta",
        output_dir=str(out_dir),
        wanted_accs=["GCF_000005845.2"],
        num_wanted=1,
        ncbi_sub_table_path=out_dir / "downloaded-assemblies-info.tsv",
        not_found_path=out_dir / "ncbi-accessions-not-found.txt",
        not_downloaded_path=out_dir / "ncbi-accessions-not-downloaded.tsv",
    )


_lineage_parse = parse_ncbi_assembly_summary


def _lineage_table(tmp_path, add_ncbi_tax=False, add_gtdb_tax=False,
                   gtdb_lineage=None, with_ranks=True):
    summary = _lineage_summary(tmp_path, with_ranks=with_ranks)
    rd = _lineage_run_data(tmp_path)
    rd.add_ncbi_tax = add_ncbi_tax
    rd.add_gtdb_tax = add_gtdb_tax
    rd.gtdb_lineage = gtdb_lineage
    _lineage_parse(summary, rd)
    lines = Path(rd.ncbi_sub_table_path).read_text().splitlines()
    header = lines[0].split("\t")
    return header, [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def test_no_lineage_columns_by_default(tmp_path):
    header, _ = _lineage_table(tmp_path)
    assert not [c for c in header if c.startswith(("ncbi_", "gtdb_"))]


def test_add_ncbi_tax_adds_prefixed_ncbi_columns(tmp_path):
    header, rows = _lineage_table(tmp_path, add_ncbi_tax=True)
    assert header[-7:] == [f"ncbi_{r}" for r in RANKS]
    assert rows[0]["ncbi_species"] == "Escherichia coli"
    assert rows[0]["ncbi_domain"] == "Bacteria"
    assert not [c for c in header if c.startswith("gtdb_")]


def test_add_gtdb_tax_adds_prefixed_gtdb_columns(tmp_path):
    mapping = {accession_core("GCF_000005845.2"): _LINEAGE}
    header, rows = _lineage_table(tmp_path, add_gtdb_tax=True, gtdb_lineage=mapping,
                                  with_ranks=False)
    assert header[-7:] == [f"gtdb_{r}" for r in RANKS]
    assert rows[0]["gtdb_species"] == "Escherichia coli"
    assert not [c for c in header if c.startswith("ncbi_")]


def test_both_taxonomies_can_be_on_at_once(tmp_path):
    """
    They're independent flags, and the prefixes are what keeps a mixed table
    unambiguous when GTDB and NCBI disagree.
    """
    gtdb = ("Bacteria", "GtdbPhylum", "GtdbClass", "GtdbOrder", "GtdbFamily",
            "GtdbGenus", "Gtdb species")
    mapping = {accession_core("GCF_000005845.2"): gtdb}
    header, rows = _lineage_table(tmp_path, add_ncbi_tax=True, add_gtdb_tax=True,
                                  gtdb_lineage=mapping)
    assert header[-14:] == ([f"ncbi_{r}" for r in RANKS] +
                            [f"gtdb_{r}" for r in RANKS])
    assert rows[0]["ncbi_phylum"] == "Pseudomonadota"
    assert rows[0]["gtdb_phylum"] == "GtdbPhylum"


def test_accession_missing_from_gtdb_gets_NA(tmp_path):
    """GTDB is bacteria/archaea only, so a miss is expected, not an error."""
    header, rows = _lineage_table(tmp_path, add_gtdb_tax=True, gtdb_lineage={},
                                  with_ranks=False)
    assert [rows[0][f"gtdb_{r}"] for r in RANKS] == list(NO_LINEAGE)


def test_lineage_columns_come_after_the_link_columns(tmp_path):
    header, _ = _lineage_table(tmp_path, add_ncbi_tax=True)
    assert header.index("local_destination") < header.index("ncbi_domain")
    assert header.index("target_link") < header.index("ncbi_domain")


def test_every_row_has_the_full_width(tmp_path):
    mapping = {accession_core("GCF_000005845.2"): _LINEAGE}
    header, rows = _lineage_table(tmp_path, add_ncbi_tax=True, add_gtdb_tax=True,
                                  gtdb_lineage=mapping)
    assert all(len(row) == len(header) for row in rows)
    assert "" not in rows[0].values()
