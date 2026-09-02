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
        ncbi_sub_table_path=out_dir / "wanted-ncbi-accessions-info.tsv",
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
