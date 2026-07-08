import gzip
import pytest # type: ignore
from pathlib import Path

from gtotree.utils.ncbi.slim_ncbi_assembly_summary import (
    KEPT_COLUMNS,
    parse_ncbi_header,
    build_slim_assembly_summary,
    write_date_retrieved,
)


# the real NCBI assembly_summary 38-column header (clean names, no leading #)
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


def _row(accession, taxid="562", org="Escherichia coli", infra="strain=K-12",
         version="latest", level="Complete Genome", asm_name="ASM584v2",
         ftp="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"):
    """Build a full 38-field NCBI data row with sentinel values in unused cols."""
    fields = [f"col{i}" for i in range(len(NCBI_COLUMNS))]
    vals = {
        "assembly_accession": accession, "taxid": taxid, "organism_name": org,
        "infraspecific_name": infra, "version_status": version,
        "assembly_level": level, "asm_name": asm_name, "ftp_path": ftp,
    }
    for name, v in vals.items():
        fields[NCBI_COLUMNS.index(name)] = v
    return "\t".join(fields)


def _write_summary(path, accessions, gzip_it=False, reorder=False):
    """
    Write a fake NCBI assembly_summary file: a '##' provenance line, the
    '#'-prefixed header, then one data row per accession. If reorder, shuffle the
    column order (and rows) to prove name-based selection is order-independent.
    """
    cols = list(NCBI_COLUMNS)
    rows = [_row(a) for a in accessions]
    if reorder:
        import random
        order = list(range(len(cols)))
        random.Random(0).shuffle(order)
        cols = [cols[i] for i in order]
        new_rows = []
        for r in rows:
            f = r.split("\t")
            new_rows.append("\t".join(f[i] for i in order))
        rows = new_rows
    text = "##  See ftp README\n#" + "\t".join(cols) + "\n" + "\n".join(rows) + "\n"
    if gzip_it:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)


def test_parse_ncbi_header_strips_hash():
    names = parse_ncbi_header("#assembly_accession\ttaxid\torganism_name\n")
    assert names == ["assembly_accession", "taxid", "organism_name"]


def test_parse_ncbi_header_rejects_non_header():
    with pytest.raises(ValueError):
        parse_ncbi_header("GCF_000005845.2\t562\tEscherichia coli\n")


def test_build_slim_keeps_only_kept_columns_in_order(tmp_path):
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    _write_summary(gb, ["GCA_000005845.2"])
    _write_summary(rs, ["GCF_000005845.2"])

    n = build_slim_assembly_summary(gb, rs, out)
    assert n == 2

    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == KEPT_COLUMNS          # clean header, no '#'
    assert not lines[0].startswith("#")
    for line in lines[1:]:
        assert len(line.split("\t")) == len(KEPT_COLUMNS)


def test_build_slim_combines_genbank_and_refseq(tmp_path):
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    _write_summary(gb, ["GCA_000005845.2", "GCA_000001405.29"])
    _write_summary(rs, ["GCF_000005845.2"])

    n = build_slim_assembly_summary(gb, rs, out)
    assert n == 3
    accs = [l.split("\t")[0] for l in out.read_text().splitlines()[1:]]
    assert "GCA_000005845.2" in accs
    assert "GCA_000001405.29" in accs
    assert "GCF_000005845.2" in accs


def test_build_slim_values_land_in_right_columns(tmp_path):
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    _write_summary(gb, ["GCA_000005845.2"])
    _write_summary(rs, [])  # empty refseq (header only)

    build_slim_assembly_summary(gb, rs, out)
    header, row = out.read_text().splitlines()[:2]
    rec = dict(zip(header.split("\t"), row.split("\t")))
    assert rec["assembly_accession"] == "GCA_000005845.2"
    assert rec["taxid"] == "562"
    assert rec["organism_name"] == "Escherichia coli"
    assert rec["ftp_path"].endswith("GCF_000005845.2_ASM584v2")
    # sentinel columns that should NOT be present
    assert "bioproject" not in rec
    assert "genome_size" not in rec


def test_build_slim_name_based_selection_survives_reordered_columns(tmp_path):
    """If NCBI reorders columns, name-based selection still picks the right
    values (positional selection would silently grab the wrong fields)."""
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    _write_summary(gb, ["GCA_000005845.2"], reorder=True)
    _write_summary(rs, [])

    build_slim_assembly_summary(gb, rs, out)
    header, row = out.read_text().splitlines()[:2]
    rec = dict(zip(header.split("\t"), row.split("\t")))
    assert rec["taxid"] == "562"
    assert rec["organism_name"] == "Escherichia coli"
    assert rec["ftp_path"].endswith("GCF_000005845.2_ASM584v2")


def test_build_slim_handles_gzipped_inputs(tmp_path):
    gb = tmp_path / "gb.txt.gz"
    rs = tmp_path / "rs.txt.gz"
    out = tmp_path / "slim.tsv"
    _write_summary(gb, ["GCA_000005845.2"], gzip_it=True)
    _write_summary(rs, ["GCF_000005845.2"], gzip_it=True)

    n = build_slim_assembly_summary(gb, rs, out)
    assert n == 2


def test_build_slim_missing_column_raises(tmp_path):
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    bad_cols = [c for c in NCBI_COLUMNS if c != "ftp_path"]
    gb.write_text("##x\n#" + "\t".join(bad_cols) + "\n")
    _write_summary(rs, [])
    with pytest.raises(ValueError, match="ftp_path"):
        build_slim_assembly_summary(gb, rs, out)


def test_build_slim_skips_short_rows(tmp_path):
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    good = _row("GCA_000005845.2")
    short = "GCA_999999999.1\t562\ttruncated"
    gb.write_text("##x\n#" + "\t".join(NCBI_COLUMNS) + "\n" + good + "\n" + short + "\n")
    _write_summary(rs, [])
    n = build_slim_assembly_summary(gb, rs, out)
    assert n == 1   # the short/malformed row is dropped


def test_build_slim_atomic_no_tmp_left_behind(tmp_path):
    gb = tmp_path / "gb.txt"
    rs = tmp_path / "rs.txt"
    out = tmp_path / "slim.tsv"
    _write_summary(gb, ["GCA_000005845.2"])
    _write_summary(rs, [])
    build_slim_assembly_summary(gb, rs, out)
    assert not (tmp_path / "slim.tsv.tmp").exists()


################################################################################
# date-retrieved stamping
################################################################################

def test_write_date_retrieved_format(tmp_path):
    from datetime import date
    p = tmp_path / "date-retrieved.txt"
    write_date_retrieved(p, when=date(2026, 4, 1))
    assert p.read_text() == "2026,04,01\n"


def test_write_date_retrieved_defaults_to_today(tmp_path):
    from datetime import date
    p = tmp_path / "date-retrieved.txt"
    write_date_retrieved(p)
    today = date.today().strftime("%Y,%m,%d")
    assert p.read_text().strip() == today
