"""
Combine the NCBI GenBank + RefSeq assembly_summary tables into a single slim
TSV holding only the columns GToTree uses, with a clean (no leading '#') header
row.

This is the same projection used by bit's scheduled GitHub Action that rebuilds
and re-hosts the slim table (the hosted asset GToTree pulls by default), and by
the local rebuild fallback here.

NCBI's assembly_summary_*.txt format:
  - line 1 is a '##' provenance comment
  - line 2 is the real header, beginning with '#assembly_accession\t...'
  - remaining lines are tab-delimited data rows in a fixed 38-column layout

KEPT_COLUMNS is the single source of truth for the slim layout and column order.
Ported from bit (bit.modules.ncbi.slim_ncbi_assembly_summary).
"""

import gzip
import os

# the columns GToTree actually reads from the assembly summary, in the order
# they are written to the slim file. Keep in sync with the reader
# (parse_ncbi_assembly_summary.NAMES), which looks these up by name.
KEPT_COLUMNS = [
    "assembly_accession",
    "taxid",
    "organism_name",
    "infraspecific_name",
    "version_status",
    "assembly_level",
    "asm_name",
    "ftp_path",
]

# NCBI's real header line starts with this token (with the leading '#')
_NCBI_HEADER_LEAD = "#assembly_accession"


def _open_maybe_gzip(path, mode="rt"):
    """open a path as text, transparently handling .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_ncbi_header(line):
    """
    turn an NCBI assembly_summary header line into a list of clean column names.
    Strips a single leading '#' from the first field if present. Raises
    ValueError if the line doesn't look like the expected header.
    """
    stripped = line.rstrip("\n")
    if not stripped.startswith("#"):
        raise ValueError("expected an NCBI assembly_summary header line "
                         "(beginning with '#assembly_accession')")
    # drop only the leading '#', keep the rest of that first field
    stripped = stripped[1:]
    names = stripped.split("\t")
    if not names or names[0] != "assembly_accession":
        raise ValueError(
            "unexpected NCBI assembly_summary header; first column is "
            f"'{names[0] if names else ''}', expected 'assembly_accession'")
    return names


def _read_summary(path):
    """
    return (header_names, data_lines) for one NCBI assembly_summary file.
    Skips the '##' provenance line, parses the '#...' header, and collects the
    remaining raw data lines (newline-stripped, non-empty).
    """
    fh = _open_maybe_gzip(path, "rt")
    header = None
    data_lines = []
    for line in fh:
        if header is None:
            if line.startswith("##"):
                continue                    # provenance comment -> skip
            if line.startswith("#"):
                header = parse_ncbi_header(line)
                continue
            # a file with no header at all -> not the expected format
            raise ValueError(f"no header found in {path} before data rows")
        s = line.rstrip("\n")
        if s:
            data_lines.append(s)
    fh.close()
    if header is None:
        raise ValueError(f"no header found in {path}")
    return header, data_lines


def _project_rows(data_lines, keep_idx):
    """
    return slimmed rows (lists of values) selecting keep_idx positions from each
    data line. Rows with too few fields to cover the needed indices are skipped
    (defensive against malformed lines).
    """
    max_needed = max(keep_idx)
    out = []
    for s in data_lines:
        fields = s.split("\t")
        if len(fields) <= max_needed:
            continue
        out.append([fields[i] for i in keep_idx])
    return out


def build_slim_assembly_summary(genbank_path, refseq_path, out_path,
                                keep_columns=None):
    """
    combine the GenBank and RefSeq assembly_summary files into a single slim TSV
    at out_path, keeping only keep_columns (default KEPT_COLUMNS) BY NAME, with a
    clean header row. Returns the number of data rows written.

    Both inputs share the kept column names (GenBank and RefSeq use the same
    assembly_summary schema). Columns are located by name in each file's own
    header, so the two files are allowed to differ in column order.
    """
    keep_columns = list(keep_columns) if keep_columns else list(KEPT_COLUMNS)

    gb_header, gb_data = _read_summary(genbank_path)
    rs_header, rs_data = _read_summary(refseq_path)

    def _indices(header, source):
        missing = [c for c in keep_columns if c not in header]
        if missing:
            raise ValueError(
                f"{source} assembly_summary is missing expected column(s): "
                f"{', '.join(missing)}")
        return [header.index(c) for c in keep_columns]

    gb_idx = _indices(gb_header, "GenBank")
    rs_idx = _indices(rs_header, "RefSeq")

    rows = _project_rows(gb_data, gb_idx)
    rows += _project_rows(rs_data, rs_idx)

    tmp_path = str(out_path) + ".tmp"
    with open(tmp_path, "w") as out:
        out.write("\t".join(keep_columns) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    os.replace(tmp_path, out_path)        # atomic swap into place

    return len(rows)


GENBANK_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

# names of the two files packaged inside the distributed tarball (kept flat at
# the archive root, mirroring the GTDB slim bundle).
TABLE_FILENAME = "ncbi-assembly-info.tsv"
DATE_FILENAME = "date-retrieved.txt"


def write_date_retrieved(path, when=None):
    """
    write a date-retrieved.txt in the stored format (YYYY,MM,DD). Defaults to
    today (the build date when called from the rebuild fallback).
    """
    from datetime import date as _date
    when = when or _date.today()
    with open(path, "w") as fh:
        fh.write(when.strftime("%Y,%m,%d") + "\n")
