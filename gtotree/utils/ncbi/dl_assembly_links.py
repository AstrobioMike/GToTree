"""
Accession -> (download link, local destination) lookup for `gtt dl-ncbi-assemblies`.

Separate from parse_ncbi_assembly_summary.parse_assembly_summary(), which does the
same lookup but for a main GToTree run: it is coupled to the run-state model (GenomeData
entries, tmp_dir, removed-genomes book-keeping) and emits only a base link, since the
main driver fetches one fixed format. The standalone subcommand has no run state,
writes into the user's output dir, and needs a per-format target link.

The link-building helpers themselves are imported from that module rather than
duplicated, so a fix to the NCBI path layout lands in both.
"""

import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.dataset as ds # type: ignore

from gtotree.utils.ncbi.parse_ncbi_assembly_summary import (build_base_link,
                                                            resolve_base_link,
                                                            _clean)


# NCBI's suffix for each format, and the extension we save it under locally
FORMAT_EXTENSIONS = {
    "genbank":     ("_genomic.gbff.gz",         ".gb.gz"),
    "fasta":       ("_genomic.fna.gz",          ".fasta.gz"),
    "protein":     ("_protein.faa.gz",          ".faa.gz"),
    "gff":         ("_genomic.gff.gz",          ".gff.gz"),
    "nt_cds":      ("_cds_from_genomic.fna.gz", "_cds_from_genomic.fna.gz"),
    "feature_tab": ("_feature_table.txt.gz",    ".tsv.gz"),
    "report":      ("_assembly_report.txt",     "_assembly_report.txt"),
    "stats":       ("_assembly_stats.txt",      "_assembly_stats.txt"),
}


_NEEDED_COLUMNS = [
    "assembly_accession", "asm_name", "taxid", "organism_name",
    "infraspecific_name", "version_status", "assembly_level", "ftp_path",
]


def _resolve_links(dl_acc, assembly_name, ftp_path):
    """
    Prefer the table's ftp_path (rewritten to https), else rebuild the path from the
    accession + assembly name, else NA.

    Returns (http_base_link, dir_basename); dir_basename is the leaf directory name,
    which is also the filename stem NCBI uses for every file in that directory.
    """
    if ftp_path and ftp_path.lower() != "na":
        http_path = resolve_base_link(ftp_path, dl_acc, assembly_name)
    elif assembly_name != "NA" and dl_acc != "NA":
        http_path, _ = build_base_link(dl_acc, assembly_name)
    else:
        return "NA", "NA"

    http_path = http_path.rstrip("/")
    dir_basename = http_path.split("/")[-1]
    return http_path, dir_basename


def parse_ncbi_assembly_summary(assembly_summary_file, run_data):
    """
    Look up run_data.wanted_accs in the hosted NCBI Parquet table, write the
    per-accession info table (run_data.ncbi_sub_table_path) with target links and
    local destinations, and write the not-found list.

    Sets run_data.num_found / num_not_found and returns run_data.
    """
    # root (version-stripped) accession -> the exact string the user asked for
    wanted_dict = {}
    for acc in run_data.wanted_accs:
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    dataset = ds.dataset(str(assembly_summary_file), format="parquet")

    root_field = pc.replace_substring_regex(ds.field("assembly_accession"),
                                            r"\..*$", "")
    predicate = pc.is_in(
        root_field,
        value_set=pa.array(sorted(set(wanted_dict)), type=pa.string()))

    with open(run_data.ncbi_sub_table_path, "w") as out_file:

        cols = ["input_accession", "found_accession", "assembly_name", "taxid",
                "organism_name", "infraspecific_name", "version_status",
                "assembly_level", "http_base_link"]
        if run_data.wanted_format:
            cols.extend(["target_link", "local_destination"])
        out_file.write("\t".join(cols) + "\n")

        scanner = dataset.scanner(columns=_NEEDED_COLUMNS, filter=predicate)

        for batch in scanner.to_batches():
            for row in batch.to_pylist():
                dl_acc = _clean(row.get("assembly_accession"))
                root = dl_acc.split(".")[0]
                if root not in wanted_dict:
                    continue

                found.add(wanted_dict[root])

                assembly_name = _clean(row.get("asm_name"))
                ftp_path = (row.get("ftp_path") or "").strip()
                http_path, dir_basename = _resolve_links(dl_acc, assembly_name,
                                                         ftp_path)

                out_fields = [
                    wanted_dict[root], dl_acc, assembly_name,
                    _clean(row.get("taxid")), _clean(row.get("organism_name")),
                    _clean(row.get("infraspecific_name")),
                    _clean(row.get("version_status")),
                    _clean(row.get("assembly_level")), http_path,
                ]
                out_line = "\t".join(out_fields)

                if run_data.wanted_format:
                    ncbi_ext, local_ext = FORMAT_EXTENSIONS[run_data.wanted_format]
                    target_link = (f"{http_path}/{dir_basename}{ncbi_ext}"
                                   if http_path.lower() != "na" else "NA")
                    local_path = f"{run_data.output_dir}/{dl_acc}{local_ext}"
                    out_line += "\t" + target_link + "\t" + local_path

                out_file.write(out_line + "\n")

    not_found = set(run_data.wanted_accs) - found

    if not_found:
        with open(run_data.not_found_path, "w") as not_found_file:
            for acc in sorted(not_found):
                not_found_file.write(acc + "\n")

    run_data.num_found = len(found)
    run_data.num_not_found = len(not_found)

    return run_data
