import re
import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.dataset as ds # type: ignore
from gtotree.utils.misc.general import REASON_NOT_FOUND_AT_NCBI
from gtotree.utils.misc.stages import GenomeRemovalStage
from gtotree.utils.misc.summary_info import write_removed_genomes_report


# columns read from the NCBI Parquet asset for the per-accession info sub-table
_NEEDED_COLUMNS = [
    "assembly_accession", "asm_name", "taxid", "organism_name",
    "infraspecific_name", "version_status", "assembly_level", "ftp_path",
]

_ACCESSION_ROOT_REGEX = r'^([^.]+)\..*$'


def sanitize_assembly_name(name):
    """mirror NCBI's directory-name sanitizing (ported from bit)."""
    sanitized = re.sub(r"[\s/,#()\[\]]", "_", name)
    sanitized = re.sub(r"_+", "_", sanitized)
    return sanitized.strip("_")


def build_base_link(dl_acc, assembly_name):
    """
    fallback directory URL builder, used only when the summary row has no usable
    ftp_path (ported from bit's parse_ncbi_assembly_summary.build_base_link).
    Returns the directory URL (with trailing slash) and its basename.
    """
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    prefix, rest = dl_acc.split("_", 1)
    number = rest.split(".")[0]
    p1, p2, p3 = number[0:3], number[3:6], number[6:9]
    dir_basename = f"{dl_acc}_{sanitize_assembly_name(assembly_name)}"
    path = f"{prefix}/{p1}/{p2}/{p3}/{dir_basename}/"
    return base_url + path, dir_basename


def resolve_base_link(ftp_path, dl_acc, assembly_name):
    """
    resolve the download directory URL for a row: prefer the summary's real
    ftp_path (normalized to https, no trailing slash); otherwise rebuild it from
    the accession + assembly name; otherwise "na". Returns the directory URL
    (no trailing slash, or "na").
    """
    ftp_path = (ftp_path or "").strip()
    if ftp_path and ftp_path.lower() != "na":
        return ftp_path.replace("ftp://", "https://").rstrip("/")
    if dl_acc and dl_acc.lower() != "na" and assembly_name and assembly_name.lower() != "na":
        built, _ = build_base_link(dl_acc, assembly_name)
        return built.rstrip("/")
    return "na"


def _clean(value):
    if value is None:
        return "NA"
    value = str(value).strip()
    return value if value else "NA"


def _matched_accession_filter(dataset, wanted_dict):
    """
    Build the scanner filter that selects the wanted accessions, or None if there
    are none
    """
    if not wanted_dict:
        return None

    acc_col = dataset.to_table(columns=["assembly_accession"]).column(0)
    roots = pc.replace_substring_regex(acc_col, _ACCESSION_ROOT_REGEX, r"\1")
    mask = pc.is_in(roots, value_set=pa.array(sorted(wanted_dict)))

    return ds.field("assembly_accession").isin(pc.filter(acc_col, mask))


def parse_assembly_summary(assembly_summary_file, run_data):
    """
    Look up the wanted NCBI accessions in the hosted NCBI Parquet table and write the
    per-accession info sub-table.

    Writes a sub-table (ncbi-accessions-info.tsv) into run_data.tmp_dir with one row
    per found wanted accession and these tab-separated columns:
        input_accession, found_accession, assembly_name, taxid, organism_name,
        infraspecific_name, version_status, assembly_level, http_base_link

    (input_accession and http_base_link are the two columns processing's
    get_base_link depends on -- their names/semantics are held stable.)

    Accessions not found have their GenomeData entries marked removed, which is what
    puts them into removed-genomes.tsv. Returns the updated run_data.
    """

    if run_data.ncbi_sub_table_path:
        return run_data

    # root (version-stripped) accession -> the exact string the user asked for
    wanted_dict = {}
    for acc in run_data.get_input_ncbi_accs():
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    gd_by_acc = {gd.id: gd for gd in run_data.ncbi_accs}

    ncbi_sub_table_path = run_data.tmp_dir + "/ncbi-accessions-info.tsv"

    dataset = ds.dataset(str(assembly_summary_file), format="parquet")
    predicate = _matched_accession_filter(dataset, wanted_dict)

    with open(ncbi_sub_table_path, "w") as out_file:

        out_file.write("input_accession\tfound_accession\tassembly_name\ttaxid\t"
                       "organism_name\tinfraspecific_name\tversion_status\t"
                       "assembly_level\thttp_base_link\n")

        # predicate=None means no wanted accessions; the scanner handles that fine
        scanner = dataset.scanner(columns=_NEEDED_COLUMNS, filter=predicate)

        for batch in scanner.to_batches():
            for row in batch.to_pylist():
                dl_acc_raw = _clean(row.get("assembly_accession"))
                root = dl_acc_raw.split(".")[0]
                if root not in wanted_dict:
                    continue

                found.add(wanted_dict[root])

                dl_acc = dl_acc_raw
                assembly_name = _clean(row.get("asm_name"))
                taxid = _clean(row.get("taxid"))
                org_name = _clean(row.get("organism_name"))
                infra_name = _clean(row.get("infraspecific_name"))
                version_status = _clean(row.get("version_status"))
                assembly_level = _clean(row.get("assembly_level"))

                matched_gd = gd_by_acc.get(wanted_dict[root])
                if matched_gd is not None:
                    matched_gd.taxid = taxid

                # prefer the real ftp_path (normalized to https), else rebuild it,
                # else "na" (this is the base_link consumed downstream).
                http_path = resolve_base_link(row.get("ftp_path"), dl_acc, assembly_name)

                out_file.write("\t".join([
                    wanted_dict[root],
                    dl_acc,
                    assembly_name,
                    taxid,
                    org_name,
                    infra_name,
                    version_status,
                    assembly_level,
                    http_path,
                ]) + "\n")

    not_found = set(run_data.get_input_ncbi_accs()) - found

    for acc_gd in run_data.ncbi_accs:
        if acc_gd.id in not_found:
            acc_gd.acc_was_found = False
            acc_gd.mark_removed(REASON_NOT_FOUND_AT_NCBI,
                                GenomeRemovalStage.NCBI_LOOKUP)
        else:
            acc_gd.acc_was_found = True

    if len(not_found) > 0:
        write_removed_genomes_report(run_data)

    run_data.ncbi_sub_table_path = ncbi_sub_table_path

    return run_data
