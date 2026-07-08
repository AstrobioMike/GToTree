import itertools
import re
from gtotree.utils.messaging import report_ncbi_accs_not_found


# logical field -> column name in the (slim or full) NCBI assembly summary
NAMES = {
    "accession":      "assembly_accession",
    "assembly_name":  "asm_name",
    "taxid":          "taxid",
    "org_name":       "organism_name",
    "infra_name":     "infraspecific_name",
    "version_status": "version_status",
    "assembly_level": "assembly_level",
    "ftp_path":       "ftp_path",
}

# fixed positions in NCBI's full assembly_summary_*.txt (headerless fallback)
LEGACY_POS = {
    "accession": 0, "assembly_name": 15, "taxid": 5, "org_name": 7,
    "infra_name": 8, "version_status": 10, "assembly_level": 11,
    "ftp_path": 19,
}


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


def resolve_positions(header_line):
    """header line (clean or '#'-prefixed) -> {logical: index}."""
    names = header_line.lstrip("#").rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(names)}
    return {logical: idx[col] for logical, col in NAMES.items()}


def _clean(value):
    value = value.strip()
    return value if value else "NA"


def parse_assembly_summary(assembly_summary_file, run_data):
    """
    parse NCBI's assembly summary file down to the provided wanted accessions.

    Writes a sub-table (ncbi-accessions-info.tsv) into run_data.tmp_dir with one
    row per found wanted accession and these tab-separated columns:
        input_accession, found_accession, assembly_name, taxid, organism_name,
        infraspecific_name, version_status, assembly_level, http_base_link

    Accessions not found are recorded (ncbi-accessions-not-found.txt) and their
    GenomeData entries are marked removed. Returns the updated run_data.
    """

    if run_data.ncbi_sub_table_path:
        return run_data

    wanted_dict = {}
    for acc in run_data.get_input_ncbi_accs():
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    ncbi_sub_table_path = run_data.tmp_dir + "/ncbi-accessions-info.tsv"

    with open(ncbi_sub_table_path, "w") as out_file, \
            open(assembly_summary_file, "r") as assemblies:

        out_file.write("input_accession\tfound_accession\tassembly_name\ttaxid\t"
                       "organism_name\tinfraspecific_name\tversion_status\t"
                       "assembly_level\thttp_base_link\n")

        # peek at the first line to decide header-based vs. legacy-positional
        first = assemblies.readline()
        first_stripped = first.lstrip("#").rstrip("\n")
        if first_stripped.split("\t", 1)[0] == "assembly_accession":
            pos = resolve_positions(first)          # header present
            remaining = assemblies                  # data starts on next line
        else:
            pos = dict(LEGACY_POS)                   # headerless -> legacy pos
            remaining = itertools.chain([first], assemblies)  # first is data

        p_acc = pos["accession"]; p_name = pos["assembly_name"]
        p_tax = pos["taxid"]; p_org = pos["org_name"]; p_infra = pos["infra_name"]
        p_ver = pos["version_status"]; p_lvl = pos["assembly_level"]
        p_ftp = pos["ftp_path"]
        max_pos = max(pos.values())

        for line in remaining:
            fields = line.rstrip("\n").split("\t")
            if not fields or len(fields) <= max_pos:
                continue

            root = fields[p_acc].split(".")[0]
            if root not in wanted_dict:
                continue

            found.add(wanted_dict[root])

            dl_acc = _clean(fields[p_acc])
            assembly_name = _clean(fields[p_name])
            taxid = _clean(fields[p_tax])
            org_name = _clean(fields[p_org])
            infra_name = _clean(fields[p_infra])
            version_status = _clean(fields[p_ver])
            assembly_level = _clean(fields[p_lvl])

            # prefer the real ftp_path (normalized to https), else rebuild it,
            # else "na" (this is the base_link consumed downstream).
            http_path = resolve_base_link(fields[p_ftp], dl_acc, assembly_name)

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

    if len(not_found) > 0:
        with open(run_data.run_files_dir + "/ncbi-accessions-not-found.txt", "w") as not_found_file:
            for acc in not_found:
                not_found_file.write(acc + "\n")

        report_ncbi_accs_not_found(len(not_found), run_data.run_files_dir_rel)

    for acc_gd in run_data.ncbi_accs:
        if acc_gd.id in not_found:
            acc_gd.acc_was_found = False
            acc_gd.mark_removed("accession not found at NCBI")
        else:
            acc_gd.acc_was_found = True

    run_data.ncbi_sub_table_path = ncbi_sub_table_path

    return run_data
