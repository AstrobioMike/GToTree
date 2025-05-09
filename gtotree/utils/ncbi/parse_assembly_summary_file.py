from gtotree.utils.messaging import report_ncbi_accs_not_found

def parse_assembly_summary(assembly_summary_file, run_data):
    """
    parse NCBI's assembly summary file down to the provided wanted accessions

    Parameters:
        assembly_summary_file (str): Path to NCBI's assembly summary file.
        wanted_accessions (list of str): A list of wanted accessions (as strings).
                                         The function will use the part before the first dot as the key.
        output_file (str): Path to write the filtered summary. Default is "Wanted.tsv".

    Returns:
        None. (The output is written to output_file.)

    The function works as follows:
        1. Builds a dictionary mapping from each wanted accession's "root" (i.e. before the first '.')
            to the full accession string.
        2. Reads through the assembly summary file, splitting each line on tab.
        3. If the first column's root accession is in the wanted dictionary,
            it extracts selected fields and writes them out (using "NA" for missing values).

    Expected fields written (tab-separated):
        - The original wanted accession (from the wanted list)
        - Downloaded accession (field 0)
        - Assembly name (field 15)
        - Taxid (field 5)
        - Organism name (field 7)
        - Infraspecific name (field 8)
        - Version status (field 10)
        - Assembly level (field 11)
        - FTP path (field 19)
        - HTTP path (built from the downloaded accession and assembly name; in case ftp path is empty, which does rarely happen)
    """

    if run_data.ncbi_sub_table_path:
        return run_data

    wanted_dict = {}
    for acc in run_data.get_input_ncbi_accs():
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    ncbi_sub_table_path = run_data.tmp_dir + "/ncbi-accessions-info.tsv"

    with open(ncbi_sub_table_path, "w") as out_file:
        out_file.write("input_accession\tfound_accession\tassembly_name\ttaxid\torganism_name\tinfraspecific_name\tversion_status\tassembly_level\thttp_base_link\n")
        with open(assembly_summary_file, "r") as assemblies:
            for line in assemblies:
                fields = line.strip().split("\t")
                if not fields:
                    continue
                root = fields[0].split(".")[0]
                if root in wanted_dict:
                    found.add(wanted_dict[root])

                    dl_acc = fields[0].strip() if fields[0].strip() else "NA"
                    assembly_name = fields[15].strip() if len(fields) > 15 and fields[15].strip() else "NA"
                    taxid = fields[5].strip() if len(fields) > 5 and fields[5].strip() else "NA"
                    org_name = fields[7].strip() if len(fields) > 7 and fields[7].strip() else "NA"
                    infra_name = fields[8].strip() if len(fields) > 8 and fields[8].strip() else "NA"
                    version_status = fields[10].strip() if len(fields) > 10 and fields[10].strip() else "NA"
                    assembly_level = fields[11].strip() if len(fields) > 11 and fields[11].strip() else "NA"

                    http_path = build_base_link(dl_acc, assembly_name)

                    out_line = "\t".join([
                        wanted_dict[root],
                        dl_acc,
                        assembly_name,
                        taxid,
                        org_name,
                        infra_name,
                        version_status,
                        assembly_level,
                        http_path
                    ]) + "\n"
                    out_file.write(out_line)

    not_found = set(run_data.get_input_ncbi_accs()) - found

    if len(not_found) > 0:
        with open(run_data.run_files_dir + "/ncbi-accessions-not-found.txt", "w") as not_found_file:
            for acc in not_found:
                not_found_file.write(acc + "\n")

        report_ncbi_accs_not_found(len(not_found), run_data.run_files_dir_rel)

    for acc_gd in run_data.ncbi_accs:
        if acc_gd.id in not_found:
            acc_gd.was_found = False
            acc_gd.removed = True
        else:
            acc_gd.was_found = True

    run_data.ncbi_sub_table_path = ncbi_sub_table_path

    return run_data


def build_base_link(dl_acc, assembly_name):
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    prefix, number = dl_acc.split("_")

    p1 = number[0:3]
    p2 = number[3:6]
    p3 = number[6:9]

    path = f"{prefix}/{p1}/{p2}/{p3}/{dl_acc}_{assembly_name}/"

    return base_url + path
