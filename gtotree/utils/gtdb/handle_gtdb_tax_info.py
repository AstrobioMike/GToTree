import os
import pandas as pd

def update_mapping_dict_with_gtdb_tax_info(args, run_data):

    sub_gtdb_tax_tab = subset_gtdb_info(run_data)
    for index, row in sub_gtdb_tax_tab.iterrows():
        curr_input_acc = row["input_acc"]
        if curr_input_acc in run_data.mapping_dict:
            continue
        else:
            new_label = ""
            for rank in args.lineage.split(","):
                new_label += row[rank.lower()] + "_"

            new_label = new_label[:-1]  # removing last underscore
            new_label = new_label.replace(" ", "_")
            run_data.mapping_dict[curr_input_acc] = new_label

    return(run_data)


def subset_gtdb_info(run_data):

    # making dict of input accessions we're searching for
    input_accs = run_data.remaining_ncbi_accs()
    wanted_dict = {}

    for acc in input_accs:
        root_acc = acc.split(".")[0][4:]
        wanted_dict[str(root_acc)] = acc


    master_gtdb_tab = os.path.join(os.environ['GTDB_dir'], "GTDB-arc-and-bac-metadata.tsv")
    sub_gtdb_tab_path = os.path.join(run_data.tmp_dir, "gtdb-subset.tsv")

    if not os.path.exists(sub_gtdb_tab_path):

        # iterate through master gtdb table and writing sub-table
        with open(master_gtdb_tab, "r") as in_file, open(sub_gtdb_tab_path, "w") as out_file:

            # writing out header
            out_file.write(in_file.readline())

            for line in in_file:
                split_line = line.strip().split("\t")

                # checking if the accession is in our wanted_dict
                acc_with_no_version = split_line[0][7:].split(".")[0]

                if acc_with_no_version in wanted_dict:
                    out_file.write(line)


    sub_gtdb_tax_path = os.path.join(run_data.tmp_dir, "gtdb-subset-tax.tsv")

    if not os.path.exists(sub_gtdb_tax_path):

        with open(sub_gtdb_tab_path, "r") as in_file, open(sub_gtdb_tax_path, "w") as out_file:
            # writing out header
            out_file.write("input_acc\tbase_gtdb_acc\tfull_gtdb_acc\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")

            # skipping header
            next(in_file)

            # iterating through the gtdb table
            for line in in_file:
                split_line = line.strip().split("\t")
                full_gtdb_acc = split_line[0]
                base_gtdb_acc = full_gtdb_acc[3:]

                input_acc = get_input_acc(base_gtdb_acc, input_accs)

                gtdb_tax_list = [split_line[1], split_line[2], split_line[3], split_line[4], split_line[5], split_line[6], split_line[7]]

                if len(gtdb_tax_list) != 7:
                    print("GTDB entry " + full_gtdb_acc + " doesn't seem to have full lineage info.")

                out_file.write(f"{input_acc}\t{base_gtdb_acc}\t{full_gtdb_acc}\t{gtdb_tax_list[0]}\t{gtdb_tax_list[1]}\t{gtdb_tax_list[2]}\t{gtdb_tax_list[3]}\t{gtdb_tax_list[4]}\t{gtdb_tax_list[5]}\t{gtdb_tax_list[6]}\n")

    sub_gtdb_tax_tab = pd.read_csv(sub_gtdb_tax_path, sep="\t", header=0)
    return(sub_gtdb_tax_tab)


def get_input_acc(base_gtdb_acc, input_accs):
    # getting the input accession that goes with a gtdb entry
    if base_gtdb_acc in input_accs:
        return base_gtdb_acc

    base_acc_no_version = base_gtdb_acc.split(".")[0]
    if base_acc_no_version in input_accs:
        return base_acc_no_version

    base_acc_no_prefix = base_acc_no_version.split("_")[1]

    # dealing with situations where we pulled the tax for a GTDB GCF_ but the input was a GCA_
    for acc in input_accs:
        input_acc_no_version_or_prefix = acc.split(".")[0].split("_")[1]
        if base_acc_no_prefix == input_acc_no_version_or_prefix:
            return acc

    return None
