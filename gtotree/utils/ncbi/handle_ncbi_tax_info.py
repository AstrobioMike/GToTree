import os
import pandas as pd
import subprocess

def update_mapping_dict_with_ncbi_tax_info(args, run_data):

    sub_ncbi_tab = pd.read_csv(os.path.join(run_data.tmp_dir, "ncbi-accessions-info.tsv"), sep="\t", header=0)

    taxid_path = os.path.join(run_data.tmp_dir, "ncbi-accession-taxids.txt")
    sub_ncbi_tab['taxid'].to_csv(taxid_path, index=False, header=False)

    sub_ncbi_tax_path = os.path.join(run_data.tmp_dir, "ncbi-accession-tax.tsv")

    success = run_taxonkit(taxid_path, sub_ncbi_tax_path)
    if not success:
        return run_data

    sub_ncbi_tax_tab = pd.read_csv(sub_ncbi_tax_path, sep="\t", header=None)
    sub_ncbi_tax_tab.columns = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

    sub_ncbi_tax_tab = reformat_ncbi_tax_tab(sub_ncbi_tax_tab)

    combined_tab = pd.concat([sub_ncbi_tab["input_accession"], sub_ncbi_tax_tab], axis=1)
    combined_tab = combined_tab.rename(columns={"input_accession": "input_acc"})
    combined_tab.to_csv(sub_ncbi_tax_path, sep="\t", index=False)

    for index, row in combined_tab.iterrows():
        curr_input_acc = row["input_acc"]
        if curr_input_acc in run_data.mapping_dict:
            continue
        elif row['domain'] == "NA":
            continue
        else:
            new_label = f"{curr_input_acc}_"
            for rank in args.lineage.split(","):
                new_label += row[rank.lower()] + "_"

            new_label = new_label[:-1]  # removing last underscore
            new_label = new_label.replace(" ", "_")
            run_data.mapping_dict[curr_input_acc] = new_label

    return run_data


def run_taxonkit(taxid_path, sub_ncbi_tax_path):

    taxonkit_reformat_pattern = '{domain|superkingdom}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}\t{strain|subspecies}'
    cmd = f'taxonkit lineage -L -r -c {taxid_path} | taxonkit reformat2 -r NA -f "{taxonkit_reformat_pattern}" | cut -f 4- > {sub_ncbi_tax_path}'
    success = subprocess.run(cmd, shell=True, check=True)

    return success


def reformat_ncbi_tax_tab(df):

    # setting all to NAs if there was a problem with taxonkit
    mask = df['domain'].str.contains(r'\{domain\|superkingdom\}', na=False)
    # set *every* column in those rows to the string "NA"
    df.loc[mask, :] = 'NA'

    # this is to remove redundancy in species/strain
    orig_species = df['species'].tolist()

    df['species'] = [
        sp.removeprefix(genus + ' ')
        for genus, sp in zip(df['genus'], df['species'])
    ]

    df['strain'] = [
        st.removeprefix(full_sp + ' ').replace(' ', '-')
        for full_sp, st in zip(orig_species, df['strain'])
    ]

    return df
