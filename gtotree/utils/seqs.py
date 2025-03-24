
from Bio import SeqIO
import os
from gtotree.utils.general import remove_file_if_exists

def filter_and_rename_fasta(prefix, run_data, in_path, max_length = 99999):

    num = 0
    infile = f"{in_path}/{prefix}_protein.faa"
    outfile = f"{run_data.all_input_genome_AA_files_dir}/{prefix}.faa"
    with open(outfile, "w") as outfile, open(infile, "r") as infile:

        for record in SeqIO.parse(infile, "fasta"):

            if len(record.seq) <= max_length:
                num += 1
                header = f">{prefix}_{num}"
                outfile.write(f"{header}\n{record.seq}\n")


def extract_filter_and_rename_cds_amino_acids_from_gb(prefix, input_gb, run_data, max_length = 99999):

    num = 0
    output_file = f"{run_data.all_input_genome_AA_files_dir}/{prefix}.faa"
    note_terms_to_exclude = ["frameshifted", "internal stop", "incomplete"]
    location_terms_to_exclude = ["join", "<", ">"]

    try:
        with open(input_gb, "r") as infile, open(output_file, "w") as outfile:
            records = list(SeqIO.parse(infile, "genbank"))
            for rec in records:
                genes = [gene for gene in rec.features if gene.type == "CDS"]
                for gene in genes:
                    location = str(gene.location)
                    if any(term in location for term in location_terms_to_exclude):
                        continue
                    if "note" in gene.qualifiers:
                        note = gene.qualifiers["note"][0]
                        if any(term in note for term in note_terms_to_exclude):
                            continue
                    if "transl_except" in gene.qualifiers or "pseudo" in gene.qualifiers:
                        continue

                    translation = gene.qualifiers.get("translation")
                    if translation:
                        if len(translation[0]) <= max_length:
                            num += 1
                            header = f">{prefix}_{num}"
                            outfile.write(f"{header}\n{translation[0]}\n")
        if num == 0:
            remove_file_if_exists(output_file)
            return False
        else:
            return True
    except:
        remove_file_if_exists(output_file)
        return False


def extract_fasta_from_gb(prefix, input_gb, run_data):

    num = 0
    output_file = f"{run_data.genbank_processing_dir}/{prefix}.fasta"

    with open(input_gb, "r") as infile, open(output_file, "w") as outfile:
        records = list(SeqIO.parse(infile, "genbank"))
        for rec in records:
            num += 1
            outfile.write(f">{prefix}_{num}\n{rec.seq}\n")
