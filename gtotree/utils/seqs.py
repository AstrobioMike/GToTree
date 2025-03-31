import statistics
from Bio import SeqIO
from gtotree.utils.general import (remove_file_if_exists,
                                   check_file_exists_and_not_empty)
from gtotree.utils.messaging import (add_border,
                                     report_notice)

def filter_and_rename_fasta(prefix, run_data, in_path, full_path = False, max_length = 99999):

    num = 0
    if full_path:
        infile = in_path
    else:
        infile = f"{in_path}/{prefix}_protein.faa"
    outpath = f"{run_data.ready_genome_AA_files_dir}/{prefix}.faa"
    with open(outpath, "w") as outfile, open(infile, "r") as infile:

        for record in SeqIO.parse(infile, "fasta"):

            if len(record.seq) <= max_length:
                num += 1
                header = f">{prefix}_{num}"
                outfile.write(f"{header}\n{record.seq}\n")

    if num == 0:
        remove_file_if_exists(outpath)
        return False, None

    else:
        return True, outpath

def extract_filter_and_rename_cds_amino_acids_from_gb(prefix, input_gb, run_data, max_length = 99999):

    num = 0
    output_file = f"{run_data.ready_genome_AA_files_dir}/{prefix}.faa"
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
            return False, None
        else:
            return True, output_file
    except:
        remove_file_if_exists(output_file)
        return False, None


def extract_fasta_from_gb(prefix, input_gb, run_data):

    num = 0
    output_file = f"{run_data.genbank_processing_dir}/{prefix}.fasta"

    with open(input_gb, "r") as infile, open(output_file, "w") as outfile:
        records = list(SeqIO.parse(infile, "genbank"))
        for rec in records:
            num += 1
            outfile.write(f">{prefix}_{num}\n{rec.seq}\n")


def check_target_SCGs_have_seqs(run_data, ext):

    SCG_targets_present_at_start = run_data.get_all_SCG_targets_remaining()
    new_SCG_targets_missing = []

    for SCG_obj in SCG_targets_present_at_start:
        SCG = SCG_obj.id
        path = run_data.found_SCG_seqs_dir + f"/{SCG}{ext}"
        present = check_file_exists_and_not_empty(path)
        if not present:
            SCG_obj.removed = True
            SCG_obj.reason_removed = "no sequences found or remaining after length-filtering"
            new_SCG_targets_missing.append(SCG)

    if len(new_SCG_targets_missing) > 0:
        add_border()
        message = f"    Some target single-copy genes were not found or were filtered out of the\n"
        message += f"    analysis.\n\n    At this point, these include:\n      {'\n      '.join(SCG_targets_missing)}"
        report_notice(message)

    return run_data


def filter_seqs_by_length(path, cutoff):

    # getting range allowed
    lengths = [len(record.seq) for record in SeqIO.parse(path, "fasta")]
    median = statistics.median(lengths)
    min_length = round(median - (median * cutoff))
    max_length = round(median + (median * cutoff))

    # filtering and writing back out
    records = list(SeqIO.parse(path, "fasta"))
    genomes_with_hits_after_filtering = []
    filtered_records = []
    for record in records:
        if min_length <= len(record.seq) <= max_length:
            filtered_records.append(record)
            genomes_with_hits_after_filtering.append(record.id)

    out_path = path.rstrip(".fasta") + "-gene-filtered.fasta"
    with open(out_path, "w") as out_handle:
        SeqIO.write(filtered_records, out_handle, "fasta")

    return genomes_with_hits_after_filtering


def filter_seqs_by_genome_ids(path, ids_to_remove):

    records = list(SeqIO.parse(path, "fasta"))
    filtered_records = [record for record in records if record.id not in ids_to_remove]
    out_path = path.rstrip("-gene-filtered.fasta") + "-genome-filtered.fasta"

    with open(out_path, "w") as out_handle:
        SeqIO.write(filtered_records, out_handle, "fasta")

