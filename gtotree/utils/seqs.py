import os
import shutil
import statistics
from Bio import SeqIO
import subprocess
from gtotree.utils.general import (remove_file_if_exists,
                                   check_file_exists_and_not_empty)

def filter_and_rename_fasta(prefix, run_data, in_path, full_path = False, max_length = 99999):
    if full_path:
        AA_infile = in_path
    else:
        AA_infile = f"{in_path}/{prefix}_protein.faa"
    AA_outpath = f"{run_data.ready_genome_files_dir}/{prefix}.faa"

    num = _filter_and_rename_fasta(AA_outpath, AA_infile, prefix, max_length)

    if num == 0:
        remove_file_if_exists(AA_outpath)
        return False, None, num, None

    else:
        if run_data.nucleotide_mode:
            nt_infile = f"{in_path}/{prefix}_cds.fasta"
            nt_outpath = f"{run_data.ready_genome_files_dir}/{prefix}.fasta"
            nt_max_length = (max_length * 3) + 3  # max length is 3x the AA max length + 3 for the stop codon
            _filter_and_rename_fasta(nt_outpath, nt_infile, prefix, nt_max_length)
        else:
            nt_outpath = None
        return True, AA_outpath, num, nt_outpath


def _filter_and_rename_fasta(outpath, infile, prefix, max_length):
    num = 0
    with open(outpath, "w") as outfile, open(infile, "r") as in_file:
        for record in SeqIO.parse(in_file, "fasta"):
            if len(record.seq) <= max_length:
                num += 1
                header = f">{prefix}_{num}"
                outfile.write(f"{header}\n{record.seq}\n")
    return num

def extract_filter_and_rename_cds_amino_acids_from_gb(prefix, input_gb, run_data, max_length = 99999):

    num = 0
    output_file = f"{run_data.ready_genome_files_dir}/{prefix}.faa"
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
            return False, None, num
        else:
            return True, output_file, num
    except:
        remove_file_if_exists(output_file)
        return False, None, num


def extract_fasta_from_gb(prefix, input_gb, run_data):

    num = 0
    output_file = f"{run_data.genbank_processing_dir}/{prefix}.fasta"

    with open(input_gb, "r") as infile, open(output_file, "w") as outfile:
        records = list(SeqIO.parse(infile, "genbank"))
        for rec in records:
            num += 1
            outfile.write(f">{prefix}_{num}\n{rec.seq}\n")


def check_target_SCGs_have_seqs(run_data, ext):

    SCG_targets_present_at_start = run_data.get_all_SCG_targets()
    SCG_targets_missing = []

    for SCG_obj in SCG_targets_present_at_start:
        SCG = SCG_obj.id
        path = run_data.found_SCG_seqs_dir + f"/{SCG}{ext}"
        present = check_file_exists_and_not_empty(path)
        if not present:
            SCG_obj.mark_removed("no seqs found or no seqs remaining after length-filtering")
            SCG_targets_missing.append(SCG)

    run_data.num_SCG_targets_removed = len(SCG_targets_missing)

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

    root, ext = os.path.splitext(path)
    out_path = f"{root}-gene-filtered{ext}"
    with open(out_path, "w") as out_handle:
        SeqIO.write(filtered_records, out_handle, "fasta")

    return genomes_with_hits_after_filtering


def filter_seqs_by_genome_ids(path, ids_to_remove, out_path):

    records = list(SeqIO.parse(path, "fasta"))
    filtered_records = [record for record in records if record.id not in ids_to_remove]

    with open(out_path, "w") as out_handle:
        SeqIO.write(filtered_records, out_handle, "fasta")


def run_muscle(id, run_data, inpath, outpath, log_path):
    if run_data.use_muscle_super5:
        align_type = "-super5"
    else:
        align_type = "-align"

    cmd = [
        "muscle",
        f"{align_type}",
        f"{inpath}",
        "-output",
        f"{outpath}",
        "-threads", f"{run_data.num_muscle_threads}"
    ]

    try:
        with open(log_path, "w") as log_out:
            subprocess.run(cmd, stdout=log_out, stderr=log_out, check=True)
        muscle_failed = False
    except Exception:
        muscle_failed = True

    return muscle_failed


def run_trimal(inpath, output, log_path):

    cmd = [
        "trimal",
        "-in", f"{inpath}",
        "-out", f"{output}",
        "-automated1"
    ]

    try:
        with open(log_path, "w") as log_out:
            subprocess.run(cmd, stdout=log_out, stderr=log_out, check=True)
        trimal_failed = False
    except:
        trimal_failed = True

    return trimal_failed


def add_needed_gap_seqs(run_data, inpath, outpath):

    all_needed_ids = run_data.get_all_remaining_input_genome_ids()
    records = list(SeqIO.parse(inpath, "fasta"))
    align_len = len(records[0].seq)
    record_dict = {record.id: record for record in records}

    with open(outpath, "w") as out_handle:
        for req_id in all_needed_ids:
            if req_id in record_dict:
                out_handle.write(f">{req_id}\n{record_dict[req_id].seq}\n")
            else:
                seq = "-" * align_len
                out_handle.write(f">{req_id}\n{seq}\n")


def concatenate_alignments(run_data):
    SCGs_ready_for_cat = run_data.get_all_SCG_targets_ready_for_concatenation()
    SCG_IDs_ready_for_cat = [SCG.id for SCG in SCGs_ready_for_cat]
    SCG_paths_to_cat = [run_data.found_SCG_seqs_dir + f"/{SCG_id}-final{run_data.general_ext}" for SCG_id in SCG_IDs_ready_for_cat]

    # initializing dictionary that will hold headers as keys and a list of all seqs to be cat'd as values
    dict_of_genomes = {genome_id: [] for genome_id in run_data.get_all_remaining_input_genome_ids()}

    # iterating through all files adding seqs to the dictionary
    for file in SCG_paths_to_cat:
        with open(file) as fasta:
            curr_header=""
            for line in fasta:
                line = line.strip()
                if line.startswith(">"):
                    curr_header=line.lstrip(">")
                else:
                    dict_of_genomes[curr_header].append(line)

    # writing out the concatenated (horizontally) alignment file
    if not run_data.nucleotide_mode:
        output_path = os.path.join(run_data.output_dir, "aligned-SCGs.faa")
        spacer = "XXXXX"
    else:
        output_path = os.path.join(run_data.output_dir, "aligned-SCGs.fasta")
        spacer = "NNNNNN"

    with open(output_path, "w") as out:
        for header, seqs in dict_of_genomes.items():
            out.write(">" + header + "\n")
            out.write(spacer.join(seqs) + "\n")

    run_data.concatenated_alignment_path = output_path

    return run_data, dict_of_genomes, SCG_IDs_ready_for_cat


def gen_partitions_file(run_data, SCG_IDs, dict_of_genomes):

    # all are same length, so just need one genome entry, then to count the bases per element in dict values list, and add for the spacers
    # getting all alignment lengths
    alignment_lengths_list = [len(x) for x in list(dict_of_genomes.values())[0]]

    curr_start = 1
    curr_stop = 0

    mol_type = "AA" if not run_data.nucleotide_mode else "DNA"
    spacer_addition = 6 if not run_data.nucleotide_mode else 7

    with open(run_data.run_files_dir + "/partitions.txt", "w") as out:

        for i in range(0,len(SCG_IDs)):
            curr_stop = curr_start + alignment_lengths_list[i] - 1

            out.write(f"{mol_type}, {SCG_IDs[i]} = {curr_start}-{curr_stop}\n")
            curr_start = curr_stop + spacer_addition


def swap_labels_in_alignment(run_data):
    orig_alignment_path = os.path.join(run_data.output_dir, f"aligned-SCGs{run_data.general_ext}")
    backup_alignment_path = os.path.join(run_data.run_files_dir, f"aligned-SCGs{run_data.general_ext}")
    new_alignment_path = os.path.join(run_data.output_dir, f"aligned-SCGs-mod-names{run_data.general_ext}")

    sequences = list(SeqIO.parse(orig_alignment_path, "fasta"))
    with open(new_alignment_path, "w") as fh:
        for seq in sequences:
            if seq.id in run_data.mapping_dict:
                seq.id = run_data.mapping_dict[seq.id]
            fh.write(f">{seq.id}\n{str(seq.seq)}\n")

    shutil.move(orig_alignment_path, backup_alignment_path)

    run_data.final_alignment_path = new_alignment_path

    return run_data


def copy_gene_alignments(run_data):

    all_SCG_targets = run_data.get_all_SCG_targets_remaining()

    for SCG in all_SCG_targets:

        source_path = os.path.join(run_data.found_SCG_seqs_dir, f"{SCG.id}-final.fasta")
        dest_path = os.path.join(run_data.individual_gene_alignments_dir, f"{SCG.id}-aligned.fasta")
        if os.path.exists(source_path):
            shutil.copyfile(source_path, dest_path)


def get_alignment_length(path):

    with open(path, "r") as f:
        next(f)
        num_sites = len(next(f).strip())

    return num_sites


def pull_out_corresponding_nt_seqs(in_AA_fasta, in_nt_fasta, out_nt_fasta):

    aa_records = SeqIO.index(in_AA_fasta, "fasta")
    nt_records = SeqIO.index(in_nt_fasta, "fasta")

    with open(out_nt_fasta, "w") as out_handle:
        for aa_id in aa_records:
            if aa_id in nt_records:
                nt_record = nt_records[aa_id]
                SeqIO.write(nt_record, out_handle, "fasta")