import os
import shutil
import statistics
from Bio import SeqIO # type: ignore
import subprocess
from gtotree.utils.misc.general import (remove_file_if_exists,
                                   file_is_usable_else_clear,
                                   atomic_write_text)
from gtotree.utils.misc.stages import SCGRemovalStage


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

    if run_data.nucleotide_mode:
        nt_infile = f"{in_path}/{prefix}_cds.fasta"
        nt_outpath = f"{run_data.ready_genome_files_dir}/{prefix}.fasta"
        nt_max_length = (max_length * 3) + 3  # max length is 3x the AA max length + 3 for the stop codon
        _filter_and_rename_fasta(nt_outpath, nt_infile, prefix, nt_max_length)
    else:
        nt_outpath = None
    return True, AA_outpath, num, nt_outpath


def fasta_has_single_record(path) -> bool:
    count = 0
    for _ in SeqIO.parse(path, "fasta"):
        count += 1
        if count > 1:
            return False
    return True


def _filter_and_rename_fasta(outpath, infile, prefix, max_length):
    num = 0
    with open(outpath, "w") as outfile, open(infile) as in_file:
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
        with open(input_gb) as infile, open(output_file, "w") as outfile:
            for rec in SeqIO.parse(infile, "genbank"):
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
        return True, output_file, num
    except:
        remove_file_if_exists(output_file)
        return False, None, num


def extract_fasta_from_gb(prefix, input_gb, run_data):

    num = 0
    output_file = f"{run_data.genbank_processing_dir}/{prefix}.fasta"

    with open(input_gb) as infile, open(output_file, "w") as outfile:
        for rec in SeqIO.parse(infile, "genbank"):
            num += 1
            outfile.write(f">{prefix}_{num}\n{rec.seq}\n")


# "no usable seqs at this path" is the same *condition* at all three of the stages
# below, but it means something different at each one, and that difference is what
# ends up in SCG-info.tsv
NO_SEQS_REASONS = {
    SCGRemovalStage.NO_HITS:       "no hits in any genome",
    SCGRemovalStage.GENE_FILTER:   "no hits remaining after gene-length filtering (`-c`)",
    SCGRemovalStage.GENOME_FILTER: "no hits remaining after genome filtering (`-G`)",
}


def check_target_SCGs_have_seqs(run_data, ext, stage, reason_fn=None):
    """
    Drop any SCG-set with no usable sequences at `ext`, recording a stage-specific reason

    `reason_fn(scg) -> str` overrides the default for the stage. The NO_HITS caller uses
    it to separate "nothing hit this at all" from "things hit it, but never as a single
    copy", which are the same empty file but very different situations for the user.
    """
    for SCG_obj in run_data.get_all_SCG_targets():

        if SCG_obj.removed:
            continue

        path = run_data.found_SCG_seqs_dir + f"/{SCG_obj.id}{ext}"
        if not file_is_usable_else_clear(path):
            reason = reason_fn(SCG_obj) if reason_fn else NO_SEQS_REASONS[stage]
            SCG_obj.mark_removed(reason, stage)

    return run_data


def _fasta_seq_lengths(path):
    """
    Sequence lengths from a FASTA, read as text without building any record objects
    """
    lengths = []
    current = 0
    started = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    lengths.append(current)
                current = 0
                started = True
            else:
                current += len(line.strip())
    if started:
        lengths.append(current)
    return lengths


def filter_seqs_by_length(path, cutoff):
    """
    Drop sequences whose length is outside `cutoff` of the set's median, writing the
    survivors to a '-gene-filtered' file and returning their genome ids
    """
    lengths = _fasta_seq_lengths(path)
    if not lengths:
        raise ValueError(f"no sequences found in {path}")

    median = statistics.median(lengths)
    min_length = round(median - (median * cutoff))
    max_length = round(median + (median * cutoff))

    root, ext = os.path.splitext(path)
    out_path = f"{root}-gene-filtered{ext}"

    genomes_with_hits_after_filtering = []
    with open(out_path, "w") as out_handle:
        for record in SeqIO.parse(path, "fasta"):
            if min_length <= len(record.seq) <= max_length:
                genomes_with_hits_after_filtering.append(record.id)
                out_handle.write(f">{record.id}\n{record.seq}\n")

    return genomes_with_hits_after_filtering


def count_fasta_records(path):
    """
    How many records are in a FASTA, without building any record objects
    """
    count = 0
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    count += 1
    except OSError:
        return 0
    return count


def filter_seqs_by_genome_ids(path, ids_to_remove, out_path):
    """
    Copy `path` to `out_path`, dropping records whose id is in `ids_to_remove`
    """
    kept = 0
    with open(out_path, "w") as out_handle:
        for record in SeqIO.parse(path, "fasta"):
            if record.id not in ids_to_remove:
                out_handle.write(f">{record.id}\n{record.seq}\n")
                kept += 1
    return kept


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
    """
    Write one record per remaining genome, gap-filling any that this SCG has no hit for
    """
    all_needed_ids = run_data.get_all_remaining_input_genome_ids()

    seq_by_id = {}
    align_len = None
    for record in SeqIO.parse(inpath, "fasta"):
        seq = str(record.seq)
        if align_len is None:
            align_len = len(seq)
        seq_by_id[record.id] = seq

    if align_len is None:
        raise ValueError(f"no sequences found in {inpath}")

    gap_row = "-" * align_len
    with open(outpath, "w") as out_handle:
        for req_id in all_needed_ids:
            out_handle.write(f">{req_id}\n{seq_by_id.get(req_id, gap_row)}\n")


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

    def _write_concatenated(out):
        for header, seqs in dict_of_genomes.items():
            out.write(">" + header + "\n")
            out.write(spacer.join(seqs) + "\n")

    atomic_write_text(output_path, _write_concatenated)

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

    txt_lines = []
    nex_lines = ["#NEXUS\n", "begin sets;\n"]

    for i in range(0, len(SCG_IDs)):

        curr_stop = curr_start + alignment_lengths_list[i] - 1

        txt_lines.append(f"{mol_type}, {SCG_IDs[i]} = {curr_start}-{curr_stop}\n")
        nex_lines.append(f"  charset {SCG_IDs[i]} = {curr_start}-{curr_stop};\n")

        curr_start = curr_stop + spacer_addition

    nex_lines.append("end;\n")

    atomic_write_text(run_data.run_files_dir + "/partitions.txt",
                      lambda f: f.writelines(txt_lines))
    atomic_write_text(run_data.run_files_dir + "/partitions.nex",
                      lambda f: f.writelines(nex_lines))


def swap_labels_in_alignment(run_data):
    ext = run_data.general_ext
    orig_alignment_path = os.path.join(run_data.output_dir, f"aligned-SCGs{ext}")
    backup_alignment_path = os.path.join(run_data.run_files_dir, f"aligned-SCGs{ext}")
    new_alignment_path = os.path.join(run_data.output_dir, f"aligned-SCGs-mod-names{ext}")

    if file_is_usable_else_clear(orig_alignment_path):
        source_path = orig_alignment_path
    elif file_is_usable_else_clear(backup_alignment_path):
        source_path = backup_alignment_path
    else:
        raise FileNotFoundError(
            "the concatenated alignment couldn't be found at either "
            f'"{orig_alignment_path}" or "{backup_alignment_path}"')

    def _write_relabeled(fh):
        for seq in SeqIO.parse(source_path, "fasta"):
            label = run_data.mapping_dict.get(seq.id, seq.id)
            fh.write(f">{label}\n{seq.seq}\n")

    atomic_write_text(new_alignment_path, _write_relabeled)

    if source_path == orig_alignment_path:
        os.replace(orig_alignment_path, backup_alignment_path)

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

    with open(path) as f:
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
