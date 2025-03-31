import shutil
import subprocess
import pandas as pd
from Bio import SeqIO
import pyhmmer.easel as easel #type: ignore
from gtotree.utils.messaging import (report_processing_stage,
                                     report_early_exit,
                                     report_hmm_search_update)
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)
from gtotree.utils.seqs import check_target_SCGs_have_seqs

def search_hmms(args, run_data):

    report_processing_stage("hmm-search")

    num_genomes_to_search = len(run_data.get_all_input_genomes_for_hmm_search())

    if num_genomes_to_search > 0:

        # writing run_data to file so it can be accessed by snakemake
        write_run_data(run_data)
        snakefile = get_snakefile_path("search-hmms.smk")
        description = "Searching HMMs"

        run_snakemake(snakefile, num_genomes_to_search, args, run_data, description)

        run_data = read_run_data(run_data.run_data_path)
        capture_hmm_search_failures(run_data)

    report_hmm_search_update(run_data)
    run_data = check_target_SCGs_have_seqs(run_data, ".fasta")

    return run_data


def run_hmm_search(id, run_data, inpath, outpath):
    cmd = [
        "hmmsearch",
        "--cut_ga",
        "--cpu", f"{run_data.num_hmm_cpus}",
        "--tblout", outpath,
        run_data.hmm_path,
        inpath
    ]

    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
        hmm_search_failed = False
    except:
        hmm_search_failed = True

    return hmm_search_failed


def read_hmmer_results(inpath):
    colnames = ["gene_id", "target_SCG", "accession", "evalue"]
    df = pd.read_csv(inpath, sep=r'\s+',
                     comment="#", header=None,
                     usecols=[0,2,3,4],
                     names=colnames)

    return df


def parse_hmmer_results(inpath, run_data):

    df = read_hmmer_results(inpath)

    remaining_SCG_targets = [SCG_target.id for SCG_target in run_data.get_all_SCG_targets_remaining()]
    dict_of_hit_counts = dict.fromkeys(remaining_SCG_targets, 0)
    for scg in remaining_SCG_targets:
        dict_of_hit_counts[scg] = df[df["target_SCG"] == scg].shape[0]

    if run_data.best_hit_mode:
        report_early_exit("Best hit mode not yet implemented!")
    else:
        pass

    dict_of_hit_gene_ids = dict.fromkeys(remaining_SCG_targets, None)

    num_SCG_hits = 0
    for scg, count in dict_of_hit_counts.items():
        ## will need additional logic for best-hit mode
        if count == 1:
            num_SCG_hits += 1
            gene_id = df.loc[df["target_SCG"] == scg, "gene_id"].iloc[0]
            dict_of_hit_gene_ids[scg] = gene_id

    return dict_of_hit_counts, dict_of_hit_gene_ids, num_SCG_hits


def get_seqs(dict_of_hit_gene_ids, AA_path):

    try:
        hit_seqs_dict = dict.fromkeys(dict_of_hit_gene_ids.keys(), None)
        reverse_lookup = {seq_id: target for target, seq_id in dict_of_hit_gene_ids.items()}
        AA_alphabet = easel.Alphabet.amino()

        with easel.SequenceFile(AA_path, digital=True, alphabet=AA_alphabet) as seq_file:
            for seq in seq_file:
                seq_id = seq.name.decode("utf8")
                if seq_id in reverse_lookup:
                    target_scg = reverse_lookup[seq_id]
                    hit_seqs_dict[target_scg] = AA_alphabet.decode(seq)
        extract_seqs_failed = False
    except:
        extract_seqs_failed = True
        hit_seqs_dict = None

    return hit_seqs_dict, extract_seqs_failed


def start_combined_SCG_hit_count_tab(run_data):
    out_file = f"{run_data.output_dir}/SCG-hit-counts.tsv"
    with open(out_file, "w") as outfile:
        outfile.write("assembly_id\t" + "\t".join(run_data.initial_SCG_targets) + "\n")


def add_to_combined_SCG_hit_count_tab(genome_id, run_data):
    table_file = f"{run_data.output_dir}/SCG-hit-counts.tsv"
    row_to_add_file = f"{run_data.hmm_results_dir}/{genome_id}/SCG-hit-counts.txt"
    with open(table_file, 'ab') as outfile:
        with open(row_to_add_file, 'rb') as infile:
            shutil.copyfileobj(infile, outfile)


def write_out_SCG_hit_seqs(genome_id, run_data):
    input_fasta = f"{run_data.hmm_results_dir}/{genome_id}/SCG-hits.fasta"
    out_dir = f"{run_data.found_SCG_seqs_dir}"
    remaining_SCG_targets = [SCG_target.id for SCG_target in run_data.get_all_SCG_targets_remaining()]
    output_paths = {target_SCG: f"{out_dir}/{target_SCG}.fasta" for target_SCG in remaining_SCG_targets}

    output_handles = {target: open(path, "a") for target, path in output_paths.items()}

    print(output_handles)

    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            target = record.id
            print(target)
            if target in output_handles:
                output_handles[target].write(f">{genome_id}\n{record.seq}\n")

    for handle in output_handles.values():
        handle.close()


def capture_hmm_search_failures(run_data):
    if len(run_data.get_failed_hmm_search_paths()) > 0:
        with open(run_data.run_files_dir + "/inputs-that-failed-at-the-hmm-search.txt", "w") as fail_file:
            for genome_id in run_data.get_failed_hmm_search_paths():
                fail_file.write(genome_id + "\n")
