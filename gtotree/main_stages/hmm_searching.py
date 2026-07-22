import os
import json
import shutil
import subprocess
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
import pyhmmer.easel as easel #type: ignore
from gtotree.utils.messaging import (report_processing_stage,
                                     report_hmm_search_update)
from gtotree.utils.general import (run_pooled_stage)
from gtotree.utils.seqs import check_target_SCGs_have_seqs


def search_hmms(args, run_data):

    report_processing_stage("hmm-search", run_data)

    genomes_to_search = run_data.get_all_input_genomes_for_hmm_search()

    if len(genomes_to_search) > 0:

        start_combined_SCG_hit_count_tab(run_data)

        run_data = run_pooled_stage(
            genomes_to_search,
            _hmm_search_worker,
            _apply_hmm_search_result,
            args, run_data,
        )

        capture_hmm_search_failures(run_data)

    report_hmm_search_update(run_data)

    run_data = check_target_SCGs_have_seqs(run_data, run_data.general_ext)

    return run_data


def _hmm_search_worker(genome, run_data):

    ID = genome.id
    AA_path = genome.final_AA_path
    out_dir = f"{run_data.hmm_results_dir}/{ID}"
    os.makedirs(out_dir, exist_ok=True)
    hmm_out_path = f"{out_dir}/SCG-hits-hmm.txt"

    hmm_search_failed = run_hmm_search(ID, run_data, AA_path, hmm_out_path)

    num_SCG_hits = 0
    num_unique_SCG_hits = 0
    extract_seqs_failed = False

    if not hmm_search_failed:
        (dict_of_hit_counts, dict_of_hit_gene_ids,
         num_SCG_hits, num_unique_SCG_hits) = parse_hmmer_results(hmm_out_path, run_data)

        with open(f"{out_dir}/SCG-hit-counts.txt", 'w') as f:
            f.write(f"{ID}")
            for count in dict_of_hit_counts.values():
                f.write(f"\t{count}")
            f.write("\n")

        AA_hit_seqs_dict, extract_seqs_failed = get_seqs(dict_of_hit_gene_ids, AA_path)

        if not extract_seqs_failed:
            with open(f"{out_dir}/SCG-hits.faa", 'w') as f:
                for gene_id, seq in AA_hit_seqs_dict.items():
                    if seq is not None:
                        f.write(f">{gene_id}\n{seq}\n")
            if run_data.nucleotide_mode:
                nt_hit_seqs_dict, _ = get_seqs(dict_of_hit_gene_ids, genome.final_nt_path)
                with open(f"{out_dir}/SCG-hits.fasta", 'w') as f:
                    for gene_id, seq in nt_hit_seqs_dict.items():
                        if seq is not None:
                            f.write(f">{gene_id}\n{seq}\n")

    return {
        "hmm_search_failed": bool(hmm_search_failed),
        "extract_seqs_failed": bool(extract_seqs_failed),
        "num_SCG_hits": int(num_SCG_hits),
        "num_unique_SCG_hits": int(num_unique_SCG_hits),
    }


def _apply_hmm_search_result(genome, status, run_data):

    hmm_search_failed = bool(status.get("hmm_search_failed", False))
    extract_seqs_failed = bool(status.get("extract_seqs_failed", False))

    if hmm_search_failed:
        genome.mark_hmm_search_failed()
        genome.mark_removed("HMM search failed")
        genome.num_SCG_hits = 0
    else:
        add_to_combined_SCG_hit_count_tab(genome.id, run_data)
        if extract_seqs_failed:
            genome.mark_extract_seqs_failed()
            genome.num_SCG_hits = 0
            genome.mark_removed("extracting sequences after HMM search failed")
        else:
            write_out_SCG_hit_seqs(genome.id, run_data)
            genome.mark_hmm_search_done()
            genome.num_SCG_hits = int(status.get("num_SCG_hits", 0))
            genome.num_unique_SCG_hits = int(status.get("num_unique_SCG_hits", 0))


def run_hmm_search(id, run_data, inpath, outpath):
    cmd = [
        "hmmsearch",
        "--cut_ga",
        "--cpu", "2",
        "--tblout", outpath,
        run_data.hmm_path,
        inpath
    ]

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=True,
        )
        hmm_search_failed = False
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode("utf-8", errors="replace") if e.stderr else ""
        print(f"[hmmsearch failed for {id}] exit {e.returncode}: {stderr.strip()}")
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

    dict_of_hit_gene_ids = dict.fromkeys(remaining_SCG_targets, None)

    num_SCG_hits = 0
    num_unique_SCG_hits = 0
    for scg, count in dict_of_hit_counts.items():

        if count == 1:
            num_SCG_hits += 1
            num_unique_SCG_hits += 1
            gene_id = df.loc[df["target_SCG"] == scg, "gene_id"].iloc[0]
            dict_of_hit_gene_ids[scg] = gene_id

        elif count > 1:
            num_SCG_hits += 1
            if run_data.best_hit_mode:
                gene_id = df.loc[df["target_SCG"] == scg, "gene_id"].iloc[0]
                dict_of_hit_gene_ids[scg] = gene_id

    return dict_of_hit_counts, dict_of_hit_gene_ids, num_SCG_hits, num_unique_SCG_hits


def get_seqs(dict_of_hit_gene_ids, path):

    try:
        hit_seqs_dict = dict.fromkeys(dict_of_hit_gene_ids.keys(), None)
        reverse_lookup = {seq_id: target for target, seq_id in dict_of_hit_gene_ids.items()}
        easel_alphabet = easel.Alphabet.amino() if path.endswith(".faa") else easel.Alphabet.dna()

        with easel.SequenceFile(path, digital=True, alphabet=easel_alphabet) as seq_file:
            for seq in seq_file:
                # pyhmmer changed seq.name from bytes to str across some versions i was testing on, but i didn't pin it down
                # so just handling either
                seq_id = seq.name.decode("utf8") if isinstance(seq.name, bytes) else seq.name
                if seq_id in reverse_lookup:
                    target_scg = reverse_lookup[seq_id]
                    hit_seqs_dict[target_scg] = easel_alphabet.decode(seq)
        extract_seqs_failed = False
    except:
        extract_seqs_failed = True
        hit_seqs_dict = None

    return hit_seqs_dict, extract_seqs_failed


def start_combined_SCG_hit_count_tab(run_data):
    out_file = f"{run_data.output_dir}/SCG-hit-counts.tsv"
    target_SCG_ids = [SCG_target.id for SCG_target in run_data.get_all_SCG_targets_remaining()]
    with open(out_file, "w") as outfile:
        outfile.write("assembly_id\t" + "\t".join(target_SCG_ids) + "\n")


def add_to_combined_SCG_hit_count_tab(genome_id, run_data):
    table_file = f"{run_data.output_dir}/SCG-hit-counts.tsv"
    row_to_add_file = f"{run_data.hmm_results_dir}/{genome_id}/SCG-hit-counts.txt"
    with open(table_file, 'ab') as outfile:
        with open(row_to_add_file, 'rb') as infile:
            shutil.copyfileobj(infile, outfile)


def write_out_SCG_hit_seqs(genome_id, run_data):
    input_fasta = f"{run_data.hmm_results_dir}/{genome_id}/SCG-hits{run_data.general_ext}"
    out_dir = f"{run_data.found_SCG_seqs_dir}"
    remaining_SCG_targets = [SCG_target.id for SCG_target in run_data.get_all_SCG_targets_remaining()]
    output_paths = {target_SCG: f"{out_dir}/{target_SCG}{run_data.general_ext}" for target_SCG in remaining_SCG_targets}

    output_handles = {target: open(path, "a") for target, path in output_paths.items()}

    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            target = record.id
            if target in output_handles:
                output_handles[target].write(f">{genome_id}\n{record.seq}\n")

    for handle in output_handles.values():
        handle.close()


def capture_hmm_search_failures(run_data):
    if len(run_data.get_failed_hmm_search_paths()) > 0:
        with open(run_data.run_files_dir + "/inputs-that-failed-at-the-hmm-search.txt", "w") as fail_file:
            for genome_id in run_data.get_failed_hmm_search_paths():
                fail_file.write(genome_id + "\n")
