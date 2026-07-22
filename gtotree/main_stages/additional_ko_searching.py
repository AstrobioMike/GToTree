import subprocess
import shutil
from collections import Counter
import os
from pathlib import Path
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
from gtotree.utils.ko.ko_handling import parse_kofamscan_targets
from gtotree.utils.messaging import (report_processing_stage,
                                     report_ko_searching_update)
from gtotree.utils.general import (run_pooled_stage)


def search_kos(args, run_data):

    report_processing_stage("additional-ko-searching", run_data)

    if run_data.additional_ko_searching_done:
        report_ko_searching_update(run_data)
        return run_data

    # generating the subset of target KOs
    run_data = parse_kofamscan_targets(run_data)

    write_out_failed_ko_targets(run_data)

    if len(run_data.found_ko_targets) > 0:

        genomes_to_search = run_data.get_all_input_genomes_for_ko_search()

        if len(genomes_to_search) > 0:

            count_rows = []

            def apply_ko(genome, status, run_data):
                if status.get("ko_search_failed"):
                    genome.mark_ko_search_failed()
                else:
                    genome.mark_ko_search_done()
                ko_results_tsv = f"{run_data.ko_results_dir}/individual-genome-results/{genome.id}/kofamscan-results.tsv"
                counts_list = get_ko_counts(run_data.found_ko_targets, ko_results_tsv)
                count_rows.append([genome.id, genome.num_genes] + counts_list)

            run_data = run_pooled_stage(
                genomes_to_search,
                _ko_search_worker,
                apply_ko,
                args, run_data,
            )

            cols = ['assembly_id', 'total_gene_count'] + run_data.found_ko_targets
            ko_counts_df = pd.DataFrame(count_rows, columns=cols)
            ko_counts_df.to_csv(
                f"{run_data.ko_results_dir}/ko-hit-counts.tsv", sep='\t', index=False)

        print("") ; print("Combining KO search results...".center(82))

        combine_all_ko_hits(run_data.found_ko_targets,
                            run_data.tmp_ko_results_dir,
                            run_data.ko_results_dir + "/ko-hit-seqs")

    run_data.additional_ko_searching_done = True
    report_ko_searching_update(run_data)

    return run_data


def _ko_search_worker(genome, run_data):
    """
    Per-genome KO search (kofamscan) in a worker thread. Writes this genome's own
    results tsv and per-KO tmp hit fastas. Returns a status dict. Must not mutate
    shared run_data state.
    """
    ID = genome.id
    AA_path = genome.final_AA_path

    base_outpath = f"{run_data.ko_results_dir}/individual-genome-results/{ID}/"
    os.makedirs(base_outpath, exist_ok=True)
    ko_search_failed = run_ko_search(run_data.target_ko_profiles_dir,
                                     run_data.target_kos_tsv,
                                     base_outpath,
                                     AA_path)

    if not ko_search_failed:
        results_tsv = f"{run_data.ko_results_dir}/individual-genome-results/{ID}/kofamscan-results.tsv"
        out_base = f"{run_data.tmp_ko_results_dir}/{ID}/"
        write_out_tmp_ko_hits(results_tsv, out_base, AA_path)

    return {"ko_search_failed": bool(ko_search_failed)}


def run_ko_search(profiles_dir, ko_file, base_outpath, AA_file):

    outpath = f"{base_outpath}/kofamscan-results.tsv"
    tmp_path = f"{base_outpath}/kofamscan-tmp"
    cmd = [
        "exec_annotation",
        "-p", profiles_dir,
        "-k", ko_file,
        "--cpu", str(2),
        "-f", "mapper",
        "--no-report-unannotated",
        "--tmp-dir", tmp_path,
        "-o", outpath,
        AA_file
    ]

    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
        kofamscan_failed = False
    except:
        kofamscan_failed = True

    shutil.rmtree(tmp_path, ignore_errors=True)

    return kofamscan_failed


def write_out_failed_ko_targets(run_data):
    if len(run_data.failed_ko_targets) > 0:
        with open(run_data.run_files_dir + "/failed-ko-targets.txt", "w") as fail_file:
            for KO in run_data.failed_ko_targets:
                fail_file.write(KO + "\n")


def get_ko_counts(ko_ids, results_tsv):

    counts = Counter()

    if not os.path.isfile(results_tsv):
        return [0 for _ in ko_ids]

    with open(results_tsv, "r") as results_file:
        for line in results_file:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            ko_id = parts[1]
            if ko_id in ko_ids:
                counts[ko_id] += 1

    counts_list = [counts[ko_id] for ko_id in ko_ids]
    return counts_list


def write_out_tmp_ko_hits(results_tsv, out_base, AA_path):
    os.makedirs(out_base, exist_ok=True)
    df = pd.read_csv(results_tsv, sep="\t", header=None, names=["gene", "ko_id"])
    seq_index = SeqIO.index(AA_path, "fasta")

    for ko_id, group in df.groupby("ko_id"):
        # get seqrecords for all gene IDs in this group
        records = []
        for gene_id in group["gene"]:
            if gene_id in seq_index:
                records.append(seq_index[gene_id])

        # writing those out for the current group
        out_path = os.path.join(out_base, f"{ko_id}.faa")
        SeqIO.write(records, out_path, "fasta")


def combine_all_ko_hits(ko_ids, tmp_ko_results_area, out_base):

    for ko in ko_ids:
        pattern = f"*/{ko}.faa"
        paths = list(Path(tmp_ko_results_area).glob(pattern))
        if paths:
            out_path = os.path.join(out_base, f"{ko}.faa")

            with open(out_path, "wb") as out_file:
                for path in paths:
                    with open(path, "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)
