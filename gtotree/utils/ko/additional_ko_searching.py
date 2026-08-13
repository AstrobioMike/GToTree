import subprocess
import shutil
from collections import Counter
import os
from pathlib import Path
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
from gtotree.utils.misc.general import (search_threads_per_genome,
                                   atomic_write_text,
                                   record_search_failure)


def _ko_search_worker(genome, run_data, aa_path=None):
    """
    Per-genome KO search (kofamscan) in a worker thread. Writes this genome's own
    results tsv and per-KO tmp hit fastas. Returns a status dict. Must not mutate
    shared run_data state.

    `aa_path` overrides the path on the GenomeData, for the fused processing stage.

    Wrapped so it cannot raise. A worker exception aborts the entire pooled stage.
    """
    try:
        return _ko_search_worker_inner(genome, run_data, aa_path)
    except BaseException as e:
        return {"ko_search_failed": True, "error": f"{type(e).__name__}: {e}"}


def _ko_search_worker_inner(genome, run_data, aa_path=None):
    ID = genome.id
    AA_path = aa_path if aa_path is not None else genome.final_AA_path

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
        "--cpu", str(search_threads_per_genome(None)),
        "-f", "mapper",
        "--no-report-unannotated",
        "--tmp-dir", tmp_path,
        "-o", outpath,
        AA_file
    ]

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=True,
        )
        kofamscan_failed = False
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode("utf-8", errors="replace") if e.stderr else ""
        record_search_failure(base_outpath, "kofamscan",
                              f"exit {e.returncode}: {stderr.strip()}")
        kofamscan_failed = True
    except OSError as e:
        # binary missing / not executable
        record_search_failure(base_outpath, "kofamscan", f"could not be run: {e}")
        kofamscan_failed = True

    shutil.rmtree(tmp_path, ignore_errors=True)

    return kofamscan_failed


def write_ko_counts_table(run_data):
    """
    Write ko-hit-counts.tsv from the per-genome result files
    """
    genomes = [gd for gd in run_data.all_input_genomes
               if gd.ko_search_done and not gd.removed]

    count_rows = []
    for gd in genomes:
        ko_results_tsv = (f"{run_data.ko_results_dir}/individual-genome-results/"
                          f"{gd.id}/kofamscan-results.tsv")
        counts_list = get_ko_counts(run_data.found_ko_targets, ko_results_tsv)
        count_rows.append([gd.id, gd.num_genes] + counts_list)

    cols = ['genome_id', 'total_gene_count'] + run_data.found_ko_targets
    ko_counts_df = pd.DataFrame(count_rows, columns=cols)
    if not ko_counts_df.empty:
        ko_counts_df = ko_counts_df.sort_values("genome_id").reset_index(drop=True)
    atomic_write_text(
        f"{run_data.ko_results_dir}/ko-hit-counts.tsv",
        lambda f: ko_counts_df.to_csv(f, sep='\t', index=False))


def write_out_failed_ko_targets(run_data):
    if len(run_data.failed_ko_targets) > 0:
        with open(run_data.run_files_dir + "/failed-ko-targets.txt", "w") as fail_file:
            for KO in run_data.failed_ko_targets:
                fail_file.write(KO + "\n")


def get_ko_counts(ko_ids, results_tsv):

    counts = Counter()

    if not os.path.isfile(results_tsv):
        return [0 for _ in ko_ids]

    with open(results_tsv) as results_file:
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
