import subprocess
import shutil
from collections import Counter
import os
from pathlib import Path
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets
from gtotree.utils.messaging import (report_processing_stage,
                                     report_pfam_searching_update)
from gtotree.utils.general import (run_pooled_stage)


def search_pfams(args, run_data):

    report_processing_stage("additional-pfam-searching", run_data)

    if run_data.additional_pfam_searching_done:
        report_pfam_searching_update(run_data)
        return run_data

    run_data = get_additional_pfam_targets(run_data)

    if len(run_data.found_pfam_targets) > 0:

        genomes_to_search = run_data.get_all_input_genomes_for_pfam_search()

        if len(genomes_to_search) > 0:

            count_rows = []

            def apply_pfam(genome, status, run_data):
                if status.get("pfam_search_failed"):
                    genome.mark_pfam_search_failed()
                else:
                    genome.mark_pfam_search_done()
                results_txt = (f"{run_data.pfam_results_dir}/individual-genome-results/"
                               f"{genome.id}/pfam-hmmsearch.txt")
                counts_list = get_pfam_counts(run_data.found_pfam_targets, results_txt)
                count_rows.append([genome.id, genome.num_genes] + counts_list)

            run_data = run_pooled_stage(
                genomes_to_search,
                _pfam_search_worker,
                apply_pfam,
                args, run_data,
            )

            cols = ['assembly_id', 'total_gene_count'] + run_data.found_pfam_targets
            pfam_counts_df = pd.DataFrame(count_rows, columns=cols)
            pfam_counts_df = pfam_counts_df.sort_values("assembly_id").reset_index(drop=True)
            pfam_counts_df.to_csv(
                f"{run_data.pfam_results_dir}/pfam-hit-counts.tsv", sep='\t', index=False)

        print("") ; print("Combining Pfam search results...".center(82))
        combine_all_pfam_hits(run_data.found_pfam_targets,
                              run_data.tmp_pfam_results_dir,
                              run_data.pfam_results_dir + "/pfam-hit-seqs")

    run_data.additional_pfam_searching_done = True
    report_pfam_searching_update(run_data)

    return run_data


def _pfam_search_worker(genome, run_data):
    """
    Runs concurrently in a thread. Touches only per-genome files: runs
    hmmsearch against the combined target-Pfam HMM, writing the tblout into this
    genome's individual-genome-results dir, then extracts the hit sequences into
    this genome's own tmp dir. Must not mutate shared run_data.
    """
    ID = genome.id
    AA_path = genome.final_AA_path

    results_dir = f"{run_data.pfam_results_dir}/individual-genome-results/{ID}"
    tmp_outpath = f"{run_data.tmp_pfam_results_dir}/{ID}"

    pfam_search_failed = run_pfam_search(run_data.all_pfam_targets_hmm_path,
                                         results_dir, AA_path)

    results_txt = f"{results_dir}/pfam-hmmsearch.txt"
    have_hits = False

    if not pfam_search_failed and os.path.isfile(results_txt):
        # any non-comment, non-empty line means at least one hit
        with open(results_txt) as f:
            have_hits = any(line.strip() and not line.startswith("#") for line in f)

        if have_hits:
            try:
                write_out_tmp_pfam_hits(results_txt, tmp_outpath, AA_path)
            except Exception:
                pfam_search_failed = True

    return {
        "pfam_search_failed": bool(pfam_search_failed),
    }


def run_pfam_search(all_pfam_targets_hmm, base_outpath, AA_file):

    os.makedirs(base_outpath, exist_ok=True)

    outpath = f"{base_outpath}/pfam-hmmsearch.txt"

    cmd = [
        "hmmsearch",
        "--cut_ga",
        "--cpu", "2",
        "--tblout", outpath,
        all_pfam_targets_hmm,
        AA_file
    ]

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=True,
        )
        pfam_search_failed = False
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode("utf-8", errors="replace") if e.stderr else ""
        print(f"[pfam hmmsearch failed for {os.path.basename(base_outpath)}] "
              f"exit {e.returncode}: {stderr.strip()}")
        pfam_search_failed = True

    return pfam_search_failed


def get_pfam_counts(pfam_ids, results_txt):

    counts = Counter()

    if not os.path.isfile(results_txt):
        return [0 for _ in pfam_ids]

    with open(results_txt, "r") as results_file:
        for line in results_file:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            pfam_id = parts[3].split(".")[0]  # Extracting the Pfam ID without version
            if pfam_id in pfam_ids:
                counts[pfam_id] += 1

    counts_list = [counts[pfam_id] for pfam_id in pfam_ids]
    return counts_list


def write_out_tmp_pfam_hits(results_txt, out_base, AA_path):

    os.makedirs(out_base, exist_ok=True)
    df = pd.read_csv(results_txt, sep=r"\s+", comment="#", header=None, engine="python")
    df = df[[0,3]]
    df.columns = ["gene", "pfam_id"]
    df["pfam_id"] = df["pfam_id"].str.split(".").str[0]
    seq_index = SeqIO.index(AA_path, "fasta")

    for pfam_id, group in df.groupby("pfam_id"):
        # get seqrecords for all gene IDs in this group
        records = []
        for gene_id in group["gene"]:
            if gene_id in seq_index:
                records.append(seq_index[gene_id])

        # writing those out for the current group
        out_path = os.path.join(out_base, f"{pfam_id}.faa")
        SeqIO.write(records, out_path, "fasta")


def combine_all_pfam_hits(pfam_ids, tmp_pfam_results_area, out_base):

    for pfam in pfam_ids:
        pattern = f"*/{pfam}.faa"
        paths = list(Path(tmp_pfam_results_area).glob(pattern))

        if paths:
            out_path = os.path.join(out_base, f"{pfam}.faa")

            with open(out_path, "wb") as out_file:
                for path in paths:
                    with open(path, "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)
