import shutil
from collections import Counter
import os
from pathlib import Path
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
from gtotree.utils.misc.general import (search_threads_per_genome,
                                   atomic_write_text,
                                   record_search_failure)
from gtotree.utils.hmms.hmm_searching_engine import search_one_genome


def _pfam_search_worker(genome, run_data, aa_path=None, pressed_base=None):
    """
    Runs concurrently in a thread. Touches only per-genome files: runs
    hmmsearch against the combined target-Pfam HMM, writing the tblout into this
    genome's individual-genome-results dir, then extracts the hit sequences into
    this genome's own tmp dir. Must not mutate shared run_data.

    `aa_path` overrides the path recorded on the GenomeData, for the fused processing
    stage which searches a genome in the same worker that just produced its FASTA.

    Wrapped so it cannot raise. A worker exception aborts the entire pooled stage.
    """
    try:
        return _pfam_search_worker_inner(genome, run_data, aa_path, pressed_base)
    except BaseException as e:
        return {"pfam_search_failed": True, "error": f"{type(e).__name__}: {e}"}


def _pfam_search_worker_inner(genome, run_data, aa_path=None, pressed_base=None):
    ID = genome.id
    AA_path = aa_path if aa_path is not None else genome.final_AA_path

    results_dir = f"{run_data.pfam_results_dir}/individual-genome-results/{ID}"
    tmp_outpath = f"{run_data.tmp_pfam_results_dir}/{ID}"

    pfam_search_failed = run_pfam_search(run_data.all_pfam_targets_hmm_path,
                                         results_dir, AA_path,
                                         pressed_base=pressed_base)

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


def write_pfam_counts_table(run_data):
    """
    Write pfam-hit-counts.tsv from the per-genome result files.

    Idempotent, and complete regardless of how much of the work happened in this
    invocation versus an earlier interrupted one.
    """
    genomes = [gd for gd in run_data.all_input_genomes
               if gd.pfam_search_done and not gd.removed]

    count_rows = []
    for gd in genomes:
        results_txt = (f"{run_data.pfam_results_dir}/individual-genome-results/"
                       f"{gd.id}/pfam-hmmsearch.txt")
        counts_list = get_pfam_counts(run_data.found_pfam_targets, results_txt)
        count_rows.append([gd.id, gd.num_genes] + counts_list)

    cols = ['genome_id', 'total_gene_count'] + run_data.found_pfam_targets
    pfam_counts_df = pd.DataFrame(count_rows, columns=cols)
    if not pfam_counts_df.empty:
        pfam_counts_df = pfam_counts_df.sort_values("genome_id").reset_index(drop=True)
    atomic_write_text(
        f"{run_data.pfam_results_dir}/pfam-hit-counts.tsv",
        lambda f: pfam_counts_df.to_csv(f, sep='\t', index=False))


def run_pfam_search(all_pfam_targets_hmm, base_outpath, AA_file, pressed_base=None):
    """
    Search one genome against the combined target-Pfam profiles, in-process via pyhmmer.

    Was a subprocess `hmmsearch --cut_ga` call. The tblout written here is
    field-for-field identical to what the CLI produced, which matters more for Pfam than
    for the SCG search: this file lands in the user-facing
    `pfam-search-results/individual-genome-results/` directory, not just a temp dir.
    """
    os.makedirs(base_outpath, exist_ok=True)

    outpath = f"{base_outpath}/pfam-hmmsearch.txt"

    try:
        search_one_genome(all_pfam_targets_hmm, AA_file, outpath,
                          pressed_base=pressed_base,
                          cpus=search_threads_per_genome(None))
        pfam_search_failed = False
    except Exception as e:
        record_search_failure(base_outpath, "Pfam search",
                              f"{type(e).__name__}: {e}")
        pfam_search_failed = True

    return pfam_search_failed


def get_pfam_counts(pfam_ids, results_txt):

    counts = Counter()

    if not os.path.isfile(results_txt):
        return [0 for _ in pfam_ids]

    with open(results_txt) as results_file:
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
