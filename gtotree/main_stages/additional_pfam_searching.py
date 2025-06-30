import subprocess
import shutil
from collections import Counter
import os
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets
from gtotree.utils.messaging import (report_processing_stage,
                                     report_pfam_searching_update)
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)


def search_pfams(args, run_data):

    report_processing_stage("additional-pfam-searching", run_data)

    if run_data.additional_pfam_searching_done:
        report_pfam_searching_update(run_data)
        return run_data

    run_data = get_additional_pfam_targets(run_data)

    if len(run_data.found_pfam_targets) > 0:

        num_genomes_to_search = len(run_data.get_all_input_genomes_for_pfam_search())

        if num_genomes_to_search > 0:
            # writing run_data to file so it can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("search-pfams.smk")
            description = "Searching Pfams"

            run_snakemake(snakefile, num_genomes_to_search, args, run_data, description)

            run_data = read_run_data(run_data.run_data_path)

        print("") ; print("Combining Pfam search results...".center(82))
        combine_all_pfam_hits(run_data.found_pfam_targets,
                              run_data.tmp_pfam_results_dir,
                              run_data.pfam_results_dir + "/pfam-hit-seqs")

    run_data.additional_pfam_searching_done = True
    report_pfam_searching_update(run_data)

    return run_data


def run_pfam_search(all_pfam_targets_hmm, base_outpath, AA_file, num_hmm_cpus):

    os.makedirs(base_outpath, exist_ok=True)

    outpath = f"{base_outpath}/pfam-results.txt"
    tmp_path = f"{base_outpath}/pfam-tmp"

    cmd = [
        "hmmsearch",
        "--cut_ga",
        "--cpu", str(num_hmm_cpus),
        "--tblout", outpath,
        all_pfam_targets_hmm,
        AA_file
    ]
    print(cmd)
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
        pfam_search_failed = False
    except:
        pfam_search_failed = True

    return pfam_search_failed


def get_pfam_counts(pfam_ids, results_txt):

    counts = Counter()

    with open(results_txt, "r") as results_file:
        for line in results_file:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 2:
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
        pattern = f"*/{pfam}*.faa"
        paths = list(Path(tmp_pfam_results_area).glob(pattern))

        if paths:
            out_path = os.path.join(out_base, f"{pfam}.faa")

            with open(out_path, "wb") as out_file:
                for path in paths:
                    with open(path, "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)
