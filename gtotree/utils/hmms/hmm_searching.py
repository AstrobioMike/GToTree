import pandas as pd
import subprocess
from gtotree.utils.messaging import (report_processing_stage,
                                     report_early_exit)
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)


def search_hmms(args, run_data):

    report_processing_stage("hmm-search")

    num_genomes_to_search = len(run_data.get_all_input_genome_for_hmm_search())

    if num_genomes_to_search > 0:

        # writing run_data to file so it can be accessed by snakemake
        write_run_data(run_data)
        snakefile = get_snakefile_path("search-hmms.smk")
        description = "Searching HMMs"

        run_snakemake(snakefile, num_genomes_to_search, args, run_data, description)

        run_data = read_run_data(run_data.run_data_path)

    return run_data


def run_hmm_search(id, run_data, inpath):
    outpath = f"{run_data.hmm_results_dir}/{id}-hmm-hits.txt"
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

    dict_of_hit_counts = dict.fromkeys(run_data.initial_SCG_targets, 0)
    for scg in run_data.initial_SCG_targets:
        dict_of_hit_counts[scg] = df[df["target_SCG"] == scg].shape[0]

    if run_data.best_hit_mode:
        report_early_exit("Best hit mode not yet implemented!")
    else:
        pass

    dict_of_hit_gene_ids = dict.fromkeys(run_data.initial_SCG_targets, None)

    for scg, count in dict_of_hit_counts.items():
        if count == 1:
            gene_id = df.loc[df["target_SCG"] == scg, "gene_id"].iloc[0]
            dict_of_hit_gene_ids[scg] = gene_id

    return dict_of_hit_counts, dict_of_hit_gene_ids
