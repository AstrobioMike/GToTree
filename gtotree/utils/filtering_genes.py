from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)
from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.seqs import check_target_SCGs_have_seqs

def filter_genes(args, run_data):

    report_processing_stage("filter-genes")
    cutoff = "{:.0f}".format(run_data.seq_length_cutoff * 100)
    print(f"\n      Keeping those with lengths within {cutoff}% of the median for each gene set.")

    if not run_data.SCG_hits_filtered:
        num_SCGs_remaining = len(run_data.remaining_SCG_targets)
        write_run_data(run_data)
        snakefile = get_snakefile_path("filter-genes.smk")
        description = "Filtering SCG hits"

        run_snakemake(snakefile, num_SCGs_remaining, args, run_data, description)

        run_data = read_run_data(run_data.run_data_path)
        run_data = check_target_SCGs_have_seqs(run_data, "-gene-filtered.fasta")

    return run_data
