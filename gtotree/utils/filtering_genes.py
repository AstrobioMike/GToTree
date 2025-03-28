from gtotree.utils.messaging import (report_processing_stage,
                                     report_early_exit)

def filter_genes(args, run_data):

    report_processing_stage("filter-genes")
    cutoff = "{:.0f}".format(args.seq_length_cutoff * 100)
    print(f"\n      Keeping those with lengths within {cutoff}% of the median for each gene set.\n")

    # for each target gene-set
        #

    return run_data