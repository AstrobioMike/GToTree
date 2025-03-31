from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)
from gtotree.utils.messaging import report_processing_stage
from gtotree.utils.seqs import check_target_SCGs_have_seqs

def align_and_prepare_gene_sets(args, run_data):

    report_processing_stage("align-and-prepare-gene-sets")


    return run_data
