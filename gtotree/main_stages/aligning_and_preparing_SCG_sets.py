from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)
from gtotree.utils.messaging import (report_processing_stage,
                                     report_no_SCGs_remaining)
from gtotree.utils.seqs import copy_gene_alignments

def align_and_prepare_SCG_sets(args, run_data):

    report_processing_stage("align-and-prepare-SCG-sets", run_data)

    if not run_data.all_SCG_sets_aligned:

        num_SCGs_to_align = len(run_data.get_all_SCG_targets_remaining())

        if num_SCGs_to_align > 0:
            write_run_data(run_data)
            snakefile = get_snakefile_path("align-and-prepare-SCG-sets.smk")
            description = "Aligning and preparing SCG sets"

            run_snakemake(snakefile, num_SCGs_to_align, args, run_data, description)

            run_data = read_run_data(run_data.run_data_path)

        else:
            report_no_SCGs_remaining(run_data)

        if args.keep_gene_alignments:
            copy_gene_alignments(run_data)

    return run_data
