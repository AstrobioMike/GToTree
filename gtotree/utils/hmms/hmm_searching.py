from gtotree.utils.messaging import (report_processing_stage)
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)

def search_hmms(args, run_data):

    report_processing_stage("hmm-search")

    num_genomes_to_search = len(run_data.get_all_remaining_input_genomes())

    if num_genomes_to_search > 0:
        # writing run_data to file so it can be accessed by snakemake
        write_run_data(run_data)
        snakefile = get_snakefile_path("search-hmms.smk")
        description = "Searching HMMs"

        run_snakemake(snakefile, num_genomes_to_search, args, run_data, description)

        run_data = read_run_data(run_data.run_data_path)

    return run_data
