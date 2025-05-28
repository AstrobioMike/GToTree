from gtotree.utils.pfam.pfam_handling import get_additional_pfam_targets
from gtotree.utils.messaging import (report_processing_stage,
                                     report_pfam_searching_update)
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)

def search_pfams(run_data):

    report_processing_stage("additional-pfam-searching", run_data)

    if run_data.additional_pfam_searching_done:
        report_pfam_searching_update(run_data)
        return run_data

    run_data = get_additional_pfam_targets(run_data)

    if len(run_data.found_pfam_targets) > 0:

        num_genomes_to_search = len(run_data.get_all_input_genomes_for_pfam_search())

        if num_genomes_to_search > 0:
            # writing run_data to file so it can be accessed by snakemake
            write_run_data()
            snakefile = get_snakefile_path("search-pfams.smk")
            description = "Searching Pfams"

            run_snakemake(snakefile, num_genomes_to_search, description)

            run_data = read_run_data()

        print("") ; print("Combining Pfam search results...".center(82))
        combine_all_pfam_hits(run_data.found_pfam_targets,
                              run_data.tmp_pfam_results_dir,
                              run_data.pfam_results_dir + "/pfam-hit-seqs")

    run_data.additional_pfam_searching_done = True
    report_pfam_searching_update(run_data)

    return run_data


def combine_all_pfam_hits(pfam_ids, tmp_pfam_results_dir, out_base):
    pass