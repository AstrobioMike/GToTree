from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import (gtotree_header,
                                     display_initial_run_info,
                                     summarize_results,
                                     copy_log_function)
from gtotree.main_stages.preprocessing_genomes import preprocess_genomes
from gtotree.main_stages.additional_pfam_searching import search_pfams
from gtotree.main_stages.additional_ko_searching import search_kos
from gtotree.main_stages.hmm_searching import search_hmms
from gtotree.main_stages.filtering_genes import filter_genes
from gtotree.main_stages.filtering_genomes import filter_genomes
from gtotree.main_stages.aligning_and_preparing_SCG_sets import align_and_prepare_SCG_sets
from gtotree.main_stages.concatenating_SCG_sets import concatenate_SCG_sets
from gtotree.main_stages.updating_headers import update_headers
from gtotree.main_stages.treeing import make_tree
from gtotree.utils.citations import generate_citations_info
from gtotree.utils.summary_info import generate_primary_summary_table

def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, run_data = preflight_checks(args)

    run_data = display_initial_run_info(args, run_data)

    run_data = preprocess_genomes(args, run_data)

    if run_data.target_pfams_file:
        run_data = search_pfams(args, run_data)

    if run_data.target_kos_file:
        run_data = search_kos(args, run_data)

    run_data = search_hmms(args, run_data)

    run_data = filter_genes(args, run_data)

    run_data = filter_genomes(args, run_data)

    run_data = align_and_prepare_SCG_sets(args, run_data)

    run_data = concatenate_SCG_sets(run_data)

    run_data = update_headers(args, run_data)

    run_data = make_tree(args, run_data)

    generate_citations_info(run_data)

    generate_primary_summary_table(args, run_data)

    # cleanup(args, run_data)
        # remove tmp dir unless debug mode (anything else?)

    summarize_results(args, run_data)

    copy_log_function(run_data)

    # print(f"\n\n{run_data}\n\n")
    # print(f"\n\n{args}\n\n")

