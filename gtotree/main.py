import sys
from gtotree.cli.helpmenu import helpmenu
from gtotree.cli.parser import parser
from gtotree.utils.misc.preflight_checks import preflight_checks
from gtotree.utils.misc.messaging import (gtotree_header,
                                     display_initial_run_info,
                                     summarize_results,
                                     copy_log_function)
from gtotree.main_stages.processing_genomes import process_genomes
from gtotree.main_stages.filtering_genes import filter_genes
from gtotree.main_stages.filtering_genomes import filter_genomes
from gtotree.main_stages.aligning_and_preparing_SCG_sets import align_and_prepare_SCG_sets
from gtotree.main_stages.concatenating_SCG_sets import concatenate_SCG_sets
from gtotree.main_stages.updating_headers import update_headers
from gtotree.main_stages.treeing import make_tree
from gtotree.utils.misc.citations import generate_citations_info
from gtotree.utils.misc.summary_info import generate_primary_summary_table
from gtotree.utils.misc.itol import generate_all_search_itol_files
from gtotree.utils.misc.general import cleanup, write_run_data
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.stages import PipelineStage

def main(args = None):
    if args is None:
        if len(sys.argv) == 1:
            sys.stdout.write(helpmenu)
            sys.exit()

        args = parser().parse_args()

    print(gtotree_header())

    run_data = None

    # phase_stats only runs if GTT_DEBUG_TIMING=1 was set as an env var (here for testing/dev/troubleshooting)
    phase_stats.begin("preflight checks")

    try:
        args, run_data = preflight_checks(args)

        phase_stats.begin("initial run info")
        run_data = display_initial_run_info(args, run_data)

        # process_genomes opens its own per-source sub-phases
        run_data = process_genomes(args, run_data)

        phase_stats.begin("filtering genes")
        run_data = filter_genes(args, run_data)

        phase_stats.begin("filtering genomes")
        run_data = filter_genomes(args, run_data)

        phase_stats.begin("aligning and preparing SCG sets")
        run_data = align_and_prepare_SCG_sets(args, run_data)

        phase_stats.begin("concatenating SCG sets")
        run_data = concatenate_SCG_sets(run_data)

        phase_stats.begin("updating headers")
        run_data = update_headers(args, run_data)

        phase_stats.begin("treeing")
        run_data = make_tree(args, run_data)

        phase_stats.begin("summaries and citations")
        generate_citations_info(run_data)

        generate_primary_summary_table(args, run_data)

        generate_all_search_itol_files(args, run_data)

        cleanup(args, run_data)

        run_data.mark_stage_complete(PipelineStage.FINALIZE)

        write_run_data(run_data)

        summarize_results(args, run_data)

        copy_log_function(run_data)

    finally:
        phase_stats.finish()
        if run_data is not None:
            phase_stats.write_tsv(run_data.run_files_dir)
        phase_stats.report()
