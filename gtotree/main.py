from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import gtotree_header
from gtotree.utils.runinfo import display_initial_run_info
from gtotree.utils.processing_genomes import process_genomes

def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, genome_data, tools_used = preflight_checks(args)

    display_initial_run_info(args, genome_data)

    process_genomes(args, genome_data)

    # filter_genes_by_length(args, genome_data)

    # filter_genomes_with_too_few_hits(args, genome_data)

    # generate_alignments(args, genome_data)

    # tree_building(args, genome_data)

    # summarize_results(args, genome_data)

    # cleanup(args, genome_data)
