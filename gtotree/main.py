from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import gtotree_header
from gtotree.utils.runinfo import display_initial_run_info
from gtotree.utils.preprocessing_genomes import preprocess_genomes

def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, run_data = preflight_checks(args)

    display_initial_run_info(args, run_data)

    preprocess_genomes(args, run_data)

    # hmm_searching(args, run_data)

    # optional ko searching

    # optional pfam searching

    # filter_genes_by_length(args, run_data)

    # filter_genomes_with_too_few_hits(args, run_data)

    # generate_alignments(args, run_data)

    # concatenate and rename

    # tree_building(args, run_data)

    # summarize_results(args, run_data)

    # cleanup(args, run_data)
