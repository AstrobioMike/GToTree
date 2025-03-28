from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import gtotree_header
from gtotree.utils.runinfo import display_initial_run_info
from gtotree.utils.preprocessing_genomes import preprocess_genomes
from gtotree.utils.hmms.hmm_searching import search_hmms
from gtotree.utils.filtering_genes import filter_genes
from gtotree.utils.filtering_genomes import filter_genomes


def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, run_data = preflight_checks(args)

    display_initial_run_info(args, run_data)

    run_data = preprocess_genomes(args, run_data)

    run_data = search_hmms(args, run_data)

    run_data = filter_genes(args, run_data)

    run_data = filter_genomes(args, run_data)
    print(f"\n\n{run_data}\n\n")

    # filter_genomes_with_too_few_hits(args, run_data)

    # generate_alignments(args, run_data)

    # concatenate (and rename if needed)

    # tree_building(args, run_data)

    # optional ko searching

    # optional pfam searching

    # summarize_results(args, run_data)

    # cleanup(args, run_data)
