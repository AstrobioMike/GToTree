from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import gtotree_header
from gtotree.utils.runinfo import display_initial_run_info

def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, input_genome_data, tools_used = preflight_checks(args)

    display_initial_run_info(args, input_genome_data)
