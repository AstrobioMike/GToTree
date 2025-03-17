from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import gtotree_header
from gtotree.utils.general import populate_input_genome_data
from gtotree.utils.runinfo import display_initial_run_info

def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, tools_used = preflight_checks(args)

    input_genome_data = populate_input_genome_data(args)

    display_initial_run_info(args, input_genome_data)
