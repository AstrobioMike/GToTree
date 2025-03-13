from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import preflight_checks
from gtotree.utils.messaging import gtotree_header

def main(args = None):
    if args is None:
        args = parser().parse_args()

    print(gtotree_header())

    args, tools_used = preflight_checks(args)


