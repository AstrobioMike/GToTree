from gtotree.cli.parser import parser
from gtotree.utils.preflight_checks import perform_preflight_checks

def main(args = None):
    if args is None:
        args = parser().parse_args()
    print(args)
    tools_used = perform_preflight_checks(args)
    print(tools_used)
