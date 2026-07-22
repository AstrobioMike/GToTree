import sys
import argparse
import importlib
from gtotree.cli.common import add_help, add_version_arg, CustomRichHelpFormatter


SUBCOMMAND_MAP = {
    "get-accs-from-gtdb":  "gtotree.utils.gtdb.get_accessions_from_gtdb",
    "get-accs-from-ncbi":  "gtotree.utils.ncbi.get_accessions_from_ncbi",
    "hmms":                "gtotree.utils.hmms.gtt_hmms",
    "midpoint-root-tree":  "gtotree.utils.helper_scripts.gtt_midpoint_root_tree",
    "data":                "gtotree.cli.data",
    "test":                "gtotree.tests.smoke",
}


PROGRAM_GROUPS = [
    {
        "title": "NCBI/GTDB-related",
        "programs": [
            {
                "name": "get-accs-from-gtdb",
                "desc": "search GTDB by taxonomy and retrieve NCBI accessions",
            },
            {
                "name": "get-accs-from-ncbi",
                "desc": "search NCBI by taxonomy or taxid and retrieve NCBI accessions",
            },
        ],
    },
    {
        "title": "SCG-HMM-related",
        "programs": [
            {
                "name": "hmms",
                "desc": "view the available pre-packaged HMM SCG-sets",
            },
        ],
    },
    {
        "title": "Misc",
        "programs": [
            {
                "name": "midpoint-root-tree",
                "desc": "midpoint-root a newick tree",
            },
        ],
    },
    {
        "title": "GToTree-data management",
        "programs": [
            {
                "name": "data",
                "desc": "",
                "subcommands": [
                    ("get",       "download/update a GToTree-utilized database"),
                    ("locations", "check or set data-location environment variables"),
                ],
            },
        ],
    },
    {
        "title": "Testing",
        "programs": [
            {
                "name": "test",
                "desc": "run a quick end-to-end test of the installed environment",
            },
        ],
    },
]


def print_overview():

    from rich.console import Console  # type: ignore
    from rich.table import Table  # type: ignore
    from gtotree.utils.messaging import get_version
    from datetime import datetime

    console = Console(highlight=False)
    ver = f"v{get_version()}"

    console.print()
    console.print(f"{'':>32}GToTree [green]{ver}[/green]")
    console.print(f"{'':>24}github.com/AstrobioMike/GToTree")
    console.print()
    console.print(f"{'':>28}OVERVIEW OF SUBCOMMANDS")
    console.print()

    name_col_width = max(
        max(
            max(len(p["name"]) for p in g["programs"]),
            max((len(s[0]) + 2 for p in g["programs"] for s in p.get("subcommands", [])), default=0),
        )
        for g in PROGRAM_GROUPS
    ) + 2

    console.rule(style="dim")

    for i, group in enumerate(PROGRAM_GROUPS):
        if i > 0:
            console.print()
            console.rule(style="dim")
        console.print()

        console.print(f"  [dark_orange]{group['title']}[/dark_orange]")

        table = Table(box=None, show_header=False, padding=(0, 2, 0, 4))
        table.add_column(no_wrap=True, min_width=name_col_width)
        table.add_column()

        for prog in group["programs"]:
            name = prog["name"]
            desc = prog["desc"]
            subcommands = prog.get("subcommands")

            table.add_row(
                f"[cyan]{name}[/cyan]",
                desc,
            )

            for sub_name, sub_desc in (subcommands or []):
                table.add_row(
                    f"  [dark_cyan]{sub_name}[/dark_cyan]",
                    sub_desc,
                )

        console.print(table)

    console.print()
    console.rule(style="dim")
    console.print()

    today = datetime.today().strftime('%A')
    signoff = f"Happy {today} :)"
    console.print(f"{'':>49}[green]{signoff}[/green]\n")


def _suppress_help_version_on_group_parsers(parser):
    """
    Suppress -h/--help/-v/--version from argcomplete on parsers that have
    subparsers, so TAB after a group shows only subcommand names.
    Leaf-level parsers are left untouched so their flags still appear
    without requiring a '-' prefix.
    """
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for a in parser._actions:
                if hasattr(a, 'option_strings') and any(
                    s in ('-h', '--help', '-v', '--version')
                    for s in a.option_strings
                ):
                    a.help = argparse.SUPPRESS
            for sub_parser in action.choices.values():
                _suppress_help_version_on_group_parsers(sub_parser)
            break


def build_parser():

    desc = """
        Print an overview of all available gtt subcommands.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    add_help(parser)
    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest='subcommand')

    # Lazy-load modules during tab completion so only the one relevant module
    # is imported per keypress instead of all of them at once.
    import os
    comp_line = os.environ.get('COMP_LINE', '')
    if comp_line:
        words = comp_line.split()
        if len(words) >= 2 and words[1] in SUBCOMMAND_MAP:
            # User is completing args/flags for a specific subcommand -- only
            # import that one module.
            module = importlib.import_module(SUBCOMMAND_MAP[words[1]])
            module.build_parser(parent_subparsers=subparsers)
        else:
            # User is still completing the subcommand name itself -- lightweight
            # stubs are all argcomplete needs to offer the names.
            for name in SUBCOMMAND_MAP:
                subparsers.add_parser(name, add_help=False)
    else:
        # Normal (non-completion) invocation -- build the full tree.
        for module_path in SUBCOMMAND_MAP.values():
            module = importlib.import_module(module_path)
            module.build_parser(parent_subparsers=subparsers)

    _suppress_help_version_on_group_parsers(parser)

    return parser


def main():

    parser = build_parser()

    try:
        import argcomplete  # type: ignore
        argcomplete.autocomplete(parser)
    except ImportError:
        pass

    # print the overview when called with no arguments or with -h/--help
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print_overview()
        sys.exit(0)

    subcommand = sys.argv[1]

    # handle -v/--version before subcommand dispatch
    if subcommand in ("-v", "--version"):
        parser.parse_args()
        return

    if subcommand not in SUBCOMMAND_MAP:
        from rich.console import Console  # type: ignore
        Console(stderr=True).print(f"\n    [yellow]Unknown subcommand:[/yellow] [cyan]{subcommand}[/cyan]\n")
        Console(stderr=True).print(f"  Run [cyan]gtt[/cyan] by itself to see available subcommands.\n")
        sys.exit(1)

    # rewrite argv so the target module's parser sees the right program name
    sys.argv = [f"gtt {subcommand}"] + sys.argv[2:]
    module = importlib.import_module(SUBCOMMAND_MAP[subcommand])
    module.main()
