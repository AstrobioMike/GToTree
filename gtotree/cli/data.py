import sys
import argparse
import importlib
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg


DATA_SOURCES = {
    "gtdb-data": {
        "help": "Download or update GTDB metadata",
        "description": "This subcommand downloads or updates the GTDB metadata.",
        "module": "gtotree.utils.gtdb.get_gtdb_data",
        "function": "get_gtdb_data",
    },
    "ncbi-assembly-data": {
        "help": "Download or update NCBI assembly-summary tables",
        "description": "This subcommand downloads or updates NCBI's assembly-summary tables.",
        "module": "gtotree.utils.ncbi.get_ncbi_assembly_data",
        "function": "get_ncbi_assembly_data",
    },
    "kofamscan-data": {
        "help": "Download or update KOFamScan data",
        "description": ("This subcommand downloads or updates the KOFamScan data "
                        "(KO profiles and list)."),
        "module": "gtotree.utils.ko.get_kofamscan_data",
        "function": "get_kofamscan_data",
    },
    "pfam-data": {
        "help": "Download or update Pfam data",
        "description": ("This subcommand downloads or updates the Pfam data "
                        "(Pfam-A HMMs and info)."),
        "module": "gtotree.utils.pfam.get_pfam_data",
        "function": "get_pfam_data",
    },
}


def build_parser(parent_subparsers=None):

    desc = ("This program manages GToTree-utilized databases and their location "
            "settings. See subcommand-specific help menus for more info.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "data",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    add_help(parser)
    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ############################################################################
    ### get subcommand ###
    ############################################################################

    get_desc = """
        This subcommand downloads or updates GToTree-utilized databases. All sub-subcommands
        accept an optional `-f/--force-update` flag.
        """

    get_parser = subparsers.add_parser(
        "get",
        help="Download or update a GToTree-utilized database",
        description=get_desc,
        epilog="Ex. usage: `gtt data get gtdb-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    add_help(get_parser)

    add_version_arg(get_parser)

    get_subparsers = get_parser.add_subparsers(dest="get_action", required=True, metavar='')
    get_parser.subparsers = get_subparsers

    def add_get_common_args(group):
        group.add_argument(
            "-f", "--force-update",
            help="Re-download the data even if it is already present",
            action="store_true",
        )

    for source_name, source_info in DATA_SOURCES.items():

        source_parser = get_subparsers.add_parser(
            source_name,
            help=source_info["help"],
            description=source_info["description"],
            epilog=f"Ex. usage: `gtt data get {source_name}`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

        source_optional = source_parser.add_argument_group("Optional Parameters")
        add_get_common_args(source_optional)
        add_help(source_optional)

        add_version_arg(source_optional)

        source_parser.set_defaults(func="get")

    ############################################################################
    ### locations subcommand ###
    ############################################################################

    locations_desc = ("This subcommand checks or sets the data-location environment "
                      "variables required by GToTree.")

    locations_parser = subparsers.add_parser(
        "locations",
        help="Check or set data-location environment variables",
        description=locations_desc,
        epilog="Ex. usage: `gtt data locations check`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    add_help(locations_parser)
    add_version_arg(locations_parser)

    locations_subparsers = locations_parser.add_subparsers(
        dest="locations_action", required=True, metavar='')
    locations_parser.subparsers = locations_subparsers

    ### locations check ###

    locations_check_desc = ("This subcommand reports the current data-location environment "
                            "variables and whether the paths are writable.")

    locations_check_parser = locations_subparsers.add_parser(
        "check",
        help="Report current data-location environment variables",
        description=locations_check_desc,
        epilog="Ex. usage: `gtt data locations check`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    locations_check_optional = locations_check_parser.add_argument_group("Optional Parameters")
    add_help(locations_check_optional)
    add_version_arg(locations_check_optional)
    locations_check_parser.set_defaults(func="locations_check")

    ### locations set ###

    locations_set_desc = ("This subcommand interactively sets the data-location environment "
                          "variables.")

    locations_set_parser = locations_subparsers.add_parser(
        "set",
        help="Interactively set data-location environment variables",
        description=locations_set_desc,
        epilog="Ex. usage: `gtt data locations set`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    locations_set_optional = locations_set_parser.add_argument_group("Optional Parameters")
    add_help(locations_set_optional)
    add_version_arg(locations_set_optional)
    locations_set_parser.set_defaults(func="locations_set")

    return parser


################################################################################


def _run_get(args):
    """Dispatch `gtt data get <source>` to the matching worker (lazy import)."""
    source_info = DATA_SOURCES[args.get_action]
    module = importlib.import_module(source_info["module"])
    getattr(module, source_info["function"])(force_update=args.force_update)


def _run_locations_check(args):
    from gtotree.utils.misc.data_locations import check_and_report_env_variables
    check_and_report_env_variables()


def _run_locations_set(args):
    from gtotree.utils.misc.data_locations import (set_env_variables,
                                              modify_conda_activate_startup_script,
                                              notify_to_reactivate_conda)
    paths_dict = set_env_variables()
    modify_conda_activate_startup_script(paths_dict)
    notify_to_reactivate_conda()


def main():

    parser = build_parser()

    try:
        import argcomplete  # type: ignore
        argcomplete.autocomplete(parser)
    except ImportError:
        pass

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args for a top-level subcommand so the appropriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
            sys.exit(0)

        if cmd in ("-v", "--version"):
            from gtotree.utils.misc.messaging import get_version
            print(f"GToTree v{get_version()}")
            sys.exit(0)

        if cmd in parser.subparsers.choices:
            parser.subparsers.choices[cmd].print_help(sys.stderr)
            sys.exit(0)
        else:
            print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)

    # handling no args for a second-level subcommand (e.g. `gtt data locations`)
    # so the appropriate help menu is printed
    if len(sys.argv) == 3:
        cmd = sys.argv[1]
        sub_cmd = sys.argv[2]

        if cmd in parser.subparsers.choices:
            top_sub_parser = parser.subparsers.choices[cmd]

            if hasattr(top_sub_parser, 'subparsers'):

                if sub_cmd in ("-h", "--help"):
                    top_sub_parser.print_help(sys.stderr)
                    sys.exit(0)

                if sub_cmd in ("-v", "--version"):
                    from gtotree.utils.misc.messaging import get_version
                    print(f"GToTree v{get_version()}")
                    sys.exit(0)

                if sub_cmd not in top_sub_parser.subparsers.choices:
                    print(f"\n  Invalid subcommand provided: '{sub_cmd}'\n\n  See help below.\n", file=sys.stderr)
                    top_sub_parser.print_help(sys.stderr)
                    sys.exit(1)
                # else: command is complete as-is (e.g. `gtt data locations check`),
                # fall through to parse_args()

    args = parser.parse_args()

    func_map = {
        "get":             _run_get,
        "locations_check": _run_locations_check,
        "locations_set":   _run_locations_set,
    }

    func_map[args.func](args)
