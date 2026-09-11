import shutil
import textwrap
import argparse
from rich_argparse import RichHelpFormatter  # type: ignore

# this is basically a port of what i came up with for bit

class CustomRichHelpFormatter(RichHelpFormatter):

    # this was added so help text appears on same line as subcommand for programs with subcommands
    def __init__(self, prog, indent_increment=2, max_help_position=24, width=None, **kwargs):
        super().__init__(prog, indent_increment=indent_increment,
                         max_help_position=max(max_help_position, 36),
                         width=width, **kwargs)

    def add_argument(self, action):
        super().add_argument(action)
        if hasattr(action, '_get_subactions'):
            self._action_max_length += 2

    def start_section(self, heading):
        if heading == "positional arguments":
            heading = "Available Subcommands"
        elif heading == "options":
            heading = "Optional Parameters"
        super().start_section(heading)
    group_name_formatter = lambda name: "Usage" if name.lower() == "usage" else name


def add_help(group):
    group.add_argument(
        "-h",
        "--help",
        action="help",
        help=wrap_help("Show this help message")
    )


class VersionAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from gtotree.utils.misc.messaging import get_version
        print(f"GToTree v{get_version()}")
        parser.exit()


def add_version_arg(group):
    group.add_argument(
        "-v",
        "--version",
        action=VersionAction,
        help=wrap_help("Show GToTree version")
    )


def wrap_help(text, margin=30):
    term_width = shutil.get_terminal_size((80, 24)).columns
    help_width = max(10, term_width - margin)
    cleaned_text = " ".join(text.split())

    return textwrap.fill(cleaned_text, width=help_width)


def run_subcommand_main(parser, driver, *program_errors):
    """
    The shared body of a `gtt` subcommand's main(): parse, run, translate exceptions.

    Every standalone subcommand wants the same shape: bare invocation prints help,
    Ctrl-C is a polite exit rather than a traceback, the taxonomy and output-dir errors
    are translated identically, and phase stats are reported whatever happened. Only
    the driver and the program's own error class differ, so those are the arguments.

    `program_errors` are the calling program's own exception types (e.g.
    GenSCGHMMsError, TargetSearchError); they're reported the same way as the shared
    ones, and are separate only because each program defines its own.

    Imports are deferred to call time: this module is imported by the argument
    definitions of every subcommand, and the taxonomy layer imports those back.
    """
    import sys

    from gtotree.utils.misc import phase_stats
    from gtotree.utils.misc.general import OutputDirExistsError
    from gtotree.utils.misc.messaging import report_very_early_exit
    from gtotree.utils.taxonomy.tax_select import (AmbiguousTaxon, CrossDomainTaxon,
                                                   TaxonNotFound)
    from gtotree.utils.taxonomy.wanted_ref_tax import WantedRefTaxError

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    try:
        driver(args)
    except KeyboardInterrupt:
        print()
        report_very_early_exit("Interrupted by user.", "yellow")
    except (TaxonNotFound, AmbiguousTaxon, CrossDomainTaxon, WantedRefTaxError) as e:
        report_very_early_exit(str(e))
    except OutputDirExistsError as e:
        report_very_early_exit(str(e), "yellow", leading_newline=False)
    except program_errors as e:
        report_very_early_exit(str(e))
    finally:
        phase_stats.report()
