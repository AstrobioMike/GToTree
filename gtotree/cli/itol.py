import sys
import argparse

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc.itol import (COLOR_MAP, SHAPE_MAP, SOURCE_CHOICES,
                                     ItolError)


COLOR_CHOICES = list(COLOR_MAP)
SHAPE_CHOICES = list(SHAPE_MAP)
SOURCE_OPTIONS = list(SOURCE_CHOICES)


################################################################################
# parser
################################################################################

def build_parser(parent_subparsers=None):

    desc = ("This program generates Interactive Tree of Life (iToL) annotation files "
            "that can be dropped onto a tree at itol.embl.de to decorate a set of "
            "target genomes. See itol.embl.de/help.cgi for more info. All subcommands "
            "take a single-column file of genome IDs that match the labels in the tree,"
            "or an input specification if including a prior gtotree-output dir.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "itol",
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

    ### shared args ###

    def add_common_selection(group):
        group.add_argument(
            "-g", "--wanted-genomes",
            metavar="<FILE>",
            help=("Single-column file of the genomes to decorate. Without `-i` these "
                  "must match the labels in the tree file; with `-i` they can be "
                  "genome IDs, labels, or the original inputs (e.g., accessions or "
                  "file paths)"),
        )

        group.add_argument(
            "-i", "--input-dir",
            metavar="<DIR>",
            help=("The output directory of a completed GToTree run. Lets `-g` take "
                  "any identifier for a genome, catches ones that match nothing, and "
                  "enables `--source`"),
        )

        group.add_argument(
            "--source",
            metavar="<STR>",
            action="append",
            choices=SOURCE_OPTIONS,
            help=("Instead of using `-g` as input, you can select genomes based on their "
                  "input source IF you are providing `-i`. "
                  "One of: " + ", ".join(SOURCE_OPTIONS)),
        )

    def add_common_optional(group):
        group.add_argument(
            "-o", "--output-file",
            metavar="<FILE>",
            help='Name of the output iToL file (default: "itol.txt")',
            default="itol.txt",
        )

        group.add_argument(
            "-c", "--color",
            metavar="<STR>",
            help=(f'Color to use. Either a name ({", ".join(COLOR_CHOICES)}) or a '
                  f'hex code like "#0000ff" (default: "blue")'),
            default="blue",
        )

    # `DATASET_LABEL` in the written file
    def add_dataset_label(group):
        group.add_argument(
            "-d", "--dataset-label",
            metavar="<STR>",
            help='Label of the dataset (default: "data")',
            default="data",
        )

    ############################################################################
    ### style subcommand ###
    ############################################################################

    style_desc = ("This subcommand creates an iToL style-dataset file "
                  "for coloring branches and/or labels.")

    style_parser = subparsers.add_parser(
        "style",
        help="Create an iToL style-dataset file for coloring branches/labels",
        description=style_desc,
        epilog="Ex. usage: `gtt itol style -g genomes.txt -d my-genomes -o my-genomes-itol.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    style_selection = style_parser.add_argument_group("Genome Selection")
    style_optional = style_parser.add_argument_group("Optional Parameters")

    add_common_selection(style_selection)
    add_common_optional(style_optional)

    add_dataset_label(style_optional)

    style_optional.add_argument(
        "--what-to-color",
        help='What to color (default: "branches")',
        choices=["branches", "labels", "both"],
        default="branches",
    )


    style_optional.add_argument(
        "-l", "--line-weight",
        metavar="<NUM>",
        help='Line weight when coloring branches (default: "3")',
        default="3",
    )

    add_help(style_optional)
    add_version_arg(style_optional)

    style_parser.set_defaults(func="style")

    ############################################################################
    ### binary-dataset subcommand ###
    ############################################################################

    binary_desc = ("This subcommand creates an iToL binary-dataset file, "
                   "which puts a filled symbol next to each target genome.")

    binary_parser = subparsers.add_parser(
        "binary-dataset",
        help="Create an iToL binary-dataset file",
        description=binary_desc,
        epilog="Ex. usage: `gtt itol binary-dataset -g genomes.txt -s star`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    binary_selection = binary_parser.add_argument_group("Genome Selection")
    binary_optional = binary_parser.add_argument_group("Optional Parameters")

    add_common_selection(binary_selection)
    add_common_optional(binary_optional)

    add_dataset_label(binary_optional)

    binary_optional.add_argument(
        "-s", "--shape",
        help='Shape to add (default: "square")',
        choices=SHAPE_CHOICES,
        default="square",
    )

    binary_optional.add_argument(
        "-H", "--height-factor",
        metavar="<NUM>",
        dest="height",
        help=('Increase or decrease symbol size. Values below 1 decrease the standard '
              'size, above 1 increase it (default: "1")'),
        default="1",
    )

    add_help(binary_optional)
    add_version_arg(binary_optional)

    binary_parser.set_defaults(func="binary-dataset")

    ############################################################################
    ### colorstrip subcommand ###
    ############################################################################

    colorstrip_desc = ("This subcommand creates an iToL colorstrip file, "
                       "which draws a colored band alongside each target genome.")

    colorstrip_parser = subparsers.add_parser(
        "colorstrip",
        help="Create an iToL colorstrip file",
        description=colorstrip_desc,
        epilog="Ex. usage: `gtt itol colorstrip -g genomes.txt -d 'my group'`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    colorstrip_selection = colorstrip_parser.add_argument_group("Genome Selection")
    colorstrip_optional = colorstrip_parser.add_argument_group("Optional Parameters")

    add_common_selection(colorstrip_selection)
    add_common_optional(colorstrip_optional)

    add_dataset_label(colorstrip_optional)

    colorstrip_optional.add_argument(
        "-W", "--strip-width",
        metavar="<NUM>",
        help='Width of the colorstrip (default: "25")',
        default="25",
    )

    colorstrip_optional.add_argument(
        "--color-branches-too",
        help="Add this flag if wanting to color branches also",
        action="store_true",
    )

    add_help(colorstrip_optional)
    add_version_arg(colorstrip_optional)

    colorstrip_parser.set_defaults(func="colorstrip")

    ############################################################################
    ### text-dataset subcommand ###
    ############################################################################

    text_desc = ("This subcommand creates an iToL text-dataset file, "
                 "which places a text label alongside each target genome.")

    text_parser = subparsers.add_parser(
        "text-dataset",
        help="Create an iToL text-dataset file",
        description=text_desc,
        epilog="Ex. usage: `gtt itol text-dataset -g genomes.txt -t 'from site A'`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    text_selection = text_parser.add_argument_group("Genome Selection")
    text_optional = text_parser.add_argument_group("Optional Parameters")

    add_common_selection(text_selection)

    text_selection.add_argument(
        "-t", "--text-to-add",
        metavar="<STR>",
        help="Text to add to the target genomes",
        required=True,
    )

    add_common_optional(text_optional)

    add_dataset_label(text_optional)

    add_help(text_optional)
    add_version_arg(text_optional)

    text_parser.set_defaults(func="text-dataset")

    return parser


################################################################################
# argument validation
################################################################################

def _as_number(value, flag, kind=float):
    """
    Coerce, then drop a pointless trailing ".0" so a default of 3 is written as
    `3` rather than `3.0` -- matching the files generated automatically after a
    Pfam/KO search, which are otherwise byte-identical.
    """
    try:
        number = kind(value)
        if kind is float and number.is_integer():
            return int(number)
        return number
    except (TypeError, ValueError):
        what = "an integer" if kind is int else "a number"
        raise ItolError(f'The value passed to `{flag}` must be {what}, but "{value}" '
                        "was given.")


def check_args(args):
    """
    Coerce the numeric options and settle how genomes are being selected.

    Everything else argparse already constrained.
    """
    if not args.wanted_genomes and not args.source:
        raise ItolError(
            "We need to know which genomes to decorate. Give a list with `-g`, or "
            "point `-i` at a finished GToTree run and select with `--source`.")

    if args.source and not args.input_dir:
        raise ItolError(
            "`--source` reads which genomes came from where out of a finished GToTree run's "
            "summary table, so it needs `-i` pointed at that run's output directory.")

    if args.func == "binary-dataset":
        args.height = _as_number(args.height, "--height-factor")

    elif args.func == "colorstrip":
        args.strip_width = _as_number(args.strip_width, "--strip-width", kind=int)

    elif args.func == "style":
        args.line_weight = _as_number(args.line_weight, "--line-weight")

    return args


################################################################################
# dispatch
################################################################################

def select_targets(args):
    """
    Work out which tree labels to annotate.

    Three shapes, in increasing amounts of help from the run directory:
      * `-g` alone    : the list is taken as tree labels verbatim (what `bit itol`
                        does, and what works when there's no run dir to consult)
      * `-i -g`       : the list may use any identifier and is translated to labels,
                        with anything that matches nothing raised as an error
      * `-i --source` : no list at all; the run's summary says which came from where
    """
    from gtotree.utils.misc.itol import (read_targets, resolve_wanted_genomes,
                                         select_by_source)

    if not args.input_dir:
        return read_targets(args.wanted_genomes), []

    if args.source:
        selection = select_by_source(args.input_dir, args.source)

        if args.wanted_genomes:
            # both given: `--source` sets the pool, `-g` narrows it
            wanted = resolve_wanted_genomes(args.input_dir,
                                            read_targets(args.wanted_genomes))
            keep = set(selection.labels) & set(wanted.labels)
            if not keep:
                raise ItolError(
                    "No genomes matched both `--source` and the list given to `-g`.")
            ordered = [l for l in wanted.labels if l in keep]
            return ordered, wanted.dropped_not_in_tree

        return selection.labels, selection.dropped_not_in_tree

    selection = resolve_wanted_genomes(args.input_dir,
                                       read_targets(args.wanted_genomes))
    return selection.labels, selection.dropped_not_in_tree


def run_itol(args):
    """
    Read the target list once, then hand off to the matching writer.
    """
    from gtotree.utils.misc.itol import (resolve_color, write_binary_dataset,
                                         write_colorstrip, write_style_dataset,
                                         write_text_dataset)

    args = check_args(args)

    targets, dropped = select_targets(args)
    color = resolve_color(args.color)

    if args.func == "style":
        write_style_dataset(args.output_file, targets, label=args.dataset_label,
                            color=color, what_to_color=args.what_to_color,
                            line_weight=args.line_weight)

    elif args.func == "binary-dataset":
        write_binary_dataset(args.output_file, targets, label=args.dataset_label,
                             color=color, shape=args.shape,
                             height_factor=args.height)

    elif args.func == "colorstrip":
        write_colorstrip(args.output_file, targets, label=args.dataset_label,
                         color=color, strip_width=args.strip_width,
                         color_branches_too=args.color_branches_too)

    elif args.func == "text-dataset":
        write_text_dataset(args.output_file, targets, args.text_to_add,
                           label=args.dataset_label, color=color)

    return targets, dropped


def report_finish(args, targets, dropped):
    from gtotree.utils.misc.messaging import color_text, report_message

    print(f"\n    Wrote annotations for {color_text(f'{len(targets):,}', 'green')} "
          f"genome(s) to:")
    print(f"      {color_text(args.output_file, 'green')}\n")

    if dropped:
        report_message(
            f"{len(dropped)} selected genome(s) aren't in the final tree, so they "
            "were left out. See the `reason_removed` column of genomes-summary-info.tsv for why.",
            "yellow", ii="    ", si="    ")


def main():

    parser = build_parser()

    try:
        import argcomplete  # type: ignore
        argcomplete.autocomplete(parser)
    except ImportError:
        pass

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args for a subcommand so the appropriate help menu is printed
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

        print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n",
              file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    from gtotree.utils.misc.messaging import report_very_early_exit

    try:
        targets, dropped = run_itol(args)
    except ItolError as e:
        report_very_early_exit(str(e))
    else:
        report_finish(args, targets, dropped)


if __name__ == "__main__":
    main()
