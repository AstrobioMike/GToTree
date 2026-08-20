"""
re-label a finished GToTree run without re-running it

This takes the finished output directory and writes a fresh tree, alignment,
and summary table carrying the new labels, leaving the original run untouched.

I added this because folks occasionally want to change their labels or add taxonomy,
and previously that would require an entire new run
"""

import os
import sys
import shutil
import argparse

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.general import OutputDirExistsError
from gtotree.utils.misc.messaging import (report_message, color_text,
                                          report_phase_header, report_very_early_exit,
                                          spinner)
from gtotree.utils.update_headers import update_headers_run as run_lib
from gtotree.utils.update_headers.update_headers_run import (
    UpdateHeadersError,
    ALIGNMENT_BASENAME,
    LABEL_MAP_FILENAME,
    SUMMARY_FILENAME,
    output_name,
)


DEFAULT_OUTPUT_DIR = "gtotree-updated-headers"
DEFAULT_LINEAGE = "domain,phylum,class,genus,species"


################################################################################
# parser
################################################################################

def build_parser(parent_subparsers=None):

    desc = ("This is a helper program that takes the output directory of a completed "
            "GToTree run and writes out a new tree, alignment, and genomes-summary "
            "table carrying updated genome headers/labels. The original output directory isn't modified.")

    example = ("Ex. usage: `gtt update-headers -i gtotree-output --add-gtdb-tax "
               "-L domain,phylum,family,species`")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "update-headers",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog=example,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters")
    labels = parser.add_argument_group("New Labels (at least one, otherwise original IDs come back)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i", "--input-dir",
        metavar="<DIR>",
        dest="input_dir",
        help="The output directory of a completed GToTree run",
        action="store",
    )

    labels.add_argument(
        "-m", "--mapping-file",
        metavar="<FILE>",
        help=("A two- or three-column tab-delimited file specifying the wanted labels, "
              "exactly as the main GToTree program's `-m` takes it"),
        action="store",
    )

    labels.add_argument(
        "-D", "--add-gtdb-tax",
        help=("Add GTDB taxonomic lineage to the labels of any genomes that came from "
              "NCBI accessions. See `--lineage-ranks`."),
        action="store_true",
    )

    labels.add_argument(
        "-t", "--add-ncbi-tax",
        help=("Add NCBI taxonomic lineage to the labels of any genomes that came from "
              "NCBI accessions. See `--lineage-ranks`."),
        action="store_true",
    )

    labels.add_argument(
        "-L", "--lineage-ranks",
        metavar="<STR>",
        dest="lineage",
        default=DEFAULT_LINEAGE,
        help=(f'Comma-delimited ranks to include when adding lineage info (default: '
              f'"{DEFAULT_LINEAGE}")'),
        action="store",
    )

    optional.add_argument(
        "-o", "--output-dir",
        metavar="<DIR>",
        dest="output_dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f'Desired output directory (default: "{DEFAULT_OUTPUT_DIR}")',
        action="store",
    )

    optional.add_argument(
        "-p", "--output-prefix",
        metavar="<STR>",
        dest="output_prefix",
        default="",
        help=("Prefix for the output file names"),
        action="store",
    )

    optional.add_argument(
        "-F", "--force-overwrite",
        action="store_true",
        help="Overwrite the output directory if it already exists",
    )

    add_help(optional)
    add_version_arg(optional)

    parser.set_defaults(func="update_headers")

    return parser


################################################################################
# argument validation
################################################################################

def check_args(args):
    """
    Everything that can be settled without reading the run
    """
    if not args.input_dir:
        raise UpdateHeadersError(
            "We need to know which finished GToTree run to re-label! Point `-i` at its "
            "output directory.")

    if args.add_gtdb_tax and args.add_ncbi_tax:
        raise UpdateHeadersError(
            "You've asked to add both GTDB and NCBI taxonomy to the labels. Please "
            "choose one or the other.")

    accepted_ranks = ["domain", "phylum", "class", "order", "family", "genus",
                      "species", "strain"]
    for rank in args.lineage.split(","):
        if rank.strip().lower() not in accepted_ranks:
            raise UpdateHeadersError(
                f'You specified "{args.lineage}" to `-L`, but "{rank}" is not an '
                "accepted taxonomic rank. Accepted ranks are any combination of: "
                f"{', '.join(accepted_ranks)}")

    if args.lineage != DEFAULT_LINEAGE and not (args.add_gtdb_tax or args.add_ncbi_tax):
        raise UpdateHeadersError(
            "You've specified a custom lineage (`-L`), but neither `-D`/`--add-gtdb-tax` "
            "nor `-t`/`--add-ncbi-tax` was given to say which taxonomy to use.")

    if os.path.abspath(args.output_dir.rstrip("/")) == \
            os.path.abspath(args.input_dir.rstrip("/")):
        raise UpdateHeadersError(
            "The output directory (`-o`) is the same as the run being re-labeled "
            "(`-i`). This program never writes into the run it reads, so please point "
            "`-o` somewhere else.")

    return args


def ensure_reference_data(args):  # pragma: no cover
    """
    Fetch whichever taxonomy asset the requested labels need
    """
    if args.add_gtdb_tax:
        from gtotree.utils.gtdb.get_gtdb_data import get_gtdb_data
        get_gtdb_data()
    if args.add_ncbi_tax:
        from gtotree.utils.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data
        get_ncbi_assembly_data()


################################################################################
# label building
################################################################################

def build_new_mapping_dict(args, run_data):
    """
    Replace the previous run's mapping dict with the one this invocation asks for

    Clearing first is what makes re-labelling work at all: both taxonomy handlers skip
    any genome already present in the dict, so carrying the old mapping in would mean
    nothing ever changed.
    """
    run_data.mapping_dict = {}
    run_data.tax_info_dict = {}
    run_data.mapping_file_path = ""

    if args.mapping_file:
        from gtotree.utils.misc.preflight_checks import check_mapping_file
        args, run_data = check_mapping_file(args, run_data)

    if args.add_gtdb_tax:
        from gtotree.utils.gtdb.handle_gtdb_tax_info import (
            update_mapping_dict_with_gtdb_tax_info)
        run_data = update_mapping_dict_with_gtdb_tax_info(args, run_data)

    if args.add_ncbi_tax:
        from gtotree.utils.ncbi.handle_ncbi_tax_info import (
            update_mapping_dict_with_ncbi_tax_info)
        run_data = update_mapping_dict_with_ncbi_tax_info(args, run_data)

    return args, run_data


################################################################################
# reporting
################################################################################

def _phase_counter():
    state = {"i": 0}

    def nxt():
        state["i"] += 1
        return state["i"]

    return nxt


def section(title):
    phase_stats.begin(title)
    report_phase_header(title)


def report_finish(out_dir, written, num_genomes, num_relabeled):
    """
    Counts come from the alignment rather than the tree, since a run made with `-N`
    has no tree but is still perfectly re-labelable.
    """
    # border = color_text("  " + "-" * 78, "green")
    # title = color_text("  " + "Headers updated!".center(78), "green")
    title = "  " + "Headers updated!".center(42)

    print()
    # print(border)
    print(title)
    # print(border)

    print(f"\n      {color_text(f'{num_relabeled:,}', 'green')} of "
          f"{color_text(f'{num_genomes:,}', 'green')} genome labels changed.\n")

    print("      Results written to:")
    print(f"        {color_text(out_dir + '/', 'green')}\n")

    print("      Including:")
    for entry in written:
        print(f"        {color_text(entry, 'green')}")
    print()


################################################################################
# driver
################################################################################

def run_update_headers(args):
    """
    The whole program, phase by phase
    """
    args = check_args(args)

    # n = _phase_counter()
    # section(f"Phase {n()}: Reading the completed GToTree run...")

    run_data = run_lib.load_completed_run(args.input_dir)
    run_data = run_lib.rehome_run_data(run_data, args.input_dir)

    previous_labels = run_lib.labels_in_completed_outputs(run_data)

    # num_in_tree = len(run_data.get_all_remaining_input_genome_ids())
    # print(f"      Found a finished run of "
    #       f"{color_text(f'{len(run_data.all_input_genomes):,}', 'green')} input "
    #       f"genome(s), {color_text(f'{num_in_tree:,}', 'green')} of which made it "
    #       f"into the final tree.")

    if not run_data.final_tree_path:
        report_message(
            "That run has no tree in it (either `-N` was used or the tree file has "
            "since been moved), so we'll re-label the alignment and summary table "
            "only.", "yellow", ii="      ", si="      ")

    # section(f"Phase {n()}: Resolving the new labels...")

    ensure_reference_data(args)

    if args.add_gtdb_tax or args.add_ncbi_tax:
        source = "GTDB" if args.add_gtdb_tax else "NCBI"
        with spinner(f"Pulling {source} lineages...", "", clear_on_done=True,
                     reclaim_line=True):
            args, run_data = build_new_mapping_dict(args, run_data)
    else:
        args, run_data = build_new_mapping_dict(args, run_data)

    new_labels = run_lib.new_labels_for(run_data)

    if not run_data.mapping_dict:
        report_message(
            "No new labels were specified (`-m`, `-D`, and `-t` were all left off), so "
            "the outputs will carry the original genome IDs.", "yellow",
            ii="      ", si="      ")
    # else:
    #     print(f"      {color_text(f'{len(run_data.mapping_dict):,}', 'green')} "
    #           f"genome(s) have a new label.")

    # section(f"Phase {n()}: Writing the updated outputs...")

    out_dir, written, num_genomes, num_relabeled = write_outputs(
        args, run_data, previous_labels, new_labels)

    print()
    report_finish(out_dir, written, num_genomes, num_relabeled)


def setup_output_dir(args):
    """
    Create the output directory, honoring `-F`

    Not `prepare_output_dir`: that carves out a working sub-directory for a run with
    intermediate state to keep, and this program has none -- everything it writes is
    a final output.
    """
    out_dir = args.output_dir.rstrip("/")

    if os.path.exists(out_dir):
        if not args.force_overwrite:
            raise OutputDirExistsError(
                f"The output directory '{out_dir}' already exists, and we don't want "
                "to overwrite anything accidentally. Use `-F` to overwrite it, or "
                "specify a different directory with `-o`.")
        shutil.rmtree(out_dir)

    os.makedirs(out_dir, exist_ok=True)

    return out_dir


def write_outputs(args, run_data, previous_labels, new_labels):
    """
    Everything that lands in the new output directory
    """
    out_dir = setup_output_dir(args)

    prefix = args.output_prefix.strip().rstrip("-")
    written = []

    alignment_name = output_name(prefix, f"{ALIGNMENT_BASENAME}{run_data.general_ext}")
    alignment_path = os.path.join(out_dir, alignment_name)
    _source, num_genomes, num_relabeled = run_lib.write_updated_alignment(
        run_data, previous_labels, new_labels, alignment_path)
    written.append(alignment_name)

    if run_data.final_tree_path:
        tree_name = output_name(prefix, os.path.basename(run_data.final_tree_path))
        tree_path = os.path.join(out_dir, tree_name)
        run_lib.write_updated_tree(run_data, previous_labels, new_labels, tree_path)
        written.append(tree_name)

    summary_name = output_name(prefix, SUMMARY_FILENAME)
    summary_path = os.path.join(out_dir, summary_name)
    run_lib.write_updated_summary_table(args, run_data, summary_path)
    written.append(summary_name)

    label_map_name = output_name(prefix, LABEL_MAP_FILENAME)
    run_lib.write_label_map(run_data, previous_labels, new_labels,
                            os.path.join(out_dir, label_map_name))
    written.append(label_map_name)

    for itol_dir in run_lib.regenerate_itol_files(run_data, summary_path, out_dir):
        written.append(os.path.relpath(itol_dir, out_dir) + "/")

    return out_dir, written, num_genomes, num_relabeled


def main():  # pragma: no cover
    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    try:
        run_update_headers(args)
    except UpdateHeadersError as e:
        report_very_early_exit(str(e))
    except OutputDirExistsError as e:
        report_very_early_exit(str(e))
    finally:
        phase_stats.finish()
        phase_stats.report()


if __name__ == "__main__":
    main()
