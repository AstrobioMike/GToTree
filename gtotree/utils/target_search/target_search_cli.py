"""
The shared CLI and driver behind `gtt search-pfams` and `gtt search-kos`.

Both subcommands call into here with their `TargetSearchSpec`; the two entry-point
modules are thin enough to read in one screen. The phased output, spinners, and
progress-bar conventions follow `gtt gen-scg-hmms` rather than the main GToTree driver.

The library-raises / CLI-translates split is the same as everywhere else: the stage
functions raise `TargetSearchError` (or let a taxonomy error propagate), and `main()`
turns those into a friendly message and a clean exit.
"""

import os
import sys
import shutil
import argparse

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.general import (read_run_data, write_run_data, CorruptRunData,
                                        OutputDirExistsError)
from gtotree.utils.misc.resume_state import (ResumeProfile, hash_strings,
                                        hash_local_genomes, hash_file_contents,
                                        STATE_VERSION)
from gtotree.utils.misc.messaging import (report_message, color_text, spinner,
                                     report_very_early_exit)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import TaxonNotFound, AmbiguousTaxon
from gtotree.utils.taxonomy.wanted_ref_tax import WantedRefTaxError
from gtotree.utils.target_search import target_search_stages as stages
from gtotree.utils.target_search import target_search_outputs as outputs
from gtotree.utils.target_search.target_search_setup import (
    TargetSearchError,
    RUN_DATA_FILENAME,
    check_args,
    check_dependencies,
    check_env_vars,
    setup_output_dir,
    build_run_data,
    adopt_genome_progress,
    ensure_processing_dirs,
    ensure_reference_data,
    validate_input_files,
    fill_in_shared_args,
)


DEFAULT_NUM_JOBS = 10


################################################################################
# resume
################################################################################

# stage names, in pipeline order
STAGE_TARGETS = "targets"
STAGE_SEARCH = "search"
STAGE_OUTPUTS = "outputs"

STAGE_ORDER = [STAGE_TARGETS, STAGE_SEARCH, STAGE_OUTPUTS]

RESUME = ResumeProfile(
    name="target-search",
    stages=STAGE_ORDER,
    field_labels={
        "state_version": "the run-state format",
        "subcommand": "the subcommand",
        "accessions_sha256": "the set of input NCBI accessions",
        "local_genomes_sha256": "the local genome files (contents, paths, or set)",
        "targets_sha256": "the list of search targets",
        "source": "--source",
        "wanted_ref_tax": "--wanted-ref-tax",
        "target_rank": "--target-rank",
        "derep_rank": "--derep-rank",
        "data_version": "the reference-database version",
    },
    # the database version isn't known until it's been fetched, so a run interrupted
    # before then legitimately has None stored
    deferred_fields=("data_version",),
)


def build_fingerprint(run_data, args, spec, data_version=None):
    """
    Everything that affects what this run produces.

    Two layers cooperate on resume, and the split is worth being explicit about:

      * this fingerprint guards the run as a whole, refusing a `-R` whose parameters
        changed in a way that would mix results from two different runs
      * the per-genome flags on GenomeData (`processing_done`, `pfam_search_done`,
        `ko_search_done`) give genome-level resume, since `genomes_needing_processing`
        filters on exactly those

    So there are no per-genome stages here: run-data.json is that record.

    Deliberately excludes `--num-jobs`, `--keep-working-dir`, `-d`, and the output
    directory name -- those change how the run executes, not what it produces.

    The target list is hashed by CONTENTS rather than path: the meaningful unit is
    which IDs were asked for, so moving or renaming the file shouldn't force a re-run,
    and editing it in place absolutely should.
    """
    accessions = run_data.get_input_ncbi_accs()
    local_genomes = (list(run_data.genbank_files) + list(run_data.fasta_files)
                     + list(run_data.amino_acid_files))

    return {
        "state_version": STATE_VERSION,
        "subcommand": spec.subcommand,
        "accessions_sha256": hash_strings(accessions),
        "num_accessions": len(set(accessions)),
        "local_genomes_sha256": hash_local_genomes(local_genomes),
        "num_local_genomes": len(local_genomes),
        "targets_sha256": hash_file_contents(spec.targets_file(args)),
        "num_targets": getattr(args, "total_targets", None),
        "source": (args.source or "").lower(),
        "wanted_ref_tax": args.wanted_ref_tax,
        "target_rank": args.target_rank,
        "derep_rank": args.derep_rank,
        "data_version": data_version,
    }


################################################################################
# parser
################################################################################

def build_parser(spec, parent_subparsers=None):
    """
    Build this subcommand's parser.

    The genome-input flags deliberately match the main GToTree program's short flags
    (`-a`, `-g`, `-f`, `-A`, `-w`), since these subcommands are a subset of it and the
    same input files should work in both without editing.
    """
    targets = spec.target_label_plural

    desc = (f"This is a helper program that searches a set of input genomes for a list "
            f"of target {targets}. It takes the same genome inputs as the main GToTree "
            f"program, preprocesses them the same way, and produces the same "
            f"{spec.target_label} search results a full GToTree run would have produced "
            f"if you'd passed `{spec.targets_flag}` to it.")

    example = (f"Ex. usage: `gtt {spec.subcommand} -f my-genomes.txt "
               f"{spec.targets_flag} my-targets.txt`")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            spec.subcommand,
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
    genomes = parser.add_argument_group("Input Genomes (at least one)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        spec.targets_flag, spec.targets_flag_long,
        metavar="<FILE>",
        dest=spec.targets_dest,
        help=(f"A single-column file of the target {spec.target_label} IDs to search "
              f"for (e.g. '{spec.example_target}')."),
        action="store",
    )

    genomes.add_argument(
        "-w", "--wanted-ref-tax",
        metavar="<STR>",
        help=("A target taxon whose reference genomes should also be searched (e.g. "
              "'Nitrospirota')."),
        action="store",
    )

    genomes.add_argument(
        "-a", "--ncbi-accessions",
        metavar="<FILE>",
        help="A single-column file of NCBI accessions to search.",
        action="store",
    )

    genomes.add_argument(
        "-f", "--fasta-files",
        metavar="<FILE>",
        help=("A single-column file listing nucleotide fasta files to search; genes "
              "will be called with prodigal."),
        action="store",
    )

    genomes.add_argument(
        "-A", "--amino-acid-files",
        metavar="<FILE>",
        help=("A single-column file listing amino-acid fasta files to search; each "
              "should hold the proteins for one genome."),
        action="store",
    )

    genomes.add_argument(
        "-g", "--genbank-files",
        metavar="<FILE>",
        help="A single-column file listing GenBank files to search.",
        action="store",
    )

    optional.add_argument(
        "--source",
        type=str.lower,
        default="gtdb",
        choices=["gtdb", "ncbi"],
        help=("Which taxonomy source to select `--wanted-ref-tax` genomes from "
              "(default: gtdb)"),
        action="store",
    )

    optional.add_argument(
        "--target-rank",
        choices=list(RANKS),
        help=("Target rank, if needed to disambiguate a taxon name that exists at "
              "multiple ranks"),
        action="store",
    )

    optional.add_argument(
        "--derep-rank",
        default="auto",
        choices=["auto", "off"] + list(RANKS),
        help=("Dereplicate `-w`-selected genomes down to one per unique value of this "
              "rank (default: auto, which uses two ranks finer than the target). Use "
              "'off' to keep them all."),
        action="store",
    )

    optional.add_argument(
        "-o", "--output-directory",
        metavar="<DIR>",
        default=spec.default_output_dir,
        dest="output_dir",
        help=f'Desired output directory (default: "{spec.default_output_dir}")',
        action="store",
    )

    optional.add_argument(
        "-j", "--num-jobs",
        metavar="<INT>",
        default=DEFAULT_NUM_JOBS,
        type=int,
        help=(f"Number of genomes to process and search concurrently (default: "
              f"{DEFAULT_NUM_JOBS})"),
        action="store",
    )

    optional.add_argument(
        "--keep-working-dir",
        action="store_true",
        help="Keep the intermediate working directory",
    )

    optional.add_argument(
        "-R", "--resume",
        action="store_true",
        help=("Resume a previous run in the same output directory, reusing any genomes "
              "that already finished. Refuses if the run parameters changed."),
    )

    optional.add_argument(
        "-F", "--force-overwrite",
        action="store_true",
        help="Overwrite the output directory if it already exists",
    )

    optional.add_argument(
        "-d", "--debug",
        action="store_true",
        help="Keep each genome's intermediate sequence files",
    )

    add_help(optional)
    add_version_arg(optional)

    parser.set_defaults(func=spec.subcommand.replace("-", "_"))

    return parser


################################################################################
# terminal styling (matching gen-scg-hmms)
################################################################################

def _phase_counter():
    """Returns a callable yielding 1, 2, 3, ... on each call (for phase labels)."""
    state = {"i": 0}

    def nxt():
        state["i"] += 1
        return state["i"]

    return nxt


def section(title):
    phase_stats.begin(title)
    print(color_text(f"\n\n  {title}\n", "yellow"))


def report_finish(out_dir, run_data, spec, summary_path, targets_with_hits):
    """
    The closing banner, mirroring gen-scg-hmms'.

    `targets_with_hits` is how many targets actually got a combined hit-seqs fasta. At
    zero there is no hit-seqs directory to point at, and saying so plainly is a lot
    more useful than listing an output that isn't there -- a run where nothing hit
    anything is a real (and usually surprising) result, not an error.
    """
    searched, removed, failed = outputs.summarize_counts(run_data, spec)
    num_targets = len(spec.found_targets(run_data))

    border = color_text("  " + "-" * 78, "green")
    title = color_text("  " + f"{spec.target_label} search complete!".center(78), "green")
    print()
    print(border)
    print(title)
    print(border)

    print(f"\n      {color_text(f'{searched:,}', 'green')} genome(s) searched for "
          f"{color_text(f'{num_targets:,}', 'green')} {spec.target_label} target(s).")

    if removed or failed:
        dropped = removed + failed
        report_message(
            f"{dropped:,} input genome(s) didn't make it through; see summary table.", "yellow", ii="      ", si="      ")
        print()

    if not targets_with_hits:
        report_message(
            f"No hits were found for any of the {num_targets:,} "
            f"{spec.target_label} target(s) in any of the genomes searched :/", "yellow", ii="      ", si="      ", newline=False)
        print()

    print("      Results written to:")
    print(f"        {color_text(out_dir + '/', 'green')}\n")

    print("      Including:")
    print(f"        {color_text(spec.counts_filename, 'green')}")
    if targets_with_hits:
        print(f"        {color_text(spec.hit_seqs_subdir + '/', 'green')}")
    print(f"        {color_text(os.path.basename(summary_path), 'green')}")
    print()


################################################################################
# driver
################################################################################

def _load_previous_run_data(work_dir):
    """
    Load a prior run's RunData, or None. A corrupt file is a hard stop rather than a
    silent fresh start: silently discarding it would re-download everything while
    looking like a resume.
    """
    path = os.path.join(work_dir, RUN_DATA_FILENAME)
    if not os.path.isfile(path):
        return None
    try:
        return read_run_data(path)
    except CorruptRunData as e:
        raise TargetSearchError(
            f"We're trying to resume a previous run, but {e}. That usually means the "
            "previous run was interrupted while saving its state. Your best bet is to "
            "start fresh with `-F`.") from e


def ensure_all_required_data(args, spec):  # pragma: no cover
    """
    Make every managed dataset this run needs is present
    """
    ensure_reference_data(args, spec)
    spec.ensure_data()

    check_env_vars(spec)

    if spec.describe_data_version is None:
        return None

    return spec.describe_data_version(None)


def run_search(args, spec):  # pragma: no cover
    """
    The whole program, phase by phase.
    """
    args = fill_in_shared_args(args)
    args = check_args(args, spec)
    check_dependencies(args, spec)

    out_dir, work_dir = setup_output_dir(args, spec)

    resuming = bool(args.resume)
    state = RESUME.load(work_dir) if resuming else None
    previous = _load_previous_run_data(work_dir) if resuming else None

    n = _phase_counter()

    args = validate_input_files(args, spec)
    data_version = ensure_all_required_data(args, spec)

    section(f"Phase {n()}: Resolving input genomes...")
    run_data = build_run_data(args, spec, out_dir, work_dir, previous=previous)

    run_data, _selection = stages.resolve_input_genomes(args, run_data)

    ensure_processing_dirs(run_data)

    adopt_genome_progress(run_data, previous)

    fingerprint = build_fingerprint(run_data, args, spec, data_version=data_version)

    if resuming and state:
        differences = RESUME.compare(state.get("fingerprint"), fingerprint)
        if differences:
            raise TargetSearchError(
                RESUME.refusal_message(differences, flag="`-R`/`--resume`"))
        print(f"\n      {color_text('Resuming from the previous run', 'green')}")
    else:
        resuming = False
        state = RESUME.new(fingerprint)
        RESUME.save(work_dir, state)

    # after the fingerprint gate, since the gate compares the accessions as given and
    # this is where some of them stop being part of the run
    run_data = stages.lookup_ncbi_accessions(run_data)

    write_run_data(run_data)

    section(f"Phase {n()}: Collecting target {spec.target_label_plural}...")
    if data_version:
        print(f"      {spec.target_label} version being used: "
              f"{color_text(data_version, 'green')}")

    state.setdefault("fingerprint", {})["data_version"] = data_version

    run_data = stages.phase_collect_targets(run_data, spec, out_dir, resuming=resuming)
    write_run_data(run_data)
    RESUME.mark_complete(state, STAGE_TARGETS, work_dir=work_dir)
    RESUME.save(work_dir, state)

    section(f"Phase {n()}: Processing and searching genomes...")
    plan = stages.build_plan(args, spec)
    run_data = stages.phase_search_genomes(args, run_data, spec, plan)
    RESUME.mark_complete(state, STAGE_SEARCH, work_dir=work_dir)
    RESUME.save(work_dir, state)

    section(f"Phase {n()}: Combining results...")
    summary_path, targets_with_hits = stages.phase_write_outputs(run_data, spec, out_dir)
    RESUME.mark_complete(state, STAGE_OUTPUTS, work_dir=work_dir)
    RESUME.save(work_dir, state)

    if not args.keep_working_dir and not args.debug:
        shutil.rmtree(work_dir, ignore_errors=True)

    # closes the final phase, so this has to happen before anything reads the table
    phase_stats.finish()
    phase_stats.write_tsv(out_dir)

    report_finish(out_dir, run_data, spec, summary_path, targets_with_hits)

    return run_data


def main(spec):  # pragma: no cover
    parser = build_parser(spec)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    try:
        run_search(args, spec)
    except KeyboardInterrupt:
        print()
        report_very_early_exit("Interrupted by user.", "yellow")
    except (TaxonNotFound, AmbiguousTaxon, WantedRefTaxError) as e:
        report_very_early_exit(str(e))
    except OutputDirExistsError as e:
        report_very_early_exit(str(e), "yellow", leading_newline=False)
    except TargetSearchError as e:
        report_very_early_exit(str(e))
    finally:
        phase_stats.report()
