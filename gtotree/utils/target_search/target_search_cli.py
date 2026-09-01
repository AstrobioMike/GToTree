"""
The driver behind `gtt search-annotations`.

`search-annotations` is the one command for searching input genomes for target
annotations -- Pfam domains (`-p`) and/or KEGG Orthologs (`-K`).

How the pieces fit together:

  * The per-type behavior lives entirely in `TargetSearchSpec` (one spec per target
    type, in `target_search_spec.py`) and the shared stage/setup/output helpers. This
    driver drives a *list* of specs; every spec-specific step -- collecting targets,
    pressing profiles, writing the counts table, combining hit seqs -- is a call
    through those shared helpers, so a one-type run and a full GToTree run with
    `-p`/`-K` produce byte-for-byte the same per-type files.
  * The genome phases (resolve inputs, look up NCBI accessions, preprocess) are
    target-type-independent, so they run once. The fused per-genome worker already
    dispatches on `plan.do_pfam` / `plan.do_ko`, so a plan with both flags on searches
    a genome for both in the worker that produced its proteins.

Output layout: each target type gets its own subdirectory under the output directory
(`<out>/pfam/`, `<out>/ko/`), holding that type's counts table, hit seqs, and a
`<type>-genomes-summary-info.tsv`. Run-level files that aren't type-specific -- the
top-level `genomes-summary-info.tsv` (genome preprocessing only), the removed-genomes
report, and phase-stats -- live at the top.
"""

import os
import sys
import shutil
import tempfile
import argparse

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.general import (read_run_data, write_run_data, CorruptRunData,
                                        OutputDirExistsError, adopt_genome_progress,
                                        run_pooled_stage,
                                        GTT_PROGRESS_BAR_FORMAT_INDENTED_6)
from gtotree.utils.misc.resume_state import (ResumeProfile, hash_strings,
                                        hash_local_genomes, hash_file_contents,
                                        STATE_VERSION)
from gtotree.utils.misc.messaging import (report_message, color_text, spinner,
                                     report_phase_header, report_very_early_exit)
from gtotree.utils.misc.summary_info import write_removed_genomes_report
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import TaxonNotFound, AmbiguousTaxon, CrossDomainTaxon
from gtotree.utils.taxonomy.wanted_ref_tax import WantedRefTaxError
from gtotree.utils.misc.general import wanted_ref_tax_list
from gtotree.utils.hmms.hmm_searching_engine import press_profiles
from gtotree.main_stages.processing_genomes import (SearchPlan, _fused,
                                                    genomes_needing_processing)
from gtotree.utils.target_search import target_search_stages as stages
from gtotree.utils.target_search import target_search_outputs as outputs
from gtotree.utils.target_search.target_search_spec import get_spec
from gtotree.utils.target_search.target_search_setup import (
    TargetSearchError,
    RUN_DATA_FILENAME,
    check_args_multi,
    requested_specs,
    check_dependencies,
    check_env_vars,
    build_run_data,
    ensure_processing_dirs,
    ensure_reference_data,
    validate_genome_input_files,
    validate_targets_file,
    fill_in_shared_args,
    wire_spec_results_dirs,
    make_spec_result_dirs,
)


SUBCOMMAND = "search-annotations"
DEFAULT_OUTPUT_DIR = "gtt-annotation-search-output"
DEFAULT_NUM_JOBS = 10

SPEC_NAMES = ["pfam", "ko"]

SUBDIRS = {"pfam": "pfam", "ko": "ko"}


def _all_specs():
    return [get_spec(name) for name in SPEC_NAMES]


################################################################################
# resume
################################################################################

STAGE_TARGETS = "targets"
STAGE_SEARCH = "search"
STAGE_OUTPUTS = "outputs"

STAGE_ORDER = [STAGE_TARGETS, STAGE_SEARCH, STAGE_OUTPUTS]

RESUME = ResumeProfile(
    name="annotation-search",
    stages=STAGE_ORDER,
    field_labels={
        "state_version": "the run-state format",
        "subcommand": "the subcommand",
        "target_types": "which annotation types are being searched (-p / -K)",
        "accessions_sha256": "the set of input NCBI accessions",
        "local_genomes_sha256": "the local genome files (contents, paths, or set)",
        "pfam_targets_sha256": "the list of target Pfams",
        "ko_targets_sha256": "the list of target KOs",
        "source": "--source",
        "ncbi_section": "--ncbi-section",
        "wanted_ref_tax": "--wanted-ref-tax",
        "target_rank": "--target-rank",
        "target_domain": "--target-domain",
        "derep_rank": "--derep-rank",
        "pfam_data_version": "the Pfam database version",
        "ko_data_version": "the KO database version",
    },
    deferred_fields=("pfam_data_version", "ko_data_version"),
)


def build_fingerprint(run_data, args, specs, data_versions=None):
    """
    Everything that affects what this combined run produces.

    Same two-layer resume as the single command (this fingerprint guards the run as a
    whole; per-genome flags give genome-level resume), extended so each target type has
    its own targets hash and data version. That independence is deliberate: editing the
    Pfam list must force the Pfam targets to be recollected without touching the KO
    side, and adding `-K` to a previously `-p`-only run is a genuine change to what the
    run produces, caught by `target_types` here.

    Target lists are hashed by CONTENTS, not path, so moving the file doesn't force a
    re-run and editing it in place does.
    """
    data_versions = data_versions or {}
    requested = {spec.subcommand for spec in requested_specs(args, specs)}

    accessions = run_data.get_input_ncbi_accs()
    local_genomes = (list(run_data.genbank_files) + list(run_data.fasta_files)
                     + list(run_data.amino_acid_files))

    pfam_spec = get_spec("pfam")
    ko_spec = get_spec("ko")

    return {
        "state_version": STATE_VERSION,
        "subcommand": SUBCOMMAND,
        "target_types": sorted(requested),
        "accessions_sha256": hash_strings(accessions),
        "num_accessions": len(set(accessions)),
        "local_genomes_sha256": hash_local_genomes(local_genomes),
        "num_local_genomes": len(local_genomes),
        "pfam_targets_sha256": hash_file_contents(pfam_spec.targets_file(args)),
        "ko_targets_sha256": hash_file_contents(ko_spec.targets_file(args)),
        "source": (args.source or "").lower(),
        "ncbi_section": (getattr(args, "ncbi_section", None) or "").lower(),
        "wanted_ref_tax": (sorted(wanted_ref_tax_list(args)) or None),
        "target_rank": args.target_rank,
        "target_domain": getattr(args, "target_domain", None),
        "derep_rank": args.derep_rank,
        "pfam_data_version": data_versions.get("pfam"),
        "ko_data_version": data_versions.get("ko"),
    }


################################################################################
# parser
################################################################################

def build_parser(parent_subparsers=None):
    """
    Build the `search-annotations` parser.

    Both target flags are optional here; `check_args_multi` enforces that at least one
    is given. The genome-input flags match the main GToTree program's short flags
    (`-a`, `-g`, `-f`, `-A`, `-w`) so the same input files work in both.
    """
    pfam_spec = get_spec("pfam")
    ko_spec = get_spec("ko")

    desc = ("This is a helper program that searches a set of input genomes for target "
            "KOs and/or Pfams. It takes the same genome inputs as the main GToTree program, "
            "preprocesses them the same way, and produces the same Pfam/KO search results "
            "a full GToTree run would produce if you passed `-p` and/or `-K` to it.")

    example = ("Ex. usage: `gtt search-annotations -w Alteromonas -f my-fasta-files.txt "
               "-p target-pfams.txt -K target-kos.txt`")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            SUBCOMMAND,
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

    targets = parser.add_argument_group("Search Targets (at least one)")
    genomes = parser.add_argument_group("Input Genomes (at least one)")
    optional = parser.add_argument_group("Optional Parameters")

    targets.add_argument(
        pfam_spec.targets_flag, pfam_spec.targets_flag_long,
        metavar="<FILE>",
        dest=pfam_spec.targets_dest,
        help=(f"A single-column file of the target Pfam IDs to search for (e.g. "
              f"'{pfam_spec.example_target}')."),
        action="store",
    )

    targets.add_argument(
        ko_spec.targets_flag, ko_spec.targets_flag_long,
        metavar="<FILE>",
        dest=ko_spec.targets_dest,
        help=(f"A single-column file of the target KO IDs to search for (e.g. "
              f"'{ko_spec.example_target}')."),
        action="store",
    )

    genomes.add_argument(
        "-w", "--wanted-ref-tax",
        metavar="<STR>",
        help=("A target taxon whose reference genomes should also be searched (e.g., "
              "'Nitrospirota'). May be given more than once to pool several taxa into "
              "one set (e.g. `-w Bacteria -w Archaea`); each is resolved and "
              "dereplicated on its own, then merged."),
        action="append",
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
        "--ncbi-section",
        dest="ncbi_section",
        type=str.lower,
        default="both",
        choices=["refseq", "genbank", "both"],
        help=("Which section of NCBI to draw `--wanted-ref-tax` genomes from "
              "(default: both). The default of 'both' is typically fine with the "
              "default `--derep-rank auto`. Ignored with `--source gtdb`."),
        action="store",
    )

    optional.add_argument(
        "--target-rank",
        type=str.lower,
        choices=list(RANKS),
        help=("Target rank, if needed to disambiguate a taxon name that exists at "
              "multiple ranks"),
        action="store",
    )

    optional.add_argument(
        "--target-domain",
        type=str,
        dest="target_domain",
        default=None,
        help=("Target domain, if needed to disambiguate a taxon name that exists in "
              "multiple domains (e.g., bacillus is both a bacterial and eukaryotic genus)"),
        action="store",
    )

    optional.add_argument(
        "--derep-rank",
        default="auto",
        type=str.lower,
        choices=["auto", "off"] + list(RANKS),
        help=("Dereplicate `-w`-selected genomes down to one per unique value of this "
              "rank (default: auto, which uses two ranks finer than the target). Use "
              "'off' to keep them all."),
        action="store",
    )

    optional.add_argument(
        "-o", "--output-dir",
        metavar="<DIR>",
        default=DEFAULT_OUTPUT_DIR,
        dest="output_dir",
        help=f'Desired output directory (default: "{DEFAULT_OUTPUT_DIR}")',
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
        help=("Keep the working directory and typically temp intermediate files"),
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

    add_help(optional)
    add_version_arg(optional)

    parser.set_defaults(func=SUBCOMMAND.replace("-", "_"))

    return parser


################################################################################
# terminal styling (matching the single-command driver)
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


def _spec_out_dir(out_dir, spec_name):
    return os.path.join(out_dir, SUBDIRS[spec_name])


def report_finish(out_dir, run_data, active, summary_path):
    """
    The closing banner for a combined run.

    `active` is the list of (spec_name, spec, targets_with_hits) actually run. One
    genome count is reported (genomes are shared), then a per-type line for what each
    search found and where it landed.
    """
    # every active spec searched the same surviving genome set, so any of them gives
    # the shared genome count
    any_spec = active[0][1]
    searched, removed, failed = outputs.summarize_counts(run_data, any_spec)

    labels = " and ".join(spec.target_label_plural for _, spec, _ in active)

    border = color_text("  " + "-" * 78, "green")
    title = color_text("  " + "Annotation search complete!".center(78), "green")
    print()
    print(border)
    print(title)
    print(border)

    print(f"\n      {color_text(f'{searched:,}', 'green')} genome(s) searched for "
          f"{labels}.\n")

    print("      Results written to:")
    print(f"        {color_text(out_dir + '/', 'green')}\n")

    for spec_name, spec, targets_with_hits in active:
        num_targets = len(spec.found_targets(run_data))
        sub = SUBDIRS[spec_name] + "/"
        print(f"      {color_text(spec.target_label, 'green')} results "
              f"({color_text(f'{num_targets:,}', 'green')} target(s)) in "
              f"{color_text(sub, 'green')}")
        if not targets_with_hits:
            report_message(
                f"No hits were found for any {spec.target_label} target in any "
                "genome searched.", "yellow", ii="        ", si="        ",
                newline=False)
        else:
            print(f"        {color_text(sub + spec.counts_filename, 'green')}")
            print(f"        {color_text(sub + spec.hit_seqs_subdir + '/', 'green')}")
        print(f"        {color_text(sub + spec.summary_filename, 'green')}")
        if spec.failed_targets(run_data):
            print(f"        {color_text(sub + spec.failed_targets_filename, 'green')}")
        print()

    print("      Genome preprocessing summary (all input genomes):")
    print(f"        {color_text(os.path.basename(summary_path), 'green')}\n")

    if removed or failed:
        report_message("Any input genomes that didn't make it through are reported in:",
                       "yellow", ii="      ", si="      ", newline=False)
        print(f"        {color_text(stages.removals_pointer(run_data), 'yellow')}\n")
    print()


################################################################################
# driver
################################################################################

def _load_previous_run_data(work_dir):
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


def ensure_all_required_data(args, specs):  # pragma: no cover
    """
    Fetch every managed dataset the requested target types need, returning a
    {spec_name: version-or-None} map for the run-info block and fingerprint.

    Only the requested types are fetched -- a `-p`-only run never touches the KO data.
    """
    ensure_reference_data(args, specs[0] if specs else get_spec("pfam"))

    versions = {}
    for name in SPEC_NAMES:
        spec = get_spec(name)
        if spec not in specs:
            continue
        spec.ensure_data()
        check_env_vars(spec)
        versions[name] = (spec.describe_data_version(None)
                          if spec.describe_data_version else None)
    return versions


def _prepare_spec(run_data, args, spec, spec_name, out_dir, tmp_dir, total):
    """
    Root one spec's results at its output subdirectory and create its result subdirs.

    Called after the shared RunData exists so each target type writes to `<out>/pfam/`
    or `<out>/ko/` instead of colliding at the top level.
    """
    spec_out = _spec_out_dir(out_dir, spec_name)
    make_spec_result_dirs(spec_out, spec)
    wire_spec_results_dirs(
        run_data, args, spec, tmp_dir,
        results_root=spec_out,
        results_root_rel=os.path.join(run_data.output_dir_rel, SUBDIRS[spec_name]),
        tmp_subdir=os.path.join("target-hit-seqs", spec_name),
        total_targets=total,
    )
    return spec_out


def _combined_plan(args, specs):
    """
    A SearchPlan turning on exactly the requested target types' searches.

    `do_scg=False` for the same reason the single command sets it: there's no SCG set
    to press or search and no tree downstream. The fused worker already runs whichever
    of Pfam/KO are on, so one plan drives a one-pass search for both.
    """
    names = {spec.subcommand for spec in specs}
    return SearchPlan(
        do_pfam=(get_spec("pfam").subcommand in names),
        do_ko=(get_spec("ko").subcommand in names),
        keep_genome_files=bool(getattr(args, "keep_working_dir", False)),
        do_scg=False,
    )


def phase_search_genomes(args, run_data, specs, plan):
    """
    Preprocess and search every not-yet-finished genome in one pool for all requested
    target types.

    This is the single-command `phase_search_genomes` generalized to press each
    profile-based type once (only Pfam presses today) and report both target counts.
    `genome_is_fully_processed` already ANDs every requested search flag, so a genome
    that finished Pfam but not KO in an interrupted run is correctly re-picked here and
    the worker skips the Pfam work it already has.
    """
    phase_stats.begin("processing and searching genomes")

    press_dirs = []
    try:
        for spec in specs:
            if not spec.presses_profiles:
                continue
            press_dir = tempfile.mkdtemp(prefix="gtt-press-")
            press_dirs.append(press_dir)
            hmm_path = getattr(run_data, spec.combined_hmm_attr)
            with spinner(f"Preparing {spec.target_label} profiles for searching...",
                         f"Prepared {spec.target_label} profiles"):
                plan.pressed_pfam_base = press_profiles(
                    hmm_path, press_dir, "target-profiles")

        to_process = genomes_needing_processing(run_data.all_input_genomes, plan)
        alive = [gd for gd in run_data.all_input_genomes if not gd.removed]
        already_done = len(alive) - len(to_process)

        labels = " and ".join(
            f"{len(spec.found_targets(run_data)):,} {spec.target_label} target(s)"
            for spec in specs)

        if not to_process:
            print(f"\n      All {len(alive):,} genome(s) were already processed and "
                  "searched in a previous run")
            return run_data

        print(f"\n      Processing and searching {len(to_process):,} genome(s) for "
              f"{labels}:")
        if already_done:
            print(f"        ({already_done:,} already done in a previous run)")

        preprocess, apply_status = stages._dispatching_worker_pair(run_data)
        worker, apply_result = _fused(preprocess, apply_status, plan)

        run_data = run_pooled_stage(to_process, worker, apply_result, args, run_data,
                                    bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED_6)
    finally:
        for press_dir in press_dirs:
            shutil.rmtree(press_dir, ignore_errors=True)

    for spec in specs:
        stages.mark_failed_searches_removed(run_data, spec)

    write_removed_genomes_report(run_data)
    write_run_data(run_data)

    # a genome is "searched" if it completed at least one requested search; per-type
    # completion is what the summary table and per-type reporting reflect
    any_searched = any(stages.count_searched_genomes(run_data, spec) for spec in specs)
    failed_here = run_data.genomes_removed_at(*stages.SEARCH_PHASE_REMOVAL_STAGES)

    if failed_here:
        report_message(
            f"{len(failed_here):,} genome(s) failed the search phase. Reported in:",
            "yellow", ii="      ", si="      ")
        print(f"        {color_text(stages.removals_pointer(run_data), 'yellow')}")

    if not any_searched:
        raise TargetSearchError(
            "None of the input genomes made it through to a completed search. The "
            f"reason for each is in {stages.removals_pointer(run_data)}.")

    return run_data


def run_search(args, specs=None):  # pragma: no cover
    """
    The whole combined program, phase by phase.

    `specs` defaults to all target types; only the ones whose flag was given are
    actually run, so passing a subset is only needed by tests.
    """
    if specs is None:
        specs = _all_specs()

    args = fill_in_shared_args(args)
    args = check_args_multi(args, specs)

    active_specs = requested_specs(args, specs)
    for spec in active_specs:
        check_dependencies(args, spec)

    out_dir, work_dir = setup_output_dir_multi(args, active_specs)

    resuming = bool(args.resume)
    state = RESUME.load(work_dir) if resuming else None
    previous = _load_previous_run_data(work_dir) if resuming else None

    n = _phase_counter()

    args = validate_genome_input_files(args)
    totals = {}
    for spec_name in SPEC_NAMES:
        spec = get_spec(spec_name)
        if spec in active_specs:
            totals[spec_name] = validate_targets_file(args, spec)

    data_versions = ensure_all_required_data(args, active_specs)

    section(f"Phase {n()}: Resolving input genomes...")
    run_data = build_run_data(args, active_specs[0], out_dir, work_dir,
                              previous=previous)

    # re-root every active spec's results at its own subdirectory (build_run_data wired
    # only the first spec, at the flat out_dir)
    for spec_name in SPEC_NAMES:
        spec = get_spec(spec_name)
        if spec in active_specs:
            _prepare_spec(run_data, args, spec, spec_name, out_dir, run_data.tmp_dir,
                          totals[spec_name])

    run_data, _selection = stages.resolve_input_genomes(args, run_data)

    ensure_processing_dirs(run_data)

    adopt_genome_progress(run_data, previous)

    fingerprint = build_fingerprint(run_data, args, specs, data_versions=data_versions)

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

    plan = _combined_plan(args, active_specs)

    if not (resuming and stages.all_genomes_already_processed(run_data, plan)):
        run_data = stages.lookup_ncbi_accessions(run_data)

    write_run_data(run_data)

    section(f"Phase {n()}: Collecting targets...")
    for spec_name in SPEC_NAMES:
        spec = get_spec(spec_name)
        if spec not in active_specs:
            continue
        version = data_versions.get(spec_name)
        if version:
            print(f"      {spec.target_label} version being used: "
                  f"{color_text(version, 'green')}")
        run_data = stages.phase_collect_targets(
            run_data, spec, _spec_out_dir(out_dir, spec_name), resuming=resuming)

    state.setdefault("fingerprint", {})["pfam_data_version"] = data_versions.get("pfam")
    state.setdefault("fingerprint", {})["ko_data_version"] = data_versions.get("ko")

    write_run_data(run_data)
    RESUME.mark_complete(state, STAGE_TARGETS, work_dir=work_dir)
    RESUME.save(work_dir, state)

    section(f"Phase {n()}: Processing and searching genomes...")
    run_data = phase_search_genomes(args, run_data, active_specs, plan)
    RESUME.mark_complete(state, STAGE_SEARCH, work_dir=work_dir)
    RESUME.save(work_dir, state)

    section(f"Phase {n()}: Combining results...")
    active = []
    summary_path = None
    for spec_name in SPEC_NAMES:
        spec = get_spec(spec_name)
        if spec not in active_specs:
            continue
        spec_out = _spec_out_dir(out_dir, spec_name)
        s_path, targets_with_hits = stages.phase_write_outputs(run_data, spec, spec_out)
        active.append((spec_name, spec, targets_with_hits))

    # one run-level genomes summary at the top, preprocessing columns only (per-type
    # hit counts and search status live in each type's own subdirectory table)
    summary_path = outputs.write_root_genomes_summary(out_dir, run_data)

    RESUME.mark_complete(state, STAGE_OUTPUTS, work_dir=work_dir)
    RESUME.save(work_dir, state)

    if not args.keep_working_dir:
        shutil.rmtree(work_dir, ignore_errors=True)

    phase_stats.finish()
    phase_stats.write_tsv(out_dir)

    report_finish(out_dir, run_data, active, summary_path)

    return run_data


def setup_output_dir_multi(args, specs):
    """
    Create the output dir and working dir, then each requested type's subdirectory.

    Mirrors the single command's `setup_output_dir` (which honors -F/-R via
    `prepare_output_dir`) but lays out `<out>/pfam/` and `<out>/ko/` instead of
    flattening result subdirs into the top level.
    """
    from gtotree.utils.misc.general import prepare_output_dir
    from gtotree.utils.target_search.target_search_setup import WORKING_DIR_NAME

    out_dir, work_dir = prepare_output_dir(args.output_dir,
                                           resume=args.resume,
                                           force_overwrite=args.force_overwrite,
                                           work_dir_name=WORKING_DIR_NAME,
                                           ii="      ", si="      ")

    for spec_name in SPEC_NAMES:
        spec = get_spec(spec_name)
        if spec in specs:
            make_spec_result_dirs(_spec_out_dir(out_dir, spec_name), spec)

    return out_dir, work_dir


def main():  # pragma: no cover
    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    try:
        run_search(args)
    except KeyboardInterrupt:
        print()
        report_very_early_exit("Interrupted by user.", "yellow")
    except (TaxonNotFound, AmbiguousTaxon, CrossDomainTaxon, WantedRefTaxError) as e:
        report_very_early_exit(str(e))
    except OutputDirExistsError as e:
        report_very_early_exit(str(e), "yellow", leading_newline=False)
    except TargetSearchError as e:
        report_very_early_exit(str(e))
    finally:
        phase_stats.report()


if __name__ == "__main__":
    main()
