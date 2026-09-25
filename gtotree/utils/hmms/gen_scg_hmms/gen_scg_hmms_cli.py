"""
`gtt gen-scg-hmms` -- generate a single-copy-gene HMM set from a set of target genomes.

This is the CLI/orchestration layer. It owns argument parsing, the phased terminal
output (named sections, spinners, progress bars styled after bit's `gen-mg`), and
translating the library layer's exceptions into friendly messages + a clean exit.
All of the actual work lives in the sibling modules:

    gen_scg_hmms.py           Pfam coverage filtering, profile extraction, single-copy calls
    gen_scg_hmms_genomes.py   target-genome resolution and amino-acid retrieval
    gen_scg_hmms_search.py    the hmmsearch stage
    gen_scg_hmms_outputs.py   output tables
"""

import os
import shutil
import argparse

from tqdm import tqdm  # type: ignore

from gtotree.cli.common import (CustomRichHelpFormatter, add_help, add_version_arg,
                                run_subcommand_main)
from gtotree.utils.misc.general import (run_pooled_stage,
                                   prepare_output_dir,
                                   adopt_genome_progress,
                                   read_run_data,
                                   resolve_input_genomes,
                                   wanted_ref_tax_list,
                                   PREPROCESSING_STAGE_BY_SOURCE,
                                   REASON_NOT_FOUND_AT_NCBI,
                                   write_run_data,
                                   GTT_PROGRESS_BAR_FORMAT_INDENTED,
                                   GTT_PROGRESS_BAR_FORMAT_NO_COUNT_INDENTED,
                                   GTT_PROGRESS_SMOOTHING)
from gtotree.utils.misc.stages import (GenomeRemovalStage,
                                       GENOME_REMOVAL_STAGE_ORDER,
                                       PREPROCESSING_REMOVAL_STAGES)
from gtotree.utils.misc.summary_info import write_removed_genomes_report
from gtotree.utils.misc import phase_stats
from gtotree.utils.misc.messaging import (report_message, color_text, spinner,
                                     report_phase_header, REMOVED_GENOMES_FILENAME)
from gtotree.utils.misc.data_locations import ensure_reference_data
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.exclusion_list import exclusion_list_help
from gtotree.utils.taxonomy.wanted_ref_tax import resolved_gtdb_section
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import (
    GenSCGHMMsError,
    DEFAULT_MIN_PFAM_COVERAGE,
    DEFAULT_MAX_SHARED_PROTEIN_PERCENT,
    DEFAULT_MAX_HMMS,
    pfam_data_paths,
    load_coverage_filtered_pfams,
    read_hmm_accessions,
    write_filtered_pfam_hmms,
    count_genomes_sharing,
)
from gtotree.utils.misc.resume_state import (ResumeProfile, hash_strings,
                                             hash_local_genomes, load_sidecar,
                                             save_sidecar,
                                             hash_file_contents, STATE_VERSION)
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import (TargetGenomeError,
                                                                  build_run_data,
                                                                  resolve_download_info,
                                                                  fetch_amino_acids_pooled,
                                                                  relabel_and_append,
                                                                  MAX_DOWNLOAD_THREADS)
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_local import process_local_genome
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_search import search_profiles
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_outputs as outputs


DEFAULT_OUTPUT_DIR = "gtt-gen-scg-hmms-output"
DEFAULT_PERCENT_SINGLE_COPY = 90

# Default thread count (only used for the hmmsearch stage in here)
DEFAULT_THREADS = 8

# low number of genomes message threshold
FEW_GENOMES_THRESHOLD = 10


################################################################################
# resume
################################################################################

# stage names, in pipeline order
STAGE_GENOMES = "genomes"
STAGE_PFAMS = "pfams"
STAGE_SEARCH = "search"

STAGE_ORDER = [STAGE_GENOMES, STAGE_PFAMS, STAGE_SEARCH]

RESUME = ResumeProfile(
    name="gen-scg-hmms",
    stages=STAGE_ORDER,
    field_labels={
        "state_version": "the run-state format",
        "accessions_sha256": "the set of target genomes",
        "local_genomes_sha256": "the local genome files (contents, paths, or set)",
        "min_pfam_coverage": "--min-pfam-coverage",
        "source": "--source",
        "ncbi_section": "--ncbi-section",
        "gtdb_section": "--gtdb-section",
        "wanted_ref_tax": "--wanted-ref-tax",
        "target_rank": "--target-rank",
        "target_domain": "--target-domain",
        "derep_rank": "--derep-rank",
        "min_completeness": "--min-completeness",
        "max_contamination": "--max-contamination",
        "exclusion_list_sha256": "the exclusion list (--exclusion-list)",
        "pfam_version": "the Pfam version",
    },
    # the Pfam version isn't resolved until the Pfam stage runs, so a run interrupted
    # before then legitimately has None stored and shouldn't be refused on that basis
    deferred_fields=("pfam_version",),
)


def build_fingerprint(run_data, args, pfam_version=None):
    """
    Everything that affects the genomes and search results a resume would reuse.

    Deliberately does NOT include `--num-jobs`, `--num-threads`, `--keep-working-dir`,
    or the output directory name: those change how the run executes, not what it
    produces, so changing them shouldn't invalidate a resume.

    Nor the selection parameters (`-p`, `--max-hmms`, `--max-shared-protein-percent`):
    they only act on the search results, in the final phase, which a resume always
    reruns anyway, so a resume can take new values for them.

    Local genomes are hashed with size and mtime as well as path, because unlike an
    NCBI accession a local file's contents can change while its path stays the same,
    and resuming across that would silently mix old results with new input.
    """
    accessions = run_data.get_input_ncbi_accs()
    local_genomes = (list(run_data.genbank_files) + list(run_data.fasta_files)
                     + list(run_data.amino_acid_files))

    return {
        "state_version": STATE_VERSION,
        "accessions_sha256": hash_strings(accessions),
        "num_accessions": len(set(accessions)),
        "local_genomes_sha256": hash_local_genomes(local_genomes),
        "num_local_genomes": len(local_genomes),
        "min_pfam_coverage": args.min_pfam_coverage,
        "source": (args.source or "").upper(),
        "ncbi_section": (getattr(args, "ncbi_section", None) or "").lower(),
        # the RESOLVED value, so a defaulted run and an explicit `--gtdb-section reps`
        # fingerprint identically -- they select identically
        "gtdb_section": resolved_gtdb_section(getattr(args, "gtdb_section", None)),
        "wanted_ref_tax": (sorted(wanted_ref_tax_list(args)) or None),
        "target_rank": args.target_rank,
        "target_domain": getattr(args, "target_domain", None),
        "derep_rank": args.derep_rank,
        "min_completeness": args.min_completeness,
        "max_contamination": args.max_contamination,
        # hashed by CONTENTS, not path: moving the file shouldn't force a re-run,
        # and editing it in place must, since it changes which genomes get selected
        "exclusion_list_sha256": hash_file_contents(
            getattr(args, "exclusion_list", None)),
        "pfam_version": pfam_version,
    }


################################################################################
# parser
################################################################################

def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to generate a new single-copy gene (SCG) HMM set "
            "for use with GToTree. It takes a set of target genomes as various possible inputs "
            "and/or selected by taxonomy with `--wanted-ref-tax`, finds which Pfam "
            "profiles are present in exactly one copy in a specified percentage of them, "
            "and writes those profiles out as a new SCG-HMM set. See the wiki for "
            "details on the process: github.com/AstrobioMike/GToTree/wiki/SCG-sets")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "gen-scg-hmms",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt gen-scg-hmms -w Nitrospirota --derep-rank genus`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters (at least one)")
    selection = parser.add_argument_group("Taxon-selection Parameters (used with `-w`)")
    optional = parser.add_argument_group("General Parameters")

    required.add_argument(
        "-w", "--wanted-ref-tax",
        metavar="<STR>",
        help=("A target taxon whose reference genomes should be used (e.g., "
              "'Nitrospirota'). May be given more than once to pool several taxa into "
              "one set (e.g. `-w Bacteria -w Archaea`); each is resolved and "
              "dereplicated on its own, then merged."),
        action="append",
    )

    required.add_argument(
        "-a", "--ncbi-accessions",
        metavar="<FILE>",
        help=("A single-column file of NCBI accessions to build the SCG set from."),
        action="store",
    )

    required.add_argument(
        "-g", "--genbank-files",
        metavar="<FILE>",
        help=("A single-column file listing GenBank files to include (the same input "
              "style as the main GToTree program)."),
        action="store",
    )

    required.add_argument(
        "-f", "--fasta-files",
        metavar="<FILE>",
        help=("A single-column file listing nucleotide fasta files to include; genes "
              "will be called with prodigal."),
        action="store",
    )

    required.add_argument(
        "-A", "--amino-acid-files",
        metavar="<FILE>",
        help=("A single-column file listing amino-acid fasta files to include; each "
              "should hold the proteins for one genome."),
        action="store",
    )

    required.add_argument(
        "--from-run",
        metavar="<DIR>",
        dest="from_run",
        default=None,
        help=("The output directory of a finished `gen-scg-hmms` run in order to build "
              "a new set with different `-p` and/or `--max-hmms`, reusing its genomes and "
              "search results rather than redoing them. The new set is written to `-o`, "
              "and the original directory isn't modified. This is incompatible with other "
              "required parameters here."),
        action="store",
    )

    selection.add_argument(
        "--source",
        type=str.lower,
        default="gtdb",
        choices=["gtdb", "ncbi"],
        help=("Which taxonomy source to select genomes from (default: gtdb)"),
        action="store",
    )

    selection.add_argument(
        "--ncbi-section",
        dest="ncbi_section",
        type=str.lower,
        default="genbank",
        choices=["refseq", "genbank", "both"],
        help=("Which section of NCBI to draw genomes from (default: genbank). Ignored with `--source gtdb`."),
        action="store",
    )

    selection.add_argument(
        "--gtdb-section",
        dest="gtdb_section",
        type=str.lower,
        default=None,
        choices=["reps", "all"],
        help=("Which section of GTDB to draw genomes from: 'reps' "
              "(GTDB species representatives) or 'all' (every genome under the "
              "requested taxon); default: reps. Ignored with `--source ncbi`."),
        action="store",
    )

    selection.add_argument(
        "--target-rank",
        type=str.lower,
        choices=list(RANKS),
        help=("Target rank, if needed to disambiguate a taxon name that exists at "
              "multiple ranks"),
        action="store",
    )

    selection.add_argument(
        "--target-domain",
        type=str,
        dest="target_domain",
        default=None,
        help=("Target domain, if needed to disambiguate a taxon name that exists in "
              "multiple domains (e.g., bacillus is both a bacterial and eukaryotic genus)"),
        action="store",
    )

    selection.add_argument(
        "--derep-rank",
        default="off",
        type=str.lower,
        choices=["auto", "off"] + list(RANKS),
        help=("Dereplicate the selected genomes down to one per unique value of this "
              "rank (default: off). Use 'auto' for two ranks finer than the target."),
        action="store",
    )

    selection.add_argument(
        "--min-completeness",
        metavar="<FLOAT>",
        default=None,
        type=float,
        help=("Minimum estimated completeness for a genome to be "
              "eligible, between 0 and 100 (default: None)."),
        action="store",
    )

    selection.add_argument(
        "--max-contamination",
        metavar="<FLOAT>",
        default=None,
        type=float,
        help=("Maximum estimated contamination for a genome to be "
              "eligible (default: None)."),
        action="store",
    )

    selection.add_argument(
        "--exclusion-list",
        metavar="<FILE>",
        dest="exclusion_list",
        default=None,
        help=exclusion_list_help("-w"),
        action="store",
    )

    optional.add_argument(
        "-o", "--output-dir",
        metavar="<DIR>",
        default=DEFAULT_OUTPUT_DIR,
        dest="output_dir",
        help=(f'Desired output directory (default: "{DEFAULT_OUTPUT_DIR}")'),
        action="store",
    )

    optional.add_argument(
        "-p", "--percent-single-copy",
        metavar="<INT>",
        default=DEFAULT_PERCENT_SINGLE_COPY,
        type=int,
        help=("The percent of target genomes that need to have exactly 1 copy of a "
              f"Pfam for it to be included, between 1 and 100 (default: "
              f"{DEFAULT_PERCENT_SINGLE_COPY})"),
        action="store",
    )

    optional.add_argument(
        "--max-hmms",
        metavar="<INT>",
        dest="max_hmms",
        default=DEFAULT_MAX_HMMS,
        type=int,
        help=("The maximum number of SCG-HMMs to keep. If more than this qualify, those "
              "with exactly 1 copy in the most genomes are kept (ties go to the Pfam with "
              "higher average coverage). Set to 0 if you want no max limit (default: "
              f"{DEFAULT_MAX_HMMS})"),
        action="store",
    )

    optional.add_argument(
        "-j", "--num-jobs",
        metavar="<INT>",
        default=10,
        type=int,
        help=("Number of concurrent jobs like downloading/processing genomes (default: 10)"),
        action="store",
    )

    optional.add_argument(
        "-t", "--num-threads",
        metavar="<INT>",
        default=DEFAULT_THREADS,
        type=int,
        help=(f"Number of threads to use during hmmsearch (default: {DEFAULT_THREADS}; this does NOT multiply with `--num-jobs`)"),
        action="store",
    )

    optional.add_argument(
        "--min-pfam-coverage",
        metavar="<FLOAT>",
        default=DEFAULT_MIN_PFAM_COVERAGE,
        type=float,
        help=("Minimum average coverage (%%) of the underlying proteins for a Pfam profile "
              f"to be considered (default: {DEFAULT_MIN_PFAM_COVERAGE:g})"),
        action="store",
    )

    optional.add_argument(
        "--max-shared-protein-percent",
        metavar="<FLOAT>",
        default=DEFAULT_MAX_SHARED_PROTEIN_PERCENT,
        type=float,
        help=argparse.SUPPRESS,
        # help=("If two Pfams that pass the single-copy cutoff hit the same protein in at "
        #       "least this percent of genomes, only one of them is kept, since they'd "
        #       "otherwise pull the same gene twice. "
        #       "Between 0 and 100 (default: "
        #       f"{DEFAULT_MAX_SHARED_PROTEIN_PERCENT:g})"),
        action="store",
    )

    optional.add_argument(
        "--keep-working-dir",
        action="store_true",
        help=("Keep the intermediate working directory"),
    )

    optional.add_argument(
        "-R",
        "--resume",
        action="store_true",
        help=("Resume a previous run in the same output directory, reusing any stages "
              "that already completed. Refuses if the run parameters changed."),
    )

    optional.add_argument(
        "-F", "--force-overwrite",
        action="store_true",
        help=("Overwrite the output directory if it already exists"),
    )

    add_help(optional)
    add_version_arg(optional)

    parser.set_defaults(func="gen_scg_hmms")

    return parser


################################################################################
# terminal styling (matching bit's gen-mg)
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
    report_phase_header(title)


################################################################################
# validation / setup
################################################################################

# what `--from-run` can't be combined with, as (dest, flag): each of these decides which
# genomes are used or what they're searched against, and a `--from-run` reuses both.
# `-j`, `-t`, and `--keep-working-dir` aren't here: they only change how a scan runs,
# so they're simply not used rather than refused
INCOMPATIBLE_WITH_FROM_RUN = (
    ("wanted_ref_tax", "-w"),
    ("ncbi_accessions", "-a"),
    ("genbank_files", "-g"),
    ("fasta_files", "-f"),
    ("amino_acid_files", "-A"),
    ("source", "--source"),
    ("ncbi_section", "--ncbi-section"),
    ("gtdb_section", "--gtdb-section"),
    ("target_rank", "--target-rank"),
    ("target_domain", "--target-domain"),
    ("derep_rank", "--derep-rank"),
    ("min_completeness", "--min-completeness"),
    ("max_contamination", "--max-contamination"),
    ("exclusion_list", "--exclusion-list"),
    ("min_pfam_coverage", "--min-pfam-coverage"),
    ("resume", "-R/--resume"),
)


def _from_run_conflicts(args):
    """
    The INCOMPATIBLE_WITH_FROM_RUN flags given a non-default value, as their flag labels.

    Compared against the parser's own defaults so there's one place they're defined;
    an incompatible flag explicitly given its default value has no effect either way.
    """
    defaults = build_parser().parse_args([])
    return [flag for dest, flag in INCOMPATIBLE_WITH_FROM_RUN
            if getattr(args, dest, None) != getattr(defaults, dest, None)]


def _check_from_run_args(args):
    conflicts = _from_run_conflicts(args)
    if conflicts:
        raise GenSCGHMMsError(
            "`--from-run` reuses the genomes and search results of the previous run, "
            f"so it can't be combined with: {', '.join(conflicts)}. With `--from-run`, "
            "`-p` and `--max-hmms` set how the new set is selected.")

    out_dir = os.path.realpath(args.output_dir.rstrip("/"))
    if out_dir == os.path.realpath(args.from_run.rstrip("/")):
        raise GenSCGHMMsError(
            "`--from-run` writes a new set to a new directory and leaves the previous "
            "run as it is, so `-o` needs to be a different directory than the "
            "`--from-run` one.")


def check_args(args):
    """Validate arguments, raising GenSCGHMMsError with a friendly message."""
    if getattr(args, "max_hmms", DEFAULT_MAX_HMMS) < 0:
        raise GenSCGHMMsError(
            "The `--max-hmms` parameter can't be negative (0 means no cap).")

    if getattr(args, "from_run", None):
        _check_from_run_args(args)

    input_flags = (args.ncbi_accessions, args.wanted_ref_tax, args.genbank_files,
                   args.fasta_files, args.amino_acid_files,
                   getattr(args, "from_run", None))
    if not any(input_flags):
        raise GenSCGHMMsError(
            "We need some target genomes to work with! Provide any combination of an "
            "accessions file (`-a`), GenBank files (`-g`), fasta files (`-f`), "
            "amino-acid files (`-A`), and/or a target taxon (`-w`). Or a previous "
            "run to build a new set from (with `--from-run`).")

    if not 0 < args.percent_single_copy <= 100:
        raise GenSCGHMMsError(
            "The `--percent-single-copy` (-p) parameter needs to be between 1 and 100.")

    max_shared = getattr(args, "max_shared_protein_percent",
                         DEFAULT_MAX_SHARED_PROTEIN_PERCENT)
    if not 0 <= max_shared <= 100:
        raise GenSCGHMMsError(
            "The `--max-shared-protein-percent` parameter needs to be between 0 and 100.")

    if args.num_threads < 1:
        raise GenSCGHMMsError("The `--num-threads` (-t) parameter needs to be at least 1.")

    if args.num_jobs < 1:
        raise GenSCGHMMsError("The `--num-jobs` (-n) parameter needs to be at least 1.")

    if args.min_completeness is not None and not 0 <= args.min_completeness <= 100:
        raise GenSCGHMMsError(
            "The `--min-completeness` parameter needs to be between 0 and 100.")

    if args.max_contamination is not None and args.max_contamination < 0:
        raise GenSCGHMMsError(
            "The `--max-contamination` parameter can't be negative.")

    if (args.min_completeness is not None or args.max_contamination is not None) \
            and not args.wanted_ref_tax:
        raise GenSCGHMMsError(
            "`--min-completeness` / `--max-contamination` filter the genomes selected "
            "by `--wanted-ref-tax` (-w), so they need `-w` to act on. Genomes given "
            "directly with `-a`, `-g`, `-f`, or `-A` are always used as provided.")

    exclusion_list = getattr(args, "exclusion_list", None)
    if exclusion_list:
        if not args.wanted_ref_tax:
            raise GenSCGHMMsError(
                "An `--exclusion-list` only has an effect alongside "
                "`--wanted-ref-tax` (-w), since it removes genomes from what `-w` "
                "selects. Genomes given directly with `-a`, `-g`, `-f`, or `-A` are "
                "always used as provided, so nothing would be excluded.")
        if not os.path.isfile(exclusion_list):
            raise GenSCGHMMsError(
                f"The specified exclusion list '{exclusion_list}' can't be found.")

    return args


def setup_output_dir(args):
    """
    Create the output dir and working dir, honoring -F and -R.

    -F and -R are mutually exclusive: one throws the previous run away, the other
    depends on it being intact.
    """
    out_dir = args.output_dir.rstrip("/")
    resume = getattr(args, "resume", False)

    if resume and args.force_overwrite:
        raise GenSCGHMMsError(
            "`--resume` and `-F`/`--force-overwrite` can't be used together, one "
            "reuses the previous run and the other deletes it.")

    if resume and os.path.isfile(os.path.join(out_dir, outputs.HMM_INFO_FILENAME)):
        report_message(
            f"The run in '{out_dir}' already finished, so there's nothing to "
            "resume. To build a set from it with a different `-p` or `--max-hmms`, "
            f"use `--from-run {out_dir}` with a new `-o`. Or use `-F` to rebuild this "
            "from scratch.\n", "yellow")
        exit(0)

    return prepare_output_dir(out_dir, resume=resume,
                              force_overwrite=args.force_overwrite)


################################################################################
# phases
################################################################################

def phase_resolve_genomes(args, run_data):
    """
    Resolve the target genomes from `-a`, `-w`, and the local-file flags
    """
    ensure_reference_data(has_ncbi_accessions=bool(args.ncbi_accessions),
                          wanted_ref_tax=args.wanted_ref_tax,
                          source=args.source)

    return resolve_input_genomes(args, run_data, GenSCGHMMsError)


def phase_get_amino_acids(run_data, work_dir, args, to_fetch, download_info):
    """
    Get amino acids for every target genome, downloading the NCBI accessions phase 1
    resolved and processing any local files, and combine them into a single fasta with
    genome-traceable headers
    """
    phase_stats.checkpoint("phase 2 start")

    combined_path = os.path.join(work_dir, "all-target-proteins.faa")

    local_genomes = [gd for gd in _local_genomes(run_data) if not gd.removed]

    if not to_fetch and not local_genomes:
        _report_removals(run_data, "None of the target genomes could be resolved to "
                                   "something usable.")

    download_workers = max(1, min(int(args.num_jobs), MAX_DOWNLOAD_THREADS,
                                  max(len(to_fetch), 1)))
    local_workers = max(1, min(int(args.num_jobs), max(len(local_genomes), 1)))

    with open(combined_path, "w") as combined:

        def absorb(gd, aa_path, error):
            """Fold one finished genome into the combined fasta and clean up."""
            if error is not None:
                gd.mark_removed(error, PREPROCESSING_STAGE_BY_SOURCE[gd.source])
            else:
                try:
                    relabel_and_append(gd.id, aa_path, combined)
                    gd.mark_processing_done()
                except TargetGenomeError as e:
                    gd.mark_removed(str(e), PREPROCESSING_STAGE_BY_SOURCE[gd.source])
            if aa_path and os.path.exists(aa_path):
                try:
                    os.remove(aa_path)
                except OSError:
                    pass

        if to_fetch:
            by_id = {gd.id: gd for gd in to_fetch}

            def on_result(acc, aa_path, used_prodigal, error):
                gd = by_id[acc]
                if used_prodigal:
                    gd.mark_prodigal_used()
                    run_data.tools_used.prodigal_used = True
                if error is not None:
                    gd.acc_was_downloaded = False
                absorb(gd, aa_path, error)

            _print_pool_label("Processing", "NCBI-derived", download_workers,
                              lead_newline=False)

            fetch_amino_acids_pooled(
                [gd.id for gd in to_fetch], download_info, work_dir, args=args,
                on_result=on_result,
                bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED,
                lead_newline=False)

        phase_stats.checkpoint("phase 2 after downloads + prodigal + relabel")

        if local_genomes:

            def local_worker(gd, rd):
                """
                Runs in a worker thread. Touches only per-genome paths and returns a
                plain status dict -- never raises, and never appends to the combined
                fasta (that happens in local_apply, on the main thread).
                """
                try:
                    aa_path, used_prodigal = process_local_genome(gd, rd["work_dir"])
                    return {"aa_path": aa_path, "prodigal": used_prodigal, "error": None}
                except BaseException as e:  # noqa: BLE001 - reported per-genome
                    return {"aa_path": None, "prodigal": False, "error": str(e)}

            def local_apply(gd, status, rd):
                if status["prodigal"]:
                    gd.mark_prodigal_used()
                    run_data.tools_used.prodigal_used = True
                absorb(gd, status["aa_path"], status["error"])

            _print_pool_label("Processing", "local", local_workers,
                              lead_newline=bool(to_fetch))

            run_pooled_stage(local_genomes, local_worker, local_apply, args,
                             {"work_dir": work_dir},
                             bar_format=GTT_PROGRESS_BAR_FORMAT_INDENTED,
                             lead_newline=False)

    write_removed_genomes_report(run_data)
    write_run_data(run_data)

    kept_ids = [gd.id for gd in run_data.all_input_genomes
                if gd.processing_done and not gd.removed]

    if not kept_ids:
        _report_removals(run_data, "No amino-acid sequences could be obtained for any "
                                   "of the target genomes.")

    _report_removals(run_data)

    print(f"\n      {color_text(f'{len(kept_ids):,} genome(s) ready to search', 'green')}")

    return combined_path, kept_ids


def _print_pool_label(action, kind, workers, lead_newline=True):
    """
    Print the header line that sits above a pooled-stage progress bar
    """
    label = f"      {action} {kind} genome(s)"
    if workers > 1:
        label += f" with {workers} concurrent job(s)"
    lead = "\n" if lead_newline else ""
    print(f"{lead}{label}:")


def _local_genomes(run_data):
    """The three local-file input pools, in input order."""
    return (list(run_data.genbank_files) + list(run_data.fasta_files)
            + list(run_data.amino_acid_files))


def _resolve_accession_downloads(run_data):
    """
    Look every accession up in the NCBI assembly table, marking the ones NCBI doesn't
    list as removed at the lookup stage.

    Returns (to_fetch, download_info): the accession GenomeData still in play, and the
    accession -> download-location map the fetch stage takes.

    Runs in phase 1, the way `gtt search-annotations` does it: an accession
    NCBI no longer lists is a fact about the *input set*, so it belongs next to the
    input-source listing rather than mixed in with the fetch stage's own failures.
    """
    accs = [gd for gd in run_data.ncbi_accs if not gd.removed]
    if not accs:
        return [], {}

    # the input-source listing above doesn't end with a blank line, so this step owns
    # the one that separates it
    print()

    with spinner("Looking up accessions in the NCBI assembly data...",
                 "Looked up accessions at NCBI"):
        info, not_found = resolve_download_info([gd.id for gd in accs])

    to_fetch = []
    for gd in accs:
        entry = info.get(gd.id)
        if entry is None:
            gd.acc_was_found = False
            gd.mark_removed(REASON_NOT_FOUND_AT_NCBI, GenomeRemovalStage.NCBI_LOOKUP)
            continue
        gd.acc_was_found = True
        gd.organism_name = entry.get("organism_name")
        to_fetch.append(gd)

    if not_found:
        report_message(f"{len(not_found):,} accession(s) not found at NCBI. Reported in:",
                       "yellow", ii="      ", si="      ")
        print(f"        {color_text(removed_genomes_path(run_data), 'yellow')}")

    return to_fetch, info


def _report_removals(run_data, fatal_message=None):
    """
    Say how many target genomes this phase lost, and point at the file explaining each
    """
    removed = run_data.genomes_removed_at(*PROCESSING_REMOVAL_STAGES)

    if removed:
        by_stage = _removals_by_stage(removed)

        if len(by_stage) == 1:
            report_message(f"{len(removed):,} genome(s) "
                           f"{_removal_stage_label(by_stage[0][0])}. Reported in:",
                           "yellow", ii="      ", si="      ")
            print(f"        {color_text(removed_genomes_path(run_data), 'yellow')}")
        else:
            report_message(f"{len(removed):,} genome(s) didn't make it through "
                           "processing:", "yellow", ii="      ", si="      ")
            for stage, count in by_stage:
                print(f"        "
                      f"{color_text(f'{count:,} {_removal_stage_label(stage)}', 'yellow')}")
            print(f"\n      {color_text('Each is reported in:', 'yellow')}")
            print(f"        {color_text(removed_genomes_path(run_data), 'yellow')}")

    if fatal_message:
        if any(gd.removed for gd in run_data.all_input_genomes):
            fatal_message += (f" The reason for each is in "
                              f"{removed_genomes_path(run_data)}.")
        raise GenSCGHMMsError(fatal_message)


PROCESSING_REMOVAL_STAGES = tuple(stage for stage in PREPROCESSING_REMOVAL_STAGES
                                  if stage != GenomeRemovalStage.NCBI_LOOKUP)


# how each removal stage this phase can reach reads in the end-of-phase summary
REMOVAL_STAGE_LABELS = {
    GenomeRemovalStage.NCBI_LOOKUP: "not found at NCBI",
    GenomeRemovalStage.NCBI_DOWNLOAD: "failed downloading or gene-calling",
    GenomeRemovalStage.GENBANK_PREP: "failed GenBank-file processing",
    GenomeRemovalStage.FASTA_PREP: "failed fasta-file processing",
    GenomeRemovalStage.AMINO_ACID_PREP: "failed amino-acid-file processing",
}


def _removal_stage_label(stage):
    """
    Reads correctly both after a count ("1 not found at NCBI") and as the tail of the
    collapsed one-cause line ("3 of 11 target genome(s) not found at NCBI")
    """
    return REMOVAL_STAGE_LABELS.get(stage, f"removed at the {stage} stage")


def _removals_by_stage(removed):
    """
    (stage, count) pairs in pipeline order, skipping stages nothing was lost at
    """
    counts = {}
    for gd in removed:
        counts[gd.removed_at] = counts.get(gd.removed_at, 0) + 1

    ordered = [(stage, counts.pop(stage)) for stage in GENOME_REMOVAL_STAGE_ORDER
               if stage in counts]

    # anything with a stage outside the known order still gets reported rather than
    # silently dropped from a total that counts it
    ordered.extend(sorted(counts.items()))

    return ordered


def removed_genomes_path(run_data):
    return f"{run_data.run_files_dir_rel}/{REMOVED_GENOMES_FILENAME}"


def phase_filter_pfams(work_dir, args, state=None, resuming=False):
    """
    Load the managed Pfam data and write out the coverage-filtered profile subset.

    Returns (filtered_hmm_path, pfam_info, filtered_accs, pfam_version).
    """
    from gtotree.utils.pfam.get_pfam_data import get_pfam_data, get_stored_pfam_version

    # ensures the managed Pfam asset is present (downloads it once if not)
    pfam_data_dir = get_pfam_data()
    master_hmm_path, info_path = pfam_data_paths(pfam_data_dir)
    pfam_version = get_stored_pfam_version(pfam_data_dir) or "NA"

    print(f"      Pfam version being used: {color_text(pfam_version, 'green')}\n")

    with spinner("Filtering Pfams by average coverage...", "Filtered Pfams by coverage"):
        pfam_info, total_profiles = load_coverage_filtered_pfams(
            info_path, min_coverage=args.min_pfam_coverage)

    print(f"        {len(pfam_info):,} profile(s) with average coverage > "
          f"{args.min_pfam_coverage}")

    filtered_hmm_path = os.path.join(work_dir, "coverage-filtered-pfams.hmm")

    # the extraction is a ~2 minute streaming pass over the 2 GB master HMM, so it's
    # well worth reusing when the previous run already produced it
    if resuming and RESUME.is_reusable(state, STAGE_PFAMS, work_dir):
        with spinner("Reusing previously filtered Pfam profiles...",
                     "Reused filtered Pfam profiles"):
            found = read_hmm_accessions(filtered_hmm_path)
    else:
        print("\n      Extracting target profiles:")
        with tqdm(total=total_profiles,
                  bar_format=GTT_PROGRESS_BAR_FORMAT_NO_COUNT_INDENTED,
                  ncols=76, smoothing=GTT_PROGRESS_SMOOTHING) as pbar:
            found = write_filtered_pfam_hmms(
                master_hmm_path, set(pfam_info), filtered_hmm_path,
                progress_callback=pbar.update)

    # keep ordering stable and aligned with what was actually written
    filtered_accs = found

    return filtered_hmm_path, pfam_info, filtered_accs, pfam_version


def phase_search(filtered_hmm_path, combined_path, num_genomes, args,
                 work_dir=None, resuming=False):
    """
    Run the hmmsearch stage with a progress bar over genomes

    Returns (hits_by_genome, shared_by_genome); see `search_profiles`
    """
    checkpoint_path = (os.path.join(work_dir, SEARCH_CHECKPOINT_FILENAME)
                       if work_dir else None)

    shared_by_genome = {}
    with tqdm(total=num_genomes, desc="      Progress", ncols=74,
              unit=" genome", smoothing=GTT_PROGRESS_SMOOTHING) as pbar:
        hits_by_genome = search_profiles(
            filtered_hmm_path, combined_path,
            threads=args.num_threads,
            progress_callback=pbar.update,
            checkpoint_path=checkpoint_path,
            resume=resuming,
            shared_proteins=shared_by_genome)
    return hits_by_genome, shared_by_genome


def phase_write_scan_record(out_dir, hits_by_genome, kept_ids, filtered_accs,
                            genomes_sharing, pfam_info, pfam_version, run_data):
    """
    Write the tables recording what the scan found, which are also everything
    `--from-run` needs to select a new set from this run later.
    """
    per_genome_counts = {genome_id: dict(hits_by_genome.get(genome_id, {}))
                         for genome_id in kept_ids}

    with spinner("Writing search-result tables...", "Wrote search-result tables"):
        outputs.write_hit_counts(out_dir, kept_ids, filtered_accs, per_genome_counts)
        outputs.write_target_genomes(out_dir, kept_ids, run_data)
        outputs.write_shared_protein_pairs(out_dir, genomes_sharing, len(kept_ids),
                                           pfam_info)
        outputs.write_pfam_version(out_dir, pfam_version)


def phase_select_and_write(out_dir, source_hmm_path, hits_by_genome, genome_ids,
                           filtered_accs, genomes_sharing, pfam_info, args,
                           from_run=None, source_total_profiles=None):
    """
    Select the SCG targets, extract their profiles, and write the selection outputs.
    Shared by a full run and `--from-run`.

    `source_hmm_path` is where the profiles are extracted from: the coverage-filtered
    subset in a full run, or the master Pfam HMM with `--from-run` (whose working dir,
    holding the subset, is usually long gone). `source_total_profiles`, if given, sizes
    a progress bar for that extraction, as streaming the master file takes a while.

    Returns (final_hmm_path, num_targets).
    """
    from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import select_scg_targets

    max_shared = getattr(args, "max_shared_protein_percent",
                         DEFAULT_MAX_SHARED_PROTEIN_PERCENT)
    max_hmms = getattr(args, "max_hmms", DEFAULT_MAX_HMMS)

    with spinner("Determining single-copy genes...", "Determined single-copy genes"):
        selection = select_scg_targets(
            hits_by_genome, genome_ids, filtered_accs, genomes_sharing, pfam_info,
            args.percent_single_copy, max_shared_percent=max_shared, max_hmms=max_hmms)

    if not selection.targets:
        raise GenSCGHMMsError(
            f"No Pfams were found in exactly one copy in >= {args.percent_single_copy}% "
            "of the target genomes, so there's no SCG set to write. You could try "
            "lowering `-p`, or check that the target genomes are as closely related as "
            "intended.")

    print(f"        {selection.num_single_copy:,} Pfam(s) present in exactly one copy in >= "
          f"{args.percent_single_copy}% of the {len(genome_ids):,} genome(s)")
    if selection.shared_exclusions:
        print(f"        {len(selection.shared_exclusions):,} of those dropped for commonly "
              "hitting the same protein as another retained Pfam")
        print(f"          (see {outputs.SHARED_PROTEIN_EXCLUSIONS_FILENAME})")
    if selection.capped_out:
        print(f"        {len(selection.capped_out):,} more left out by the `--max-hmms` cap "
              f"of {max_hmms:,}, keeping those single-copy in the most genomes")
    if selection.shared_exclusions or selection.capped_out:
        print(f"        {len(selection.targets):,} Pfam(s) retained as SCG targets")

    hmm_filename = outputs.default_hmm_filename(out_dir, len(selection.targets))
    final_hmm_path = os.path.join(out_dir, hmm_filename)

    if source_total_profiles:
        print("\n      Extracting target profiles:")
        with tqdm(total=source_total_profiles,
                  bar_format=GTT_PROGRESS_BAR_FORMAT_NO_COUNT_INDENTED,
                  ncols=76, smoothing=GTT_PROGRESS_SMOOTHING) as pbar:
            found = write_filtered_pfam_hmms(source_hmm_path, selection.targets,
                                             final_hmm_path, progress_callback=pbar.update)
    else:
        with spinner("Writing the new SCG-HMM set...", "Wrote the new SCG-HMM set"):
            found = write_filtered_pfam_hmms(source_hmm_path, selection.targets,
                                             final_hmm_path)

    missing = set(selection.targets) - set(found)
    if missing:
        _remove_quietly(final_hmm_path)
        raise GenSCGHMMsError(
            f"{len(missing):,} of the selected profiles weren't found in "
            f"'{source_hmm_path}' (e.g. {sorted(missing)[0]}), so the set couldn't be "
            "written.")

    with spinner("Writing summary tables...", "Wrote summary tables"):
        outputs.write_shared_protein_exclusions(out_dir, selection.shared_exclusions,
                                                pfam_info)
        outputs.write_selection_params(out_dir, [
            ("percent_single_copy", args.percent_single_copy),
            ("max_hmms", max_hmms),
            ("max_shared_protein_percent", f"{max_shared:g}"),
            ("from_run", os.path.abspath(from_run) if from_run else None),
        ])
        # last, since its presence is what marks the run as finished
        outputs.write_scg_targets_info(
            out_dir, selection.targets, pfam_info, selection.single_copy_counts,
            selection.multi_copy_counts, selection.num_genomes)

    return final_hmm_path, len(selection.targets)


def report_finish(out_dir, final_hmm_path, num_targets, num_genomes,
                  removed_report_path=None):
    border = color_text("  " + "-" * 78, "green")
    title = color_text("  " + "SCG-HMM set complete!".center(78), "green")
    print()
    print(border)
    print(title)
    print(border)

    print(f"\n      New SCG-HMM set holding {color_text(f'{num_targets:,}', 'green')} target gene(s), "
          f"built from {color_text(f'{num_genomes:,}', 'green')} genome(s),")
    print(f"      written to:\n        {color_text(final_hmm_path, 'green')}\n")

    print("      Supporting tables written to:")
    print(f"        {color_text(out_dir + '/', 'green')}\n")

    if removed_report_path:
        report_message("Any input genomes that didn't make it through are reported in:",
                       "yellow", ii="      ", si="      ", newline=False)
        print(f"        {color_text(removed_report_path, 'yellow')}\n")
    print()

    # if os.environ.get("GToTree_HMM_dir"):
    #     hmm_path = os.environ["GToTree_HMM_dir"]
    #     wprint("If you'd like to add this new SCG-HMM set to the stored GToTree ones, "
    #            f"you can copy it into:\n  {hmm_path} :)",
    #            ii="      ", si="      ")
    #     print()


SEARCH_STAGE_SIDECAR = "search-hits.json"

# the per-genome shared-protein pair counts from the search; kept as its own sidecar
# so a work dir from before these were recorded is simply seen as lacking them (and
# the search redone) rather than misread
SEARCH_SHARED_SIDECAR = "search-shared-proteins.json"

SEARCH_CHECKPOINT_FILENAME = "search-checkpoint.jsonl"


def _remove_quietly(path):
    """Delete a work-dir file if it's there; its absence is never a problem."""
    try:
        os.remove(path)
    except OSError:
        pass


################################################################################
# driver
################################################################################

def gen_scg_hmms(args):  # pragma: no cover
    args = check_args(args)

    if args.from_run:
        return gen_scg_hmms_from_run(args)

    out_dir, work_dir = setup_output_dir(args)

    resuming = bool(getattr(args, "resume", False))
    state = RESUME.load(work_dir) if resuming else None

    n = _phase_counter()

    run_data = build_run_data(args, out_dir, work_dir)
    previous = read_run_data(run_data.run_data_path) if resuming else None

    section(f"Phase {n()}: Resolving target genomes...")
    run_data, _selections = phase_resolve_genomes(args, run_data)

    adopt_genome_progress(run_data, previous)

    total_genomes = len(run_data.all_input_genomes)
    if total_genomes < FEW_GENOMES_THRESHOLD:
        report_message(
            f"Just so you know, {total_genomes} genomes is on "
            "the low side for this. "
            "The single-copy percentage becomes a weak signal with few genomes, and the "
            "resulting SCG set may be less reliable.", "orange",
            ii="      ", si="      ")

    fingerprint = build_fingerprint(run_data, args, pfam_version=None)

    if resuming and state:
        differences = RESUME.compare(state.get("fingerprint"), fingerprint)
        if differences:
            raise GenSCGHMMsError(RESUME.refusal_message(differences))
        print(f"\n      {color_text('Resuming from the previous run', 'green')}")
    else:
        resuming = False
        state = RESUME.new(fingerprint)
        RESUME.save(work_dir, state)

    reusing_genomes = bool(resuming
                           and RESUME.is_reusable(state, STAGE_GENOMES, work_dir))

    if reusing_genomes:
        to_fetch, download_info = [], {}
    else:
        to_fetch, download_info = _resolve_accession_downloads(run_data)
        phase_stats.checkpoint("phase 1 after resolve_download_info (NCBI parquet)")

    write_run_data(run_data)

    section(f"Phase {n()}: Getting target-genome amino acids...")
    combined_path = os.path.join(work_dir, "all-target-proteins.faa")

    if reusing_genomes:
        kept_ids = [gd.id for gd in run_data.all_input_genomes
                    if gd.processing_done and not gd.removed]
        with spinner("Reusing previously downloaded amino acids...",
                     f"Reused amino acids for {len(kept_ids):,} genome(s)"):
            pass
    else:
        combined_path, kept_ids = phase_get_amino_acids(run_data, work_dir, args,
                                                        to_fetch, download_info)
        RESUME.mark_complete(
            state, STAGE_GENOMES,
            [combined_path, run_data.run_data_path],
            work_dir=work_dir)
        RESUME.save(work_dir, state)

    section(f"Phase {n()}: Preparing Pfam profiles...")
    filtered_hmm_path, pfam_info, filtered_accs, pfam_version = phase_filter_pfams(
        work_dir, args, state=state, resuming=resuming)

    # now that the Pfam version is known, make sure it matches what a resumed run used
    if resuming and state.get("fingerprint", {}).get("pfam_version") not in (None, pfam_version):
        raise GenSCGHMMsError(
            "`--resume` was specified, but the Pfam version changed since the previous "
            f"run ({state['fingerprint']['pfam_version']} -> {pfam_version}). Use a new "
            "output directory with `-o`, or start over with `-F`.")
    state.setdefault("fingerprint", {})["pfam_version"] = pfam_version
    RESUME.mark_complete(state, STAGE_PFAMS, [filtered_hmm_path],
                         work_dir=work_dir)
    RESUME.save(work_dir, state)

    section(f"Phase {n()}: Searching genomes with Pfam profiles...")
    cached_hits = load_sidecar(work_dir, SEARCH_STAGE_SIDECAR) if resuming else None
    # an empty dict is a legitimate "nothing shared", so only None means missing
    cached_shared = load_sidecar(work_dir, SEARCH_SHARED_SIDECAR) if resuming else None
    if resuming and RESUME.is_reusable(state, STAGE_SEARCH, work_dir) \
            and cached_hits and cached_shared is not None:
        with spinner("Reusing previous search results...", "Reused previous search results"):
            hits_by_genome = cached_hits
            shared_by_genome = cached_shared
    else:
        hits_by_genome, shared_by_genome = phase_search(
            filtered_hmm_path, combined_path, len(kept_ids),
            args, work_dir=work_dir, resuming=resuming)
        save_sidecar(work_dir, SEARCH_STAGE_SIDECAR, hits_by_genome)
        save_sidecar(work_dir, SEARCH_SHARED_SIDECAR, shared_by_genome)
        RESUME.mark_complete(
            state, STAGE_SEARCH,
            [os.path.join(work_dir, SEARCH_STAGE_SIDECAR),
             os.path.join(work_dir, SEARCH_SHARED_SIDECAR)],
            work_dir=work_dir)
        RESUME.save(work_dir, state)
        _remove_quietly(os.path.join(work_dir, SEARCH_CHECKPOINT_FILENAME))

    section(f"Phase {n()}: Determining single-copy genes and writing outputs...")
    genomes_sharing = count_genomes_sharing(shared_by_genome, kept_ids)
    phase_write_scan_record(out_dir, hits_by_genome, kept_ids, filtered_accs,
                            genomes_sharing, pfam_info, pfam_version, run_data)
    final_hmm_path, num_targets = phase_select_and_write(
        out_dir, filtered_hmm_path, hits_by_genome, kept_ids, filtered_accs,
        genomes_sharing, pfam_info, args)

    if not args.keep_working_dir:
        shutil.rmtree(work_dir, ignore_errors=True)

    # closes the final phase, so this has to happen before anything reads the table
    phase_stats.finish()
    phase_stats.write_tsv(out_dir)

    removed_report = (removed_genomes_path(run_data)
                      if any(gd.removed for gd in run_data.all_input_genomes) else None)
    report_finish(out_dir, final_hmm_path, num_targets, len(kept_ids), removed_report)


def gen_scg_hmms_from_run(args):  # pragma: no cover
    """
    Build a new SCG-HMM set from a finished run's output tables, with new selection
    parameters, into a new output dir. Nothing is searched: the previous run's hit
    counts and shared-protein pairs are all the selection needs, and the chosen
    profiles are extracted from the managed master Pfam HMM.
    """
    from gtotree.utils.pfam.get_pfam_data import get_pfam_data, get_stored_pfam_version

    from_run = args.from_run.rstrip("/")
    outputs.check_previous_run(from_run)

    # checked here rather than left to prepare_output_dir, whose message suggests `-R`
    if os.path.exists(args.output_dir) and not args.force_overwrite:
        raise GenSCGHMMsError(
            f"The output directory '{args.output_dir}' already exists, and we don't "
            "want to overwrite anything accidentally. Use `-F` to overwrite it, or "
            "specify a different directory with `-o`.")

    n = _phase_counter()

    section(f"Phase {n()}: Reading the previous run...")
    print(f"      Building a new set from: {color_text(from_run + '/', 'green')}\n")
    with spinner("Reading search results...", "Read search results"):
        genome_ids, filtered_accs, hits_by_genome = outputs.read_hit_counts(
            os.path.join(from_run, outputs.HIT_COUNTS_FILENAME))
        genomes_sharing = outputs.read_shared_protein_pairs(
            os.path.join(from_run, outputs.SHARED_PROTEIN_PAIRS_FILENAME))
        with open(os.path.join(from_run, outputs.PFAM_VERSION_FILENAME)) as f:
            previous_pfam_version = f.readline().strip()
    print(f"        {len(genome_ids):,} genome(s) searched with "
          f"{len(filtered_accs):,} Pfam profile(s)")

    section(f"Phase {n()}: Preparing Pfam profiles...")
    pfam_data_dir = get_pfam_data()
    master_hmm_path, info_path = pfam_data_paths(pfam_data_dir)
    pfam_version = get_stored_pfam_version(pfam_data_dir) or "NA"
    if pfam_version != previous_pfam_version:
        raise GenSCGHMMsError(
            f"The previous run was built with Pfam {previous_pfam_version}, but the "
            f"Pfam data available now is {pfam_version}, so its profiles can't be "
            "pulled from here to match. A new set would need a full run.")
    print(f"      Pfam version being used: {color_text(pfam_version, 'green')}\n")

    with spinner("Loading Pfam info...", "Loaded Pfam info"):
        pfam_info, total_profiles = load_coverage_filtered_pfams(info_path,
                                                                 min_coverage=None)

    section(f"Phase {n()}: Determining single-copy genes and writing outputs...")

    # only created once everything above has checked out, so a refusal leaves nothing
    out_dir, work_dir = prepare_output_dir(args.output_dir, resume=False,
                                           force_overwrite=args.force_overwrite)
    # nothing here needs a working dir
    shutil.rmtree(work_dir, ignore_errors=True)

    with spinner("Copying search-result tables...", "Copied search-result tables"):
        copied = outputs.copy_scan_record(from_run, out_dir)
    final_hmm_path, num_targets = phase_select_and_write(
        out_dir, master_hmm_path, hits_by_genome, genome_ids, filtered_accs,
        genomes_sharing, pfam_info, args, from_run=from_run,
        source_total_profiles=total_profiles)

    phase_stats.finish()
    phase_stats.write_tsv(out_dir)

    removed_report = (os.path.join(args.output_dir.rstrip("/"), REMOVED_GENOMES_FILENAME)
                      if REMOVED_GENOMES_FILENAME in copied else None)
    report_finish(out_dir, final_hmm_path, num_targets, len(genome_ids), removed_report)


def main():  # pragma: no cover
    run_subcommand_main(build_parser(), gen_scg_hmms, GenSCGHMMsError)


if __name__ == "__main__":
    main()
