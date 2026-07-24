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
import sys
import json
import shutil
import argparse

from tqdm import tqdm  # type: ignore

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.general import run_pooled_stage
from gtotree.utils.messaging import (report_message, color_text, spinner,
                                     report_very_early_exit, wprint)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import TaxonNotFound, AmbiguousTaxon
from gtotree.utils.taxonomy.wanted_ref_tax import (resolve_wanted_ref_tax_accessions,
                                                   WantedRefTaxError)
from gtotree.utils.hmms.gen_scg_hmms import (GenSCGHMMsError, DEFAULT_MIN_PFAM_COVERAGE,
                                             pfam_data_paths,
                                             load_coverage_filtered_pfams,
                                             read_hmm_accessions,
                                             write_filtered_pfam_hmms)
from gtotree.utils.hmms import gen_scg_hmms_resume as resume
from gtotree.utils.hmms.gen_scg_hmms_genomes import (TargetGenomeError,
                                                     read_accessions_file,
                                                     resolve_download_info,
                                                     fetch_amino_acids_pooled,
                                                     relabel_and_append,
                                                     MAX_DOWNLOAD_THREADS,
                                                     MISSED_NOT_FOUND)
from gtotree.utils.hmms.gen_scg_hmms_local import (build_local_genomes,
                                                   process_local_genome,
                                                   SOURCE_GENBANK, SOURCE_FASTA,
                                                   SOURCE_AMINO_ACID)
from gtotree.utils.hmms.gen_scg_hmms_search import search_profiles
from gtotree.utils.hmms import gen_scg_hmms_outputs as outputs


DEFAULT_OUTPUT_DIR = "gtt-gen-scg-hmms-output"
DEFAULT_PERCENT_SINGLE_COPY = 90

# Default thread count (only used for the hmmsearch stage in here)
DEFAULT_THREADS = 8

# low number of genomes message threshold
FEW_GENOMES_THRESHOLD = 10


################################################################################
# parser
################################################################################

def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to generate a new single-copy gene (SCG) HMM set "
            "for use with GToTree. It takes a set of target genomes as various possible inputs "
            "and/or selected by taxonomy with `--wanted-ref-tax`, finds "
            "which Pfam profiles are present in exactly one copy in most of them, and "
            "writes those profiles out as a new SCG-HMM set. See the wiki for "
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
            epilog="Ex. usage: `gtt gen-scg-hmms -W Nitrospirota --derep-rank genus`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters (at least one)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-a", "--target-accessions",
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
        "-W", "--wanted-ref-tax",
        metavar="<STR>",
        help=("A target taxon whose reference genomes should be used (e.g. "
              "'Nitrospirota')."),
        action="store",
    )

    optional.add_argument(
        "-S", "--source",
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
        default="off",
        choices=["auto", "off"] + list(RANKS),
        help=("Dereplicate the selected genomes down to one per unique value of this "
              "rank (default: off). Use 'auto' for two ranks finer than the target."),
        action="store",
    )

    optional.add_argument(
        "-o", "--output-directory",
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
        help=("Minimum average coverage of the underlying proteins for a Pfam profile "
              f"to be considered (default: {DEFAULT_MIN_PFAM_COVERAGE})"),
        action="store",
    )

    optional.add_argument(
        "--keep-working-dir",
        action="store_true",
        help=("Keep the intermediate working directory"),
    )

    optional.add_argument(
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
    print(color_text(f"\n\n  {title}\n", "yellow"))


def section_border():
    print(color_text("      " + "- " * 34, "yellow"))


################################################################################
# validation / setup
################################################################################

def check_args(args):
    """Validate arguments, raising GenSCGHMMsError with a friendly message."""
    input_flags = (args.target_accessions, args.wanted_ref_tax, args.genbank_files,
                   args.fasta_files, args.amino_acid_files)
    if not any(input_flags):
        raise GenSCGHMMsError(
            "We need some target genomes to work with! Provide any combination of an "
            "accessions file (`-a`), GenBank files (`-g`), fasta files (`-f`), "
            "amino-acid files (`-A`), and/or a target taxon (`-W`).")

    if not 0 < args.percent_single_copy <= 100:
        raise GenSCGHMMsError(
            "The `--percent-single-copy` (-p) parameter needs to be between 1 and 100.")

    if args.num_threads < 1:
        raise GenSCGHMMsError("The `--num-threads` (-t) parameter needs to be at least 1.")

    if args.num_jobs < 1:
        raise GenSCGHMMsError("The `--num-jobs` (-n) parameter needs to be at least 1.")

    return args


def setup_output_dir(args):
    """
    Create the output dir and working dir, honoring -F and --resume.

    -F and --resume are mutually exclusive: one throws the previous run away, the other
    depends on it being intact.
    """
    out_dir = args.output_dir.rstrip("/")
    work_dir = os.path.join(out_dir, "working-dir")

    if getattr(args, "resume", False) and args.force_overwrite:
        raise GenSCGHMMsError(
            "`--resume` and `-F`/`--force-overwrite` can't be used together -- one "
            "reuses the previous run and the other deletes it.")

    if os.path.exists(out_dir):
        if getattr(args, "resume", False):
            if not os.path.isdir(work_dir):
                report_message(
                    "There's no working directory from a previous run to resume from, "
                    "so we'll start fresh.", "yellow")
            os.makedirs(work_dir, exist_ok=True)
            return out_dir, work_dir

        if not args.force_overwrite:
            raise GenSCGHMMsError(
                f"The output directory '{out_dir}' already exists, and we don't want to "
                "overwrite anything accidentally. Use `--resume` to continue that run, "
                "`-F` to overwrite it, or specify a different directory with `-o`.")
        shutil.rmtree(out_dir)

    elif getattr(args, "resume", False):
        report_message(
            f"`--resume` was specified, but '{out_dir}' doesn't exist yet, so we'll "
            "just start fresh.", "yellow")

    os.makedirs(work_dir, exist_ok=True)

    return out_dir, work_dir


################################################################################
# phases
################################################################################

def phase_resolve_genomes(args):
    """
    Resolve the target genomes from `-a` and/or `-W`.

    Returns (accessions, sources) where `sources` maps accession -> where it came
    from, for the target-genomes table.
    """
    accessions = []
    sources = {}

    if args.target_accessions:
        with spinner("Reading target accessions...", "Read target accessions"):
            from_file = read_accessions_file(args.target_accessions)
        for acc in from_file:
            if acc not in sources:
                accessions.append(acc)
                sources[acc] = "input-accessions"
        print(f"        {len(from_file):,} accession(s) read from "
              f"{args.target_accessions}")

    if args.wanted_ref_tax:
        with spinner(f"Selecting reference genomes for '{args.wanted_ref_tax}'...",
                     "Selected reference genomes"):
            selected, selection = resolve_wanted_ref_tax_accessions(
                args.source, args.wanted_ref_tax,
                target_rank=args.target_rank,
                derep_rank=args.derep_rank)

        added = 0
        for acc in selected:
            if acc not in sources:
                accessions.append(acc)
                sources[acc] = f"{args.source.upper()}:{selection.canonical}"
                added += 1

        detail = f"        {len(selected):,} genome(s) selected"
        if selection.resolved_rank:
            detail += f" ({selection.canonical} at rank {selection.resolved_rank})"
        print(detail)
        if selection.effective_derep_rank:
            print(f"        dereplicated to one genome per {selection.effective_derep_rank}")
        if added != len(selected):
            print(f"        {len(selected) - added:,} already present from `-a`")

        for warning in selection.warnings:
            report_message(warning, "orange", ii="        ", si="        ")

    local_genomes, local_missing = build_local_genomes(args)
    if local_genomes or local_missing:
        with spinner("Reading local genome files...", "Read local genome files"):
            pass
        by_source = {}
        for gd in local_genomes:
            by_source[gd.source] = by_source.get(gd.source, 0) + 1
        for source in (SOURCE_GENBANK, SOURCE_FASTA, SOURCE_AMINO_ACID):
            if by_source.get(source):
                print(f"        {by_source[source]:,} {source} file(s)")
        if local_missing:
            report_message(
                f"{len(local_missing):,} listed file(s) couldn't be found and will be "
                "skipped.", "yellow", ii="        ", si="        ")

    total = len(accessions) + len(local_genomes)
    if not total:
        raise GenSCGHMMsError("No target genomes were resolved to work with.")

    print(f"\n      {color_text(f'{total:,} total target genome(s)', 'green')}")

    return accessions, sources, local_genomes, local_missing


def phase_get_amino_acids(accessions, local_genomes, local_missing, work_dir, args):
    """
    Get amino acids for every target genome -- downloading NCBI accessions and
    preprocessing any local files -- and combine them into a single fasta with
    genome-traceable headers.

    Returns (combined_path, kept_ids, missed, organism_names, sources_extra).
    """
    missed = list(local_missing)
    organism_names = {}
    sources_extra = {}
    kept_ids = []
    prodigal_used = 0

    info = {}
    to_fetch = []
    if accessions:
        with spinner("Resolving download locations...", "Resolved download locations"):
            info, not_found = resolve_download_info(accessions)

        missed.extend((acc, MISSED_NOT_FOUND) for acc in not_found)
        if not_found:
            report_message(f"{len(not_found):,} accession(s) were not found in the NCBI "
                           "assembly data.", "yellow", ii="      ", si="      ")

        to_fetch = [acc for acc in accessions if acc in info]
        for acc in to_fetch:
            organism_names[acc] = info[acc].get("organism_name")

    if not to_fetch and not local_genomes:
        raise GenSCGHMMsError(
            "None of the target genomes could be resolved to something usable.")

    for gd in local_genomes:
        sources_extra[gd.id] = gd.source

    combined_path = os.path.join(work_dir, "all-target-proteins.faa")

    workers = max(1, min(int(args.num_jobs), MAX_DOWNLOAD_THREADS, max(len(to_fetch), 1)))
    if to_fetch and workers > 1:
        print(f"      Downloading with {workers} concurrent job(s)")

    print()
    with open(combined_path, "w") as combined:

        def absorb(genome_id, aa_path, used_prodigal, error):
            """Fold one finished genome into the combined fasta and clean up."""
            nonlocal prodigal_used
            if error is not None:
                missed.append((genome_id, error))
            else:
                try:
                    relabel_and_append(genome_id, aa_path, combined)
                    kept_ids.append(genome_id)
                    if used_prodigal:
                        prodigal_used += 1
                except TargetGenomeError as e:
                    missed.append((genome_id, str(e)))
            if aa_path and os.path.exists(aa_path):
                try:
                    os.remove(aa_path)
                except OSError:
                    pass

        if to_fetch:
            fetch_amino_acids_pooled(
                to_fetch, info, work_dir, args=args, on_result=absorb)

        if local_genomes:
            if to_fetch:
                print()

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
                absorb(gd.id, status["aa_path"], status["prodigal"], status["error"])

            run_pooled_stage(local_genomes, local_worker, local_apply, args,
                             {"work_dir": work_dir})

    # keep output ordering tied to the inputs rather than download timing
    order = {acc: i for i, acc in enumerate(to_fetch)}
    for i, gd in enumerate(local_genomes):
        order[gd.id] = len(to_fetch) + i
    kept_ids.sort(key=lambda a: order.get(a, 0))

    print()
    if prodigal_used:
        print(f"      {prodigal_used:,} genome(s) needed gene-calling with prodigal")

    if not kept_ids:
        raise GenSCGHMMsError(
            "No amino-acid sequences could be obtained for any of the target genomes.")

    total_wanted = len(to_fetch) + len(local_genomes)
    if len(kept_ids) != total_wanted:
        report_message(
            f"{total_wanted - len(kept_ids):,} of {total_wanted:,} target genome(s) "
            "didn't make it through.", "yellow", ii="      ", si="      ")

    print(f"\n      {color_text(f'{len(kept_ids):,} genome(s) ready to search', 'green')}")

    return combined_path, kept_ids, missed, organism_names, sources_extra


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

    print(f"      Pfam version being used: {color_text(pfam_version, 'green')}")

    with spinner("Filtering Pfams by average coverage...", "Filtered Pfams by coverage"):
        pfam_info = load_coverage_filtered_pfams(
            info_path, min_coverage=args.min_pfam_coverage)

    print(f"        {len(pfam_info):,} profile(s) with average coverage > "
          f"{args.min_pfam_coverage}")

    filtered_hmm_path = os.path.join(work_dir, "coverage-filtered-pfams.hmm")

    # the extraction is a ~2 minute streaming pass over the 2 GB master HMM, so it's
    # well worth reusing when the previous run already produced it
    if resuming and resume.stage_is_reusable(state, resume.STAGE_PFAMS, work_dir):
        with spinner("Reusing previously filtered Pfam profiles...",
                     "Reused filtered Pfam profiles"):
            found = read_hmm_accessions(filtered_hmm_path)
    else:
        print()
        with tqdm(desc="    Extracting profiles", ncols=78, unit=" profile") as pbar:
            found = write_filtered_pfam_hmms(
                master_hmm_path, set(pfam_info), filtered_hmm_path,
                progress_callback=pbar.update)
        print()

    print(f"      {color_text(f'{len(found):,} profile(s) to search against', 'green')}")

    # keep ordering stable and aligned with what was actually written
    filtered_accs = found

    return filtered_hmm_path, pfam_info, filtered_accs, pfam_version


def phase_search(filtered_hmm_path, combined_path, filtered_accs, args):
    """Run the hmmsearch stage with a progress bar over profiles."""
    print()
    with tqdm(total=len(filtered_accs), desc="    Progress", ncols=78,
              unit=" profile") as pbar:
        hits_by_genome = search_profiles(
            filtered_hmm_path, combined_path,
            threads=args.num_threads,
            progress_callback=pbar.update)
    print()
    return hits_by_genome


def phase_determine_and_write(out_dir, filtered_hmm_path, hits_by_genome, kept_ids,
                              filtered_accs, pfam_info, pfam_version, sources,
                              organism_names, missed, args):
    """Determine the single-copy set, extract it, and write all outputs."""
    from gtotree.utils.hmms.gen_scg_hmms import count_single_copy_hits

    with spinner("Determining single-copy genes...", "Determined single-copy genes"):
        wanted_accs, per_genome_counts = count_single_copy_hits(
            hits_by_genome, kept_ids, filtered_accs, args.percent_single_copy)

    if not wanted_accs:
        raise GenSCGHMMsError(
            f"No Pfams were found in exactly one copy in >= {args.percent_single_copy}% "
            "of the target genomes, so there's no SCG set to write. You could try "
            "lowering `-p`, or check that the target genomes are as closely related as "
            "intended.")

    print(f"        {len(wanted_accs):,} Pfam(s) present in exactly one copy in >= "
          f"{args.percent_single_copy}% of the {len(kept_ids):,} genome(s)")

    hmm_filename = outputs.default_hmm_filename(out_dir, len(wanted_accs))
    final_hmm_path = os.path.join(out_dir, hmm_filename)

    with spinner("Writing the new SCG-HMM set...", "Wrote the new SCG-HMM set"):
        # re-extract from the already-filtered subset rather than the 2 GB master
        write_filtered_pfam_hmms(filtered_hmm_path, wanted_accs, final_hmm_path)

    with spinner("Writing summary tables...", "Wrote summary tables"):
        outputs.write_scg_targets_info(out_dir, wanted_accs, pfam_info)
        outputs.write_hit_counts(out_dir, kept_ids, filtered_accs, per_genome_counts)
        outputs.write_target_genomes(out_dir, kept_ids, sources, organism_names)
        outputs.write_pfam_version(out_dir, pfam_version)
        missed_path = outputs.write_missed_accessions(out_dir, missed)

    return final_hmm_path, len(wanted_accs), missed_path


def report_finish(out_dir, final_hmm_path, num_targets, num_genomes, pfam_version,
                  missed_path, args):
    """Final summary block."""
    print()
    report_message("-" * 78, "green", ii="  ", newline=False)
    report_message(f"SCG-HMM set complete!".center(78), "green", ii="  ",
                   newline=False)
    report_message("-" * 78, "green", ii="  ", newline=False)

    print(f"\n      New SCG-HMM set holding {color_text(f'{num_targets:,}', 'green')} "
          f"target genes (from Pfam v{pfam_version}),")
    print(f"      built from {color_text(f'{num_genomes:,}', 'green')} genome(s), written to:")
    print(f"        {color_text(final_hmm_path, 'green')}\n")

    print(f"      Supporting tables written to:")
    print(f"        {color_text(out_dir + '/', 'green')}\n")

    if missed_path:
        report_message("Any input accessions that didn't make it through are reported in:",
                       "yellow", ii="      ", si="      ")
        print(f"        {color_text(missed_path, 'yellow')}\n")

    if os.environ.get("GToTree_HMM_dir"):
        wprint("If you'd like to add this new SCG-HMM set to the stored GToTree ones, "
               "you can do so with the `gtt-store-SCG-HMMs` program :)",
               ii="      ", si="      ")
        print()


GENOME_STAGE_SIDECAR = "genome-stage.json"
SEARCH_STAGE_SIDECAR = "search-hits.json"


def _save_json(work_dir, filename, payload):
    """Write a stage sidecar atomically."""
    path = os.path.join(work_dir, filename)
    tmp_path = path + ".part"
    try:
        with open(tmp_path, "w") as f:
            json.dump(payload, f)
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise
    return path


def _load_json(work_dir, filename):
    path = os.path.join(work_dir, filename)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


################################################################################
# driver
################################################################################

def gen_scg_hmms(args):  # pragma: no cover
    args = check_args(args)
    out_dir, work_dir = setup_output_dir(args)

    resuming = bool(getattr(args, "resume", False))
    state = resume.load_state(work_dir) if resuming else None

    n = _phase_counter()

    section(f"Phase {n()}: Resolving target genomes...")
    accessions, sources, local_genomes, local_missing = phase_resolve_genomes(args)

    if len(accessions) + len(local_genomes) < FEW_GENOMES_THRESHOLD:
        report_message(
            f"Just so you know, {len(accessions) + len(local_genomes)} genomes is on "
            "the low side for this. "
            "The single-copy percentage becomes a weak signal with few genomes, and the "
            "resulting SCG set may be less reliable.", "orange",
            ii="      ", si="      ")

    # the fingerprint can only be built once the genome set is known; the Pfam version
    # isn't known until later, so it starts as None and is filled in below
    fingerprint = resume.build_fingerprint(
        accessions, args, pfam_version=None, local_genomes=local_genomes)

    if resuming and state:
        differences = resume.compare_fingerprints(state.get("fingerprint"), fingerprint)
        if differences:
            raise GenSCGHMMsError(
                "`--resume` was specified, but this run doesn't match the previous one "
                "in that directory:\n        - " + "\n        - ".join(differences) +
                "\n\n      Resuming would mix results from two different runs. Use a new "
                "output directory with `-o`, or start over with `-F`.")
        print(f"\n      {color_text('Resuming from the previous run', 'green')}")
    else:
        state = resume.new_state(fingerprint)
        resume.save_state(work_dir, state)

    section(f"Phase {n()}: Getting target-genome amino acids...")
    combined_path = os.path.join(work_dir, "all-target-proteins.faa")
    sidecar = _load_json(work_dir, GENOME_STAGE_SIDECAR) if resuming else None

    if resuming and resume.stage_is_reusable(state, resume.STAGE_GENOMES, work_dir) \
            and sidecar:
        kept_ids = sidecar["kept_ids"]
        missed = [tuple(m) for m in sidecar["missed"]]
        organism_names = sidecar["organism_names"]
        sources.update(sidecar.get("sources_extra") or {})
        with spinner("Reusing previously downloaded amino acids...",
                     f"Reused amino acids for {len(kept_ids):,} genome(s)"):
            pass
    else:
        (combined_path, kept_ids, missed, organism_names,
         sources_extra) = phase_get_amino_acids(
            accessions, local_genomes, local_missing, work_dir, args)
        sources.update(sources_extra)
        _save_json(work_dir, GENOME_STAGE_SIDECAR, {
            "kept_ids": kept_ids,
            "missed": [list(m) for m in missed],
            "organism_names": organism_names,
            "sources_extra": sources_extra,
        })
        resume.mark_stage_complete(
            state, resume.STAGE_GENOMES,
            [combined_path, os.path.join(work_dir, GENOME_STAGE_SIDECAR)],
            work_dir=work_dir)
        resume.save_state(work_dir, state)

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
    resume.mark_stage_complete(state, resume.STAGE_PFAMS, [filtered_hmm_path],
                               work_dir=work_dir)
    resume.save_state(work_dir, state)

    section(f"Phase {n()}: Searching genomes for Pfam profiles...")
    cached_hits = _load_json(work_dir, SEARCH_STAGE_SIDECAR) if resuming else None
    if resuming and resume.stage_is_reusable(state, resume.STAGE_SEARCH, work_dir) \
            and cached_hits:
        with spinner("Reusing previous search results...", "Reused previous search results"):
            hits_by_genome = cached_hits
    else:
        hits_by_genome = phase_search(filtered_hmm_path, combined_path, filtered_accs, args)
        _save_json(work_dir, SEARCH_STAGE_SIDECAR, hits_by_genome)
        resume.mark_stage_complete(
            state, resume.STAGE_SEARCH,
            [os.path.join(work_dir, SEARCH_STAGE_SIDECAR)],
            work_dir=work_dir)
        resume.save_state(work_dir, state)

    section(f"Phase {n()}: Determining single-copy genes and writing outputs...")
    final_hmm_path, num_targets, missed_path = phase_determine_and_write(
        out_dir, filtered_hmm_path, hits_by_genome, kept_ids, filtered_accs, pfam_info,
        pfam_version, sources, organism_names, missed, args)

    if not args.keep_working_dir:
        shutil.rmtree(work_dir, ignore_errors=True)

    report_finish(out_dir, final_hmm_path, num_targets, len(kept_ids), pfam_version,
                  missed_path, args)


def main():  # pragma: no cover
    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    try:
        gen_scg_hmms(args)
    except KeyboardInterrupt:
        print()
        report_very_early_exit("Interrupted by user.", "yellow")
    except (TaxonNotFound, AmbiguousTaxon, WantedRefTaxError) as e:
        report_very_early_exit(str(e))
    except GenSCGHMMsError as e:
        report_very_early_exit(str(e))


if __name__ == "__main__":
    main()
