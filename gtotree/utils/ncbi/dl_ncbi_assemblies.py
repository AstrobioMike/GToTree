"""
Targets are given either as assembly accessions (`-a`) or pulled by taxonomy (`-w`)
from the NCBI or GTDB tables, optionally dereplicated to one best genome per rank.

NOTE ON SEMI-CODE-DUPLICATION, MIKE: GToTree already downloads NCBI assemblies during a run,
in utils/misc/processing_genomes.py. That one downloads-and-processes a single format as
part of building a tree; this one is a standalone fetcher for any of eight formats and
does no processing. They are deliberately separate rather than one parameterized
routine. The retry/backoff engine below is shared in spirit with bit's
`bit dl-ncbi-assemblies` -- if you fix a download bug here, check there too.
"""

import os
import sys
import time
import shutil
import argparse
import urllib.request
import urllib.error
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from tqdm import tqdm # type: ignore

from gtotree.cli.common import CustomRichHelpFormatter, add_version_arg
# The retry/backoff POLICY is shared with the in-run downloader rather than
# reimplemented: same throttle split, same sawtooth, same ceilings. Only the transfer
# differs (this one guards against NCBI error pages and writes atomically).
from gtotree.utils.misc.processing_genomes import (_sleep_backoff,
                                                   NCBI_DOWNLOAD_MAX_RETRIES,
                                                   NCBI_DOWNLOAD_TIMEOUT,
                                                   NCBI_THROTTLE_STATUS,
                                                   NCBI_TRANSIENT_STATUS)
from gtotree.utils.misc.messaging import report_message, wprint, color_text
from gtotree.utils.ncbi.get_ncbi_assembly_data import (get_ncbi_assembly_data,
                                                       ncbi_data_table_path)
from gtotree.utils.ncbi.dl_assembly_links import parse_ncbi_assembly_summary
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, AmbiguousTaxon,
                                               CrossDomainTaxon)
from gtotree.utils.taxonomy.get_accs_shared import (ASSEMBLY_LEVELS,
                                                    resolved_derep_rank,
                                                    parse_assembly_levels)
from gtotree.utils.taxonomy.wanted_ref_tax import (resolve_wanted_ref_tax_accessions,
                                                   expand_wanted_ref_tax,
                                                   describe_all_expansion,
                                                   WantedRefTaxError)
from gtotree.utils.taxonomy.exclusion_list import (load_exclusion_cores,
                                                   exclusion_list_help)


# NCBI_TRANSIENT_STATUS (imported above) is the set worth another attempt;
# NCBI_THROTTLE_STATUS is the subset meaning we are specifically being rate-limited
# and so earns true exponential backoff. Both come from the in-run downloader so the
# two GToTree downloaders can't drift on what counts as retryable.

max_threads = 20
max_retries = NCBI_DOWNLOAD_MAX_RETRIES

FORMAT_CHOICES = ["fasta", "protein", "genbank", "gff", "nt_cds",
                  "feature_tab", "report", "stats"]
SOURCE_CHOICES = ["gtdb", "ncbi"]
SECTION_CHOICES = ["refseq", "genbank", "both"]


def build_parser(parent_subparsers=None, show_detailed=False):

    desc = """
        This program downloads assembly files for NCBI genomes. Targets can be pulled
        based on taxonomy (`-w`) from either GTDB or NCBI (also see the --derep-rank parameter),
        and/or they can be explicitly specified as assembly accessions in a file (passed to `-a`).
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "dl-ncbi-assemblies",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt dl-ncbi-assemblies -w Nitrospirota`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    def h(text):
        return argparse.SUPPRESS if not show_detailed else text

    required = parser.add_argument_group("Required Parameters (one or both)")
    selection = parser.add_argument_group("Taxon-selection Parameters (used with `-w`)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-w", "--wanted-tax",
        metavar="<STR>",
        dest="wanted_ref_tax",
        action="append",
        default=None,
        help=("Target taxon to pull genomes for (a name, or 'all'). Can be given "
              "multiple times, and can be combined with `-a`."),
    )

    required.add_argument(
        "-a", "--ncbi-accessions",
        metavar="<FILE>",
        dest="ncbi_accessions",
        default=None,
        help="Single-column file of wanted NCBI assembly accessions",
    )


    selection.add_argument(
        "--source",
        type=str.lower,
        choices=SOURCE_CHOICES,
        default="gtdb",
        help="Which taxonomy to pull `-w` targets from (default: gtdb)",
    )

    selection.add_argument(
        "--ncbi-section",
        type=str.lower,
        dest="ncbi_section",
        choices=SECTION_CHOICES,
        default="both",
        help=("Which part of NCBI to draw from (default: both). You probably only "
              "need to worry about changing this from 'both' if you are setting "
              "`--derep-rank off` and/or targeting a single species. Ignored with "
              "`--source gtdb`."),
    )

    selection.add_argument(
        "--derep-rank",
        type=str.lower,
        choices=["auto", "off"] + list(RANKS),
        default="auto",
        help=("Dereplicate down to one genome per unique value of this rank "
              "(default: auto). 'auto' is two ranks finer than the "
              "target (one rank finer for eukaryotes). 'off' downloads every "
              "genome under the target taxon, so use with care :)"),
    )

    selection.add_argument(
        "--target-rank",
        type=str.lower,
        dest="target_rank",
        choices=list(RANKS),
        default=None,
        help="Target rank (if needed to disambiguate a taxon name that exists at multiple ranks)",
    )

    selection.add_argument(
        "--target-domain",
        metavar="<STR>",
        dest="target_domain",
        default=None,
        help="Target domain (if needed to disambiguate a taxon name that exists in multiple domains)",
    )

    selection.add_argument(
        "--exclusion-list",
        metavar="<FILE>",
        dest="exclusion_list",
        default=None,
        help=h(exclusion_list_help("-w")),
    )

    selection.add_argument(
        "--representatives-only",
        dest="representatives_only",
        action="store_true",
        help=h("With `--source gtdb`, only pull GTDB representative genomes; "
               "with `--source ncbi`, only pull NCBI reference genomes. If the goal is "
               "removing redundancy, the `--derep-rank` parameter can handle that while "
               "ensuring the breadth of available diversity is maintained."),
    )

    selection.add_argument(
        "--assembly-level",
        action="append",
        choices=list(ASSEMBLY_LEVELS),
        default=None,
        help=h("Only include genomes (from `-w`) at this assembly level. "
               "Can be provided multiple times."),
    )

    selection.add_argument(
        "--min-completeness",
        metavar="<FLOAT>",
        type=float,
        dest="min_completeness",
        default=None,
        help=h("Don't include any genomes (from `-w`) below this checkm completeness "
               "(default: None). If set, genomes with no recorded "
               "value are also excluded."),
    )

    selection.add_argument(
        "--max-contamination",
        metavar="<FLOAT>",
        type=float,
        dest="max_contamination",
        default=None,
        help=h("Don't include any genomes (from `-w`) above this checkm contamination "
               "(default: None). If set, genomes with no recorded "
               "value are also excluded."),
    )

    selection.add_argument(
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help=("Run the selection (based on `-w`), but just report how many genomes would be downloaded"),
    )

    optional.add_argument(
        "-f",
        "--format",
        type=str.lower,
        choices=FORMAT_CHOICES,
        default="fasta",
        help='Format to download (default: "fasta")',
    )

    optional.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=10,
        help="Number of concurrent downloads (default: 10; capped at 20 to keep NCBI happy)",
    )

    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        dest="output_dir",
        default=".",
        help="Directory for output files (default: current directory)",
    )

    optional.add_argument(
        "-s",
        "--show-detailed-help",
        dest="show_detailed_help",
        action="store_true",
        help="Show detailed help, including additional taxon-selection parameters",
    )

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show basic help",
    )
    add_version_arg(optional)

    return parser


def main():

    argv = sys.argv[1:]
    if "-s" in argv or "--show-detailed-help" in argv:
        build_parser(show_detailed=True).print_help(sys.stderr)
        sys.exit(0)

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dl_ncbi_assemblies(args)


def dl_ncbi_assemblies(args):

    preflight_checks(args)

    wanted_accs, selections, expansion_note = resolve_targets(args)

    if getattr(args, "dry_run", False):
        report_selection(wanted_accs, selections, expansion_note, args)
        return

    run_data = setup(args, wanted_accs)

    if selections:
        report_selection(wanted_accs, selections, expansion_note, args)

    run_data = parse_main_assembly_table(run_data)

    summarize_search(run_data)

    run_data = download_assemblies(run_data)

    report_finish(run_data)


def preflight_checks(args):

    accessions_file = getattr(args, "ncbi_accessions", None)
    wanted = getattr(args, "wanted_ref_tax", None)

    if not accessions_file and not wanted:
        report_message("Nothing to download. Provide an accessions file with `-a`, a "
                       "taxon with `-w`, or both.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(1)

    if accessions_file and not os.path.isfile(accessions_file):
        report_message(f"The specified accessions file, '{accessions_file}', "
                       "was not found.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(1)

    if (getattr(args, "assembly_level", None)
            and getattr(args, "source", "gtdb") == "gtdb"):
        report_message("`--assembly-level` is only applicable with `--source ncbi`. "
                       "Re-run with `--source ncbi`, or drop `--assembly-level`.",
                       "yellow", ii="    ", si="    ", width=100,
                       trailing_newline=True)
        sys.exit(1)

    exclusion_list = getattr(args, "exclusion_list", None)
    if exclusion_list:
        if not wanted:
            report_message("An `--exclusion-list` was provided, but that only has an "
                           "effect alongside `-w` (it removes genomes from what `-w` "
                           "pulls in). Accessions given directly with `-a` are always "
                           "downloaded as provided, so nothing would be excluded. Stopping "
                           "here to let you adjust input parameters we we're all on the same page :)",
                           "yellow", ii="    ", si="    ", width=100,
                           trailing_newline=True)
            sys.exit(1)
        if not os.path.isfile(exclusion_list):
            report_message(f"The specified exclusion list, '{exclusion_list}', wasn't "
                           "found :(", "yellow",
                           ii="    ", si="    ", width=100, trailing_newline=True)
            sys.exit(1)

    # a dry run reports and stops, so it must not create anything
    if not getattr(args, "dry_run", False):
        if args.output_dir and not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)

    get_ncbi_assembly_data()


def _selection_kwargs(args):
    """
    The selection knobs shared by every `-w` target in a run.
    """
    # Deliberately NOT the source's own default (`None`), which would make a GTDB
    # pull species-representatives-only. Two reasons:
    #   1. --derep-rank is `auto` here, so volume is already controlled by
    #      dereplication. Reps-only was a volume brake, and on top of derep it only
    #      narrows what each group can be represented BY -- quietly excluding a
    #      possibly-better genome for no benefit.
    #   2. `gtt get-accs-from-gtdb` already defaults to all genomes, so deferring to
    #      the source default here would make this subcommand the odd one out.
    # With this, --representatives-only means the same thing for both sources, and
    # matches my `bit dl-ncbi-assemblies`.
    reps_only = bool(getattr(args, "representatives_only", False))

    return {
        "target_rank": getattr(args, "target_rank", None),
        "derep_rank": resolved_derep_rank(args),
        "target_domain": getattr(args, "target_domain", None),
        "ncbi_section": getattr(args, "ncbi_section", "refseq"),
        "reps_only": reps_only,
        "assembly_levels": parse_assembly_levels(getattr(args, "assembly_level", None)),
        "min_completeness": getattr(args, "min_completeness", None),
        "max_contamination": getattr(args, "max_contamination", None),
        "exclude_cores": load_exclusion_cores(getattr(args, "exclusion_list", None)),
        "include_rows": False,
    }


def resolve_targets(args):
    """
    Turn `-a` and/or `-w` into ONE deduplicated accession list.

    Returns (accessions, selections, (expansion_note, num_from_file)). `selections` is
    one record per resolved `-w` target carrying what it selected and how many of those
    were NEW after deduplication -- with a repeatable `-w`, overlapping taxa mean the
    per-taxon counts will not sum to the total, and that gap is the thing worth showing.
    """
    accessions = []
    seen = set()

    accessions_file = getattr(args, "ncbi_accessions", None)
    if accessions_file:
        with open(accessions_file, "r") as f:
            for line in f:
                acc = line.strip()
                if acc and acc not in seen:
                    seen.add(acc)
                    accessions.append(acc)

    num_from_file = len(accessions)

    wanted = getattr(args, "wanted_ref_tax", None)
    if not wanted:
        return accessions, [], (None, num_from_file)

    source = getattr(args, "source", "gtdb")

    try:
        expanded, domains = expand_wanted_ref_tax(source, wanted)
    except WantedRefTaxError as err:
        _exit_with(str(err))
    expansion_note = describe_all_expansion(source, domains)

    kwargs = _selection_kwargs(args)
    selections = []

    for taxon in expanded:
        try:
            taxon_accs, selection = resolve_wanted_ref_tax_accessions(
                source, taxon, **kwargs)
        except AmbiguousTaxon:
            _exit_with(f"Since the `-w` taxon '{taxon}' occurs at more than one rank, "
                       "you'll need to specify which rank is wanted with "
                       "`--target-rank`.")
        except CrossDomainTaxon as err:
            _exit_with(f"The `-w` taxon '{err.taxon}' occurs in more than one domain "
                       f"({', '.join(err.domains_found)}), so pulling on the name alone "
                       "would mix genomes from different domains. Specify which domain "
                       "is wanted with `--target-domain` "
                       f"(e.g. `--target-domain {err.domains_found[0]}`).")
        except TaxonNotFound:
            _exit_with(f"The `-w` taxon '{taxon}' doesn't seem to exist at any rank in "
                       f"the {source} taxonomy :(")
        except (WantedRefTaxError, ValueError) as err:
            _exit_with(str(err))

        num_new = 0
        for acc in taxon_accs:
            if acc not in seen:
                seen.add(acc)
                accessions.append(acc)
                num_new += 1

        selections.append(TaxonSelection(
            taxon=taxon,
            canonical=selection.canonical,
            resolved_rank=selection.resolved_rank,
            effective_derep_rank=selection.effective_derep_rank,
            num_selected=len(taxon_accs),
            num_new=num_new,
            num_excluded=getattr(selection, "num_excluded", 0),
            warnings=list(selection.warnings)))

    return accessions, selections, (expansion_note, num_from_file)


def _exit_with(message):
    report_message(message, "yellow", ii="    ", si="    ", width=100,
                   trailing_newline=True)
    sys.exit(1)


@dataclass
class TaxonSelection:
    """What one `-w` target resolved to, for the reporting layer."""
    taxon: str = None
    canonical: str = None
    resolved_rank: str = None
    effective_derep_rank: str = None
    num_selected: int = 0
    num_new: int = 0
    num_excluded: int = 0
    warnings: list = None


def report_selection(accessions, selections, expansion_note, args):
    """
    What each `-w` resolved to, and the combined total.

    `num_selected` is what the taxonomy core returned for that taxon; `num_new` is what
    survived deduplication against `-a` and against earlier `-w` targets. They differ
    whenever two taxa overlap, which is exactly when a bare per-taxon count would
    mislead. Shared with `--dry-run` so the numbers can't drift apart.
    """
    note, num_from_file = (expansion_note if expansion_note else (None, 0))

    print("")
    if note:
        report_message(note, "yellow", ii="    ", si="    ", width=100)

    if num_from_file:
        wprint(f"{num_from_file:,} accession(s) read from "
               f"{color_text(args.ncbi_accessions)}", ii="    ", si="    ")
        print("")

    for sel in selections:
        derep = (f"dereplicated to one genome per {sel.effective_derep_rank}"
                 if sel.effective_derep_rank else "dereplication off")
        wprint(f"- {color_text(sel.canonical)}", ii="    ", si="    ")
        wprint(f"- resolved to {sel.resolved_rank}; {derep}", ii="        ", si="        ")

        line = f"- {sel.num_selected:,} genome(s) selected"
        overlap = sel.num_selected - sel.num_new
        if overlap:
            line += f", {sel.num_new:,} new ({overlap:,} already counted)"
        wprint(line, ii="        ", si="        ")

        for warning in (sel.warnings or []):
            report_message(warning, "yellow", ii="        ", si="        ", width=100)
        print("")

    total = len(accessions)
    if getattr(args, "dry_run", False):
        wprint(color_text(f"Total that would be downloaded: {total:,} genome(s) in "
                          f"{args.format} format.", "green"), ii="    ", si="    ")
        print("")
    else:
        wprint(color_text(f"Total to download: {total:,} genome(s)", "green"),
               ii="    ", si="    ")


def setup(args, wanted_accs=None):

    if wanted_accs is None:
        with open(args.ncbi_accessions, "r") as f:
            wanted_accs = [line.strip() for line in f if line.strip()]

    return RunData(
        wanted_format=args.format,
        num_jobs=args.jobs,
        output_dir=args.output_dir,
        wanted_accs=wanted_accs,
        num_wanted=len(wanted_accs),
        ncbi_sub_table_path=Path(args.output_dir) / "wanted-ncbi-accessions-info.tsv",
        not_found_path=Path(args.output_dir) / "ncbi-accessions-not-found.txt",
        not_downloaded_path=Path(args.output_dir) / "ncbi-accessions-not-downloaded.tsv",
        from_taxon=bool(getattr(args, "wanted_ref_tax", None)),
    )


def parse_main_assembly_table(run_data):

    report_message(f"Targeting {run_data.num_wanted:,} accession(s) in "
                   f"{run_data.wanted_format} format...", "yellow",
                   ii="    ", si="    ", width=100, trailing_newline=False)

    return parse_ncbi_assembly_summary(ncbi_data_table_path(), run_data)


@dataclass
class RunData:
    wanted_format: str = None
    num_jobs: int = 10
    output_dir: str = None
    wanted_accs: list = None
    num_wanted: int = 0
    num_found: int = 0
    num_not_found: int = 0
    num_downloaded: int = 0
    num_skipped: int = 0
    num_not_downloaded: int = 0
    ncbi_sub_table_path: str = None
    not_found_path: str = None
    not_downloaded_path: str = None
    quiet: bool = False
    from_taxon: bool = False

    @property
    def not_found_reason(self):
        """
        Why accessions might be missing depends on where they came from. Telling
        someone who passed `-w` that their input "may be invalid" is wrong -- they
        never typed an accession -- so that path just states the count.
        """
        if self.from_taxon:
            return ""
        return "They may be invalid or suppressed."

    @property
    def none_found_hint(self):
        if self.from_taxon:
            return "None of the selected genomes could be matched at NCBI."
        return "Are the inputs assembly accessions?"


def summarize_search(run_data):

    if run_data.num_found == run_data.num_wanted:
        return

    if run_data.num_found > 0:
        reason = f" {run_data.not_found_reason}" if run_data.not_found_reason else ""
        report_message(f"{run_data.num_not_found:,} accession(s) not found at NCBI."
                       f"{reason} Reported in:", "yellow",
                       ii="    ", si="    ", width=100)
        wprint(color_text(str(run_data.not_found_path)), ii="        ", si="        ")
        print("")
    else:
        report_message(f"None of the {run_data.num_wanted:,} target accession(s) were "
                       f"found at NCBI. {run_data.none_found_hint}", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        os.remove(run_data.not_found_path)
        sys.exit(1)


def _is_error_page(resp):
    """
    NCBI sometimes answers a missing file with an HTML/XML error page and a 200, so a
    status check alone isn't enough. Checked before writing anything so we never
    replace a good file with an error page.
    """
    content_type = (resp.headers.get("Content-Type") or "").lower()
    return "xml" in content_type or "html" in content_type


def download_one(target_link, local_dest, retries=max_retries):
    """
    Fetch one file, atomically, with the shared retry/backoff policy.

    Returns (local_dest, error_or_None, status) where status is one of
    "downloaded" / "skipped" / "failed" (permanent) / "failed_transient" (worth
    another pass).
    """
    local_path = Path(local_dest)
    tmp_path = local_path.with_name(local_path.name + ".tmp")

    # atomic writes (below) guarantee any 'final file' completed this run
    if local_path.exists() and local_path.stat().st_size > 0:
        return (local_dest, None, "skipped")

    local_path.parent.mkdir(parents=True, exist_ok=True)

    if not os.access(local_path.parent, os.W_OK):
        return (local_dest,
                f"download directory '{local_path.parent}' is not writable", "failed")

    # clear any stale tmp left behind by a prior interrupted run
    tmp_path.unlink(missing_ok=True)

    req = urllib.request.Request(target_link, headers={"User-Agent": "curl/8.0"})

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=NCBI_DOWNLOAD_TIMEOUT) as resp:

                if _is_error_page(resp):
                    if attempt == retries:
                        return (local_dest,
                                f"NCBI returned an error page after {retries} attempts",
                                "failed_transient")
                    _sleep_backoff(attempt)
                    continue

                with open(tmp_path, "wb") as out:
                    shutil.copyfileobj(resp, out, length=1024 * 256)

            if tmp_path.stat().st_size == 0:
                tmp_path.unlink(missing_ok=True)
                if attempt == retries:
                    return (local_dest, "Downloaded file was empty", "failed_transient")
                _sleep_backoff(attempt)
                continue

            os.replace(tmp_path, local_path)
            return (local_dest, None, "downloaded")

        except urllib.error.HTTPError as err:
            tmp_path.unlink(missing_ok=True)

            if err.code == 404:
                return (local_dest,
                        "Not available in requested format (404)", "failed")

            if err.code not in NCBI_TRANSIENT_STATUS:
                return (local_dest, f"HTTP {err.code}", "failed")

            if attempt == retries:
                return (local_dest, f"HTTP {err.code} after {retries} attempts",
                        "failed_transient")

            # a 429, or any response that bothered to tell us when to come back, is a
            # throttle -- back off properly. A bare 5xx is just a hiccup.
            throttled = (err.code in NCBI_THROTTLE_STATUS
                         or err.headers.get("Retry-After") is not None)
            _sleep_backoff(attempt, err, throttled=throttled)

        except (urllib.error.URLError, OSError) as err:
            tmp_path.unlink(missing_ok=True)
            if attempt == retries:
                return (local_dest, str(err), "failed_transient")
            _sleep_backoff(attempt)


def run_download_pass(targets, run_data, desc="Progress"):
    """
    runs one pooled download pass over a list of (target_link, local_dest) tuples.
    returns (permanent_failures, transient_failures, num_skipped) where each
    failure list holds (target_link, local_dest, error) so transient ones can be retried
    """
    permanent = []
    transient = []
    num_skipped = 0

    link_by_dest = {dest: link for link, dest in targets}

    pool = ThreadPoolExecutor(max_workers=min(run_data.num_jobs, max_threads))
    try:
        futures = {
            pool.submit(download_one, link, dest): dest
            for link, dest in targets
        }

        desc_buffer = "      "
        ncols = 70
        with tqdm(total=len(targets), desc=f"{desc_buffer}{desc}", unit=" file", ncols=ncols, smoothing=0.05) as pbar:
            for future in as_completed(futures):
                dest, error, status = future.result()
                if status == "failed":
                    permanent.append((link_by_dest[dest], dest, error))
                elif status == "failed_transient":
                    transient.append((link_by_dest[dest], dest, error))
                elif status == "skipped":
                    num_skipped += 1
                pbar.update(1)
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        pool.shutdown(wait=True)

    return permanent, transient, num_skipped


def download_assemblies(run_data):

    report_message("Downloading assemblies...", "yellow",
                   ii="    ", si="    ", width=100, trailing_newline=True)

    targets = []
    with open(run_data.ncbi_sub_table_path, "r") as f:
        header = f.readline().strip().split("\t")
        link_idx = header.index("target_link")
        dest_idx = header.index("local_destination")
        for line in f:
            fields = line.strip().split("\t")
            targets.append((fields[link_idx], fields[dest_idx]))

    permanent, transient, num_skipped = run_download_pass(targets, run_data)

    # second pass on transient-only failures
    if transient:
        retry_targets = [(link, dest) for link, dest, _ in transient]
        print("")
        report_message(f"{len(transient):,} file(s) failed with transient errors, "
                       "doing another pass to see if we can grab them...", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)

        time.sleep(3)
        retry_permanent, retry_transient, retry_skipped = run_download_pass(
            retry_targets, run_data, desc="Progress"
        )
        num_skipped += retry_skipped
        # anything still failing after the retry is final, regardless of category
        permanent.extend(retry_permanent)
        permanent.extend(retry_transient)

    failed = [(dest, error) for _, dest, error in permanent]

    if failed:
        with open(run_data.not_downloaded_path, "w") as fh:
            fh.write("accession\terror\n")
            for dest, error in failed:
                acc = Path(dest).stem.split(".")[0]
                fh.write(f"{acc}\t{error}\n")
        run_data.num_not_downloaded = len(failed)
    else:
        Path(run_data.not_downloaded_path).unlink(missing_ok=True)

    run_data.num_skipped = num_skipped
    run_data.num_downloaded = len(targets) - len(failed)

    return run_data


def report_finish(run_data):

    skipped_note = ""
    if run_data.num_skipped > 0:
        skipped_note = f" ({run_data.num_skipped} already present, skipped)"

    if run_data.num_downloaded == run_data.num_wanted:
        report_message(f"All {run_data.num_wanted:,} file(s) downloaded successfully!"
                       f"{skipped_note}", "green",
                       ii="    ", si="    ", width=100, trailing_newline=True)

    elif run_data.num_downloaded == run_data.num_found:
        report_message(f"All {run_data.num_found:,} found file(s) downloaded "
                       f"successfully!{skipped_note}", "green",
                       ii="    ", si="    ", width=100, trailing_newline=True)

    elif run_data.num_not_downloaded > 0:
        report_message(f"{run_data.num_not_downloaded:,} file(s) failed to download "
                       "from NCBI. They may not be available in the requested format, "
                       "or it may have been a transient problem. Reported in:",
                       "yellow", ii="    ", si="    ", width=100)
        wprint(color_text(str(run_data.not_downloaded_path)),
               ii="        ", si="        ")
        print("")

        if run_data.num_downloaded > 0:
            report_message(f"The remaining {run_data.num_downloaded:,} found file(s) "
                           f"downloaded successfully.{skipped_note}", "green",
                           ii="    ", si="    ", width=100, trailing_newline=True)
        else:
            report_message(f"No files were successfully downloaded...{skipped_note}",
                           "yellow", ii="    ", si="    ", width=100,
                           trailing_newline=True)
            sys.exit(1)
