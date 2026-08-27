import os
import sys
import shutil
import time
import pandas as pd # type: ignore
import tempfile
from collections import Counter
from gtotree.utils.misc.messaging import (color_text,
                                     bullet_list,
                                     report_message,
                                     report_very_early_exit,
                                     report_run_already_complete,
                                     report_missing_input_genomes_file,
                                     report_missing_pfam_targets_file,
                                     report_missing_ko_targets_file,
                                     report_missing_mapping_file,
                                     report_missing_exclusion_list_file,
                                     report_problem_with_mapping_file,
                                     report_notice,
                                     report_run_info_banner,
                                     RUN_INFO_BANNER,
                                     RUN_INFO_BANNER_TRAILING_BLANKS,
                                     many_genomes_notice,
                                     few_genomes_notice,
                                     absurd_number_of_genomes_notice,
                                     gtotree_header,
                                     stdout_and_log,
                                     spinner)
from gtotree.utils.hmms.scg_hmms_setup import (autopick_scg_set,
                                          check_viruses_have_their_own_hmms,
                                          ViralTaxonNeedsOwnHMMs,
                                          resolve_hmm_arg,
                                          resolve_hmm_source,
                                          populate_SCG_targets)
from gtotree.utils.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data
from gtotree.utils.gtdb.get_gtdb_data import get_gtdb_data
from gtotree.utils.ko.get_kofamscan_data import get_kofamscan_data
from gtotree.utils.pfam.get_pfam_data import get_pfam_data
from gtotree.utils.taxonomy.tax_ranks import RANKS, accession_core
from gtotree.utils.taxonomy.tax_select import AmbiguousTaxon, TaxonNotFound, CrossDomainTaxon
from gtotree.utils.taxonomy.tax_targets import is_all_target
from gtotree.utils.taxonomy.wanted_ref_tax import (resolve_wanted_ref_tax_accessions,
                                                   expand_wanted_ref_tax,
                                                   describe_all_expansion,
                                                   WantedRefTaxError)
from gtotree.utils.misc.general import (ToolsUsed,
                                   read_single_column_file,
                                   duplicate_entries,
                                   CorruptRunData,
                                   SOURCE_ACCESSION, SOURCE_GENBANK, SOURCE_FASTA,
                                   SOURCE_AMINO_ACID,
                                   populate_run_data,
                                   read_run_data,
                                   wanted_ref_tax_list)
from gtotree.utils.misc.target_id_checks import (check_target_id_file,
                                                 TargetIDFormatError)
from gtotree.utils.misc.stages import PipelineStage
from gtotree.utils.misc.resume_state import (ResumeProfile, hash_file_contents,
                                        STATE_VERSION)
from gtotree.utils.misc.context import log_file_var


def preflight_checks(args):
    """
    Everything that has to be settled before the run starts
    """
    check_for_essential_deps()
    args = primary_args_validation(args)
    previous_run_data = load_previous_run_data(args)
    check_for_required_dbs(args, previous_run_data)
    report_run_info_banner()
    selections = select_wanted_ref_tax(args, previous_run_data)
    args = resolve_hmm(args, selections, previous_run_data)
    args, run_data = setup_run_data(args, previous_run_data)
    run_data = merge_wanted_ref_tax(run_data, selections,
                                    exclusion_list=getattr(args, "exclusion_list", None))
    check_for_genome_id_collisions(run_data)
    check_for_min_input_genomes(run_data)
    run_data = track_tools_used(args, run_data)
    args, run_data = final_setups(args, run_data)

    return args, run_data


RESUME = ResumeProfile(
    name="GToTree",
    field_labels={
        "state_version": "the run-state format",
        "ncbi_accessions_sha256": "the NCBI accessions file (-a)",
        "genbank_files_sha256": "the GenBank files list (-g)",
        "fasta_files_sha256": "the fasta files list (-f)",
        "amino_acid_files_sha256": "the amino-acid files list (-A)",
        "mapping_file_sha256": "the mapping file (-m)",
        "target_pfams_sha256": "the target Pfams file (-p)",
        "target_kos_sha256": "the target KOs file (-K)",
        "exclusion_list_sha256": "the exclusion list (--exclusion-list)",
        "hmm": "-H/--hmm",
        "wanted_ref_tax": "-w/--wanted-ref-tax",
        "target_rank": "--target-rank",
        "target_domain": "--target-domain",
        "derep_rank": "--derep-rank",
        "source": "--source",
        "add_gtdb_tax": "-D/--add-gtdb-tax",
        "add_ncbi_tax": "-t/--add-ncbi-tax",
        "lineage": "-L/--lineage-ranks",
        "seq_length_cutoff": "-c/--seq-length-cutoff",
        "gene_representation_cutoff": "-r/--gene-rep-cutoff",
        "genome_hits_cutoff": "-G/--genome-hits-cutoff",
        "best_hit_mode": "-B/--best-hit-mode",
        "no_super5": "-X/--no-super5",
        "no_tree": "-N/--no-tree",
        "tree_program": "-T/--tree-program",
        "nucleotide_mode": "-z/--nucleotide-mode",
        "keep_gene_alignments": "-k/--keep-gene-alignments",
    },
)

# input files are fingerprinted by contents
_INPUT_FILE_FIELDS = (
    ("ncbi_accessions_sha256", "ncbi_accessions"),
    ("genbank_files_sha256", "genbank_files"),
    ("fasta_files_sha256", "fasta_files"),
    ("amino_acid_files_sha256", "amino_acid_files"),
    ("mapping_file_sha256", "mapping_file"),
    ("target_pfams_sha256", "target_pfams_file"),
    ("target_kos_sha256", "target_kos_file"),
    ("exclusion_list_sha256", "exclusion_list"),
)

# args that change how the run executes, not what it produces. Listed explicitly
# so adding any new flags is a deliberate decision about if they impact a resume attempt or not
_EXECUTION_ONLY_ARGS = ("num_jobs", "num_muscle_threads", "keep_working_dir",
                        "tmp_dir", "force_overwrite", "resume", "output_dir",
                        "run_files_dir", "run_files_dir_rel",
                        "output_already_existed")


def build_fingerprint(args):
    """
    Everything about this run that affects what it produces
    """
    fingerprint = {"state_version": STATE_VERSION}

    for field, dest in _INPUT_FILE_FIELDS:
        fingerprint[field] = hash_file_contents(getattr(args, dest, None))

    for field in RESUME.field_labels:
        if field in fingerprint:
            continue
        fingerprint[field] = getattr(args, field, None)

    fingerprint["wanted_ref_tax"] = sorted(wanted_ref_tax_list(args)) or None

    return fingerprint

def check_for_essential_deps():
    commands = ["muscle", "trimal"]
    for cmd in commands:
        program_check(cmd, essential = True)


def program_check(cmd, essential = False):
    path = shutil.which(cmd)
    if not path:
        if essential:
            report_message(f"{cmd} is an essential dependency, but it's not in your PATH :(")
            report_very_early_exit(None, color = "yellow")
        else:
            report_message(f"You specified to use {cmd}, but it's not in your PATH :(", color = "yellow")
            report_very_early_exit(None, color = "yellow")


def primary_args_validation(args):
    check_for_minimum_args(args)
    check_optional_deps(args)
    check_set_values(args)
    check_lineage(args)
    check_tree_program(args)
    checks_for_nucleotide_mode(args)
    check_wanted_ref_tax_args(args)
    check_target_id_formats(args)
    args = check_input_genome_files(args)
    args = check_output_dir(args)
    return args


TARGET_ID_FILES = (
    ("target_pfams_file", "-p", "pfam"),
    ("target_kos_file", "-K", "ko"),
)


def check_target_id_formats(args):
    """
    Confirm any `-p`/`-K` file actually holds IDs of the type its flag expects
    """
    for dest, flag, key in TARGET_ID_FILES:
        path = getattr(args, dest, None)
        if not path:
            continue

        check_path(path, flag)
        check_for_whitespace(path, flag)

        try:
            check_target_id_file(path, key)
        except TargetIDFormatError as err:
            report_message(str(err))
            report_very_early_exit(suggest_help=True)


def check_for_minimum_args(args):
    if (not args.ncbi_accessions and not args.genbank_files and not args.fasta_files
            and not args.amino_acid_files and not args.wanted_ref_tax):
        report_message("You need to provide at least one input-genome source!")
        report_very_early_exit(suggest_help=True)
    # `-H` is only mandatory without `-w`; with `-w` provided, gtotree can auto-select the used SCG-set
    if not args.hmm and not args.wanted_ref_tax:
        report_message("You need to specify the target-SCG HMMs (`-H`) you want to use. You "
                       "can view the available gene-sets packaged with GToTree by running `gtt hmms`. "
                       "Or, if you add reference genomes by taxonomy with `-w`/`--wanted-ref-tax`, "
                       "GToTree will select a suitable pre-packaged set for you.")
        report_very_early_exit(suggest_help=True)


def check_optional_deps(args):
    pass


def check_set_values(args):

    if args.seq_length_cutoff < 0 or args.seq_length_cutoff > 1:
        report_message("The sequence-length-cutoff value (passed to `-c`) must be between 0 and 1 (inclusive).")
        report_very_early_exit()

    if args.gene_representation_cutoff <0 or args.gene_representation_cutoff > 1:
        report_message("The gene-representation-cutoff value (passed to `-g`) must be between 0 and 1 (inclusive).")
        report_very_early_exit()

    if args.genome_hits_cutoff < 0 or args.genome_hits_cutoff > 1:
        report_message("The genome-hits-cutoff value (passed to `-G`) must be between 0 and 1 (inclusive).")
        report_very_early_exit()

    if args.num_jobs < 1:
        report_message("The number of jobs to run in parallel (passed to `-j`) must be at least 1.")
        report_very_early_exit()

    if args.num_muscle_threads < 1:
        report_message("The number of threads to use for MUSCLE (passed to `-M`) must be at least 1.")
        report_very_early_exit()


def check_lineage(args):
    accepted_ranks = ["domain", "phylum", "class", "order", "family", "genus", "species", "strain"]
    lineage_list = args.lineage.split(",")

    for rank in lineage_list:
        if rank.lower() not in accepted_ranks:
            report_message(f'You specified "{args.lineage}" to the `-L` argument, but "{rank}" is not an accepted taxonomic rank.')
            print(f"\n  Accepted ranks are any combination of the below entered as a comma-delimited list:\n\n        {'\n        '.join(accepted_ranks)}")
            report_very_early_exit()

    if args.lineage != "domain,phylum,class,genus,species" and not args.add_ncbi_tax and not args.add_gtdb_tax:
        report_message("You've specified a custom lineage (`-L`), but neither the "
                       "`--add-gtdb-tax` or `--add-ncbi-tax` flags were provided to indicate which taxonomy to use.")
        report_very_early_exit(suggest_help=True)

    if args.add_ncbi_tax and args.add_gtdb_tax:
        report_message("You've specified add taxonomic info based on GTDB and NCBI taxonomies. "
                       "Please choose one or the other.")
        report_very_early_exit(suggest_help=True)


TREE_PROGRAM_EXECUTABLES = {
    "FastTree":     "FastTree",
    "FastTreeMP":   "FastTreeMP",
    "VeryFastTree": "VeryFastTree",
    "IQTREE":       "iqtree",
}


def check_tree_program(args):
    accepted_programs = list(TREE_PROGRAM_EXECUTABLES)
    if args.tree_program not in accepted_programs:
        report_message(f'You specified "{args.tree_program}" to the `-T` argument, but that\'s not an available treeing program.')
        print(f"\n  Available treeing programs are:\n\n        {'\n        '.join(accepted_programs)}")
        report_very_early_exit(suggest_help=True)
    program_check(TREE_PROGRAM_EXECUTABLES[args.tree_program])


def checks_for_nucleotide_mode(args):
    if args.nucleotide_mode:
        if args.amino_acid_files:
            report_message("You've specified wanting to work with nucleotide sequences (by passing the `-z` parameter), "
                           "but also provided some input genomes as amino-acid files (passed to `-A`). We can't confidently reverse-translate "
                           "amino-acid seqs to nucleotide seqs, so we can't take both of those options.")
            report_very_early_exit(suggest_help=True)
        if args.genbank_files:
            report_message("You've specified wanting to work with nucleotide sequences (by passing the `-z` parameter), "
                           "but also provided some input genomes as genbank files (passed to `-g`). Input genbank files are currently "
                           "not supported with nucleotide mode.")
            report_very_early_exit(suggest_help=True)


def check_wanted_ref_tax_args(args):
    """
    Cross-argument validation for the --wanted-ref-tax (-w) family, done up front
    (before any asset is fetched or taxon resolved) so bad combinations fail fast:

      * --target-rank / --derep-rank only mean something alongside -w.
      * --target-rank, when given, must be one of the 7 taxonomic ranks.
      * --derep-rank must be 'auto', 'off', or one of the 7 ranks.
    """
    ranks = list(RANKS)

    if not args.wanted_ref_tax:
        dangling = []
        if args.target_rank is not None:
            dangling.append("`--target-rank`")
        # --derep-rank has a non-None default ("auto"); only a user-changed value is dangling
        if args.derep_rank not in (None, "auto"):
            dangling.append("`--derep-rank`")
        if dangling:
            joined = " and ".join(dangling)
            report_message(f"You provided {joined}, but that only applies when adding "
                           "reference genomes by taxonomy with `-w`/`--wanted-ref-tax`.")
            report_very_early_exit(suggest_help=True)
        return

    if args.target_rank is not None and any(is_all_target(t)
                                            for t in wanted_ref_tax_list(args)):
        report_message("You passed `-w all` along with `--target-rank`, but 'all' is "
                       "expanded to the source's domains, so there is no name left "
                       "for `--target-rank` to disambiguate. Drop `--target-rank`, or "
                       "name the taxon you want instead of 'all'.")
        report_very_early_exit(suggest_help=True)

    if args.target_rank is not None and args.target_rank.strip().lower() not in ranks:
        report_message(f'You specified "{args.target_rank}" to `--target-rank`, but that '
                       "is not an accepted taxonomic rank.")
        print(f"\n  Accepted ranks are:\n\n        {'\n        '.join(ranks)}\n")
        report_very_early_exit(suggest_help=True)

    if args.derep_rank is not None:
        dr = args.derep_rank.strip().lower()
        if dr not in ("auto", "off") and dr not in ranks:
            report_message(f'You specified "{args.derep_rank}" to `--derep-rank`, but that '
                           "is not 'auto', 'off', or an accepted taxonomic rank.")
            print(f"\n  Accepted ranks are:\n\n        {'\n        '.join(ranks)}\n")
            report_very_early_exit(suggest_help=True)


def wanted_ref_tax_already_resolved(previous_run_data):
    """
    True when a resume can reuse the `-w` accessions the previous run selected
    """
    return (previous_run_data is not None
            and bool(previous_run_data.get_wanted_ref_tax_accs()))


def select_wanted_ref_tax(args, previous_run_data=None):
    """
    Resolve `-w` to reference accessions, or reuse the previous run's

    `-w` may be given more than once. Each taxon is resolved and dereplicated on its
    own, then merged

    Returns a list of RefGenomeSelection (empty when there's nothing to resolve).
    """
    if not args.wanted_ref_tax:
        return []

    if wanted_ref_tax_already_resolved(previous_run_data):
        return []

    selections = []

    try:
        wanted, domains = expand_wanted_ref_tax(args.source,
                                                wanted_ref_tax_list(args))
    except WantedRefTaxError as err:
        report_message(str(err))
        report_very_early_exit(suggest_help=True)

    note = describe_all_expansion(args.source, domains)
    if note:
        report_message(note, color=None, ii="    ", si="    ", newline=False)
        print("")

    for taxon in wanted:
        try:
            with spinner(f"Gathering references for '{taxon}'...", "",
                         clear_on_done=True):
                _accessions, selection = resolve_wanted_ref_tax_accessions(
                    args.source, taxon,
                    target_rank=args.target_rank, derep_rank=args.derep_rank,
                    target_domain=getattr(args, "target_domain", None),
                    building_tree=True)
        except AmbiguousTaxon:
            report_message(f"Since the `-w` taxon '{taxon}' occurs at more than "
                           "one rank, you'll need to specify which rank is wanted with "
                           "`--target-rank`.")
            report_very_early_exit(suggest_help=True)
        except CrossDomainTaxon as err:
            report_message(f"The `-w` taxon '{taxon}' occurs in more than one domain "
                           f"({', '.join(err.domains_found)}). Specify which domain is "
                           "wanted with `--target-domain` "
                           f"(e.g. `--target-domain {err.domains_found[0]}`).")
            report_very_early_exit(suggest_help=True)
        except TaxonNotFound as err:
            # a bad --target-domain resolves to TaxonNotFound with a specific message;
            # surface it, else fall back to the generic "doesn't exist" wording
            detail = str(err)
            if "domain" in detail:
                report_message(detail)
            else:
                report_message(f"The `-w` taxon '{taxon}' doesn't seem to exist at any "
                               f"rank in the {args.source} taxonomy :(")
            report_very_early_exit(suggest_help=True)
        except (WantedRefTaxError, ValueError) as err:
            report_message(str(err))
            report_very_early_exit(suggest_help=True)

        selections.append(selection)

    return selections


def merge_wanted_ref_tax(run_data, selections, exclusion_list=None):
    """
    Fold each `-w` selection's accessions into run_data's NCBI-accession input pool
    (deduping against user-provided accessions, and against each other). The run_data
    half of the job select_wanted_ref_tax() starts.

    `exclusion_list` is an optional path to a single-column file of accessions. Any
    accession it names is dropped from every `-w` selection BEFORE the selection is
    merged, so an excluded genome never reaches processing. The exclusion set is
    applied only to `-w`-selected accessions, not to anything the user provided
    directly through `-a`. Matching is on the accession core only (ignores GCA/GCF
    and version suffix)
    """
    # a resume passes no selections (they were resolved and merged by the original
    # run); leave the recorded exclusion count from that run untouched
    if not selections:
        return run_data

    # keyed by core, empty cores dropped so junk can't match an unrecognized accession
    excluded_cores = set()
    if exclusion_list:
        for entry in read_single_column_file(exclusion_list):
            core = accession_core(entry)
            if core:
                excluded_cores.add(core)

    total_excluded = 0

    for selection in selections or []:
        if excluded_cores:
            kept = [acc for acc in selection.accessions
                    if accession_core(acc) not in excluded_cores]
            num_excluded = len(selection.accessions) - len(kept)
        else:
            kept = selection.accessions
            num_excluded = 0
        total_excluded += num_excluded

        added = run_data.merge_wanted_ref_tax_accessions(kept,
                                                         taxon=selection.canonical)
        run_data.record_wanted_ref_tax_selection(selection, taxon=selection.canonical,
                                                 num_added=added,
                                                 num_excluded=num_excluded)

    run_data.wanted_ref_tax_num_excluded = total_excluded

    return run_data


def resolve_hmm(args, selections=None, previous_run_data=None):
    """
    Figure out `-H` being used

    Canonicalizing afterwards, whichever way the name arrived, is what keeps
    `-H universal` and an auto-selected `Universal-Hug-et-al` the same string in the
    fingerprint, in the reporting, and as a filename under GToTree_HMM_dir.
    """
    # ahead of both routes: a viral `-w` target can't take a pre-built set by either one
    try:
        check_viruses_have_their_own_hmms(args, selections)
    except ViralTaxonNeedsOwnHMMs as err:
        report_message(str(err))
        report_very_early_exit(suggest_help=True)

    if not args.hmm:
        carried_over = (previous_run_data.fingerprint.get("hmm")
                        if previous_run_data is not None else None)

        if carried_over:
            args.hmm = carried_over
            args.hmm_auto_selected = "carried over from the run being resumed"
        else:
            picked = autopick_scg_set(args.source, selections)
            args.hmm = picked.name
            args.hmm_auto_selected = picked.reason

    args = resolve_hmm_arg(args)

    return args


def check_input_genome_files(args):
    """
    Validate the single-column input-genome lists, up front so a malformed one fails
    before any asset is downloaded or any taxon resolved.
    """
    if args.ncbi_accessions:
        args.ncbi_accessions = check_expected_single_column_input(args.ncbi_accessions, "-a")
    if args.genbank_files:
        args.genbank_files = check_expected_single_column_input(args.genbank_files, "-g")
    if args.fasta_files:
        args.fasta_files = check_expected_single_column_input(args.fasta_files, "-f")
    if args.amino_acid_files:
        args.amino_acid_files = check_expected_single_column_input(args.amino_acid_files, "-A")

    if getattr(args, "exclusion_list", None):
        if not args.wanted_ref_tax:
            report_message(
                "An `--exclusion-list` was provided, but it only has an effect alongside "
                "`-w`/`--wanted-ref-tax` (it removes accessions from the references `-w` "
                "pulls in). Nothing would be excluded, so this is almost certainly a "
                "mistake."
            )
            report_very_early_exit(suggest_help=True)
        args.exclusion_list = check_expected_single_column_input(
            args.exclusion_list, "--exclusion-list")

    return args


def load_previous_run_data(args):
    """
    The previous run's RunData if we're resuming one and it's on disk, else None
    """
    if not args.resume:
        return None

    run_data_path = os.path.join(args.run_files_dir, "run-data.json")
    if not os.path.exists(run_data_path):
        return None

    try:
        return read_run_data(run_data_path)
    except CorruptRunData as e:
        report_message(
            "We are trying to resume a previous run (specified by the `-R` or `--resume` "
            f"flag), but {e}. That usually means the previous run was interrupted while "
            "saving its state. Your best bet is to start a fresh run by adding the `-F` "
            "flag to force-overwrite the previous outputs.")
        report_very_early_exit()


def setup_run_data(args, previous_run_data=None):
    run_data = previous_run_data

    if run_data is not None:
        differences = RESUME.compare(run_data.fingerprint, build_fingerprint(args))
        if differences:
            report_message(
                "We are trying to resume a previous run (specified by the `-R` or "
                "`--resume` flag), but this run doesn't match the previous one:\n\n"
                f"{bullet_list(differences)}\n\n"
                "Resuming would mix results from two different runs. Your best "
                "bet is to start a fresh run by adding the `-F` flag to "
                "force-overwrite the previous outputs or specify a new output dir.")
            report_very_early_exit()

        if run_data.stage_is_complete(PipelineStage.FINALIZE):
            report_run_already_complete(args.output_dir)

    fresh_run = run_data is None

    if fresh_run:
        run_data = populate_run_data(args)
        run_data.fingerprint = build_fingerprint(args)

    run_data = resolve_hmm_source(args, run_data)

    if fresh_run:
        run_data = populate_SCG_targets(run_data)

    if args.mapping_file:
        args, run_data = check_mapping_file(args, run_data)

    if args.target_pfams_file:
        args.target_pfams_file, total_pfam_targets = check_expected_single_column_input(args.target_pfams_file, "-p", get_count=True)
        run_data.target_pfams_file = args.target_pfams_file
        run_data.total_pfam_targets = total_pfam_targets

    if args.target_kos_file:
        args.target_kos_file, total_ko_targets = check_expected_single_column_input(args.target_kos_file, "-K", get_count=True)
        run_data.target_kos_file = args.target_kos_file
        run_data.total_ko_targets = total_ko_targets

    return args, run_data


GENOME_SOURCE_FLAGS = {
    SOURCE_ACCESSION: "-a",
    SOURCE_GENBANK: "-g",
    SOURCE_FASTA: "-f",
    SOURCE_AMINO_ACID: "-A",
}

# how many colliding ids to spell out before summarizing the rest
MAX_REPORTED_ID_COLLISIONS = 50


def collect_genome_id_collisions(run_data):
    """
    {genome id -> the GenomeData entries that share it}, for ids claimed more than once

    Split out from the reporting so the condition can be tested without going through a
    terminal-and-exit path.
    """
    by_id = {}
    for gd in run_data.all_input_genomes:
        by_id.setdefault(gd.id, []).append(gd)

    return {id: gds for id, gds in by_id.items() if len(gds) > 1}


def check_for_genome_id_collisions(run_data):
    """
    No two input genomes from any source can be allowed to resolve to the same genome id
    """
    collisions = collect_genome_id_collisions(run_data)
    if not collisions:
        return

    plural = "" if len(collisions) == 1 else "s"
    report_message(
        f"{len(collisions)} genome ID{plural} would be claimed by more than one input "
        "genome. Since the ID is used as both a filename and a sequence-header prefix, "
        "those genomes would clash and cause problems.\n\n"
        "Problematic ones include:"
    )

    shown = list(collisions.items())[:MAX_REPORTED_ID_COLLISIONS]
    lines = []
    for id, gds in shown:
        lines.append(id)
        for gd in gds:
            flag = GENOME_SOURCE_FLAGS.get(gd.source, "?")
            lines.append(f"  {flag}  {gd.provided_path or gd.id}")
    report_message("\n".join(lines), ii="    ", si="    ")

    remaining = len(collisions) - len(shown)
    if remaining:
        report_message(f"...and {remaining} more.", ii="    ", si="    ")

    report_message(
        "Each input genome needs to resolve to a unique ID. Renaming the input files so "
        "their names differ (ignoring any `.gz` and the sequence extension) will sort "
        "this out."
    )
    report_very_early_exit()


def check_for_min_input_genomes(run_data):
    total_input_genomes = len(run_data.get_all_input_genome_ids())
    if total_input_genomes < 4:
        word = "was" if total_input_genomes == 1 else "were"
        message = f"\n  {color_text(f"At least 4 input genomes are required, but only {total_input_genomes} {word} detected.\n\n", "yellow")}"
        message += "  See `GToTree -h` for more info.\n\nExiting for now :(\n"
        print(message)
        exit(1)


def check_output_dir(args):
    if os.path.exists(args.output_dir):
        args.output_already_existed = True

        if not args.force_overwrite and not args.resume:
            report_message(f'The specified output directory "{args.output_dir}" already exists. '
                           'Please either remove it, provide the `-F` flag to overwrite it, or '
                           'provide the `-R` flag to attempt to resume a previous run.')
            report_very_early_exit()

        if args.force_overwrite:
            shutil.rmtree(args.output_dir)
            os.makedirs(args.output_dir)
    else:
        args.output_already_existed = False

    args.run_files_dir_rel = os.path.join(args.output_dir, "run-files")
    args.run_files_dir = os.path.abspath(args.run_files_dir_rel)

    return args


def check_path(path, flag):

    if not os.path.isfile(path):
        if flag == "-a" or flag == "-g" or flag == "-f" or flag == "-A":
            report_missing_input_genomes_file(path, flag)
        if flag == "-p":
            report_missing_pfam_targets_file(path, flag)
        if flag == "-K":
            report_missing_ko_targets_file(path, flag)
        if flag == "-m":
            report_missing_mapping_file(path, flag)
        if flag == "--exclusion-list":
            report_missing_exclusion_list_file(path, flag)


def check_expected_single_column_input(path, flag, get_count=False):
    """
    Validate one single-column input file, returning the path unchanged
    """
    check_path(path, flag)
    check_for_whitespace(path, flag)
    report_duplicate_entries(path, flag)
    check_inputs_exist(path, flag)

    if get_count:
        return path, len(read_single_column_file(path))

    return path


def check_for_whitespace(path, flag):

    # there should be no tabs or spaces in these input files
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    for lineno, line in enumerate(lines, start=1):
        for label, char in (("spaces", " "), ("tabs", "\t")):
            if char in line:
                message = (
                    f'The specified input file "{path}" (passed to `{flag}`) contains '
                    f'{label} in one or more entries, which is not expected and will break '
                    f'things downstream.\n\nThe first one is on line {lineno}:\n\n'
                    f'    {line}\n\n'
                    'Please remove the whitespace from that file and try again.')
                report_very_early_exit(message, suggest_help=True)


def check_line_endings(path, flag):
    """
    Rewrite a CRLF file as LF, returning the path to use

    Only the mapping file (`-m`) needs this. The single-column inputs are read
    through `read_single_column_file()`, whose `.strip()` absorbs a trailing `\r`; the
    mapping file goes through `pd.read_csv(sep="\t")`, where the `\r` ends up inside
    the last field of every row and silently corrupts the output labels.
    """

    with open(path, 'rb') as f:
        content = f.read()

    if b'\r\n' in content:
        new_filename = f"{path}-unix"
        new_content = content.replace(b'\r\n', b'\n')
        with open(new_filename, 'wb') as f:
            f.write(new_content)

        report_message(f'Input file "{path}" (passed to `{flag}`) had Windows-formatting that would have caused problems. '
                       f'A modified version was created, "{new_filename}", which will be used.')

        path = new_filename
        time.sleep(2)

    return path


# how many repeated entries to name before summarizing the rest
MAX_REPORTED_DUPLICATES = 10


def report_duplicate_entries(path, flag):
    """
    Tell the user about repeated entries
    (`read_single_column_file()` drops them, so no action needed and no new file created anymore)
    """
    duplicates = duplicate_entries(path)
    if not duplicates:
        return

    plural = "y" if len(duplicates) == 1 else "ies"
    shown = duplicates[:MAX_REPORTED_DUPLICATES]
    lines = list(shown)
    remaining = len(duplicates) - len(shown)
    if remaining:
        lines.append(f"...and {remaining} more.")

    report_message(f'Input file "{path}" (passed to `{flag}`) has {len(duplicates):,} '
                   f'repeated entr{plural}, which will only be used once:')
    report_message("\n".join(lines), ii="    ", si="    ")


def check_inputs_exist(path, flag):
    if flag in ["-g", "-f", "-A"]:
        # checking that all the input files actually exist
        for entry in read_single_column_file(path):
            if not os.path.exists(entry):
                report_message(f'The specified input-genome file "{entry}" (passed to `{flag}`) does not exist.')
                report_very_early_exit()


def check_mapping_file(args, run_data, flag = "-m"):

    check_path(args.mapping_file, flag)
    mapping_file_problems = check_mapping_file_problem_chars_and_fields(args.mapping_file)
    if mapping_file_problems:
        report_problem_with_mapping_file(mapping_file_problems, args.mapping_file, flag = "-m")

    args.mapping_file = check_line_endings(args.mapping_file, flag)

    run_data.mapping_file_path = args.mapping_file

    # this is to handle when resuming runs, so that we don't overwrite the mapping_dict
    if len(run_data.mapping_dict) == 0 and len(run_data.suffix_dict) == 0:
        mapping_dict, suffix_dict = make_mapping_dict(run_data.mapping_file_path)

        check_mapping_dict_for_duplicate_labels(mapping_dict, run_data.mapping_file_path)

        # every key (full-replacement and suffix-only) must map to a real input genome
        check_all_mapping_file_entries_are_in_input_genomes(
            {**mapping_dict, **suffix_dict}, run_data)

        # rekeyed to genome ids here, at the one place the dicts are built, so every
        # consumer downstream can just look up by id
        run_data.mapping_dict = normalize_mapping_dict_keys(mapping_dict, run_data)
        run_data.suffix_dict = normalize_mapping_dict_keys(suffix_dict, run_data)

    return args, run_data


def build_mapping_key_lookup(run_data):
    """
    {form the user may legitimately write in column 1 -> genome id}

    Both the validator and the key normalization work off this one dict, so what
    gets accepted and what can be resolved can't drift apart.

    Accepted forms are the path as provided (which is what a user copying their
    `-g`/`-f`/`-A` listing file will naturally write) and the genome id itself
    (which is what an NCBI accession from `-a` already is, and is also the form
    that shows up in the alignment headers and the summary table).
    """
    lookup = {}
    for gd in run_data.all_input_genomes:
        if gd.provided_path:
            lookup.setdefault(gd.provided_path, gd.id)
        lookup.setdefault(gd.id, gd.id)
    return lookup


def normalize_mapping_dict_keys(mapping_dict, run_data):
    """
    Rekey the user's mapping dict from whatever they wrote in column 1 to the
    genome id.

    Without this, `mapping_dict` ends up with two key spaces in it: the user's
    entries keyed by the raw column-1 string (a path like "genome-1.faa" for
    file-based input) and the entries the GTDB/NCBI taxonomy handlers add later,
    keyed by genome id. Nothing that consumes the dict can then look anything up
    with one form of key.

    Key order is preserved, so a later duplicate-label check or report still sees
    the user's file order.
    """
    lookup = build_mapping_key_lookup(run_data)
    return {lookup.get(key, key): value for key, value in mapping_dict.items()}


def check_mapping_file_problem_chars_and_fields(path):

    problematic_chars='()*&^#$@!/\\|[]:;'
    errors = []
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            columns = line.split('\t')

            # Check if the number of columns is either 2 or 3.
            if len(columns) not in (2, 3):
                errors.append(f"Line {lineno} has {len(columns)} column(s) (expected 2 or 3): {line}")

            # Check each field for any problematic characters.
            for i, field in enumerate(columns, start=1):
                for char in problematic_chars:
                    if char not in field:
                        continue
                    # "/" is the one allowed exception, and only in the first column
                    # (input genome paths legitimately contain it)
                    if i == 1 and char == "/":
                        continue
                    errors.append(
                        f"Line {lineno}, column {i} contains at least one problematic character ('{char}'):\n  {field}"
                    )
                    break
    return errors


def make_mapping_dict(path):

    # makes a dictionary mapping input genomes to the wanted labels, plus a
    # separate dictionary of suffix-only entries to be appended later
    # (after any taxonomy info is added)

    #   1. Checks that every first-column entry is unique
    #   2. Uses the first column as the key
    #   3. If there are 2 columns, the value is the second column (a full
    #      replacement label)
    #   4. If there are 3 columns and the second column is non-empty, the value
    #      is "second-column_third-column" (a full replacement label with the
    #      suffix baked in)
    #   5. If there are 3 columns and the second column is empty, the third
    #      column is a suffix only. It is NOT baked into the mapping value here,
    #      because the base label may still be built from taxonomy info that
    #      isn't known yet. The suffix is stored in `suffix_dict` and appended
    #      later (see append_pending_suffixes); if no taxonomy is added, it gets
    #      appended to the original input label.

    df = pd.read_csv(path, sep="\t", header=None, dtype=str, names=[0, 1, 2]).fillna("")

    # checking none of the first column entries are duplicates
    if df[0].duplicated().any():
        report_message(
            f'The mapping file "{path}" (passed to `-m`) has some duplicate entries in the first column.'
            ' Please address that and try again.'
        )
        report_very_early_exit()

    mapping_dict = {}
    suffix_dict = {}
    for _idx, row in df.iterrows():
        key = row[0].strip()
        col2 = row[1].strip() if len(row) > 1 else ""
        col3 = row[2].strip() if len(row) > 2 else ""

        if col2:
            # full replacement label; suffix (if any) baked in now
            value = f"{col2}_{col3}" if col3 else col2
            mapping_dict[key] = value.replace(" ", "_")
        elif col3:
            # suffix only; deferred so a taxonomy-built base can be used
            suffix_dict[key] = col3.replace(" ", "_")
        # else: empty row for this key, nothing to do

    return mapping_dict, suffix_dict


def check_mapping_dict_for_duplicate_labels(mapping_dict, path):
    # checking none of the desired labels are duplicates
    counts = Counter(mapping_dict.values())
    duplicates = [val for val, count in counts.items() if count > 1]
    if duplicates:
        report_message(
            f'The mapping file "{path}" (passed to `-m`) specifies to have duplicate output genome labels.\n\n'
            f'Problematic ones include:'
        )
        report_message(f'{"\n".join(duplicates)}', ii="    ", si="    ")
        report_message(
            'Each input genome must map to a unique label. Please address that and try again.'
        )
        report_very_early_exit()


def check_all_mapping_file_entries_are_in_input_genomes(mapping_dict, run_data):
    entries_in_mapping_file = set(mapping_dict.keys())
    accepted_entries = set(build_mapping_key_lookup(run_data))
    missing_keys = entries_in_mapping_file - accepted_entries
    if missing_keys:
        report_message(
            f'The mapping file "{run_data.mapping_file_path}" (passed to `-m`) specifies some input genomes that are not found in any of the the input-genome sources.\n\n'
            f'Problematic ones include:'
        )
        report_message(f'{"\n".join(missing_keys)}', ii="    ", si="    ")
        report_message(
            "Each input genome in the mapping file (column 1) must be present in one of the input-genome sources. Please address this and try again."
        )
        report_very_early_exit()


def check_for_required_dbs(args, previous_run_data=None):

    reusing_wanted_ref_tax = wanted_ref_tax_already_resolved(previous_run_data)

    if args.wanted_ref_tax and not reusing_wanted_ref_tax:
        wanted_ref_tax_source = args.source.strip().lower()
    else:
        wanted_ref_tax_source = None

    auto_selecting_hmm = (bool(args.wanted_ref_tax) and not args.hmm
                          and not reusing_wanted_ref_tax)

    if (args.ncbi_accessions or args.add_ncbi_tax
            or wanted_ref_tax_source in ("ncbi", "gtdb")):
        get_ncbi_assembly_data()
    if args.add_gtdb_tax or wanted_ref_tax_source == "gtdb" or auto_selecting_hmm:
        get_gtdb_data()
    if args.target_kos_file:
        get_kofamscan_data()
    if args.target_pfams_file:
        get_pfam_data()


def track_tools_used(args, run_data):

    tools_used = ToolsUsed()

    if args.fasta_files:
        tools_used.prodigal_used = True
    if args.add_gtdb_tax:
        tools_used.gtdb_used = True
    if args.tree_program == "FastTreeMP" or args.tree_program == "FastTree":
        tools_used.fasttree_used = True
    if args.tree_program == "VeryFastTree":
        tools_used.veryfasttree_used = True
        tools_used.fasttree_used = True
    if args.tree_program == "IQTREE":
        tools_used.iqtree_used = True
    if args.hmm == "Universal" or args.hmm == "Universal-Hug-et-al" or args.hmm == "Universal-Hug-et-al.hmm":
        tools_used.universal_SCGs_used = True
    if args.target_pfams_file:
        tools_used.pfam_db_used = True
    if args.target_kos_file:
        tools_used.kofamscan_used = True

    run_data.tools_used = tools_used

    return run_data


def check_input_genomes_amount(total_input_genomes, args):
    if total_input_genomes >= 500 and total_input_genomes < 12500 and not args.no_super5:
        message = many_genomes_notice(total_input_genomes)
        report_notice(message)
        time.sleep(30)
    if total_input_genomes <= 20:
        message = few_genomes_notice(total_input_genomes, args)
        report_notice(message)
        # time.sleep(5)
    if total_input_genomes >= 12500:
        message = absurd_number_of_genomes_notice(total_input_genomes)
        report_notice(message)
        time.sleep(60)


def final_setups(args, run_data):

    log_file = os.path.abspath(os.path.join(args.output_dir, "gtotree-runlog.txt"))
    args.log_file = log_file
    run_data.log_file = log_file
    log_file_var.set(log_file)

    os.makedirs(args.run_files_dir, exist_ok=True)

    full_execution_command = f"{' '.join(sys.argv)}"
    stdout_and_log(gtotree_header(), log_file=args.log_file, log_only=True, restart_log=True)
    stdout_and_log("    Command executed:\n       ", full_execution_command, "\n", log_file=args.log_file, log_only=True)
    stdout_and_log(RUN_INFO_BANNER, log_file=args.log_file, log_only=True,
                   end="\n" * (RUN_INFO_BANNER_TRAILING_BLANKS + 1))

    args, run_data = setup_tmp_dir(args, run_data)

    logs_dir_rel = os.path.join(args.output_dir, "logs")
    logs_dir = os.path.abspath(logs_dir_rel)
    os.makedirs(logs_dir, exist_ok=True)
    run_data.logs_dir_rel = logs_dir_rel
    run_data.logs_dir = logs_dir
    run_data.output_dir_rel = args.output_dir
    run_data.gtotree_logs_dir = os.path.join(logs_dir, "gtotree-logs")
    os.makedirs(run_data.gtotree_logs_dir, exist_ok=True)

    if run_data.ncbi_accs:
        ncbi_processing_dir_rel = os.path.join(args.tmp_dir, "ncbi-acc-processing")
        ncbi_processing_dir = os.path.abspath(ncbi_processing_dir_rel)
        os.makedirs(ncbi_processing_dir, exist_ok=True)
        run_data.ncbi_processing_dir_rel = ncbi_processing_dir_rel
        run_data.ncbi_processing_dir = ncbi_processing_dir

    if args.genbank_files:
        genbank_processing_dir = os.path.join(args.tmp_dir, "genbank-processing")
        genbank_processing_dir = os.path.abspath(genbank_processing_dir)
        os.makedirs(genbank_processing_dir, exist_ok=True)
        run_data.genbank_processing_dir = genbank_processing_dir

    if args.fasta_files:
        fasta_processing_dir = os.path.join(args.tmp_dir, "fasta-processing")
        fasta_processing_dir = os.path.abspath(fasta_processing_dir)
        os.makedirs(fasta_processing_dir, exist_ok=True)
        run_data.fasta_processing_dir = fasta_processing_dir

    if args.amino_acid_files:
        AA_processing_dir = os.path.join(args.tmp_dir, "amino-acid-processing")
        AA_processing_dir = os.path.abspath(AA_processing_dir)
        os.makedirs(AA_processing_dir, exist_ok=True)
        run_data.AA_processing_dir = AA_processing_dir

    run_data.ready_genome_files_dir = os.path.join(args.tmp_dir, "ready-genome-files")
    os.makedirs(run_data.ready_genome_files_dir, exist_ok=True)
    run_data.hmm_results_dir = os.path.join(args.tmp_dir, "hmm-results")
    os.makedirs(run_data.hmm_results_dir, exist_ok=True)
    run_data.found_SCG_seqs_dir = os.path.join(args.tmp_dir, "found-SCG-seqs")
    os.makedirs(run_data.found_SCG_seqs_dir, exist_ok=True)

    run_data.best_hit_mode = args.best_hit_mode
    run_data.output_dir = os.path.abspath(args.output_dir)
    run_data.seq_length_cutoff = args.seq_length_cutoff
    run_data.gene_representation_cutoff = args.gene_representation_cutoff

    num_all_input_genomes = len(run_data.get_all_input_genome_ids())
    if num_all_input_genomes > 500 and not args.no_super5:
        run_data.use_muscle_super5 = True

    run_data.num_muscle_threads = args.num_muscle_threads
    run_data.nucleotide_mode = args.nucleotide_mode
    run_data.general_ext = ".fasta" if run_data.nucleotide_mode else ".faa"

    if len(run_data.mapping_dict) > 0 or len(run_data.suffix_dict) > 0 or args.add_ncbi_tax or args.add_gtdb_tax:
        run_data.updating_headers = True
    else:
        run_data.updating_headers = False

    if args.keep_gene_alignments:
        run_data.individual_gene_alignments_dir_rel = os.path.join(args.run_files_dir_rel, "individual-alignments")
        run_data.individual_gene_alignments_dir = os.path.abspath(run_data.individual_gene_alignments_dir_rel)
        os.makedirs(run_data.individual_gene_alignments_dir, exist_ok=True)

    if args.target_pfams_file:
        run_data = setup_pfam_dirs(run_data)

    if args.target_kos_file:
        run_data = setup_ko_dirs(run_data)

    return args, run_data


def setup_pfam_dirs(run_data):
    run_data.pfam_results_dir = os.path.join(run_data.output_dir, "pfam-search-results")
    run_data.pfam_results_dir_rel = os.path.join(run_data.output_dir_rel, "pfam-search-results")
    run_data.tmp_pfam_results_dir = os.path.join(run_data.tmp_dir, "pfam-search-results")


    dirs = [run_data.pfam_results_dir,
            run_data.tmp_pfam_results_dir,
            os.path.join(run_data.pfam_results_dir, "info"),
            os.path.join(run_data.pfam_results_dir, "individual-genome-results"),
            os.path.join(run_data.pfam_results_dir, "itol-files"),
            os.path.join(run_data.pfam_results_dir, "pfam-hit-seqs"),
            os.path.join(run_data.pfam_results_dir, "target-pfam-profiles")]

    for d in dirs:
        os.makedirs(d, exist_ok=True)

    return run_data


def setup_ko_dirs(run_data):
    run_data.ko_results_dir = os.path.join(run_data.output_dir, "ko-search-results")
    run_data.ko_results_dir_rel = os.path.join(run_data.output_dir_rel, "ko-search-results")
    run_data.tmp_ko_results_dir = os.path.join(run_data.tmp_dir, "ko-search-results")

    dirs = [run_data.ko_results_dir,
            run_data.tmp_ko_results_dir,
            os.path.join(run_data.ko_results_dir, "individual-genome-results"),
            os.path.join(run_data.ko_results_dir, "itol-files"),
            os.path.join(run_data.ko_results_dir, "ko-hit-seqs"),
            os.path.join(run_data.ko_results_dir, "target-ko-profiles")]

    for d in dirs:
        os.makedirs(d, exist_ok=True)

    return run_data


def setup_tmp_dir(args, run_data):

    if args.resume:
        if run_data.tmp_dir:
            args.tmp_dir = run_data.tmp_dir
        else:
            tmp_dir = tempfile.mkdtemp(prefix = "gtt-tmp-", dir = args.output_dir)
            args.tmp_dir = tmp_dir
            run_data.tmp_dir = tmp_dir

        return args, run_data

    if args.tmp_dir:
        try:
            os.makedirs(args.tmp_dir, exist_ok=True)
            tmp_dir = tempfile.mkdtemp(dir = args.tmp_dir)
            args.tmp_dir = tmp_dir
            run_data.tmp_dir = tmp_dir
        except OSError:
            report_message(f"We could not create a temporary directory in the location you specified: {args.tmp_dir}")
            report_message("Maybe you don't have write permissions there?")
            report_very_early_exit()
    else:
        tmp_dir = tempfile.mkdtemp(prefix = "gtt-tmp-", dir = args.output_dir)
        args.tmp_dir = tmp_dir
        run_data.tmp_dir = tmp_dir

    return args, run_data
