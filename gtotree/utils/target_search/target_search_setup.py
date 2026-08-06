"""
Setup and validation for `gtt search-pfams` / `gtt search-kos`.

The main driver's `preflight_checks()` can't be reused wholesale here: it requires an
SCG-HMM set (`-H`), insists muscle and trimal are on PATH, validates a treeing program,
and enforces a four-genome floor. None of that applies when all we're doing is
searching genomes for a list of targets -- searching a single genome for one Pfam is a
perfectly reasonable thing to ask for.

So this follows the `gen_scg_hmms` precedent: a small `check_args` plus setup that
calls the *shared* preflight helpers directly (`check_expected_single_column_input`,
`setup_tmp_dir`) rather than threading a mode flag through the main preflight. The
input-file validation, genome parsing, and directory conventions all still come from
the same code the main run uses.

The output layout is flattened by one level relative to a full GToTree run: what would
land in `<out>/pfam-search-results/` lands directly in `<out>/`, since the whole output
directory is already the search results.
"""

import os
import shutil
import argparse

from gtotree.utils.general import GenomeData, RunData, ToolsUsed
from gtotree.utils.messaging import report_message


class TargetSearchError(Exception):
    """
    A user-facing problem with the run, translated to a friendly message + clean exit
    by the CLI layer (the same library-raises / CLI-translates split used everywhere
    else in the package).
    """


WORKING_DIR_NAME = "working-dir"
RUN_DATA_FILENAME = "run-data.json"
SUMMARY_FILENAME = "genomes-summary-info.tsv"


################################################################################
# argument validation
################################################################################

def check_args(args, spec):
    """
    Validate arguments, raising TargetSearchError with a friendly message.
    """
    input_flags = (args.ncbi_accessions, args.genbank_files, args.fasta_files,
                   args.amino_acid_files, args.wanted_ref_tax)
    if not any(input_flags):
        raise TargetSearchError(
            "We need some input genomes to work with! Provide any combination of an "
            "NCBI accessions file (`-a`), GenBank files (`-g`), fasta files (`-f`), "
            "amino-acid files (`-A`), and/or a target taxon (`-w`).")

    if not spec.targets_file(args):
        raise TargetSearchError(
            f"We need to know which {spec.target_label_plural} to search for! Provide a "
            f"single-column file of {spec.target_label} IDs (e.g. "
            f"'{spec.example_target}') with `{spec.targets_flag}`.")

    if args.num_jobs < 1:
        raise TargetSearchError("The `--num-jobs` (-j) parameter needs to be at least 1.")

    if args.resume and args.force_overwrite:
        raise TargetSearchError(
            "`-R`/`--resume` and `-F`/`--force-overwrite` can't be used together, one "
            "reuses the previous run and the other deletes it.")

    _check_dangling_ref_tax_args(args)

    return args


def _check_dangling_ref_tax_args(args):
    """
    `--target-rank` and `--derep-rank` only mean something alongside `-w`. Caught up
    front so a misunderstanding fails before any dataset is fetched.
    """
    if args.wanted_ref_tax:
        return

    dangling = []
    if args.target_rank is not None:
        dangling.append("`--target-rank`")
    if args.derep_rank not in (None, "auto"):
        dangling.append("`--derep-rank`")

    if dangling:
        joined = " and ".join(dangling)
        raise TargetSearchError(
            f"You provided {joined}, but that only applies when adding reference "
            "genomes by taxonomy with `-w`/`--wanted-ref-tax`.")


################################################################################
# dependency checks
################################################################################

def check_dependencies(args, spec):
    """
    Check only what this run actually needs.

    Deliberately separate from the main driver's `check_for_essential_deps()`, which
    demands muscle and trimal -- neither subcommand aligns or trims anything, so
    requiring them would turn a working setup into a spurious failure.

    prodigal is conditional: it's only invoked for nucleotide fastas (`-f`), for NCBI
    accessions with no protein file, and as the GenBank fallback. Amino-acid input
    alone never needs it.
    """
    missing = []

    for binary in spec.required_binaries:
        if not shutil.which(binary):
            missing.append(binary)

    if _might_need_prodigal(args) and not shutil.which("prodigal"):
        missing.append("prodigal")

    if missing:
        joined = ", ".join(missing)
        was = "is" if len(missing) == 1 else "are"
        raise TargetSearchError(
            f"{joined} {was} needed for this run, but not in your PATH :(")


def _might_need_prodigal(args):
    """
    True when any input source could route through gene-calling.

    NCBI accessions are included because an assembly with no protein fasta falls back
    to prodigal, and GenBank files because a file with no usable CDS translations does
    the same. Only an amino-acid-only run is guaranteed not to need it.
    """
    return bool(args.fasta_files or args.genbank_files
                or args.ncbi_accessions or args.wanted_ref_tax)


def check_env_vars(spec):
    """
    Confirm the data-location environment variables this target type reads directly.

    `spec.ensure_data()` sets these up as a side effect of fetching the managed
    dataset, so this runs after it as a guard rather than instead of it -- the search
    helpers read the variables straight out of `os.environ` and would otherwise fail
    deep inside a worker thread with a bare KeyError.
    """
    for var in spec.required_env_vars:
        if not os.environ.get(var):
            raise TargetSearchError(
                f"The environment variable '{var}' doesn't seem to be set :(\n\n"
                "      This shouldn't happen -- check on things with "
                "`gtt data locations check`.")


################################################################################
# output directory
################################################################################

def setup_output_dir(args, spec):
    """
    Create the output dir and its working dir, honoring -F and -R.

    Returns (out_dir, work_dir). The working dir holds the resume state, the saved
    run-data, and all intermediates; it's removed on success unless
    `--keep-working-dir` was given, which is also what makes a resumed run possible
    after a failure.
    """
    out_dir = args.output_dir.rstrip("/")
    work_dir = os.path.join(out_dir, WORKING_DIR_NAME)

    if os.path.exists(out_dir):
        if args.resume:
            if not os.path.isdir(work_dir):
                # an output dir with no working dir is a run that never got far enough
                # to leave anything resumable, so start fresh rather than refusing
                report_message(
                    f"There's no working directory in '{out_dir}' from a previous run "
                    "to resume from, so we'll start fresh.", "yellow",
                    ii="      ", si="      ")
            os.makedirs(work_dir, exist_ok=True)
            _make_result_dirs(out_dir, spec)
            return out_dir, work_dir

        if not args.force_overwrite:
            raise TargetSearchError(
                f"The output directory '{out_dir}' already exists, and we don't want "
                "to overwrite anything accidentally. Use `-R` to continue that run, "
                "`-F` to overwrite it, or specify a different directory with `-o`.")
        shutil.rmtree(out_dir)

    elif args.resume:
        report_message(
            f"`-R`/`--resume` was specified, but '{out_dir}' doesn't exist yet, so "
            "we'll just start fresh.", "yellow", ii="      ", si="      ")

    os.makedirs(work_dir, exist_ok=True)
    _make_result_dirs(out_dir, spec)

    return out_dir, work_dir


def _make_result_dirs(out_dir, spec):
    for sub in spec.result_subdirs:
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)


################################################################################
# run_data
################################################################################

def build_run_data(args, spec, out_dir, work_dir, previous=None):
    """
    Build (or adopt) the RunData this run will thread through the shared machinery.

    The directory attributes are deliberately pointed at the flattened layout:

        results_dir     -> <out_dir>            (not <out_dir>/pfam-search-results)
        run_files_dir   -> <out_dir>            (so failed-targets files land on top)
        tmp_dir         -> <out_dir>/working-dir/tmp
        run_data_path   -> <out_dir>/working-dir/run-data.json

    `run_files_dir` and the results dir being the same directory is intentional. In a
    full GToTree run they're separate because there's a whole tree pipeline's worth of
    other run files to keep apart; here the entire output *is* the search results, so
    splitting them would just add a level back.

    The genome set is ALWAYS re-parsed from the input files, even when resuming. This
    matters: adopting the previous run's genome list wholesale would mean edits to the
    listing files were silently ignored -- a genome added between runs would never be
    processed, and the resume fingerprint (built from run_data) would compare the old
    set against itself and see no change to refuse over.

    What `previous` contributes instead is run-*level* state: the resolved target set
    and the assets it produced. Per-genome progress is folded in separately by
    `adopt_genome_progress`, once the genome set is final (i.e. after any `-w`
    accessions have been merged in).
    """
    run_data = RunData()
    _populate_input_genomes(args, run_data)
    run_data.tools_used = ToolsUsed()

    if previous is not None:
        _adopt_run_level_state(run_data, previous, spec)

    tmp_dir = os.path.join(work_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    run_data.output_dir = os.path.abspath(out_dir)
    run_data.output_dir_rel = out_dir
    run_data.run_files_dir = os.path.abspath(out_dir)
    run_data.run_files_dir_rel = out_dir
    run_data.run_data_path = os.path.join(work_dir, RUN_DATA_FILENAME)
    run_data.tmp_dir = tmp_dir
    run_data.log_file = os.path.join(work_dir, "gtotree-runlog.txt")

    # searching is amino-acid only here; nucleotide mode belongs to the tree pipeline
    run_data.nucleotide_mode = False
    run_data.general_ext = ".faa"

    setattr(run_data, spec.results_dir_attr, os.path.abspath(out_dir))
    setattr(run_data, spec.results_dir_rel_attr, out_dir)
    setattr(run_data, spec.tmp_results_dir_attr,
            os.path.join(tmp_dir, "target-hit-seqs"))
    setattr(run_data, spec.targets_file_attr, spec.targets_file(args))
    # set unconditionally, including on a resumed run: it's how the reporting says
    # "N of M requested targets found", and a resumed RunData predates this run's
    # (possibly re-validated) target file
    setattr(run_data, spec.total_targets_attr, getattr(args, "total_targets", 0) or 0)

    run_data.ready_genome_files_dir = os.path.join(tmp_dir, "ready-genome-files")

    dirs = [getattr(run_data, spec.tmp_results_dir_attr),
            run_data.ready_genome_files_dir]

    if run_data.ncbi_accs:
        run_data.ncbi_processing_dir = os.path.join(tmp_dir, "ncbi-acc-processing")
        run_data.ncbi_processing_dir_rel = run_data.ncbi_processing_dir
        dirs.append(run_data.ncbi_processing_dir)
    if run_data.genbank_files:
        run_data.genbank_processing_dir = os.path.join(tmp_dir, "genbank-processing")
        dirs.append(run_data.genbank_processing_dir)
    if run_data.fasta_files:
        run_data.fasta_processing_dir = os.path.join(tmp_dir, "fasta-processing")
        dirs.append(run_data.fasta_processing_dir)
    if run_data.amino_acid_files:
        run_data.AA_processing_dir = os.path.join(tmp_dir, "amino-acid-processing")
        dirs.append(run_data.AA_processing_dir)

    for d in dirs:
        os.makedirs(d, exist_ok=True)

    return run_data


# run-level fields worth carrying across a resume, beyond the spec-named ones.
# Everything else is either re-derived each run or a path this run relays out itself.
_SHARED_CARRIED_FIELDS = ("pfam_dict", "all_pfam_targets_hmm_path", "target_kos_tsv",
                          "target_ko_profiles_dir", "wanted_ko_targets",
                          "wanted_pfam_targets", "tax_info_dict")

# GenomeData fields that identify a genome rather than describe its progress; these
# come from the current input files and must never be overwritten by a previous run
_GENOME_IDENTITY_FIELDS = frozenset(
    {"id", "source", "full_path", "provided_path", "basename"})


def _adopt_run_level_state(run_data, previous, spec):
    """
    Carry the previous run's resolved target set and its assets onto a fresh RunData.

    Whether any of it is actually trusted is decided later by
    `target_stage_is_reusable`, which checks the artifacts are still on disk -- this
    just makes the state available to check.
    """
    for attr in (spec.found_targets_attr, spec.failed_targets_attr,
                 spec.searching_done_attr, *_SHARED_CARRIED_FIELDS):
        value = getattr(previous, attr, None)
        if value:
            setattr(run_data, attr, value)

    # the NCBI sub-table lets parse_assembly_summary early-return, but only if the file
    # it points at survived; a stale path would make it skip the lookup and then fail
    # to find any download links
    sub_table = getattr(previous, "ncbi_sub_table_path", "")
    if sub_table and os.path.isfile(sub_table):
        run_data.ncbi_sub_table_path = sub_table

    if getattr(previous, "tools_used", None):
        run_data.tools_used = previous.tools_used

    return run_data


def adopt_genome_progress(run_data, previous):
    """
    Copy per-genome progress flags from a resumed run onto the current genome set.

    Matched by genome ID, which is stable across runs because it's derived by the same
    `GenomeData.from_path` / `from_acc` factories from the same inputs. A genome present
    now but not before simply keeps its fresh (all-false) flags and gets processed.

    Called after the genome set is final rather than inside `build_run_data`, because
    `-w` accessions are merged in during phase 1 and would otherwise come back with
    their progress wiped.

    Returns the number of genomes whose progress was carried over.
    """
    if previous is None:
        return 0

    from dataclasses import fields

    by_id = {gd.id: gd for gd in previous.all_input_genomes}
    carried = 0

    for gd in run_data.all_input_genomes:
        old = by_id.get(gd.id)
        # source is checked too: the same basename arriving as a genbank file in one
        # run and a fasta in the next is a different genome, and reusing its flags
        # would skip work that needs redoing
        if old is None or old.source != gd.source:
            continue
        for field_info in fields(gd):
            if field_info.name in _GENOME_IDENTITY_FIELDS:
                continue
            setattr(gd, field_info.name, getattr(old, field_info.name))
        carried += 1

    return carried


def _populate_input_genomes(args, run_data):
    """
    Turn the four input flags into GenomeData lists.

    Mirrors `general.populate_run_data`, which can't be called directly because it
    reads `run_files_dir` off args and sets run-data paths this program lays out
    differently. The GenomeData construction itself goes through the same
    `from_acc`/`from_path` factories, so a given file yields the same genome ID here as
    it would in a full GToTree run or in gen-scg-hmms.
    """
    if args.ncbi_accessions:
        run_data.ncbi_accs = [GenomeData.from_acc(entry)
                              for entry in _read_entries(args.ncbi_accessions)]

    if args.genbank_files:
        run_data.genbank_files = [GenomeData.from_path(entry, "genbank-file")
                                  for entry in _read_entries(args.genbank_files)]

    if args.fasta_files:
        run_data.fasta_files = [GenomeData.from_path(entry, "nt-fasta-file")
                                for entry in _read_entries(args.fasta_files)]
        for gd in run_data.fasta_files:
            gd.prodigal_used = True

    if args.amino_acid_files:
        run_data.amino_acid_files = [GenomeData.from_path(entry, "aa-fasta-file")
                                     for entry in _read_entries(args.amino_acid_files)]

    run_data.update_all_input_genomes()


def _read_entries(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


################################################################################
# input-file validation
################################################################################

def validate_input_files(args, spec):
    """
    Run the shared single-column input checks over every provided input file.

    These are the same checks the main driver applies -- whitespace, CRLF line endings,
    duplicate entries, and (for the file-listing flags) that every listed genome file
    actually exists -- so an input list that works with `GToTree` works here.
    """
    from gtotree.utils.preflight_checks import check_expected_single_column_input

    for dest, flag in (("ncbi_accessions", "-a"),
                       ("genbank_files", "-g"),
                       ("fasta_files", "-f"),
                       ("amino_acid_files", "-A")):
        value = getattr(args, dest, None)
        if value:
            setattr(args, dest, check_expected_single_column_input(value, flag))

    targets_path = spec.targets_file(args)
    if targets_path:
        checked, total = check_expected_single_column_input(
            targets_path, spec.targets_flag, get_count=True)
        setattr(args, spec.targets_dest, checked)
        args.total_targets = total

    return args


################################################################################
# args shim
################################################################################

def fill_in_shared_args(args):
    """
    Add the attributes the shared machinery reads off `args` but this CLI doesn't ask
    for.

    `run_pooled_stage` reads `num_jobs`; `build_search_plan` reads `debug`; the
    preprocessing helpers reach `nucleotide_mode` through run_data but a couple of
    shared paths still check args. Setting them explicitly here (rather than letting an
    AttributeError surface from inside a worker thread) keeps the shim in one visible
    place.
    """
    defaults = {
        "nucleotide_mode": False,
        "debug": False,
        "target_pfams_file": None,
        "target_kos_file": None,
        "add_ncbi_tax": False,
        "add_gtdb_tax": False,
        "mapping_file": None,
        "best_hit_mode": False,
    }
    for key, value in defaults.items():
        if not hasattr(args, key):
            setattr(args, key, value)
    return args


def make_args(**overrides):
    """
    Build a fully populated args Namespace, for tests and programmatic use.
    """
    base = {
        "ncbi_accessions": None,
        "genbank_files": None,
        "fasta_files": None,
        "amino_acid_files": None,
        "wanted_ref_tax": None,
        "source": "gtdb",
        "target_rank": None,
        "derep_rank": "auto",
        "output_dir": "gtt-search-output",
        "num_jobs": 4,
        "resume": False,
        "force_overwrite": False,
        "keep_working_dir": False,
        "debug": False,
    }
    base.update(overrides)
    return fill_in_shared_args(argparse.Namespace(**base))
