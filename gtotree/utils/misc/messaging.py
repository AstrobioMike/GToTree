import os
import textwrap
import sys
import time
import re
import contextlib
import shutil
import itertools
import threading
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from gtotree.utils.misc.context import log_file_var
from gtotree.utils.misc.stages import (GenomeRemovalStage,
                                  SCGRemovalStage,
                                  SCG_GENE_FILTERING_STAGES)


tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'orange': '\033[38;5;208m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m',
    'teal' : '\033[0;36m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m'
}

color_codes = {
    'green': '32',
    'yellow': '33',
    'red': '31',
    'teal': '36',
    'orange': '38;5;208'
}


def color_text(text, color = 'green', bold = False):
    if not sys.stdout.isatty():
        return text

    code = color_codes.get(color, color_codes['green'])
    prefix = f"\033[1;{code}m" if bold else f"\033[0;{code}m"

    return f"{prefix}{text}\033[0m"


def get_version():
    """
    The installed GToTree version, or "(dev)" if distribution metadata isn't present
    """
    try:
        return version("GToTree")
    except PackageNotFoundError:
        return "(dev)"


def gtotree_header():
    header = f"""

                                  {color_text(f"GToTree v{get_version()}", "green")}
                          github.com/AstrobioMike/GToTree
    """
    return header


def wprint(text, width = 80, ii = "  ", si = "  "):
    wrapper = textwrap.TextWrapper(width=width,
                                   initial_indent=ii,
                                   subsequent_indent=si,
                                   break_on_hyphens=False)
    paragraphs = text.splitlines()
    wrapped_paragraphs = [wrapper.fill(par) if par.strip() else par for par in paragraphs]
    print("\n".join(wrapped_paragraphs), flush=True)


class Tee:
    def __init__(self, *files):
        self.files = files

    ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')

    def write(self, s):
        for f in self.files:
            if hasattr(f, "isatty") and f.isatty():
                f.write(s)
            else:
                f.write(self.ansi_escape.sub("", s))

    def flush(self):
        for f in self.files:
            f.flush()

    def isatty(self):
        if self.files and hasattr(self.files[0], "isatty"):
            return self.files[0].isatty()
        return False


def _real_terminal_stream():
    """
    Return the underlying real terminal stream to animate the spinner on, or
    None if there isn't one.

    During a GToTree run stdout/stderr are `Tee`s that mirror everything into
    `gtotree-runlog.txt`. We do not want the per-frame carriage-return redraws
    written to that log file, so the animation is sent to the real terminal
    directly (unwrapping the Tee to find it). If nothing in the chain is an
    interactive tty, we animate nowhere and the final line still prints normally.
    """
    stream = sys.stderr
    if isinstance(stream, Tee):
        for f in stream.files:
            if hasattr(f, "isatty") and f.isatty():
                return f
        return None
    if hasattr(stream, "isatty") and stream.isatty():
        return stream
    return None


@contextlib.contextmanager
def spinner(in_progress_msg, complete_msg, indent="    ", clear_on_done=False, reclaim_line=False):  # pragma: no cover
    """
    Show a spinner while a block runs, then print a single completion line.

    Ported from bit's `spinner`, adapted for GToTree's logging: the animated
    frames are written straight to the real terminal (never the run log), while
    the final `✔ complete_msg` line goes through the normal stream so it lands
    in the log too. Elapsed time is appended only if the block took >= 60 s.

    Usage:
        with spinner("Gathering references...", "References gathered"):
            ...work...
    """
    done = threading.Event()
    term = _real_terminal_stream()

    def spin():
        if term is None:
            done.wait()
            return
        for char in itertools.cycle("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"):
            if done.is_set():
                break
            term.write(f"\r{indent}{char} {in_progress_msg}")
            term.flush()
            time.sleep(0.1)
        # clear the animated line so the completion line (or an error) starts clean
        term.write("\r" + " " * (len(indent) + len(in_progress_msg) + 4) + "\r")
        if clear_on_done and reclaim_line:
            term.write("\033[F")
        term.flush()

    t = threading.Thread(target=spin)
    t.start()
    start_time = time.monotonic()
    failed = False
    try:
        yield
    except BaseException:
        failed = True
        raise
    finally:
        done.set()
        t.join()
        # only announce completion on success; on failure the caller's own
        # error/exit message takes over on a freshly-cleared line. When
        # clear_on_done is set, the line has already been wiped and we print
        # nothing further, leaving the cursor on a clean line.
        if not failed and not clear_on_done:
            elapsed = time.monotonic() - start_time
            if elapsed >= 60:
                mins, secs = divmod(int(elapsed), 60)
                time_str = f" (took ~{mins} min and {secs} sec)"
            else:
                time_str = ""
            # final line goes through the normal stream so the run log records it
            print(f"{indent}✔ {complete_msg}{time_str}", flush=True)


def capture_stdout_to_log(log_file):
    def decorator(func):
        def wrapper(*args, **kwargs):
            if isinstance(sys.stdout, Tee):
                return func(*args, **kwargs)
            file_path = log_file() if callable(log_file) else log_file
            try:
                handle = open(file_path, 'a', buffering=1)
            except OSError:
                # The log being unavailable (directory removed, disk full, read-only
                # mount) must not take down the run: these wrappers sit on the
                # *reporting* functions, so raising here would turn a cosmetic problem
                # into a crash -- and often into a crash that hides the message the
                # user actually needed to see. Report to the terminal only.
                return func(*args, **kwargs)
            try:
                tee = Tee(sys.stdout, handle)
                with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
                    return func(*args, **kwargs)
            finally:
                handle.close()
        return wrapper
    return decorator


def format_runtime(start: datetime, now: datetime) -> str:
    elapsed = now - start
    hours = elapsed.seconds // 3600
    minutes = (elapsed.seconds % 3600) // 60
    return f"{hours} hours and {minutes} minutes"


def report_time_status(start: datetime):
    now = datetime.now()
    print("")
    print(f"It is currently {now:%I:%M %p}; the process started at {start:%I:%M %p}.".center(82))
    print(f"Current process runtime: {format_runtime(start, now)}.".center(82))


def report_final_time_status(start: datetime):
    now = datetime.now()
    run_time = format_runtime(start, now)
    print("")
    print(f"Total program runtime: {run_time}".center(82))
    goodbye_message = f"Happy {now.strftime("%A")} :)".center(82)
    print(f"{color_text(goodbye_message, 'green')}")
    print("")


def report_message(message, color = "yellow", width = 80, ii = "  ", si = "  ", newline = True, trailing_newline = False):
    if newline:
        print("", flush=True)
    if color:
        wprint(color_text(message, color), width = width, ii = ii, si = si)
    else:
        wprint(message, width = width, ii = ii, si = si)
    if trailing_newline:
        print("", flush=True)


def report_missing_input_genomes_file(path, flag):
    report_message(f'You specified "{path}" as a source of input genomes to use (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit(None, copy_log = False)


def print_suggest_help():
    print("\n  See `GToTree -h` for more info.")

@capture_stdout_to_log(lambda: log_file_var.get())
def report_early_exit(run_data, message = None, color = "red", suggest_help = False, copy_log = True):
    if message:
        print("")
        wprint(color_text(message, color))
    if suggest_help:
        print_suggest_help()
    print("\nExiting for now :(\n")
    if copy_log:
        copy_log_function(run_data)
    sys.exit(1)


def report_run_already_complete(output_dir):
    print("")
    wprint(color_text(
        "This run has already finished! Its "
        f'outputs are all in "{output_dir}".', "green"))
    wprint(color_text(
        "If you meant to run it over again from scratch, use `-F` to overwrite "
        "those outputs, or point `-o` at a new directory.", "yellow"))
    print("")
    sys.exit(0)


def report_very_early_exit(message = None, color = "red", suggest_help = False, leading_newline = True):
    if leading_newline:
        print("")
    if message:
        print("")
        wprint(color_text(message, color))
    if suggest_help:
        print_suggest_help()
    print("\nExiting for now :(\n")
    sys.exit(1)


def report_missing_pfam_targets_file(path, flag):
    report_message(f'You specified "{path}" as a source of Pfam targets to search each genome for (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit(None, copy_log = False)


def report_missing_ko_targets_file(path, flag):
    report_message(f'You specified "{path}" as a source of KO targets to search each genome for (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit(None, copy_log = False)


def report_missing_mapping_file(path, flag):
    report_message(f'You specified "{path}" as a mapping file to use (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit(None, copy_log = False)


def report_problem_with_mapping_file(mapping_file_problems, path, flag = "-m"):
    report_message(f'Unfortunately, there were problems detected in the mapping file "{path}" (passed to `{flag}`):\n')
    for problem in mapping_file_problems:
        report_message(problem, width = 100, ii = "    ", si = "    ", newline=False)
    report_message("Please correct these issues and try again!")
    report_early_exit(None, copy_log = False)


def stdout_and_log(*args, log_file="gtotree-runlog.txt", sep=" ", end="\n\n", flush=False, log_only=False, restart_log=False):
    message = sep.join(str(arg) for arg in args) + end
    if not log_only:
        print(message, end="", flush=flush)
    if restart_log:
        with open(log_file, "w") as f:
            f.write(message)
    else:
        with open(log_file, "a") as f:
            f.write(message)


################################################################################
# input-genome sources (shared by the main driver and the search subcommands)
################################################################################

# bullets sit at column 6 in both surfaces, and the total sits under them at 29 --
# `gtt search-pfams` / `gtt search-kos` are a subset of a GToTree run, so the block
# describing their inputs should be indistinguishable from GToTree's
INPUT_SOURCE_INDENT = " " * 6
INPUT_SOURCE_SUB_INDENT = " " * 10
TOTAL_LINE_INDENT = " " * 29


def _wanted_ref_tax_detail_lines(selection, indent):
    """
    The sub-bullets under one `-w` line: how that selection was narrowed, and how much
    of it had already been asked for.

    Takes a recorded selection dict off run_data rather than a live RefGenomeSelection,
    so a resumed run (which doesn't re-resolve `-w`) describes its selection exactly
    the way the original run did.
    """
    lines = []

    rank = selection.get("rank")
    derep = selection.get("derep_rank")

    if rank and derep:
        lines.append(f'{indent}- input-rank {rank} was dereplicated to one genome '
                     f'per {derep}')
    elif rank:
        lines.append(f'{indent}- all genomes under the input-rank {rank} were kept '
                     f'(--derep-rank off)')

    already_had = selection.get("num_selected", 0) - selection.get("num_added", 0)
    if already_had > 0:
        # with a repeated `-w`, the overlap can be with a previously-resolved taxon as
        # well as with `-a`, so this stays deliberately vague about which
        lines.append(f"{indent}- {already_had} more were selected, but were already "
                     f"counted")

    return lines


def input_genome_source_lines(args, run_data, indent=INPUT_SOURCE_INDENT):
    sub_indent = indent + "    "
    lines = []

    # driven off the recorded selections rather than `args.wanted_ref_tax`, because
    # `-w` is repeatable and each taxon gets its own line with its own detail
    for selection in getattr(run_data, "wanted_ref_tax_selections", None) or []:
        lines.append(f"{indent}- NCBI accessions selected for --wanted-ref-tax "
                     f"'{selection.get('taxon')}' "
                     f"({selection.get('num_added', 0):,} genomes)")
        lines.extend(_wanted_ref_tax_detail_lines(selection, sub_indent))

    if args.ncbi_accessions:
        lines.append(f"{indent}- NCBI accessions listed in {args.ncbi_accessions} "
                     f"({len(run_data.get_user_provided_ncbi_accs()):,} genomes)")
    if args.fasta_files:
        lines.append(f"{indent}- Fasta files listed in {args.fasta_files} "
                     f"({len(run_data.get_input_fasta_ids()):,} genomes)")
    if args.amino_acid_files:
        lines.append(f"{indent}- Amino-acid files listed in {args.amino_acid_files} "
                     f"({len(run_data.get_input_amino_acid_ids()):,} genomes)")
    if args.genbank_files:
        lines.append(f"{indent}- Genbank files listed in {args.genbank_files} "
                     f"({len(run_data.get_input_genbank_ids()):,} genomes)")

    return lines


def total_input_genomes_line(run_data, indent=TOTAL_LINE_INDENT):
    return f"{indent}Total input genomes: {len(run_data.all_input_genomes):,}"


RUN_INFO_BANNER = " ----------------------------------- RUN INFO ----------------------------------- "


def report_run_info_banner():
    print(f"\n{RUN_INFO_BANNER}\n", flush=True)


def report_wanted_ref_tax_info(args, run_data):
    """
    The `-w` provenance block that opens RUN INFO: which taxonomy asset the reference
    genomes came from, and anything worth flagging about how they were selected
    """
    if not getattr(args, "wanted_ref_tax", None):
        return

    from gtotree.utils.misc.general import wanted_ref_tax_source_line

    source_line = wanted_ref_tax_source_line(args.source)
    if source_line:
        print(source_line, flush=True)

    for selection in run_data.wanted_ref_tax_selections or []:
        for warning in selection.get("warnings") or []:
            report_message(warning, "yellow", ii="    ", si="    ")

    print(flush=True)


@capture_stdout_to_log(lambda: log_file_var.get())
def display_initial_run_info(args, run_data):

    # time.sleep(1)

    # this is here instead of above to prevent circular import problems (in other words, i suck at this)
    from gtotree.utils.misc.preflight_checks import check_input_genomes_amount

    report_wanted_ref_tax_info(args, run_data)

    report_message("  Input-genome sources include:")

    for line in input_genome_source_lines(args, run_data):
        print(line, flush=True)

    report_message(total_input_genomes_line(run_data, " " * 27), "green")
    # time.sleep(1)

    report_message("  Single-copy gene HMM source to be used:")
    print(f"      - {args.hmm} ({len(run_data.get_all_SCG_targets())} targets)", flush=True)
    auto_selected = getattr(args, "hmm_auto_selected", None)
    if auto_selected is not None:
        if auto_selected:
            wprint(f"(auto-selected; {auto_selected})",
                   ii=" " * 10, si=" " * 10)
        else:
            print("          (auto-selected)", flush=True)
    # time.sleep(1)

    check_and_report_any_changed_default_behavior(args, run_data)
    # time.sleep(3)

    with contextlib.redirect_stdout(sys.__stdout__):
        check_input_genomes_amount(len(run_data.all_input_genomes), args)

    run_data.start_time = datetime.now()

    return run_data


def copy_log_function(run_data):
    run_log_path = run_data.output_dir + "/gtotree-runlog.txt"
    timestamp = run_data.start_time.strftime("%Y-%m-%d_%I-%M-%S")
    new_log_path = run_data.logs_dir + f"/gtotree-logs/gtotree-runlog-{timestamp}.txt"
    shutil.copy(run_log_path, new_log_path)


@capture_stdout_to_log(lambda: log_file_var.get())
def add_border(extra_line=True):
    if extra_line:
        print()
    print("\n -------------------------------------------------------------------------------- ")


def report_phase_header(title, color = "yellow"):
    """
    The phase banner used by `gtt gen-scg-hmms`, `gtt search-pfams`, and `gtt search-kos`
    """
    rule = "  " + "-" * len(title)
    print(color_text(f"\n\n{rule}\n  {title}\n{rule}\n", color))


REMOVED_GENOMES_FILENAME = "removed-genomes.tsv"
SCG_INFO_FILENAME = "SCG-info.tsv"


def removed_genomes_pointer(run_data, stage):
    """
    Point at the one removals file, and say which rows this report is about
    """
    return (f"        {run_data.run_files_dir_rel}/{REMOVED_GENOMES_FILENAME}\n"
            f"        (stage_removed = {stage})\n\n")


@capture_stdout_to_log(lambda: log_file_var.get())
def report_notice(message, color = "yellow"):
    print("")
    print(f"{color_text("  ***********************************", color)} NOTICE {color_text("***********************************  ", color)}")
    print(message)
    print(color_text("  ******************************************************************************  ", color))


@capture_stdout_to_log(lambda: log_file_var.get())
def report_update(message, color = "green"):
    print("\n")
    print(f"{color_text("  ***********************************", color)} UPDATE {color_text("***********************************  ", color)}")
    print(message)
    print(color_text("  ******************************************************************************  ", color))


@capture_stdout_to_log(lambda: log_file_var.get())
def report_section_info(message):
    """
    Print an informational block with no banner around it
    """
    print("")
    print(message)


### specific notices
def check_and_report_any_changed_default_behavior(args, run_data):

    lines = []

    # output-dir / resume / overwrite: mutually exclusive messaging in one place
    if args.resume and args.output_already_existed:
        lines.append(f"Attempting to resume a previous run with outputs in \"{args.output_dir}\"")
    elif args.force_overwrite and args.output_already_existed:
        lines.append(f"Due to the `-F` flag, this output directory is being overwritten:\n          \"{args.output_dir}\"")
    elif args.output_dir != "gtotree-output":
        lines.append(f"The output directory has been set to: \"{args.output_dir}\"")

    if args.derep_rank != "auto":
        lines.append(f"The dereplication rank has been set to: \"{args.derep_rank}\"")

    if args.mapping_file:
        lines.append(f"Labels of the specified input genomes will be modified based on: \"{args.mapping_file}\"")

    if args.nucleotide_mode:
        lines.append("Working towards nucleotide alignments, as the `-z` flag was provided\n"
                     "          (amino-acid seqs are still used for HMM-searching of target genes)")

    if args.no_tree:
        lines.append("Only generating alignment, and no tree, as the `-N` flag was provided")

    if args.add_gtdb_tax:
        lines.append("GTDB taxonomic info will be added to labels where possible")
        if args.add_ncbi_tax:
            lines.append("NCBI taxonomic info will be added where possible when GTDB is not")
    elif args.add_ncbi_tax:
        lines.append("NCBI taxonomic info will be added to labels where possible")

    if args.lineage != "domain,phylum,class,genus,species":
        lines.append(f"Lineage info added to labels will be: \"{args.lineage}\"")

    if args.tree_program != "FastTreeMP":
        lines.append(f"The treeing program used will be: \"{args.tree_program}\"")

    if args.best_hit_mode:
        lines.append("Running in \"best-hit\" mode")

    if args.seq_length_cutoff != 0.2:
        lines.append(f"Gene-length filtering cutoff threshold (`-c`) has been set to: {args.seq_length_cutoff}")

    if args.gene_representation_cutoff != 0.1:
        lines.append(f"Gene-representation cutoff threshold (`-r`) has been set to: {args.gene_representation_cutoff}")

    if args.genome_hits_cutoff != 0.5:
        lines.append(f"Genome minimum target-gene-copy threshold (`-G`) has been set to: {args.genome_hits_cutoff}")

    if args.num_jobs != 4:
        lines.append(f"The number of jobs to run during parallelizable steps has been set to: {args.num_jobs}")

    if args.num_muscle_threads != 5:
        lines.append(f"The number of threads used for `muscle` calls will be: {args.num_muscle_threads}")

    if args.no_super5:
        lines.append("The 'super5' muscle algorithm will not be used even with greater than 1,000 input genomes")

    if args.keep_gene_alignments:
        lines.append("Individual protein-alignment files will retained, due to the `-k` flag being provided")

    if args.keep_working_dir:
        lines.append("The working directory and its intermediate files will be kept, due to the `-d` flag being provided")

    if args.target_pfams_file:
        lines.append(f"Genomes will be searched for Pfams listed in: {args.target_pfams_file} ({run_data.total_pfam_targets} targets)")

    if args.target_kos_file:
        lines.append(f"Genomes will be searched for KOs listed in: {args.target_kos_file} ({run_data.total_ko_targets} targets)")

    if not lines:
        return

    report_message("  Other options set:")
    for line in lines:
        print(f"      - {line}")


def many_genomes_notice(total_input_genomes):
    return (
    f"""    We seem to be aiming to work with {total_input_genomes} genomes. This is quite a bit, and
    the time individual gene alignments can take can quickly become prohibitive
    with many genomes like this.

    By default, GToTree is going to use the 'super5' muscle algorithm to help speed
    up the alignments for this run. If you don't want this to happen, you should
    cancel this run with `ctrl + c` now, and then add the `-X` flag to the GToTree
    command and re-run it.

    More info can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider

    And while we're chatting, you may also want to consider using \"representative\" genomes
    if you're not already. More info on that can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider

    We will wait 30 seconds before continuing with our regularly scheduled program :)"""
)


def few_genomes_notice(total_input_genomes, args):
    return (
    f"""    We seem to be aiming to work with {total_input_genomes} genomes. This is just a note that
    filtering by gene-length using the median length of a gene-set becomes
    less reliable with fewer genomes. The length-filtering is controlled by
    the `-c` parameter. If a lot of sequences are being dropped, you may want
    to consider increasing it.

    More info can be found here:
      github.com/AstrobioMike/GToTree/wiki/Things-to-consider

                   Moving forward with `-c` set to {args.seq_length_cutoff} this run."""
    )


def absurd_number_of_genomes_notice(total_input_genomes):
    return (
    f"""    The alignment and treeing steps, particularly the alignments, can become
    prohibitively computationally heavy with many input genomes. With {total_input_genomes:,} genomes,
    this job may not be feasible :(

    Often it is useful to slim down how many genomes of closely related organisms
    we are including when looking across a broad level of diversity, as having many
    closely related organisms may not add much to the final tree.

    Have you considered using a finer `--derep-rank`? It can help with reducing some
    redundancy when working at a broad level with many genomes.

    More info on that can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider

    You can consider cancelling this run now with `ctrl + c`, otherwise we will
    continue with our regularly scheduled program 60 seconds after this message
    was displayed :)"""
    )


@capture_stdout_to_log(lambda: log_file_var.get())
def report_processing_stage(stage, run_data):

    stages_dict = {
        "ncbi":                        "PROCESSING THE GENOMES PROVIDED AS NCBI ACCESSIONS",
        "genbank":                     "PROCESSING THE GENOMES PROVIDED AS GENBANK FILES",
        "fasta":                       "PROCESSING THE GENOMES PROVIDED AS FASTA FILES",
        "amino-acid":                  "PROCESSING THE GENOMES PROVIDED AS AMINO-ACID FILES",
        "processing-update":           "OVERALL SUMMARY OF INPUT-GENOME PROCESSING",
        "hmm-search":                  "SUMMARY OF THE TARGET SINGLE-COPY-GENE SEARCH",
        "filter-genes":                "FILTERING GENES BY LENGTH AND REPRESENTATION",
        "filter-genomes":              "FILTERING GENOMES WITH TOO FEW HITS",
        "align-and-prepare-SCG-sets":  "ALIGNING, TRIMMING, AND PREPARING SCG-SETS",
        "concatenate-SCG-sets":        "CONCATENATING ALL SCG-SET ALIGNMENTS TOGETHER",
        "updating-headers":            "ADDING MORE INFORMATIVE HEADERS",
        "treeing":                     "MAKING THE PHYLOGENOMIC TREE",
        "done":                        "DONE!!",
    }

    try:
        desc = stages_dict[stage]
    except KeyError as e:
        allowed = ", ".join(stages_dict)
        raise ValueError(f"Invalid stage: {stage!r}. Must be one of: {allowed!r}") from e

    add_border()

    if stage == "done":
        width = 82
        border = "#" * width
        bumper = "#" * 8
        inner_width = width - 2 * len(bumper)
        print()
        print(f"\n{border}")
        print(f"{border}")
        print(f"{bumper}{color_text(f"{desc.center(inner_width)}", "green")}{bumper}")
        print(f"{border}")
        print(f"{border}")

    else:
        width = 78
        border = "#" * width
        inner_width = width - 2 * len("####")
        print()
        print(f"\n  {border}")
        print(f"  ####{desc.center(inner_width)}####")
        print(f"  {border}")
        report_time_status(run_data.start_time)
    # time.sleep(1)


def report_ncbi_accs_not_found(num_accs, path):
    report_notice(f"    {num_accs} accession(s) not successfully found at NCBI.\n\n"
                f"    Reported in {path}/{REMOVED_GENOMES_FILENAME}")
    time.sleep(1)


def report_ncbi_update(run_data):
    num_input = len(run_data.ncbi_accs)
    num_not_found_at_ncbi = len(run_data.get_ncbi_accs_not_found())
    num_not_downloaded = len(run_data.get_ncbi_accs_not_downloaded())
    num_removed = len(run_data.get_removed_ncbi_accs())
    num_prepared = num_input - num_removed

    if num_removed == 0:
        message = (f"{color_text(f"All {num_input} input accessions were successfully downloaded and prepared!".center(82), "green")}")
    else:
        message = f"    Of the input genomes provided as {color_text("NCBI accessions", "yellow")}:\n\n"
        if num_not_found_at_ncbi > 0:
            message += (f"      {color_text(f"{num_not_found_at_ncbi} not found at NCBI", "yellow")}, reported in:\n"
                        + removed_genomes_pointer(run_data, GenomeRemovalStage.NCBI_LOOKUP))
        if num_not_downloaded > 0:
            message += (f"      {color_text(f"{num_not_downloaded} found but not successfully downloaded", "yellow")}, reported in:\n"
                        + removed_genomes_pointer(run_data, GenomeRemovalStage.NCBI_DOWNLOAD))

        message += (f"    {color_text(f"Overall, {num_prepared} of the input {num_input} accessions were successfully downloaded and\n    prepared.", "yellow")}")

    report_update(message)


def report_genbank_update(run_data):
    num_input = len(run_data.genbank_files)
    num_failed = len(run_data.get_failed_genbank_ids())
    num_prodigal_used = len(run_data.get_prodigal_used_genbank_ids())

    prodigal_note = (f"{color_text(f"{num_prodigal_used} had no CDS entries", "yellow")}, "
                     "so prodigal was used on the nucleotide sequences.")

    if num_failed == 0:
        message = color_text(f"All {num_input} input genbank files were successfully parsed and prepared!".center(82), "green")
        if num_prodigal_used > 0:
            message += f"\n\n    {prodigal_note}"
    else:
        message = f"    Of the input genomes provided as {color_text("genbank files", "yellow")}:\n\n"
        if num_prodigal_used > 0:
            message += f"      {prodigal_note}\n\n"
        message += (f"      {color_text(f"{num_failed} failed to be successfully parsed", "yellow")}, reported in:\n"
                    + removed_genomes_pointer(run_data, GenomeRemovalStage.GENBANK_PREP))
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} genbank files were successfully parsed and\n    prepared.", "yellow")}")

    report_update(message)


def report_fasta_update(run_data):
    num_input = len(run_data.fasta_files)
    num_failed = len(run_data.get_failed_fasta_ids())

    if num_failed == 0:
        message = (f"{color_text(f"All {num_input} input fasta files were successfully prepared!".center(82), "green")}")
    else:
        message = f"    Of the input genomes provided as {color_text("fasta files", "yellow")}:\n\n"
        message += (f"      {color_text(f"{num_failed} failed to be successfully processed", "yellow")}, reported in:\n"
                    + removed_genomes_pointer(run_data, GenomeRemovalStage.FASTA_PREP))
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} fasta files were successfully processed.", "yellow")}")

    report_update(message)

def report_AA_update(run_data):
    num_input = len(run_data.amino_acid_files)
    num_failed = len(run_data.get_failed_amino_acid_ids())

    if num_failed == 0:
        message = (f"{color_text(f"All {num_input} input amino-acid files were successfully prepared!".center(82), "green")}")
    else:
        message = f"    Of the input genomes provided as {color_text("amino-acid files", "yellow")}:\n\n"
        message += (f"      {color_text(f"{num_failed} failed to be successfully processed", "yellow")}, reported in:\n"
                    + removed_genomes_pointer(run_data, GenomeRemovalStage.AMINO_ACID_PREP))
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} amino-acid files were successfully processed.", "yellow")}")

    report_update(message)


def report_genome_processing_update(run_data, searched=False):
    """
    The one place the input-genome funnel is reported
    """
    num_input = len(run_data.all_input_genomes)
    num_failed_processing = len(run_data.get_genomes_removed_during_processing())
    num_failed_search = len(run_data.genomes_removed_at(GenomeRemovalStage.HMM_SEARCH)) \
        if searched else 0

    # counted as of the last stage folded in here, so removals further downstream
    # (e.g. the SCG-hit filter, already recorded on a resume) don't land in this count
    last_stage = (GenomeRemovalStage.HMM_SEARCH if searched
                  else GenomeRemovalStage.AMINO_ACID_PREP)
    num_remaining = len(run_data.genomes_alive_through(last_stage))

    verb = "processed and searched" if searched else "processed"

    if num_input == num_remaining:
        message = f"{color_text(f"All {num_input} input genomes were successfully {verb}.".center(82), "green")}"
    else:
        message = "    Of the input genomes provided:\n"
        if num_failed_processing > 0:
            message += f"      {color_text(f"{num_failed_processing} failed processing", "yellow")} as described above.\n\n"
        if num_failed_search > 0:
            message += (f"      {color_text(f"{num_failed_search} failed the target-gene search", "yellow")}, reported in:\n"
                        + removed_genomes_pointer(run_data, GenomeRemovalStage.HMM_SEARCH))
        message += f"    {color_text(f"Overall, {num_remaining} of the input {num_input} genomes were successfully {verb}.", "yellow")}"

        if num_remaining >= 4:
            message += "\n\n    Moving forward with those :)"

    report_section_info(message)

    if num_remaining < 4:
        report_too_few_genomes(run_data)


def report_pfam_searching_update(run_data):

    num_pfam_targets = run_data.total_pfam_targets
    num_pfams_found = len(run_data.found_pfam_targets)
    num_pfams_failed = len(run_data.failed_pfam_targets)

    if num_pfams_found == num_pfam_targets:
        message = f"    {color_text(f"Genomes were searched for all {num_pfam_targets} input Pfam targets.", 'green')}"
    elif num_pfams_found == 0:
        message = f"    {color_text("None of the input Pfam targets were found in the Pfam database", 'yellow')}, reported in:\n"
        message += f"      {run_data.run_files_dir_rel}/failed-pfam-targets.txt\n\n"
        message += f"    {color_text("So the input genomes were not searched for any Pfams :(", 'yellow')}"
    else:
        message = f"    {color_text(f"{num_pfams_failed} target Pfam(s) failed to be found in the Pfam database", "yellow")}, reported in:\n"
        message += f"      {run_data.run_files_dir_rel}/failed-pfam-targets.txt\n\n"
        message += f"    {color_text(f"Genomes were searched for the remaining {num_pfams_found} specified Pfams.", 'yellow')}"

    report_section_info(message)


def report_ko_searching_update(run_data):

    num_ko_targets = run_data.total_ko_targets
    num_kos_found = len(run_data.found_ko_targets)
    num_kos_failed = len(run_data.failed_ko_targets)

    if num_kos_found == num_ko_targets:
        message = f"    {color_text(f"Genomes were searched for all {num_ko_targets} input KO targets.", 'green')}"
    elif num_kos_found == 0:
        message = f"    {color_text("None of the input KO targets were found in the KO database", 'yellow')}, reported in:\n"
        message += f"      {run_data.run_files_dir_rel}/failed-ko-targets.txt\n\n"
        message += f"    {color_text("So the input genomes were not searched for any KOs :(", 'yellow')}"
    else:
        message = f"    {color_text(f"{num_kos_failed} target KO(s) failed to be found in the KO database", "yellow")}, reported in:\n"
        message += f"      {run_data.run_files_dir_rel}/failed-ko-targets.txt\n\n"
        message += f"    {color_text(f"Genomes were searched for the remaining {num_kos_found} specified KOs.", 'yellow')}"

    report_section_info(message)


def report_too_few_genomes(run_data):
    message = f"\n    {color_text("Unfortunately, there aren't enough genomes remaining to proceed...", 'red')}"
    print(message)
    report_early_exit(run_data)


def report_SCG_alignment_update(run_data):
    """
    Report SCG-sets lost during alignment
    """
    num_failed = len(run_data.SCG_targets_removed_at(SCGRemovalStage.ALIGNMENT))
    num_remaining = len(run_data.get_all_SCG_targets_remaining())

    if num_failed == 0:
        message = (f"{color_text(f"All {num_remaining} SCG-sets were successfully aligned and trimmed!".center(82), "green")}")
        report_update(message)
        return

    plural = "" if num_failed == 1 else "s"
    message = f"    {color_text(f"{num_failed} SCG-set{plural} failed to align or trim", "yellow")}, reported in:\n"
    message += f"      {run_data.run_files_dir_rel}/SCG-info.tsv\n\n"
    message += ("    This means muscle or trimal returned an error on those sets, not that\n"
                "    they were filtered out. Their logs were kept in:\n")
    message += f"      {run_data.logs_dir_rel}/failed-SCG-alignments/\n\n"
    message += (f"    {color_text(f"Moving forward with the remaining {num_remaining} target gene(s).", "yellow")}")

    report_update(message, "yellow")

    if num_remaining == 0:
        report_no_SCGs_remaining(run_data)


def report_genome_filtering_update(run_data):
    num_removed_due_to_hit_cutoff = len(
        run_data.genomes_removed_at(GenomeRemovalStage.SCG_HIT_FILTER))
    num_remaining = len(
        run_data.genomes_alive_through(GenomeRemovalStage.SCG_HIT_FILTER))
    num_input = len(run_data.all_input_genomes)

    if num_removed_due_to_hit_cutoff == 0:
        message = (f"{color_text("No genomes were removed due to having too few SCG hits!".center(82), "green")}")
    else:
        message = "    Of the input genomes remaining:\n\n"
        if not run_data.best_hit_mode:
            message += (f"      {color_text(f"{num_removed_due_to_hit_cutoff} genome(s) removed due to having too few unique SCG hits", "yellow")}, reported in:\n")
        else:
            message += (f"      {color_text(f"{num_removed_due_to_hit_cutoff} genome(s) removed due to having too few SCG hits", "yellow")}, reported in:\n")
        message += removed_genomes_pointer(run_data, GenomeRemovalStage.SCG_HIT_FILTER)

        if not run_data.best_hit_mode:
            message += ("    If this is a problem for the genomes you're working with, you could\n")
            message += ("    consider running GToTree in \"best-hit\" mode or adjusting the `-G`\n")
            message += ("    parameter. See the help menu for more info.\n\n")
        else:
            message += ("    If this is a problem for the genomes you're working with, you could\n")
            message += ("    consider adjusting the `-G` parameter. See the help menu for more info.\n\n")
        message += (f"    {color_text(f"Overall, {num_remaining} of the input {num_input} genomes are moving onto the treeing stage.", "yellow")}")

        if num_remaining >= 4:
            message += "\n\n"
            message += "Moving onto alignments with those :)".center(82)

    report_update(message)

    if num_remaining < 4:
        report_too_few_genomes(run_data)


def report_SCG_genome_filtering_update(run_data):
    """
    Report what genome filtering (`-G`) did to the SCG-sets.

    Two separate things, both previously invisible: sets that lost every one of their
    genomes and left the run entirely, and sets that survived but whose representation
    dropped below the `-r` the user asked for. The second isn't filtered on -- `-r` is
    evaluated once, before `-G` runs, and re-applying it afterwards would change the
    gene set out from under the user -- so it's surfaced here instead and left as their
    call. This is also where a total wipeout has to exit, since by the time the align
    stage notices, nothing has said why.
    """
    from gtotree.utils.misc.summary_info import (SCG_sets_below_representation_cutoff,
                                                 scg_info_denominators)

    num_removed = len(run_data.SCG_targets_removed_at(SCGRemovalStage.GENOME_FILTER))
    num_remaining = len(run_data.get_all_SCG_targets_remaining())
    below_cutoff = SCG_sets_below_representation_cutoff(run_data)

    if num_removed == 0 and not below_cutoff:
        return

    lines = []

    if num_removed:
        plural = "" if num_removed == 1 else "s"
        lines.append(
            f"      {color_text(f"{num_removed} SCG-set{plural} lost every genome with a hit", "yellow")}"
            " and left the run.")

    if below_cutoff:
        cutoff = f"{run_data.gene_representation_cutoff * 100:.0f}"
        _searched, retained = scg_info_denominators(run_data)
        num_below = len(below_cutoff)
        # the noun agrees with the total it's being drawn from, the verb with the subset
        verb = "is" if num_below == 1 else "are"
        plural = "" if num_remaining == 1 else "s"
        lines.append(
            f"      {color_text(f"{num_below} of the {num_remaining} remaining SCG-set{plural}", "yellow")} "
            f"{verb} now represented in fewer than\n      {cutoff}% of the {retained} retained "
            f"genomes (the `-r` cutoff, which is applied before\n      `-G` and not re-checked "
            "after it). They're kept, but contribute mostly gaps.")

    message = "    Of the SCG-sets that made it this far:\n\n" + "\n\n".join(lines)
    message += f"\n\n    Per-SCG details are in:\n      {run_data.run_files_dir_rel}/SCG-info.tsv"

    report_section_info(message)

    if num_remaining == 0:
        report_no_SCGs_remaining(run_data)


def report_SCG_set_filtering_update(run_data):
    total_SCG_targets = len(run_data.get_all_SCG_targets())
    num_SCG_targets_dropped = len(
        run_data.SCG_targets_removed_at(*SCG_GENE_FILTERING_STAGES))
    num_SCG_targets_remaining = total_SCG_targets - num_SCG_targets_dropped

    if num_SCG_targets_dropped == 0:
        message = f"{color_text(f"All {total_SCG_targets} SCG-targets were successfully aligned and prepared!".center(82), 'green')}"
    else:
        message = (f"    Of the initial {total_SCG_targets} SCG-targets:\n\n")
        message += (f"        {color_text(f"{num_SCG_targets_dropped} had no hits or were filtered out", 'yellow')}, reported in:\n")
        message += (f"          {run_data.run_files_dir_rel}/SCG-info.tsv")

        if num_SCG_targets_remaining != 0:
            message += "\n\n"
            message += (f"Moving forward with the remaining {num_SCG_targets_remaining} :)".center(82))

    report_update(message)

    if num_SCG_targets_remaining == 0:
        report_no_SCGs_remaining(run_data)


def report_no_SCGs_remaining(run_data):
    message = f"\n    {color_text("Unfortunately, there are no remaining SCG-targets to proceed with...", 'red')}"
    print(message)
    report_early_exit(run_data)


@capture_stdout_to_log(lambda: log_file_var.get())
def summarize_results(args, run_data):

    report_processing_stage("done", run_data)

    num_initial_genomes = len(run_data.get_all_input_genome_ids())
    num_remaining_genomes = len(run_data.get_all_remaining_input_genome_ids())

    print(f"\n  Overall, {num_remaining_genomes:,} of the initial {num_initial_genomes:,} genomes were retained (see notes below).\n")

    num_genes = len(run_data.get_all_SCG_targets_remaining())
    num_sites = run_data.final_alignment_length
    print(f"  The final alignment utilized {num_genes:,} target genes and contains {num_sites:,} total sites.\n")

    if not args.no_tree:
        final_tree_path = get_path_rel_to_outdir(run_data.final_tree_path, args)
        print(f"    Tree written to:\n        {color_text(final_tree_path, "green")}\n")

    final_alignment_path = get_path_rel_to_outdir(run_data.final_alignment_path, args)
    print(f"    Alignment written to:\n        {color_text(final_alignment_path, 'green')}\n")

    if run_data.header_update_error:
        wprint(color_text("Note: more informative headers failed swapping in for some reason, "
                          "so the alignment and tree above hold the original genome "
                          "IDs. See the log for why.", "yellow"),
               ii="    ", si="    ")
        print("")

    if args.keep_gene_alignments:
        gene_alignments_path = run_data.individual_gene_alignments_dir_rel
        print(f"    Individual target-gene alignments written to:\n        {color_text(gene_alignments_path, 'green')}\n")

    genome_summary_path = args.output_dir + "/genomes-summary-info.tsv"
    print(f"    Input-genomes summary table written to:\n        {color_text(genome_summary_path, 'green')}\n")

    SCG_hits_path = args.output_dir + "/SCG-hit-counts.tsv"
    print(f"    Summary table with hits per target-gene per genome written to:\n        {color_text(SCG_hits_path, 'green')}\n")

    if run_data.target_pfams_file:
        print(f"    Outputs from Pfam searching written to:\n        {color_text(run_data.pfam_results_dir_rel + "/", 'green')}\n")

    if run_data.target_kos_file:
        print(f"    Outputs from KO searching written to:\n        {color_text(run_data.ko_results_dir_rel + "/", 'green')}\n")

    partitions_files = f"{args.output_dir}/run-files/partitions.txt\n        {args.output_dir}/run-files/partitions.nex"
    print(f"    Partitions files (for downstream use with mixed-model treeing) written to:\n        {color_text(partitions_files, 'green')}")

    add_border(extra_line=False)

    num_genes_removed = len(run_data.get_all_SCG_targets()) - len(run_data.get_all_SCG_targets_remaining())

    if num_remaining_genomes < num_initial_genomes or num_genes_removed > 0:

        print("\n  Notes:\n")

        num_accs_not_found = len(run_data.get_ncbi_accs_not_found())
        if num_accs_not_found > 0:
            print(f"        {num_accs_not_found} accession(s) not successfully found at NCBI")

        num_accs_not_downloaded = len(run_data.get_ncbi_accs_not_downloaded())
        if num_accs_not_downloaded > 0:
            print(f"        {num_accs_not_downloaded} accession(s) did not download properly")

        num_failed_fasta_files = len(run_data.get_failed_fasta_ids())
        if num_failed_fasta_files > 0:
            print(f"        {num_failed_fasta_files} fasta file(s) failed to be processed")
        num_failed_genbank_files = len(run_data.get_failed_genbank_ids())
        if num_failed_genbank_files > 0:
            print(f"        {num_failed_genbank_files} genbank file(s) failed to be processed")
        num_failed_amino_acid_files = len(run_data.get_failed_amino_acid_ids())
        if num_failed_amino_acid_files > 0:
            print(f"        {num_failed_amino_acid_files} amino-acid file(s) failed to be processed")

        num_genomes_filtered_for_too_few_hits = len(run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering())
        if num_genomes_filtered_for_too_few_hits > 0:
            if not args.best_hit_mode:
                print(f"        {num_genomes_filtered_for_too_few_hits} genome(s) removed for having too few unique hits to the targeted SCGs")
            else:
                print(f"        {num_genomes_filtered_for_too_few_hits} genome(s) removed for having too few hits to the targeted SCGs")

        if num_genes_removed > 0:
            print(f"        {num_genes_removed} target gene(s) not found or filtered out of the analysis")

        print(f"\n    Reported along with additional informative files in:\n        {color_text(f"{run_data.run_files_dir_rel}/", 'green')}")

        add_border(extra_line=False)

    run_log_relative_path = run_data.output_dir_rel + "/gtotree-runlog.txt"
    print("\n    Log file written to:")
    print(f"        {color_text(run_log_relative_path, 'green')}")

    add_border(extra_line=False)

    citations_relative_path = args.output_dir + "/citations.txt"
    print("\n    Programs used and their citations have been written to:")
    print(f"        {color_text(citations_relative_path, 'green')}")

    add_border(extra_line=False)

    report_final_time_status(run_data.start_time)


def get_path_rel_to_outdir(path, args):
    """
    Re-express an absolute path as a relative path including the output dir
    """
    key_dir = args.output_dir
    rel = os.path.relpath(os.path.abspath(path), os.path.abspath(key_dir))
    if rel == os.pardir or rel.startswith(os.pardir + os.sep):
        raise ValueError(f"Directory {key_dir!r} not found in {path!r}")
    return os.path.join(key_dir, rel)
