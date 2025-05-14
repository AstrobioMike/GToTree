import textwrap
import sys
import time
import re
import contextlib
from datetime import datetime
from importlib.metadata import version
from gtotree.utils.context import log_file_var


tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


def gtotree_header():
    header = f"""

                                   GToTree v{version('GToTree')}
                         (github.com/AstrobioMike/GToTree)
    """
    return header


def get_version():
    return version('GToTree')


def color_text(text, color = 'green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


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


def capture_stdout_to_log(log_file):
    def decorator(func):
        def wrapper(*args, **kwargs):
            file_path = log_file() if callable(log_file) else log_file
            with open(file_path, 'a') as f:
                tee = Tee(sys.stdout, f)
                with contextlib.redirect_stdout(tee):
                    return func(*args, **kwargs)
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


def report_message(message, color = "yellow", width = 80, ii = "  ", si = "  ", newline = True):
    if newline:
        print("", flush=True)
    if color:
        wprint(color_text(message, color), width = width, ii = ii, si = si)
    else:
        wprint(message, width = width, ii = ii, si = si)


def report_missing_input_genomes_file(path, flag):
    report_message(f'You specified "{path}" as a source of input genomes to use (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_missing_pfam_targets_file(path, flag):
    report_message(f'You specified "{path}" as a source of Pfam targets to search each genome for (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_missing_ko_targets_file(path, flag):
    report_message(f'You specified "{path}" as a source of KO targets to search each genome for (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_missing_mapping_file(path, flag):
    report_message(f'You specified "{path}" as a mapping file to use (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_problem_with_mapping_file(mapping_file_problems, path, flag = "-m"):
    report_message(f'Unfortunately, there were problems detected in the mapping file "{path}" (passed to `{flag}`):\n')
    for problem in mapping_file_problems:
        report_message(problem, width = 100, ii = "    ", si = "    ", newline=False)
    report_message("Please correct these issues and try again!")
    report_early_exit()


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


@capture_stdout_to_log(lambda: log_file_var.get())
def display_initial_run_info(args, run_data):

    # this is here instead of above to prevent circular import problems (in other words, i suck at this)
    from gtotree.utils.preflight_checks import check_input_genomes_amount

    # time.sleep(1)

    print("\n ----------------------------------- RUN INFO ----------------------------------- \n", flush=True)

    # time.sleep(1)

    report_message("  Input-genome sources include:")

    if args.ncbi_accessions:
        print(f"      - NCBI accessions listed in {args.ncbi_accessions} ({len(run_data.get_input_ncbi_accs())} genomes)", flush=True)
    if args.genbank_files:
        print(f"      - Genbank files listed in {args.genbank_files} ({len(run_data.get_input_genbank_ids())} genomes)", flush=True)
    if args.fasta_files:
        print(f"      - Fasta files listed in {args.fasta_files} ({len(run_data.get_input_fasta_ids())} genomes)", flush=True)
    if args.amino_acid_files:
        print(f"      - Amino-acid files listed in {args.amino_acid_files} ({len(run_data.get_input_amino_acid_ids())} genomes)", flush=True)

    report_message(f"                           Total input genomes: {len(run_data.all_input_genomes)}", "green")
    # time.sleep(1)

    report_message("  Single-copy gene HMM source to be used:")
    print(f"      - {args.hmm} ({len(run_data.get_all_SCG_targets())} targets)", flush=True)
    # time.sleep(1)

    check_and_report_any_changed_default_behavior(args, run_data)
    # time.sleep(3)

    with contextlib.redirect_stdout(sys.__stdout__):
        check_input_genomes_amount(len(run_data.all_input_genomes), args)

    run_data.start_time = datetime.now()

    return run_data


@capture_stdout_to_log(lambda: log_file_var.get())
def report_early_exit(message = None, color = "red", suggest_help = False):
    if message:
        print("")
        wprint(color_text(message, color))
    if suggest_help:
        print("\n  See `GToTree -h` for more info.")
    print("\nExiting for now :(\n")
    sys.exit(1)

@capture_stdout_to_log(lambda: log_file_var.get())
def add_border():
    print("\n -------------------------------------------------------------------------------- ")

@capture_stdout_to_log(lambda: log_file_var.get())
def report_notice(message, color = "yellow"):
    print("")
    print(f"{color_text("  ***********************************", color)} NOTICE {color_text("***********************************  ", color)}")
    print(message)
    print(color_text("  ******************************************************************************  ", color))


@capture_stdout_to_log(lambda: log_file_var.get())
def report_update(message, color = "green"):
    print("")
    print(f"{color_text("  ***********************************", color)} UPDATE {color_text("***********************************  ", color)}")
    print(message)
    print(color_text("  ******************************************************************************  ", color))


### specific notices
def check_and_report_any_changed_default_behavior(args, run_data):

    conditions = [
        args.output_dir != "gtotree-output",
        args.mapping_file,
        args.nucleotide_mode,
        args.no_tree,
        args.add_gtdb_tax,
        args.add_ncbi_tax,
        args.lineage != "Domain,Phylum,Class,Genus,Species",
        args.tree_program != "FastTreeMP",
        args.best_hit_mode,
        args.seq_length_cutoff != 0.2,
        args.genome_hits_cutoff != 0.5,
        args.num_jobs != 1,
        args.num_hmm_cpus != 2,
        args.num_muscle_threads != 5,
        args.no_super5,
        args.keep_gene_alignments,
        args.resume and args.output_already_existed,
        args.force_overwrite and args.output_already_existed,
        args.debug,
        args.target_pfam_file,
        args.target_ko_file,
    ]

    if any(conditions):
        report_message("  Other options set:")

    if args.resume and args.output_already_existed:
        print(f"      - Attempting to resume a previous run with outputs in \"{args.output_dir}\"")

    if args.force_overwrite:
        if args.output_already_existed:
            print(f"      - The `-F` flag was provided, so this output directory is being overwritten: \"{args.output_dir}\"")

    if args.output_dir != "gtotree-output" and not args.resume:
        print(f"      - The output directory has been set to: \"{args.output_dir}\"")

    if args.mapping_file:
        print(f"      - Labels of the specified input genomes will be modified based on: \"{args.mapping_file}\"")

    if args.nucleotide_mode:
        print("      - Working towards nucleotie alignments, as the `-z` flag was provided\n"
              "          (amino-acid seqs are still used for HMM-searching of target genes)")

    if args.no_tree:
        print("      - Only generating alignment, and no tree, as the `-N` flag was provided")

    if args.add_gtdb_tax:
        print("      - GTDB taxonomic info will be added to labels where possible")
        if args.add_ncbi_tax:
            print("      - NCBI taxonomic info will be added where possible when GTDB is not")

    if args.add_ncbi_tax and not args.add_gtdb_tax:
        print("      - NCBI taxonomic info will be added to labels where possible")

    if args.lineage != "Domain,Phylum,Class,Genus,Species":
        print(f"      - Lineage info added to labels will be: \"{args.lineage}\"")

    if args.tree_program != "FastTreeMP":
        print(f"      - The treeing program used will be: \"{args.tree_program}\"")

    if args.best_hit_mode:
        print("      - Running in \"best-hit\" mode")

    if args.seq_length_cutoff != 0.2:
        print(f"      - Gene-length filtering cutoff threshold (`-c`) has been set to: {args.seq_length_cutoff}")

    if args.genome_hits_cutoff != 0.5:
        print(f"      - Genome minimum target-gene-copy threshold (`-G`) has been set to: {args.genome_hits_cutoff}")

    if args.num_jobs != 1:
        print(f"      - The number of jobs to run during parallelizable steps has been set to: {args.num_jobs}")

    if args.num_hmm_cpus != 2:
        print(f"      - The number of CPUs used for `hmmsearch` calls will be: {args.num_hmm_cpus}")

    if args.num_muscle_threads != 5:
        print(f"      - The number of threads used for `muscle` calls will be: {args.muscle_threads}")

    if args.no_super5:
        print("      - The 'super5' muscle algorithm will not be used even with greater than 1,000 input genomes")

    if args.keep_gene_alignments:
        print("      - Individual protein-alignment files will retained, due to the `-k` flag being provided")

    if args.debug:
        print("      - Debug mode is enabled")

    if args.target_pfam_file:
        print(f"      - Genomes will be searched for Pfams listed in: {args.target_pfam_file} ({run_data.total_pfam_targets} targets)")

    if args.target_ko_file:
        print(f"      - Genomes will be searched for KOs listed in: {args.target_ko_file} ({run_data.total_ko_targets} targets)")


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
      github.com/AstrobioMike/GToTree/wiki/things-to-consider#working-with-many-genomes

    And while we're chatting, you may also want to consider using \"representative\" genomes
    if you're not already. More info on that can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider#consider-using-representative-genomes

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
      github.com/AstrobioMike/GToTree/wiki/Things-to-consider#filtering-hits-by-gene-length

                   Moving forward with `-c` set to {args.seq_length_cutoff} this run."""
    )


def absurd_number_of_genomes_notice(total_input_genomes):
    return (
    f"""    The alignment and treeing steps, particularly the alignments, can become
    prohibitively memory-intensive with many input genomes. With {total_input_genomes} genomes,
    this job may not be feasible :(

    Often it is useful to slim down how many genomes of closely related organisms
    we are including when looking across a broad level of diversity, as having many
    closely related organisms may not add much to the final tree.

    Have you considered using "representative" genomes only (either from NCBI or
    GTDB)? Those both provide helpful systems for reducing some redundancy when
    working at a broad level with many genomes.

    More info on that can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider#consider-using-representative-genomes

    You can consider cancelling this run now with `ctrl + c`, otherwise we will
    continue with our regularly scheduled program 60 seconds after this message
    was displayed :)"""
    )


@capture_stdout_to_log(lambda: log_file_var.get())
def report_processing_stage(stage, run_data):

    stages_dict = {
        "ncbi":                        "PREPROCESSING THE GENOMES PROVIDED AS NCBI ACCESSIONS",
        "genbank":                     "PREPROCESSING THE GENOMES PROVIDED AS GENBANK FILES",
        "fasta":                       "PREPROCESSING THE GENOMES PROVIDED AS FASTA FILES",
        "amino-acid":                  "PREPROCESSING THE GENOMES PROVIDED AS AMINO-ACID FILES",
        "preprocessing-update":        "OVERALL SUMMARY OF INPUT-GENOME PREPROCESSING",
        "hmm-search":                  "SEARCHING GENOMES FOR TARGET SINGLE-COPY GENES",
        "filter-genes":                "FILTERING GENES BY LENGTH",
        "filter-genomes":              "FILTERING GENOMES WITH TOO FEW HITS",
        "align-and-prepare-SCG-sets":  "ALIGNING, TRIMMING, AND PREPARING SCG-SETS",
        "concatenate-SCG-sets":        "CONCATENATING ALL SCG-SET ALIGNMENTS TOGETHER",
        "updating-headers":            "ADDING MORE INFORMATIVE HEADERS",
        "treeing":                     "MAKING THE PHYLOGENOMIC TREE",
        "done":                        "DONE!!",
    }

    try:
        desc = stages_dict[stage]
    except KeyError:
        allowed = ", ".join(stages_dict)
        raise ValueError(f"Invalid stage: {stage!r}. Must be one of: {allowed!r}")

    add_border()

    if stage == "done":
        width = 82
        border = "#" * width
        bumper = "#" * 8
        inner_width = width - 2 * len(bumper)
        print()
        print(f"{border}")
        print(f"{border}")
        print(f"{bumper}{color_text(f"{desc.center(inner_width)}", "green")}{bumper}")
        print(f"{border}")
        print(f"{border}")

    else:
        width = 78
        border = "#" * width
        inner_width = width - 2 * len("####")
        print()
        print(f"  {border}")
        print(f"  ####{desc.center(inner_width)}####")
        print(f"  {border}")
        report_time_status(run_data.start_time)
    # time.sleep(1)


def report_snakemake_failure(description, log):
    time.sleep(1)
    report_message(f"\nSnakemake failed while running the \"{description}\" workflow.\n", width = 90, color = "red")
    print(color_text(f"  Check the log at:\n    {log}", "yellow"))
    report_early_exit()


def report_ncbi_accs_not_found(num_accs, path):
    report_notice(f"    {num_accs} accession(s) not successfully found at NCBI.\n\n"
                f"    Reported in {path}/ncbi-accessions-not-found.txt")
    time.sleep(1)


def report_ncbi_update(run_data):
    num_input = len(run_data.ncbi_accs)
    num_not_found_at_ncbi = len(run_data.get_ncbi_accs_not_found())
    num_not_downloaded = len(run_data.get_ncbi_accs_not_downloaded())
    num_prepared = len(run_data.get_done_ncbi_accs())
    num_removed = len(run_data.get_removed_ncbi_accs())

    if num_removed == 0:
        message = (f"{color_text(f"All {num_input} input accessions were successfully downloaded and prepared!".center(82), "green")}")
    else:
        message = f"    Of the input genomes provided as {color_text("NCBI accessions", "yellow")}:\n\n"
        if num_not_found_at_ncbi > 0:
            message += (f"      {color_text(f"{num_not_found_at_ncbi} not found at NCBI", "yellow")}, reported in:\n"
                        f"        {run_data.run_files_dir_rel}/ncbi-accessions-not-found.txt\n\n")
        if num_not_downloaded > 0:
            message += (f"      {color_text(f"{num_not_downloaded} found but not successfully downloaded", "yellow")}, reported in:\n"
                        f"        {run_data.run_files_dir_rel}/ncbi-accessions-not-downloaded.txt\n\n")
        if num_removed > 0:
            message += (f"    {color_text(f"Overall, {num_prepared} of the input {num_input} accessions were successfully downloaded and\n    prepared.", "yellow")}")

    report_update(message)


def report_genbank_update(run_data):
    num_input = len(run_data.genbank_files)
    num_failed = len(run_data.get_failed_genbank_ids())
    num_prodigal_used = len(run_data.get_prodigal_used_genbank_ids())

    if num_failed == 0 and num_prodigal_used == 0:
        message = (f"{color_text(f"All {num_input} input genbank files were successfully parsed and prepared!".center(82), "green")}")
    elif num_failed == 0:
        message = (f"{color_text(f"All {num_input} input genbank files were successfully parsed and prepared!".center(82), "green")}\n\n")
    else:
        message = f"    Of the input genomes provided as {color_text("genbank files", "yellow")}:\n\n"
    if num_prodigal_used > 0 and num_failed == 0:
        message += (f"      {color_text(f"{num_prodigal_used} had no CDS entries", "yellow")}, so prodigal was used on the nucleotide sequences.")
    elif num_prodigal_used > 0:
        message += (f"      {color_text(f"{num_prodigal_used} had no CDS entries", "yellow")}, so prodigal was used on the nucleotide sequences.\n\n")
    if num_failed > 0:
        message += (f"      {color_text(f"{num_failed} failed to be successfully parsed", "yellow")}, reported in:\n"
                    f"        {run_data.run_files_dir_rel}/failed-genbank-files.txt\n\n")
    if num_failed > 0:
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} genbank files were successfully parsed and\n    prepared.", "yellow")}")

    report_update(message)


def report_fasta_update(run_data):
    num_input = len(run_data.fasta_files)
    num_failed = len(run_data.get_failed_fasta_ids())

    if num_failed == 0:
        message = (f"{color_text(f"All {num_input} input fasta files were successfully prepared!".center(82), "green")}")
    else:
        message = f"    Of the input genomes provided as {color_text("fasta files", "yellow")}:\n\n"
        message += (f"      {color_text(f"{num_failed} failed to be successfully preprocessed", "yellow")}, reported in:\n"
                    f"        {run_data.run_files_dir_rel}/failed-fasta-files.txt\n\n")
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} fasta files were successfully preprocessed.", "yellow")}")

    report_update(message)

def report_AA_update(run_data):
    num_input = len(run_data.amino_acid_files)
    num_failed = len(run_data.get_failed_amino_acid_ids())

    if num_failed == 0:
        message = (f"{color_text(f"All {num_input} input amino-acid files were successfully prepared!".center(82), "green")}")
    else:
        message = f"    Of the input genomes provided as {color_text("amino-acid files", "yellow")}:\n\n"
        message += (f"      {color_text(f"{num_failed} failed to be successfully preprocessed", "yellow")}, reported in:\n"
                    f"        {run_data.run_files_dir_rel}/failed-amino-acid-files.txt\n\n")
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} amino-acid files were successfully preprocessed.", "yellow")}")

    report_update(message)


def report_genome_preprocessing_update(run_data):
    num_input = len(run_data.all_input_genomes)
    num_removed = len(run_data.get_all_removed_input_genomes())
    num_remaining = num_input - num_removed

    if num_input == num_remaining:
        message = (f"{color_text(f"All {num_input} input genomes were successfully preprocessed!".center(82), "green")}")
    else:
        message = f"    Of all the input genomes provided:\n\n"
        message += (f"      {color_text(f"{num_removed} failed preprocessing", "yellow")} as described above.\n\n")
        message += (f"    {color_text(f"Overall, {num_remaining} of the input {num_input} genomes were successfully preprocessed.", "yellow")}")

        if num_remaining >= 4:
            message += "\n\n"
            message += "Moving forward with those :)".center(82)

    report_update(message)

    if num_remaining < 4:
        report_too_few_genomes()


def report_too_few_genomes():
    message = f"\n    {color_text("Unfortunately, there aren't enough genomes remaining to proceed...", 'red')}"
    print(message)
    report_early_exit()

def report_hmm_search_update(run_data):
    num_searched = len(run_data.get_all_preprocessed_genomes())
    num_failed = len(run_data.get_failed_hmm_search_paths())
    num_successful = num_searched - num_failed

    if num_failed == 0:
        message = (f"{color_text(f"All {num_searched} genomes were successfully searched for target genes!".center(82), "green")}")
    else:
        message = f"    Of the {num_searched} input genomes, {num_successful} were successfully searched for the target genes.\n\n"

    report_update(message)

    if num_successful < 4:
        report_too_few_genomes()


def report_genome_filtering_update(run_data):
    num_removed_due_to_hit_cutoff = len(run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering())
    num_remaining = len(run_data.get_all_remaining_input_genomes())
    num_input = len(run_data.all_input_genomes)

    if num_removed_due_to_hit_cutoff == 0:
        message = (f"{color_text(f"No genomes were removed due to having too few SCG hits!".center(82), "green")}")
    else:
        message = f"    Of the input genomes remaining:\n\n"
        if not run_data.best_hit_mode:
            message += (f"      {color_text(f"{num_removed_due_to_hit_cutoff} genome(s) removed due to having too few unique SCG hits", "yellow")}, reported in:\n")
        else:
            message += (f"      {color_text(f"{num_removed_due_to_hit_cutoff} genome(s) removed due to having too few SCG hits", "yellow")}, reported in:\n")
        message += (f"        {run_data.run_files_dir_rel}/genomes-removed-for-too-few-SCG-hits.txt\n\n")

        if not run_data.best_hit_mode:
            message += (f"    If this is a problem for the genomes you're working with, you could\n")
            message += (f"    consider running GToTree in \"best-hit\" mode or adjusting the `-G`\n")
            message += (f"    parameter. See the help menu for more info.\n\n")
        else:
            message += (f"    If this is a problem for the genomes you're working with, you could\n")
            message += (f"    consider adjusting the `-G` parameter. See the help menu for more info.\n\n")
        message += (f"    {color_text(f"Overall, {num_remaining} of the input {num_input} made it through the preprocessing gauntlet.", "yellow")}")

        if num_remaining >= 4:
            message += "\n\n"
            message += "Moving onto alignments with those :)".center(82)

    report_update(message)

    if num_remaining < 4:
        report_too_few_genomes()


def report_SCG_set_alignment_update(run_data):
    total_SCG_targets = len(run_data.get_all_SCG_targets())
    num_SCG_targets_remaining = len(run_data.get_all_SCG_targets_remaining())
    num_SCG_targets_dropped = total_SCG_targets - num_SCG_targets_remaining

    if num_SCG_targets_dropped == 0:
        message = f"{color_text(f"All {total_SCG_targets} SCG-targets were successfully aligned and prepared!".center(82), 'green')}"
    else:
        message = (f"    Of the initial {total_SCG_targets} SCG-targets:\n\n")
        message += (f"        {color_text(f"{num_SCG_targets_dropped} had no hits or were filtered out", 'yellow')}, reported in:\n")
        message += (f"          {run_data.run_files_dir_rel}/target-SCGs-filtered-out-or-not-found.txt")

        if num_SCG_targets_remaining != 0:
            message += "\n\n"
            message += (f"Moving forward with the remaining {num_SCG_targets_remaining} :)".center(82))

    report_update(message)

    if num_SCG_targets_remaining == 0:
        report_no_SCGs_remaining()

def report_no_SCGs_remaining():
    message = f"\n    {color_text("Unfortunately, there are no remaining SCG-targets to proceed with...", 'red')}"
    print(message)
    report_early_exit()

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

    if args.keep_gene_alignments:
        gene_alignments_path = run_data.individual_gene_alignments_dir_rel
        print(f"    Individual target-gene alignments written to:\n        {color_text(gene_alignments_path, 'green')}\n")

    genome_summary_path = args.output_dir + "/genomes-summary-info.tsv"
    print(f"    Input-genomes summary table written to:\n        {color_text(genome_summary_path, 'green')}\n")

    SCG_hits_path = args.output_dir + "/SCG-hit-counts.tsv"
    print(f"    Summary table with hits per target-gene per genome written to:\n        {color_text(SCG_hits_path, 'green')}\n")

    partitions_file = args.output_dir + "/run-files/partitions.txt"
    print(f"    Partitions file (for downstream use with mixed-model treeing) written to:\n        {color_text(partitions_file, 'green')}")

    add_border()

    if num_remaining_genomes < num_initial_genomes:
        print(f"\n  Notes:\n")

        num_accs_not_found = len(run_data.get_ncbi_accs_not_found())
        if num_accs_not_found > 0:
            print(f"        {num_accs_not_found} accession(s) not successfully found at NCBI")

        num_accs_not_downloaded = len(run_data.get_ncbi_accs_not_downloaded())
        if num_accs_not_downloaded > 0:
            print(f"        {num_accs_not_downloaded} accession(s) did not download properly")

        num_failed_fasta_files = len(run_data.get_failed_fasta_ids())
        if num_failed_fasta_files > 0:
            print(f"        {num_failed_fasta_files} fasta file(s) failed to be preprocessed")
        num_failed_genbank_files = len(run_data.get_failed_genbank_ids())
        if num_failed_genbank_files > 0:
            print(f"        {num_failed_genbank_files} genbank file(s) failed to be preprocessed")
        num_failed_amino_acid_files = len(run_data.get_failed_amino_acid_ids())
        if num_failed_amino_acid_files > 0:
            print(f"        {num_failed_amino_acid_files} amino-acid file(s) failed to be preprocessed")

        num_genomes_filtered_for_too_few_hits = len(run_data.get_all_input_genomes_due_for_SCG_min_hit_filtering())
        if num_genomes_filtered_for_too_few_hits > 0:
            if not args.best_hit_mode:
                print(f"        {num_genomes_filtered_for_too_few_hits} genome(s) removed for having too few unique hits to the targeted SCGs")
            else:
                print(f"        {num_genomes_filtered_for_too_few_hits} genome(s) removed for having too few hits to the targeted SCGs")

        num_genes_removed = len(run_data.get_all_SCG_targets()) - len(run_data.get_all_SCG_targets_remaining())
        if num_genes_removed > 0:
            if not args.best_hit_mode:
                print(f"        {num_genes_removed} gene(s) excluded for having no hits or only multiple hits in each genome")
            else:
                print(f"        {num_genes_removed} gene(s) excluded for having no hits in the input genomes")

        print(f"\n    Reported along with additional informative files in:\n        {color_text(f"{run_data.run_files_dir_rel}/", 'green')}")

        add_border()

    run_log_relative_path = run_data.output_dir_rel + "/gtotree-runlog.txt"
    print(f"\n  Log file written to:\n      {color_text(run_log_relative_path, 'green')}")

    add_border()

    citations_relative_path = args.output_dir + "/citations.txt"
    print(f"\n  {color_text("Programs used and their citations have been written to:", 'yellow')}")
    print(f"      {color_text(citations_relative_path, 'green')}")

    add_border()

    report_final_time_status(run_data.start_time)


def get_path_rel_to_outdir(path, args):

    key_dir = args.output_dir
    idx = path.find(key_dir)
    if idx != -1:
        sub_path = path[idx:]
        return(sub_path)
    else:
        raise ValueError(f"Directory {key_dir!r} not found in {path!r}")
