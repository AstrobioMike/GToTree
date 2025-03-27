import textwrap
import sys
import time
import re
import contextlib
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
def check_and_report_any_changed_default_behavior(args):

    conditions = [
        args.output_dir != "gtotree-output",
        args.mapping_file,
        args.nucleotide_mode,
        args.no_tree,
        args.add_gtdb_tax,
        args.add_ncbi_tax,
        args.lineage != "Domain,Phylum,Class,Species",
        args.tree_program != "FastTreeMP",
        args.best_hit_mode,
        args.seq_length_cutoff != 0.2,
        args.genome_hits_cutoff != 0.5,
        args.num_jobs != 1,
        args.num_hmm_cpus != 2,
        args.muscle_threads != 5,
        args.no_super5,
        args.keep_gene_alignments,
        args.resume and args.output_already_existed,
        args.force_overwrite and args.output_already_existed,
        args.debug,
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

    if args.lineage != "Domain,Phylum,Class,Species":
        print(f"      - Lineage info added to labels will be: \"{args.lineage}\"")

    if args.tree_program != "FastTreeMP":
        print(f"      - The treeing program used will be: \"{args.tree_program}\"")

    if args.best_hit_mode:
        print("      - Running in \"best-hit\" mode")

    if args.seq_length_cutoff != 0.2:
        print(f"      - Gene-length filtering cutoff threshold (`-c`) has been set to: {args.seq_length_cutoff}")

    if args.genome_hits_cutoff != 0.5:
        print(f"      - Genome minimum gene-copy threshold (`-G`) has been set to: {args.genome_hits_cutoff}")

    if args.num_jobs != 1:
        print(f"      - The number of jobs to run during parallelizable steps has been set to: {args.num_jobs}")

    if args.num_hmm_cpus != 2:
        print(f"      - The number of CPUs used for `hmmsearch` calls will be: {args.num_hmm_cpus}")

    if args.muscle_threads != 5:
        print(f"      - The number of threads used for `muscle` calls will be: {args.muscle_threads}")

    if args.no_super5:
        print("      - The 'super5' muscle algorithm will not be used even with greater than 1,000 input genomes")

    if args.keep_gene_alignments:
        print("      - Individual protein-alignment files will retained, due to the `-k` flag being provided")

    if args.debug:
        print("      - Debug mode is enabled")

    if args.target_pfam_file:
        print(f"      - Genomes will be searched for Pfams listed in: {args.target_pfam_file} ({args.total_pfam_targets} targets)")

    if args.target_ko_file:
        print(f"      - Genomes will be searched for KOs listed in: {args.target_ko_file} ({args.total_ko_targets} targets)")


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

    And while we're chatting, you may also want to consider using \"prepresentative\" genomes
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
def report_processing_stage(stage):
    allowed_stages = ["ncbi", "genbank", "fasta", "amino-acid",
                      "preprocessing-update", "hmm-search", "filter-seqs",
                      "filter-genomes", "align", "tree"]

    if stage not in allowed_stages:
        raise ValueError(f"Invalid stage: {stage}. Must be one of: {', '.join(allowed_stages)}")

    if stage == "ncbi":
        message = ("\n  ##############################################################################\n"
                    "  ####        Preprocessing the genomes provided as NCBI accessions         ####\n"
                    "  ##############################################################################")
    elif stage == "genbank":
        message = ("\n  ##############################################################################\n"
                    "  ####         Preprocessing the genomes provided as genbank files          ####\n"
                    "  ##############################################################################")
    elif stage == "fasta":
        message = ("\n  ##############################################################################\n"
                    "  ####          Preprocessing the genomes provided as fasta files           ####\n"
                    "  ##############################################################################")
    elif stage == "amino-acid":
        message = ("\n  ##############################################################################\n"
                    "  ####        Preprocessing the genomes provided as amino-acid files        ####\n"
                    "  ##############################################################################")
    elif stage == "preprocessing-update":
        message = ("\n  ##############################################################################\n"
                    "  ####                Summary of input-genome preprocessing                 ####\n"
                    "  ##############################################################################")
    elif stage == "hmm-search":
        message = ("\n  ##############################################################################\n"
                    "  ####            Searching genomes for target single-copy genes            ####\n"
                    "  ##############################################################################")
    else:
        report_early_exit(f"Invalid stage ('{stage}'provided to `report_processing_stage`")

    add_border()
    print(message)
    time.sleep(1)


def report_snakemake_failure(description, snakemake_log):
    time.sleep(1)
    report_message(f"\nSnakemake failed while running the \"{description}\" workflow.\n", width = 90, color = "red")
    print(color_text(f"  Check the log at:\n    {snakemake_log}", "yellow"))
    report_early_exit()


def report_ncbi_accs_not_found(num_accs, path):
    report_notice(f"    {num_accs} accession(s) not successfully found at NCBI.\n"
                f"    Reported in {path}/ncbi-accessions-not-found.txt")
    time.sleep(1)


def report_ncbi_update(run_data):
    num_input = len(run_data.ncbi_accs)
    num_not_found_at_ncbi = len(run_data.get_ncbi_accs_not_found())
    num_not_downloaded = len(run_data.get_ncbi_accs_not_downloaded())
    num_prepared = len(run_data.get_done_ncbi_accs())
    num_removed = len(run_data.get_removed_ncbi_accs())

    if num_removed == 0:
        message = (f"    {color_text(f"All {num_input} input accessions were successfully downloaded and prepared!", "green")}")
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
        message = (f"    {color_text(f"All {num_input} input genbank files were successfully parsed and prepared!", "green")}")
    elif num_failed == 0:
        message = (f"    {color_text(f"All {num_input} input genbank files were successfully parsed and prepared!", "green")}\n\n")
    else:
        message = f"    Of the input genomes provided as {color_text("genbank files", "yellow")}:\n\n"
    if num_prodigal_used > 0 and num_failed == 0:
        message += (f"      {color_text(f"{num_prodigal_used} had no CDS entries", "yellow")}, so prodigal was used on the nucleotide sequences.")
    elif num_prodigal_used > 0:
        message += (f"      {color_text(f"{num_prodigal_used} had no CDS entries", "yellow")}, so prodigal was used on the nucleotide sequences.\n\n")
    if num_failed > 0:
        message += (f"      {color_text(f"{num_failed} failed to be successfully parsed", "yellow")}, reported in:\n"
                    f"        {run_data.run_files_dir_rel}/genbank-files-not-parsed.txt\n\n")
    if num_failed > 0:
        message += (f"    {color_text(f"Overall, {num_input - num_failed} of the input {num_input} genbank files were successfully parsed and\n    prepared.", "yellow")}")

    report_update(message)


def report_fasta_update(run_data):
    num_input = len(run_data.fasta_files)
    num_failed = len(run_data.get_failed_fasta_ids())

    if num_failed == 0:
        message = (f"    {color_text(f"All {num_input} input fasta files were successfully prepared!", "green")}")
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
        message = (f"    {color_text(f"All {num_input} input amino-acid files were successfully prepared!", "green")}")
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
        message = (f"    {color_text(f"All {num_input} input genomes were successfully preprocessed!\n\n", "green")}")
    else:
        message = f"    Of all the input genomes provided:\n\n"
        message += (f"      {color_text(f"{num_removed} failed preprocessing", "yellow")} as described above.\n\n")
        message += (f"    {color_text(f"Overall, {num_remaining} of the input {num_input} genomes were successfully preprocessed.\n\n", "yellow")}")

    message += f"                           Moving forward with those :)"
    report_update(message)
