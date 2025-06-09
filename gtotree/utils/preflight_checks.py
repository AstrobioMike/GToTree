import os
import sys
import shutil
import time
import pandas as pd
import tempfile
from collections import Counter
import json
import hashlib
from gtotree.utils.messaging import (color_text,
                                     report_message,
                                     report_early_exit,
                                     report_missing_input_genomes_file,
                                     report_missing_pfam_targets_file,
                                     report_missing_ko_targets_file,
                                     report_missing_mapping_file,
                                     report_problem_with_mapping_file,
                                     report_notice,
                                     many_genomes_notice,
                                     few_genomes_notice,
                                     absurd_number_of_genomes_notice,
                                     gtotree_header,
                                     stdout_and_log
                                     )
from gtotree.utils.hmms.scg_hmm_setup import check_hmm_file
from gtotree.utils.ncbi.get_ncbi_assembly_tables import get_ncbi_assembly_data
from gtotree.utils.ncbi.get_ncbi_tax_data import get_ncbi_tax_data
from gtotree.utils.gtdb.get_gtdb_data import get_gtdb_data
from gtotree.utils.ko.get_kofamscan_data import get_kofamscan_data
from gtotree.utils.general import (ToolsUsed,
                                   populate_run_data,
                                   read_run_data)
from gtotree.utils.context import log_file_var


def preflight_checks(args):
    check_for_essential_deps()
    args, run_data = primary_args_validation(args)
    check_for_required_dbs(args)
    run_data = track_tools_used(args, run_data)
    args, run_data = final_setups(args, run_data)

    return args, run_data


def args_hash_function(args):
    args_dict = vars(args).copy()
    args_dict.pop("force_overwrite", False)
    args_dict.pop("resume", False)
    args_dict.pop("run_files_dir_rel", None)
    args_dict.pop("run_files_dir", None)
    args_dict.pop("output_already_existed", None)
    args_json = json.dumps(args_dict, sort_keys=True)
    return hashlib.sha256(args_json.encode("utf-8")).hexdigest()

def check_for_essential_deps():
    commands = ["muscle", "hmmsearch", "trimal"]
    for cmd in commands:
        program_check(cmd, essential = True)


def program_check(cmd, essential = False):
    path = shutil.which(cmd)
    if not path:
        if essential:
            report_early_exit(None, message = f"{cmd} is an essential dependency, but it's not in your PATH :(", copy_log = False)
        else:
            report_early_exit(None, message = f"You specified to use {cmd}, but it's not in your PATH :(", copy_log = False)


def primary_args_validation(args):
    check_for_minimum_args(args)
    check_optional_deps(args)
    check_lineage(args)
    check_tree_program(args)
    checks_for_nucleotide_mode(args)
    args = check_output_dir(args)
    args, run_data = check_input_files(args)
    check_for_min_input_genomes(run_data)
    return args, run_data


def check_for_minimum_args(args):
    if not args.ncbi_accessions and not args.genbank_files and not args.fasta_files and not args.amino_acid_files:
        report_message("You need to provide at least one input-genome source!")
        report_early_exit(None, suggest_help=True, copy_log = False)
    if not args.hmm:
        report_message("You need to specify the HMM file of the target-SCGs you want to tree! "
                       "You can view the available gene-sets packaged with GToTree by running `gtt-hmms`.")
        report_early_exit(None, suggest_help=True, copy_log = False)


def check_optional_deps(args):
    if args.add_ncbi_tax:
        program_check("taxonkit")


def check_lineage(args):
    accepted_ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
    lineage_list = args.lineage.split(",")

    for rank in lineage_list:
        if rank.capitalize() not in accepted_ranks:
            report_message(f'You specified "{args.lineage}" to the `-L` argument, but "{rank}" is not an accepted taxonomic rank.')
            print(f"\n  Accepted ranks are any combination of the below entered as a comma-delimited list:\n        {'\n        '.join(accepted_ranks)}")
            report_early_exit(None, copy_log = False)

    if args.lineage != "Domain,Phylum,Class,Genus,Species" and not args.add_ncbi_tax and not args.add_gtdb_tax:
        report_message("You've specified a custom lineage (`-L`), but neither the "
                       "`-t` nor `-D` flags were provided to indicate which taxonomy to use.")
        report_early_exit(None, suggest_help=True, copy_log = False)

    if args.add_ncbi_tax and args.add_gtdb_tax:
        report_message("You've specified to use both the NCBI and GTDB taxonomies. "
                       "Please choose one or the other.")
        report_early_exit(None, suggest_help=True, copy_log = False)


def check_tree_program(args):
    accepted_programs = ["FastTree", "FastTreeMP", "VeryFastTree", "IQTREE"]
    if args.tree_program not in accepted_programs:
        report_message(f'You specified "{args.tree_program}" to the `-T` argument, but that\'s not an available treeing program.')
        print(f"\n  Available programs are:\n        {'\n        '.join(accepted_programs)}")
        report_early_exit(None, copy_log = False)
    program_check(args.tree_program)


def checks_for_nucleotide_mode(args):
    if args.nucleotide_mode:
        if args.amino_acid_files:
            report_message("You've specified wanting to work with nucleotide sequences (by passing the `-z` parameter), "
                           "but also provided some input genomes as amino-acid files (passed to `-A`). We can't confidently reverse-translate "
                           "amino-acid seqs to nucleotide seqs, so we can't take both of those options.")
            report_early_exit(None, copy_log = False)
        if args.genbank_files:
            report_message("You've specified wanting to work with nucleotide sequences (by passing the `-z` parameter), "
                           "but also provided some input genomes as genbank files (passed to `-g`). Input genbank files are currently "
                           "not supported with nucleotide mode.")
            report_early_exit(None, copy_log = False)


def check_input_files(args):
    if args.ncbi_accessions:
        args.ncbi_accessions = check_expected_single_column_input(args.ncbi_accessions, "-a")
    if args.genbank_files:
        args.genbank_files = check_expected_single_column_input(args.genbank_files, "-g")
    if args.fasta_files:
        args.fasta_files = check_expected_single_column_input(args.fasta_files, "-f")
    if args.amino_acid_files:
        args.amino_acid_files = check_expected_single_column_input(args.amino_acid_files, "-A")

    if args.resume:
        if os.path.exists(args.run_files_dir + "/run-data.json"):
            run_data = read_run_data(args.run_files_dir + "/run-data.json")
            if run_data.args_hash != args_hash_function(args):
                report_message("We are trying to resume a previous run (specified by the `-R` or `--resume` flag), "
                               "but it looks like the current arguments don't match the intial run's arguments. "
                               "Your best bet may be to just start a completely fresh run by adding the `-F` flag to "
                               "force-overwrite the previous outputs.")
                report_early_exit(None, copy_log = False)

    if "run_data" not in locals():
        run_data = populate_run_data(args)
        run_data.args_hash = args_hash_function(args)
        run_data = check_hmm_file(args, run_data)

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

    run_data.num_hmm_cpus = args.num_hmm_cpus

    return args, run_data


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
            report_early_exit(None, copy_log = False)

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


def check_expected_single_column_input(path, flag, get_count=False):

    check_path(path, flag)
    check_for_whitespace(path, flag)
    path = check_line_endings(path, flag)
    path = check_for_duplicates(path, flag)
    check_inputs_exist(path, flag)

    if get_count:
        with open(path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        return path, len(lines)

    return path


def check_for_whitespace(path, flag):

    # there should be no tabs or spaces in these input files
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    for line in lines:
        if " " in line:
            message = f'The specified input genomes file "{path}" (passed to `{flag}`) contains spaces in one or more entries. '
            message += f'That is not expected and will break things. Please double-check things and remove any whitespace from the file.'
            report_early_exit(None, suggest_help = True, copy_log = False)
        if "\t" in line:
            message = f'The specified input genomes file "{path}" (passed to `{flag}`) contains tabs in one or more entries. '
            message += f'That is not expected and will break things. Please double-check things and remove any whitespace from the file.'
            report_early_exit(None, suggest_help = True, copy_log = False)


def check_line_endings(path, flag):

    # checking for any CLRF line endings and creating a new file if so
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


def check_for_duplicates(path, flag):

    # checking for duplicates in the input file
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) != len(set(lines)):
        new_filename = f"{path}-unique"
        new_lines = list(set(lines))
        with open(new_filename, 'w') as f:
            f.write("\n".join(new_lines))

        report_message(f'Input file "{path}" (passed to `{flag}`) had duplicate entries that would have caused problems. '
                       f'A modified version was created, "{new_filename}", which will be used.')

        path = new_filename
        time.sleep(2)

    return path


def check_inputs_exist(path, flag):
    if flag in ["-g", "-f", "-A"]:
        # checking that all the input files actually exist
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if not os.path.exists(line):
                    report_message(f'The specified input-genome file "{line}" (passed to `{flag}`) does not exist.')
                    report_early_exit(None, copy_log = False)


def check_mapping_file(args, run_data, flag = "-m"):

    check_path(args.mapping_file, flag)
    mapping_file_problems = check_mapping_file_problem_chars_and_fields(args.mapping_file)
    if mapping_file_problems:
        report_problem_with_mapping_file(mapping_file_problems, args.mapping_file, flag = "-m")

    args.mapping_file = check_line_endings(args.mapping_file, flag)

    run_data.mapping_file_path = args.mapping_file

    # this is to handle when resuming runs, so that we don't overwrite the mapping_dict
    if len(run_data.mapping_dict) == 0:
        mapping_dict = make_mapping_dict(run_data.mapping_file_path)

        check_all_mapping_file_entries_are_in_input_genomes(mapping_dict, run_data)

        run_data.mapping_dict = mapping_dict

        # storing this so we are sure not to change any of these if we are adding taxonomic info later
        run_data.initial_mapping_IDs_from_user = list(mapping_dict.keys())

    return args, run_data


def check_mapping_file_problem_chars_and_fields(path):

    problematic_chars='()*&^#$@!/\\|[]:;'
    errors = []
    with open(path, 'r') as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            columns = line.split('\t')

            # Check if the number of columns is either 2 or 3.
            if len(columns) not in (2, 3):
                errors.append(f"Line {lineno} has {len(columns)} column(s) (expected 2 or 3): {line}")

            # Check each field for any problematic characters.
            for i, field in enumerate(columns, start=1):
                for char in problematic_chars:
                    if char in field:
                        # this is a special case for the "/" character, which is allowed in the first column
                        if not i == 1 and char == "/":
                            errors.append(
                                f"Line {lineno}, column {i} contains at least one problematic character ('{char}'):\n  {field}"
                            )
                            break
    return errors


def make_mapping_dict(path):
    """
    makes a dictionary mapping input genomes to the wanted labels
    key = original-input-genome-label, value = wanted-label

      1. Checks that every first-column entry is unique
      2. Uses the first column as the key
      3. If there are 2 columns, the value is the second column
      4. If there are 3 columns and the second column is non-empty,
         the value is "second-column_third-column" (joined with an underscore)
      5. If there are 3 columns and the second column is empty,
         the value is "first-column_third-column" (joined with an underscore)
    """

    df = pd.read_csv(path, sep="\t", header=None, dtype=str).fillna("")

    # checking none of the first column entries are duplicates
    if df[0].duplicated().any():
        report_message(
            f'The mapping file "{path}" (passed to `-m`) has some duplicate entries in the first column.'
            ' Please address that and try again.'
        )
        report_early_exit(None, copy_log = False)

    mapping_dict = {}
    for idx, row in df.iterrows():
        key = row[0].strip()
        col2 = row[1].strip() if len(row) > 1 else ""
        col3 = row[2].strip() if len(row) > 2 else ""

        if col3:
            if col2:
                value = f"{col2}_{col3}"
            else:
                value = f"{key}_{col3}"
        else:
            value = col2

        value = value.replace(" ", "_")
        mapping_dict[key] = value

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
            f'Each input genome must map to a unique label. Please address that and try again.'
        )
        report_early_exit(None, copy_log = False)

    return mapping_dict


def check_all_mapping_file_entries_are_in_input_genomes(mapping_dict, run_data):
    entries_in_mapping_file = set(mapping_dict.keys())
    # taking the basenames here because some inputs might have full/rel paths, but the mapping file shouldn't
    entries_in_input_genomes = set(run_data.get_all_input_genome_provided_paths()) | set(run_data.get_input_ncbi_accs())
    missing_keys = entries_in_mapping_file - entries_in_input_genomes
    if missing_keys:
        report_message(
            f'The mapping file "{run_data.mapping_file_path}" (passed to `-m`) specifies some input genomes that are not found in any of the the input-genome sources.\n\n'
            f'Problematic ones include:'
        )
        report_message(f'{"\n".join(missing_keys)}', ii="    ", si="    ")
        report_message(
            f"Each input genome in the mapping file (column 1) must be present in one of the input-genome sources. Please address this and try again."
        )
        report_early_exit(None, copy_log = False)


def check_for_required_dbs(args):
    if args.ncbi_accessions or args.add_ncbi_tax:
        get_ncbi_assembly_data()
    if args.add_ncbi_tax:
        get_ncbi_tax_data()
    if args.add_gtdb_tax:
        get_gtdb_data()
    if args.target_kos_file:
        get_kofamscan_data()


def track_tools_used(args, run_data):

    tools_used = ToolsUsed()

    if args.fasta_files:
        tools_used.prodigal_used = True
    if args.add_ncbi_tax:
        tools_used.taxonkit_used = True
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
    if total_input_genomes >= 1000 and total_input_genomes < 12500 and not args.no_super5:
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
    stdout_and_log("    Command entered:\n       ", full_execution_command, log_file=args.log_file, log_only=True)

    args, run_data = setup_tmp_dir(args, run_data)

    logs_dir_rel = os.path.join(args.output_dir, "logs")
    logs_dir = os.path.abspath(logs_dir_rel)
    os.makedirs(logs_dir, exist_ok=True)
    run_data.logs_dir_rel = logs_dir_rel
    run_data.logs_dir = logs_dir
    run_data.output_dir_rel = args.output_dir
    run_data.gtotree_logs_dir = os.path.join(logs_dir, "gtotree-logs")
    os.makedirs(run_data.gtotree_logs_dir, exist_ok=True)

    if args.ncbi_accessions:
        ncbi_downloads_dir_rel = os.path.join(args.tmp_dir, "ncbi-downloads")
        ncbi_downloads_dir = os.path.abspath(ncbi_downloads_dir_rel)
        os.makedirs(ncbi_downloads_dir, exist_ok=True)
        run_data.ncbi_downloads_dir_rel = ncbi_downloads_dir_rel
        run_data.ncbi_downloads_dir = ncbi_downloads_dir

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

    run_data.ready_genome_AA_files_dir = os.path.join(args.tmp_dir, "ready-genome-AA-files")
    os.makedirs(run_data.ready_genome_AA_files_dir, exist_ok=True)
    run_data.hmm_results_dir = os.path.join(args.tmp_dir, "hmm-results")
    os.makedirs(run_data.hmm_results_dir, exist_ok=True)
    run_data.found_SCG_seqs_dir = os.path.join(args.tmp_dir, "found-SCG-seqs")
    os.makedirs(run_data.found_SCG_seqs_dir, exist_ok=True)

    run_data.best_hit_mode = args.best_hit_mode
    run_data.output_dir = os.path.abspath(args.output_dir)
    run_data.seq_length_cutoff = args.seq_length_cutoff

    num_all_input_genomes = len(run_data.get_all_input_genome_ids())
    if num_all_input_genomes > 1000 and not args.no_super5:
        run_data.use_muscle_super5 = True

    run_data.num_muscle_threads = args.num_muscle_threads
    run_data.nucleotide_mode = args.nucleotide_mode

    if len(run_data.mapping_dict) > 0 or args.add_ncbi_tax or args.add_gtdb_tax:
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


# def check_resume_acceptable(run_data):
#     # this checks the run-data.json that exists has the same core components as the current run-data
#     pass


def setup_pfam_dirs(run_data):
    run_data.pfam_results_dir = os.path.join(run_data.output_dir, "pfam-search-results")
    run_data.pfam_results_dir_rel = os.path.join(run_data.output_dir_rel, "pfam-search-results")
    run_data.tmp_pfam_results_dir = os.path.join(run_data.tmp_dir, "pfam-search-results")


    dirs = [run_data.pfam_results_dir,
            run_data.tmp_pfam_results_dir,
            os.path.join(run_data.pfam_results_dir, "info"),
            os.path.join(run_data.pfam_results_dir, "individual-genome-results"),
            os.path.join(run_data.pfam_results_dir, "iToL-files"),
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
            os.path.join(run_data.ko_results_dir, "iToL-files"),
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
            report_early_exit(None, copy_log = False)
    else:
        tmp_dir = tempfile.mkdtemp(prefix = "gtt-tmp-", dir = args.output_dir)
        args.tmp_dir = tmp_dir
        run_data.tmp_dir = tmp_dir

    return args, run_data
