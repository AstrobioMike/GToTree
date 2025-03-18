import os
import sys
import shutil
import time
from gtotree.utils.messaging import (report_message,
                                     report_early_exit,
                                     report_missing_input_genomes_file,
                                     report_missing_pfam_targets_file,
                                     report_missing_ko_targets_file,
                                     report_notice,
                                     many_genomes_notice,
                                     few_genomes_notice,
                                     absurd_number_of_genomes_notice,
                                     gtotree_header,
                                     stdout_and_log
                                     )
from gtotree.utils.hmm_handling import check_hmm_file
from gtotree.utils.ncbi.get_ncbi_assembly_tables import get_ncbi_assembly_data
from gtotree.utils.ncbi.get_ncbi_tax_data import get_ncbi_tax_data
from gtotree.utils.gtdb.get_gtdb_data import get_gtdb_data
from gtotree.utils.kos.get_kofamscan_data import get_kofamscan_data
from gtotree.utils.general import ToolsUsed, log_file_var


def preflight_checks(args):
    check_for_essential_deps()
    args = primary_args_validation(args)
    check_for_required_dbs(args)
    tools_used = track_tools_used(args)
    args = setup_outputs_and_tmp_dir(args)
    return args, tools_used


def check_for_essential_deps():
    commands = ["muscle", "hmmsearch", "trimal"]
    for cmd in commands:
        program_check(cmd, essential = True)


def program_check(cmd, essential = False):
    path = shutil.which(cmd)
    if not path:
        if essential:
            report_early_exit(f"{cmd} is an essential dependency, but it's not in your PATH :(")
        else:
            report_early_exit(f"You specified to use {cmd}, but it's not in your PATH :(")


def primary_args_validation(args):
    check_for_minimum_args(args)
    check_optional_deps(args)
    check_lineage(args)
    check_tree_program(args)
    checks_for_nucleotide_mode(args)
    args = check_output_dir(args)
    args = check_input_files(args)
    return args


def check_for_minimum_args(args):
    if not args.ncbi_accessions and not args.genbank_files and not args.fasta_files and not args.amino_acid_files:
        report_message("You need to provide at least one input-genome source!")
        report_early_exit(suggest_help=True)
    if not args.hmm:
        report_message("You need to specify the HMM file of the target-SCGs you want to tree! "
                       "You can view the available gene-sets packaged with GToTree by running `gtt-hmms`.")
        report_early_exit(suggest_help=True)


def check_optional_deps(args):
    if args.add_ncbi_tax:
        program_check("taxonkit")


def check_lineage(args):
    accepted_ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
    lineage_list = args.lineage.split(",")

    for rank in lineage_list:
        if rank not in accepted_ranks:
            report_message(f'You specified "{args.lineage}" to the `-L` argument, but "{rank}" is not an accepted taxonomic rank.')
            print(f"\n  Accepted ranks are any combination of the below entered as a comma-delimited list:\n        {'\n        '.join(accepted_ranks)}")
            report_early_exit()

    if args.lineage != "Domain,Phylum,Class,Species" and not args.add_ncbi_tax and not args.add_gtdb_tax:
        report_message("You've specified a custom lineage (`-L`), but neither the "
                       "`-t` nor `-D` flags were provided to indicate which taxonomy to use.")
        report_early_exit(suggest_help=True)


def check_tree_program(args):
    accepted_programs = ["FastTree", "FastTreeMP", "VeryFastTree", "IQTREE"]
    if args.tree_program not in accepted_programs:
        report_message(f'You specified "{args.tree_program}" to the `-T` argument, but that\'s not an available treeing program.')
        print(f"\n  Available programs are:\n        {'\n        '.join(accepted_programs)}")
        report_early_exit()
    program_check(args.tree_program)


def checks_for_nucleotide_mode(args):
    if args.nucleotide_mode:
        if args.amino_acid_files:
            report_message("You've specified wanting to work with nucleotide sequences (by passing the `-z` parameter), "
                           "but also provided some input genomes as amino-acid files (passed to `-A`). We can't confidently reverse-translate "
                           "amino-acid seqs to nucleotide seqs, so we can't take both of those options.")
            report_early_exit()
        if args.genbank_files:
            report_message("You've specified wanting to work with nucleotide sequences (by passing the `-z` parameter), "
                           "but also provided some input genomes as genbank files (passed to `-g`). Input genbank files are currently "
                           "not supported with nucleotide mode.")
            report_early_exit()


def check_input_files(args):
    if args.ncbi_accessions:
        args.ncbi_accessions = check_expected_single_column_input(args.ncbi_accessions, "-a")
    if args.genbank_files:
        args.genbank_files = check_expected_single_column_input(args.genbank_files, "-g")
    if args.fasta_files:
        args.fasta_files = check_expected_single_column_input(args.fasta_files, "-f")
    if args.amino_acid_files:
        args.amino_acid_files = check_expected_single_column_input(args.amino_acid_files, "-A")

    args = check_hmm_file(args)

    if args.mapping_file:
        check_mapping_file(args)

    if args.target_pfam_file:
        args.target_pfam_file, total_pfam_targets = check_expected_single_column_input(args.target_pfam_file, "-p", get_count=True)
        args.total_pfam_targets = total_pfam_targets

    if args.target_ko_file:
        args.target_ko_file, total_ko_targets = check_expected_single_column_input(args.target_ko_file, "-K", get_count=True)
        args.total_ko_targets = total_ko_targets

    return args


def check_output_dir(args):
    if os.path.exists(args.output):
        if not args.force_overwrite:
            report_message(f'The output directory "{args.output}" already exists. '
                           'Please choose a different directory or add the `-F` flag to try to continue a prior run.')
            report_early_exit()
        else:
            args.output_already_existed = True
    else:
        args.output_already_existed = False
    return args


def check_path(path, flag):

    if not os.path.isfile(path):
        if flag == "-a" or flag == "-g" or flag == "-f" or flag == "-A":
            report_missing_input_genomes_file(path, flag)
        if flag == "-p":
            report_missing_pfam_targets_file(path, flag)
        if flag == "-K":
            report_missing_ko_targets_file(path, flag)


def check_expected_single_column_input(path, flag, get_count=False):

    check_path(path, flag)
    check_for_whitespace(path, flag)
    path = check_line_endings(path, flag)
    path = check_for_duplicates(path, flag)

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
            report_early_exit(f'The specified input genomes file "{path}" (passed to `{flag}`) contains spaces in one or more entries. '
                              f'That is not expected and will break things. Please double-check things and remove any whitespace from the file.', suggest_help=True)
        if "\t" in line:
            report_early_exit(f'The specified input genomes file "{path}" (passed to `{flag}`) contains tabs in one or more entries. '
                              f'That is not expected and will break things. Please double-check things and remove any whitespace from the file.', suggest_help=True)


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


def check_mapping_file(args):
    ## NEEDED ##
    # check for problem_chars='()*&^#$@!/\|[]'
    # check and fix windows line-endings
    # check for duplicate desired labels (might need to exit if this happens)
        # don't forget to check the related github issue for this
    pass


def check_for_required_dbs(args):
    if args.ncbi_accessions or args.add_ncbi_tax:
        get_ncbi_assembly_data()
    if args.add_ncbi_tax:
        get_ncbi_tax_data()
    if args.add_gtdb_tax:
        get_gtdb_data()
    if args.target_ko_file:
        get_kofamscan_data()


def track_tools_used(args):

    tools_used = ToolsUsed()

    if args.num_jobs > 1:
        tools_used.parallel_used = True
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
    if args.target_pfam_file:
        tools_used.pfam_db_used = True
    if args.target_ko_file:
        tools_used.kofamscan_used = True

    return tools_used


def check_input_genomes_amount(total_input_genomes, args):
    if total_input_genomes >= 1000 and total_input_genomes < 12500 and not args.no_super5:
        message = many_genomes_notice(total_input_genomes)
        report_notice(message)
        time.sleep(30)
    if total_input_genomes <= 20:
        message = few_genomes_notice(total_input_genomes, args)
        report_notice(message)
        time.sleep(5)
    if total_input_genomes >= 12500:
        message = absurd_number_of_genomes_notice(total_input_genomes)
        report_notice(message)
        time.sleep(60)


def check_and_report_any_changed_default_behavior(args):

    conditions = [
        args.output != "gtotree-output",
        args.mapping_file,
        args.nucleotide_mode,
        args.no_tree,
        args.add_gtdb_tax,
        args.add_ncbi_tax,
        args.lineage != "Domain,Phylum,Class,Species",
        args.tree_program != "FastTreeMP",
        args.best_hit,
        args.seq_length_cutoff != 0.2,
        args.genome_hits_cutoff != 0.5,
        args.num_jobs != 1,
        args.num_hmm_cpus != 2,
        args.muscle_threads != 5,
        args.no_super5,
        args.keep_gene_alignments,
        args.force_overwrite,
        args.debug,
    ]

    if any(conditions):
        report_message("  Other options set:")

    if args.output != "gtotree-output":
        print(f"      - The output directory has been set to: \"{args.output}\"")

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

    if args.best_hit:
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

    if args.force_overwrite:
        if args.output_already_existed:
            print("      - A previously generated output directory is being used, as the `-F` flag was provided")

    if args.debug:
        print("      - Debug mode is enabled")

    if args.target_pfam_file:
        print(f"      - Genomes will be searched for Pfams listed in: {args.target_pfam_file} ({args.total_pfam_targets} targets)")

    if args.target_ko_file:
        print(f"      - Genomes will be searched for KOs listed in: {args.target_ko_file} ({args.total_ko_targets} targets)")

    time.sleep(3)


def setup_outputs_and_tmp_dir(args):
    run_files_dir = os.path.join(args.output, "run-files")
    os.makedirs(run_files_dir, exist_ok=True)
    args.run_files_dir = run_files_dir
    log_file = os.path.join(args.output, "gtotree-runlog.txt")
    args.log_file = log_file
    log_file_var.set(log_file)
    full_execution_command = f"{' '.join(sys.argv)}"
    stdout_and_log(gtotree_header(), log_file=args.log_file, log_only=True, restart_log=True)
    stdout_and_log("    Command entered:\n       ", full_execution_command, log_file=args.log_file, log_only=True)

    # os.makedirs("gtotree-tmp")
    return args

