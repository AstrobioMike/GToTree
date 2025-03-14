import os
import shutil
import time
from gtotree.utils.messaging import (report_message,
                                     report_early_exit,
                                     report_missing_input_genomes_file,
                                     )
from gtotree.utils.ncbi.get_ncbi_assembly_tables import get_ncbi_assembly_data
from gtotree.utils.ncbi.get_ncbi_tax_data import get_ncbi_tax_data
from gtotree.utils.kos.get_kofamscan_data import get_kofamscan_data
from gtotree.utils.general import ToolsUsed


def preflight_checks(args):
    check_for_essential_deps()
    args = primary_args_validation(args)
    check_for_required_dbs(args)
    tools_used = track_tools_used(args)
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
    args = checks_for_input_files(args)
    check_mapping_file(args)
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


def checks_for_input_files(args):
    if args.ncbi_accessions:
        check_path(args.ncbi_accessions, "-a")
        path = check_expected_single_column_input(args.ncbi_accessions, "-a")
        args.ncbi_accessions = path
    if args.genbank_files:
        check_path(args.genbank_files, "-g")
        path = check_expected_single_column_input(args.genbank_files, "-g")
        args.genbank_files = path
    if args.fasta_files:
        check_path(args.fasta_files, "-f")
        path = check_expected_single_column_input(args.fasta_files, "-f")
        args.fasta_files = path
    if args.amino_acid_files:
        check_path(args.amino_acid_files, "-A")
        path = check_expected_single_column_input(args.amino_acid_files, "-A")
        args.amino_acid_files = path
    return args


def check_path(path, flag):
    if not os.path.isfile(path):
        report_missing_input_genomes_file(path, flag)


def check_expected_single_column_input(path, flag):

    check_for_whitespace(path, flag)
    path = check_line_endings(path, flag)
    path = check_for_duplicates(path, flag)

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
        get_ncbi_assembly_data(use_http=args.use_http)
    if args.add_ncbi_tax:
        get_ncbi_tax_data()
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
    if args.hmm == "Universal" or args.hmm == "Universal-Hug-et-al":
        tools_used.universal_SCGs_used = True
    if args.target_pfam_file:
        tools_used.pfam_db_used = True
    if args.target_ko_file:
        tools_used.kofamscan_used = True

    return tools_used
