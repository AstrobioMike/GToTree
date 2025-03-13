import shutil
from gtotree.utils.messaging import report_message, report_early_exit
from gtotree.utils.general import ToolsUsed

def perform_preflight_checks(args):
    check_for_essential_deps()
    check_optional_deps(args)
    check_lineage(args)
    check_tree_program(args)
    tools_used = track_tools_used(args)
    return tools_used


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

