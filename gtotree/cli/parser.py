import argparse
import sys
from gtotree.cli.helpmenu import helpmenu
from gtotree.utils.messaging import gtotree_header, get_version
from gtotree.utils.preflight_checks import check_for_essential_deps

class CustomHelpAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings,
                         dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        sys.stdout.write(helpmenu)
        parser.exit()

def parser():

    parser = argparse.ArgumentParser(
        add_help=False,
        usage=argparse.SUPPRESS,
        description="",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-h", "--help", action=CustomHelpAction,
        help="Show this help message and exit"
    )
    parser.add_argument("-v", "--version", action="version", version=f"GToTree v{get_version()}")


    # --- Primary Inputs ---
    primary = parser.add_argument_group("Primary Inputs")

    primary.add_argument("-a", "--ncbi-accessions", metavar="<file>", type=str)
    primary.add_argument("-g", "--genbank-files", metavar="<file>", type=str)
    primary.add_argument("-f", "--fasta-files", metavar="<file>", type=str)
    primary.add_argument("-A", "--amino-acid-files", metavar="<file>", type=str)
    primary.add_argument("-H", "--hmm", metavar="<file>", type=str)

    # --- Optional Inputs ---
    opt = parser.add_argument_group("Optional Inputs")
    opt.add_argument("-o", "--output", metavar="<dir>", type=str, default="gtotree-output")
    opt.add_argument("-m", "--mapping-file", metavar="<file>", type=str)

    # --- Taxonomy Options ---
    opt.add_argument("-t", "--add-ncbi-tax", action="store_true")
    opt.add_argument("-D", "--add-gtdb-tax", action="store_true")
    opt.add_argument("-L", "--lineage", metavar="<str>", type=str, default="Domain,Phylum,Class,Species")

    # --- Filtering Settings ---
    opt.add_argument("-c", "--seq-length-cutoff", metavar="<float>", type=float, default=0.2)
    opt.add_argument("-G", "--genome-hits-cutoff", metavar="<float>", type=float, default=0.5)
    opt.add_argument("-B", "--best-hit", action="store_true")

    # --- Addtional Target Searching ---
    opt.add_argument("-K", "--target-ko-file", metavar="<file>", type=str)
    opt.add_argument("-p", "--target-pfam-file", metavar="<file>", type=str)

    # --- General Run Settings ---
    opt.add_argument("-z", "--nucleotide-mode", action="store_true")
    opt.add_argument("-N", "--no-tree", action="store_true")
    opt.add_argument("-k", "--keep-gene-alignments", action="store_true")
    opt.add_argument("-T", "--tree-program", metavar="<str>", type=str, default="FastTreeMP")
    opt.add_argument("-j", "--num-jobs", metavar="<int>", type=int, default=1)
    opt.add_argument("-n", "--num-hmm-cpus", metavar="<int>", type=int, default=2)
    opt.add_argument("-M", "--muscle-threads", metavar="<int>", type=int, default=5)
    opt.add_argument("-X", "--no-super5", action="store_true")
    opt.add_argument("-P", "--use-http", action="store_true")
    opt.add_argument("-F", "--force-overwrite", action="store_true")
    opt.add_argument("--tmp-dir", default=None)
    opt.add_argument("-d", "--debug", action="store_true")

    if len(sys.argv) == 1:
        sys.stdout.write(helpmenu)
        sys.exit()

    return parser
