import argparse
import sys
from gtotree.cli.helpmenu import print_helpmenu
from gtotree.utils.misc.messaging import get_version


class CustomHelpAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings,
                         dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print_helpmenu(detailed=False)
        parser.exit()


class DetailedHelpAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings,
                         dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print_helpmenu(detailed=True)
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
        help="Show the condensed help menu and exit"
    )
    parser.add_argument(
        "-s", "--show-detailed-help", action=DetailedHelpAction,
        help="Show the detailed help menu and exit"
    )
    parser.add_argument("-v", "--version", action="version", version=f"GToTree v{get_version()}")


    # --- Primary Inputs ---
    primary = parser.add_argument_group("Primary Inputs")

    primary.add_argument("-w", "--wanted-ref-tax", metavar="<str>", type=str, default=None)
    primary.add_argument("-a", "--ncbi-accessions", metavar="<file>", type=str)
    primary.add_argument("-g", "--genbank-files", metavar="<file>", type=str)
    primary.add_argument("-f", "--fasta-files", metavar="<file>", type=str)
    primary.add_argument("-A", "--amino-acid-files", metavar="<file>", type=str)
    primary.add_argument("-H", "--hmm", metavar="<file>", type=str)

    # not a flag: set by preflight when `-H` was left off and a set was chosen from `-w`,
    # so the reporting can say so. Declared here to keep it on args from the start
    primary.set_defaults(hmm_auto_selected=None)

    # --- Optional Inputs ---
    opt = parser.add_argument_group("Optional Inputs")

    opt.add_argument("-o", "--output-dir", metavar="<dir>", type=str, default="gtotree-output")
    opt.add_argument("-m", "--mapping-file", metavar="<file>", type=str)

    # --- Adding Reference Genomes by Taxonomy ---
    opt.add_argument("--source", metavar="<str>", type=str.lower,
                     choices=["gtdb", "ncbi"], default="gtdb")
    opt.add_argument("--target-rank", metavar="<str>", type=str, default=None)
    opt.add_argument("--derep-rank", metavar="<str>", type=str, default="auto")

    # --- Taxonomy Options ---
    opt.add_argument("-D", "--add-gtdb-tax", action="store_true")
    opt.add_argument("-t", "--add-ncbi-tax", action="store_true")
    opt.add_argument("-L", "--lineage-ranks", metavar="<str>", type=str, dest="lineage", default="domain,phylum,class,genus,species")

    # --- Filtering Settings ---
    opt.add_argument("-c", "--seq-length-cutoff", metavar="<float>", type=float, default=0.2)
    opt.add_argument("-r", "--gene-rep-cutoff", metavar="<float>", type=float, dest="gene_representation_cutoff", default=0.1)
    opt.add_argument("-G", "--genome-hits-cutoff", metavar="<float>", type=float, default=0.5)
    opt.add_argument("-B", "--best-hit-mode", action="store_true")

    # --- Addtional Target Searching ---
    opt.add_argument("-K", "--target-kos", metavar="<file>", type=str, dest="target_kos_file")
    opt.add_argument("-p", "--target-pfams", metavar="<file>", type=str, dest="target_pfams_file")

    # --- General Run Settings ---
    opt.add_argument("-j", "--num-jobs", metavar="<int>", type=int, default=4)
    opt.add_argument("-M", "--muscle-threads", metavar="<int>", type=int, default=5, dest="num_muscle_threads")
    opt.add_argument("-X", "--no-super5", action="store_true")
    opt.add_argument("-N", "--no-tree", action="store_true")
    opt.add_argument("-T", "--tree-program", metavar="<str>", type=str, default="FastTreeMP")
    opt.add_argument("-z", "--nucleotide-mode", action="store_true")
    opt.add_argument("-k", "--keep-gene-alignments", action="store_true")
    opt.add_argument("-R", "--resume", action="store_true")
    opt.add_argument("-F", "--force-overwrite", action="store_true")
    opt.add_argument("--tmp-dir", default=None)
    opt.add_argument("-d", "--keep-working-dir", action="store_true")

    return parser
