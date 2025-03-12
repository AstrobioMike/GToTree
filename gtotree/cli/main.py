import argparse
import sys
from gtotree.cli.helpmenu import helpmenu
from gtotree.utils.messaging import (wprint,
                             color_text,
                             report_failure)


def main():

    print(helpmenu)
    exit()
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        prog="GToTree",
        description="This program takes input genomes from various sources and ultimately produces a phylogenomic tree.\n"
                    "You can find detailed usage information at: github.com/AstrobioMike/GToTree/wiki",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- Required Inputs ---
    req = parser.add_argument_group("REQUIRED INPUTS")

    req.add_argument("-a", "--ncbi-accessions", metavar="<file>", type=str,
                     help="Single-column file of NCBI assembly accessions")
    req.add_argument("-g", "--genbank", metavar="<file>", type=str,
                     help="Single-column file with the paths to each GenBank file")
    req.add_argument("-f", "--fasta", metavar="<file>", type=str,
                     help="Single-column file with the paths to each fasta file")
    req.add_argument("-A", "--aminoacid", metavar="<file>", type=str,
                     help="Single-column file with the paths to each amino acid file (each file should hold the coding sequences for just one genome)")
    req.add_argument("-H", "--hmm", metavar="<file>", type=str, required=True,
                     help=("Location of the uncompressed target SCGs HMM file being used, or just the HMM name if the\n"
                           "'GToTree_HMM_dir' environment variable is set (set by conda install). Run 'gtt-hmms' by itself to view available gene-sets."))

    # --- Optional Inputs ---
    opt = parser.add_argument_group("OPTIONAL INPUTS")
    opt.add_argument("-o", "--output", metavar="<str>", type=str, default="GToTree_output",
                     help="Specify the desired output directory (default: GToTree_output)")
    opt.add_argument("-m", "--mapping", metavar="<file>", type=str,
                     help=("Mapping file specifying desired genome labels.\n"
                           "A two- or three-column tab-delimited file: column 1 holds the file name or NCBI accession,\n"
                           "column 2 holds the desired new genome label, and column 3 (optional) holds a string to append."))

    # --- Taxonomy Options ---
    tax = parser.add_argument_group("OPTIONS FOR ADDING TAXONOMY INFORMATION")
    tax.add_argument("-t", "--ncbi-tax", action="store_true",
                     help="Add NCBI taxonomy info to the sequence headers (effective for NCBI accessions and GenBank files)")
    tax.add_argument("-D", "--gtdb-tax", action="store_true",
                     help="Add GTDB taxonomy (only effective for input genomes provided as NCBI accessions)")
    tax.add_argument("-L", "--lineage", metavar="<str>", type=str, default="Domain,Phylum,Class,Species,Strain",
                     help=("Comma-separated list of taxonomic ranks to add (default: Domain,Phylum,Class,Species,Strain).\n"
                           "For example: -L Domain,Phylum,Class,Order,Family,Genus,Species,Strain"))

    # --- Filtering Settings ---
    filt = parser.add_argument_group("FILTERING SETTINGS")
    filt.add_argument("-c", "--seq-length-cutoff", metavar="<float>", type=float, default=0.2,
                      help="Sequence length cutoff as a float between 0-1 (default: 0.2)")
    filt.add_argument("-G", "--genome-hits-cutoff", metavar="<float>", type=float, default=0.5,
                      help="Minimum fraction of SCG hits required per genome (default: 0.5)")
    filt.add_argument("-B", "--best-hit", action="store_true",
                      help="Run GToTree in best-hit mode (by default, multiple hits per SCG are excluded)")

    # --- KO Searching ---
    ko = parser.add_argument_group("KO SEARCHING")
    ko.add_argument("-K", "--ko-file", metavar="<file>", type=str,
                    help="Single-column file of KO targets to search each genome for")

    # --- Pfam Searching ---
    pfam = parser.add_argument_group("PFAM SEARCHING")
    pfam.add_argument("-p", "--pfam-file", metavar="<file>", type=str,
                      help="Single-column file of Pfam targets to search each genome for")

    # --- General Run Settings ---
    run = parser.add_argument_group("GENERAL RUN SETTINGS")
    run.add_argument("-z", "--nucleotide-mode", action="store_true",
                     help="Use nucleotide sequences for alignment/tree (default uses amino acid sequences)\n"
                          "Note: Only works for NCBI accessions (-a) and fasta files (-f)")
    run.add_argument("-N", "--no-tree", action="store_true",
                     help="Generate alignment only, do not create a tree")
    run.add_argument("-k", "--keep-alignments", action="store_true",
                     help="Keep individual target gene alignment files")
    run.add_argument("-T", "--tree-program", metavar="<str>", type=str, default="FastTreeMP",
                     help=("Tree program to use for tree generation. Options include: FastTree, FastTreeMP,\n"
                           "VeryFastTree, and IQ-TREE. (Default: FastTreeMP if available, otherwise FastTree)"))
    run.add_argument("-n", "--num-cpus", metavar="<int>", type=int, default=2,
                     help="Number of CPUs for the HMM search (default: 2)")
    run.add_argument("-M", "--muscle-threads", metavar="<int>", type=int, default=5,
                     help="Number of threads for muscle alignment (default: 5)")
    run.add_argument("-j", "--num-jobs", metavar="<int>", type=int, default=1,
                     help="Number of parallel jobs (default: 1)")
    run.add_argument("-X", "--no-super5", action="store_true",
                     help="Do not use the 'super5' muscle alignment algorithm (useful when working with >1,000 genomes)")
    run.add_argument("-P", "--http", action="store_true",
                     help="Use http instead of ftp for downloads")
    run.add_argument("-F", "--force-overwrite", action="store_true",
                     help="Force overwrite of the output directory if it exists")
    run.add_argument("-d", "--debug", action="store_true",
                     help="Enable debug mode (keeps temporary directory for debugging)")

    # --- Parse Arguments ---
    args = parser.parse_args()

    # For now, we'll just print the parsed arguments.
    # In your final implementation, you would pass these to your workflow functions.
    print(args)

if __name__ == "__main__":
    main()
