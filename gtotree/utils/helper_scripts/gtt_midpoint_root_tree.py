#!/usr/bin/env python

import argparse
import sys
import dendropy # type: ignore


def main():
    args = get_args()
    midpoint_root_tree(args.input_tree, args.output_tree)


def get_args():
    parser = argparse.ArgumentParser(description = "This script takes a newick tree as input and midpoint roots it.", add_help = False)

    required = parser.add_argument_group('REQUIRED PARAMETERS')
    optional = parser.add_argument_group('OPTIONAL PARAMETERS')

    required.add_argument("-i", "--input-tree", help = "input tree file", action = "store")
    optional.add_argument("-o", "--output-tree", help = 'output midpoint-rooted tree file (default: "midpoint-rooted.tre")', action = "store", default = "midpoint-rooted.tre")

    if len(sys.argv)==1 or sys.argv[1] in ["-h", "--help"]:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    return args


def midpoint_root_tree(input_tree, output_tree):

    tree = dendropy.Tree.get(path = input_tree, schema = "newick", preserve_underscores = True)
    tree.reroot_at_midpoint()
    tree.encode_bipartitions()
    tree.write(path = output_tree, schema = "newick", suppress_rooting = True, unquoted_underscores = True)


if __name__ == "__main__":
    main()
