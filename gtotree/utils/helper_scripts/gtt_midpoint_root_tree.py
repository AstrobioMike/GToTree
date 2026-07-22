#!/usr/bin/env python

import argparse
import sys
import dendropy # type: ignore
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg


def build_parser(parent_subparsers=None):

    desc = "This program takes a newick tree as input and midpoint roots it."

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "midpoint-root-tree",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt midpoint-root-tree -i input.tre -o midpoint-rooted.tre`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument("-i", "--input-tree", metavar="<FILE>",
                          help="input tree file", action="store")
    optional.add_argument("-o", "--output-tree", metavar="<FILE>",
                          help='output midpoint-rooted tree file (default: "midpoint-rooted.tre")',
                          action="store", default="midpoint-rooted.tre")

    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    midpoint_root_tree(args.input_tree, args.output_tree)


def midpoint_root_tree(input_tree, output_tree):

    tree = dendropy.Tree.get(path = input_tree, schema = "newick", preserve_underscores = True)
    tree.reroot_at_midpoint()
    tree.encode_bipartitions()
    tree.write(path = output_tree, schema = "newick", suppress_rooting = True, unquoted_underscores = True)


if __name__ == "__main__":
    main()
