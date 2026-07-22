#!/usr/bin/env python

import os
import sys
import argparse
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.messaging import color_text, report_message
from gtotree.utils.hmms.scg_hmm_setup import (check_hmm_location_var_is_set,
                                              read_in_hmm_summary_table)


def build_parser(parent_subparsers=None):

    desc = ("This program lists the pre-packaged HMM SCG-sets available in GToTree. "
            "See github.com/AstrobioMike/GToTree/wiki/SCG-sets for more info.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "hmms",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt hmms`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    optional = parser.add_argument_group("Optional Parameters")
    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()
    parser.parse_args()

    print(color_text("\n                   GToTree pre-packaged HMM SCG-sets", "yellow"))
    print("   See github.com/AstrobioMike/GToTree/wiki/SCG-sets for more info\n")

    hmm_data_dir = get_writable_hmm_dir()
    df = read_in_hmm_summary_table()

    report_available_scg_sets(hmm_data_dir, df)


def get_writable_hmm_dir():
    hmm_data_dir = check_hmm_location_var_is_set()

    if not os.path.isdir(hmm_data_dir):
        try:
            os.makedirs(hmm_data_dir, exist_ok=True)
        except OSError:
            report_message(
                "The 'GToTree_HMM_dir' location does not exist and can't be "
                "created :( Use `gtt data locations check` and `gtt data locations set` to check and configure.",
                "yellow",
            )
            sys.exit(1)

    if not os.access(hmm_data_dir, os.W_OK):
        report_message(
            "The 'GToTree_HMM_dir' location is not writable for you :( "
            "Use `gtt data locations check` and `gtt data locations set` to check and configure.",
            "yellow",
        )
        sys.exit(1)

    return hmm_data_dir


def report_available_scg_sets(hmm_data_dir, df):
    print(f"   The environment variable {color_text('GToTree_HMM_dir', 'green')} is set to:")
    print(f"     {hmm_data_dir}\n")

    print(f"   The {len(df)} available pre-packaged HMM SCG-sets include:\n")

    for _, row in df.iterrows():
        gene_set = str(row["file"]).replace(".hmm", "")
        num_genes = row["num_genes"]
        print("\t   {:<30} {:>14}".format(gene_set, f"({num_genes} genes)"))

    table_path = os.path.join(hmm_data_dir, "hmm-sources-and-info.tsv")
    print(f"\n   Details can be found in: \n     {table_path}\n")


if __name__ == "__main__":
    main()
