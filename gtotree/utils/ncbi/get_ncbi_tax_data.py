#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
to download NCBI tax data for using TaxonKit (https://bioinf.shenwei.me/taxonkit/) to add NCBI taxonomy.

For examples, please visit the GToTree wiki here: https://github.com/AstrobioMike/GToTree/wiki/example-usage
"""

import sys
import os
import argparse
import tarfile
from gtotree.utils.messaging import wprint, color_text, report_message, report_early_exit
from gtotree.utils.general import download_with_tqdm

################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Helper program to setup NCBI taxonomy data for TaxonKit to add taxonomy info.",
        epilog="Example usage: gtt-get-ncbi-tax-data"
    )

    get_ncbi_tax_data()

################################################################################

def check_tax_location_var_is_set():
    try:
        ncbi_data_dir = os.environ['TAXONKIT_DB']
    except KeyError:
        wprint(color_text("The environment variable 'TAXONKIT_DB' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return ncbi_data_dir


def check_if_data_present(location):

    names_path = os.path.join(location, "names.dmp")
    nodes_path = os.path.join(location, "nodes.dmp")

    if not os.path.isfile(names_path) or not os.path.isfile(nodes_path):
        if os.path.exists(names_path):
            os.remove(names_path)
        if os.path.exists(nodes_path):
            os.remove(nodes_path)
        return False
    return True


def download_ncbi_tax_data(location):

    taxdump_path = os.path.join(location, "taxdump.tar.gz")

    taxdump_link = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

    try:
        download_with_tqdm(taxdump_link, "        NCBI Taxonomy data", taxdump_path)
    except Exception as e:
        report_message(f"Downloading the NCBI taxonomy data failed with the following error:\n{e}", "red")
        report_early_exit()

    with tarfile.open(taxdump_path) as tarball:
        tarball.extractall(location)

    os.remove(taxdump_path)


def get_ncbi_tax_data():

    ncbi_data_dir = check_tax_location_var_is_set()
    if check_if_data_present(ncbi_data_dir):
        return
    else:
        print(color_text("    Downloading required NCBI taxonomy data (only needs to be done once)...\n", "yellow"))
        download_ncbi_tax_data(ncbi_data_dir)


################################################################################


if __name__ == "__main__":
    main()