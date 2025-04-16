#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
for setting up reference files for the glorious Genome Taxonomy Database (gtdb.ecogenomic.org/).

For examples, please visit the GToTree wiki here: https://github.com/AstrobioMike/GToTree/wiki/example-usage
"""

import sys
import os
import pandas as pd
import urllib
import argparse
from gtotree.utils.messaging import wprint, color_text
from gtotree.utils.general import download_with_tqdm



################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="This is a helper program to facilitate setting up the reference files for the \
                     glorious Genome Taxonomy Database (gtdb.ecogenomic.org). It's really meant for internal \
                     use only by the main GToTree program."
    )

    get_gtdb_data()

################################################################################


def check_gtdb_location_var_is_set():
    try:
        gtdb_data_dir = os.environ['GTDB_dir']
    except KeyError:
        wprint(color_text("The environment variable 'GTDB_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return gtdb_data_dir

    # # making sure path is writable for the user
    # path_writable = os.access(path, os.W_OK)

    # if not path_writable:
    #     print()
    #     wprint(color_text("The environment variable '" + str(variable) + "' does not seem to be writable, and we need it to be if wanting to add GTDB taxonomic lineages :(", "red"))
    #     print()
    #     wprint("Try to set it somewhere else with `gtt-data-locations set`, then run GToTree again.")
    #     print("\nExiting for now.\n")
    #     sys.exit(1)

    # return()


def gen_gtdb_tab(location):
    """ downloads and parses the GTDB info tables """

    # getting archaea
    arc_link = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz"
    arc_tsv_gz = download_with_tqdm(arc_link, "        GTDB archaeal data")
    arc_tab = pd.read_csv(arc_tsv_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
    arc_tab.rename(columns={arc_tab.columns[0]:"accession"}, inplace=True)
    arc_tab.dropna(inplace=True, how="all")

    # getting bacteria
    bac_link = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz"
    bac_tsv_gz = download_with_tqdm(bac_link, "        GTDB bacterial data")
    bac_tab = pd.read_csv(bac_tsv_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
    bac_tab.rename(columns={bac_tab.columns[0]:"accession"}, inplace=True)
    bac_tab.dropna(inplace=True, how="all")

    # combining
    gtdb_tab = pd.concat([arc_tab, bac_tab])

    # splitting gtdb taxonomy column into 7 and dropping the single column
    domain, phylum, rclass, order, family, genus, species = [], [], [], [], [], [], []

    for index, row in gtdb_tab.iterrows():
        curr_acc = row["accession"]
        tax_list = row["gtdb_taxonomy"].split(";")

        if len(tax_list) != 7:
            wprint(color_text("GTDB entry " + curr_acc + " doesn't seem to have 7-column lineage info. Something is likely wrong :(", "yellow"))
            print("")
            wprint("If this continues to happen, please file an issue at github.com/AstrobioMike/GToTree/issues")
            print("")
            wprint("Aborting for now.")
            print("")
            sys.exit(0)

        else:
            domain.append(tax_list[0][3:])
            phylum.append(tax_list[1][3:])
            rclass.append(tax_list[2][3:])
            order.append(tax_list[3][3:])
            family.append(tax_list[4][3:])
            genus.append(tax_list[5][3:])
            species.append(tax_list[6][3:])

    gtdb_tab.insert(1, "species", species)
    gtdb_tab.insert(1, "genus", genus)
    gtdb_tab.insert(1, "family", family)
    gtdb_tab.insert(1, "order", order)
    gtdb_tab.insert(1, "class", rclass)
    gtdb_tab.insert(1, "phylum", phylum)
    gtdb_tab.insert(1, "domain", domain)

    # writing out
    gtdb_tab.to_csv(location + "GTDB-arc-and-bac-metadata.tsv", index=False, sep="\t")

    gtdb_version_info = urllib.request.urlretrieve("https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION.txt", location + "GTDB-version-info.txt")


def check_and_or_get_gtdb_files(GTDB_dir):
    """ checks for and sets up ref GTDB files if needed """

    if os.path.exists(GTDB_dir + "GTDB-arc-and-bac-metadata.tsv") and os.path.exists(GTDB_dir + "GTDB-version-info.txt"):

        return None

    # generating when table doesn't exist yet
    else:
        wprint(color_text("Downloading and parsing archaeal and bacterial metadata tables from GTDB (only needs to be done once)...", "yellow"))
        print("")

        gen_gtdb_tab(GTDB_dir)


def get_gtdb_data():
    GTDB_dir = check_gtdb_location_var_is_set()
    check_and_or_get_gtdb_files(GTDB_dir)


if __name__ == "__main__":
    main()
