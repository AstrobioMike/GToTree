#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
for setting up reference files for the glorious Genome Taxonomy Database
(gtdb.ecogenomic.org/).

It ensures the (slim) GTDB metadata table is present, downloading bit's
pre-slimmed tarball if it isn't.
"""

import sys
import os
import socket
import shutil
import tarfile
import time
import argparse
import urllib
import urllib.error
import pandas as pd # type: ignore

from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm


GTDB_BASE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"

# pre-slimmed table + version-info, packaged as a single .tar.gz and hosted on a
# bit GitHub release for a fast default download. The archive contains exactly
# two files at its root:
#   GTDB-arc-and-bac-metadata.tsv   (already slimmed to GTDB_KEPT_COLUMNS)
#   GTDB-version-info.txt
# Used by default; `-f/--force-update` re-pulls this same tarball. If the
# download fails, the default path falls back to rebuilding from GTDB_BASE_URL.
GTDB_SLIM_TARBALL_URL = "https://github.com/AstrobioMike/bit/releases/download/gtdb-metadata-latest/GTDB-arc-and-bac-metadata.tar.gz"

# the stored GTDB-arc-and-bac-metadata.tsv is slimmed to only the columns
# GToTree uses. The upstream GTDB metadata tables carry ~120 columns; keeping
# only these shrinks the stored table substantially and speeds every downstream
# read. The order here is the on-disk column order. Any column absent from a
# given GTDB release is silently skipped at write time.
GTDB_KEPT_COLUMNS = [
    "accession", "ncbi_genbank_assembly_accession", "ncbi_taxid",
    "gtdb_representative", "ncbi_refseq_category",
    "domain", "phylum", "class", "order", "family", "genus", "species",
    "checkm2_completeness", "checkm2_contamination",
    "checkm_completeness", "checkm_contamination",
    "genome_size", "contig_count", "gc_count", "gc_percentage", "ambiguous_bases",
    "coding_bases", "coding_density",
]

# names of the two files packaged inside the distributed tarball
GTDB_TABLE_FILENAME = "GTDB-arc-and-bac-metadata.tsv"
GTDB_VERSION_FILENAME = "GTDB-version-info.txt"


################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Setup the GTDB data")

    parser.add_argument("-f", "--force-update",
                        help="Re-download the slim GTDB table even if it is "
                             "already present",
                        action="store_true")

    args = parser.parse_args()

    get_gtdb_data(force_update=args.force_update)

################################################################################


def get_gtdb_data(force_update=False):
    GTDB_dir = check_gtdb_location_var_is_set()
    data_present = check_if_gtdb_data_present(GTDB_dir)

    if data_present and not force_update:
        return GTDB_dir

    get_slim_gtdb_tab(GTDB_dir)
    return GTDB_dir


def check_gtdb_location_var_is_set():
    try:
        gtdb_data_dir = os.environ['GTDB_dir']
    except KeyError:
        wprint(color_text("The environment variable 'GTDB_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return gtdb_data_dir


def check_if_gtdb_data_present(location):
    """
    True if both the slim metadata table and version-info file are present and
    non-empty. If either is missing/empty, any stray copy is cleaned up and we
    return False so a fresh copy is pulled.
    """
    metadata_path = os.path.join(str(location), GTDB_TABLE_FILENAME)
    version_info_path = os.path.join(str(location), GTDB_VERSION_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(metadata_path) or not is_nonempty_file(version_info_path):
        for p in (metadata_path, version_info_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False

    return True


def report_gtdb_unreachable(err):

    print("")
    wprint(color_text("  The GTDB data server could not be reached :(", "red"))
    report_message(f"While trying to download from: {GTDB_BASE_URL}", color=None,
                   ii="    ", si="      ")
    report_message("This is typically a network/connectivity problem on the system running "
           "GToTree (e.g., no internet access, a firewall or proxy blocking the "
           "connection, or the GTDB server being temporarily unavailable), rather "
           "than a problem with GToTree itself.", color=None, ii="    ",
           si="    ")
    report_message(("Things to try:"), color=None, ii="    ", si="    ")
    print("          - confirm the system has internet access")
    print(f"          - confirm {GTDB_BASE_URL} is reachable (e.g., in a browser or with curl)")
    print("          - if behind a proxy/firewall, check that outbound HTTPS is allowed")
    print("          - wait and try again later in case the server is temporarily down")
    report_message(f"Underlying error: {err}", ii="    ", si="    ")
    print("")
    report_early_exit(None, copy_log=False)


def get_slim_gtdb_tab(location):
    """
    default path: download bit's pre-slimmed GTDB tarball and extract its two
    files (the slimmed metadata table and the version-info file) into `location`.
    Falls back to rebuilding from the upstream GTDB release (gen_gtdb_tab) if no
    tarball URL is configured or the download/extract fails, so the user always
    ends up with a usable table.
    """
    metadata_path = os.path.join(location, GTDB_TABLE_FILENAME)
    version_info_path = os.path.join(location, GTDB_VERSION_FILENAME)

    if not GTDB_SLIM_TARBALL_URL:
        gen_gtdb_tab(location)
        return

    print(color_text("\n    Downloading the prepared GTDB table (only needs to be done once)...\n", "yellow"))

    tarball_path = os.path.join(location, "GTDB-slim.tar.gz")
    expected = {GTDB_TABLE_FILENAME, GTDB_VERSION_FILENAME}

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        download_with_tqdm(GTDB_SLIM_TARBALL_URL, "        GTDB prepared data", tarball_path, speed_gate=True)

        # a truncated/corrupt download trips here (full-stream read via
        # getmembers), which triggers the fallback below rather than writing a
        # corrupt table.
        with tarfile.open(tarball_path, "r:gz") as tar:
            members = tar.getmembers()
            # take only the two expected files, matched by basename, and guard
            # against path traversal / nested dirs by extracting to flat names.
            wanted = {}
            for m in members:
                base = os.path.basename(m.name)
                if base in expected and m.isfile():
                    wanted[base] = m
            missing = expected - set(wanted)
            if missing:
                raise ValueError(
                    f"prepared GTDB archive is missing expected file(s): {', '.join(sorted(missing))}")

            for base, m in wanted.items():
                src = tar.extractfile(m)
                if src is None:
                    raise ValueError(f"could not read '{base}' from prepared GTDB archive")
                dest = os.path.join(location, base)
                with src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            tarfile.TarError, ValueError, OSError) as err:
        # clean up partial artifacts, then fall back to the upstream rebuild
        for p in (tarball_path, metadata_path, version_info_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        print("")
        wprint(color_text("  Couldn't get the prepared GTDB table; "
                          "rebuilding from the upstream GTDB release instead.", "yellow"))
        report_message(f"Underlying issue: {err}", ii="    ",
                       si="    ")
        print("")
        gen_gtdb_tab(location)
        return
    finally:
        socket.setdefaulttimeout(default_timeout)
        if os.path.exists(tarball_path):
            try:
                os.remove(tarball_path)
            except OSError:
                pass
    print("")


def gen_gtdb_tab(location):
    """
    rebuild fallback: download and parse the GTDB metadata tables directly from
    the upstream release, split the taxonomy string into 7 rank columns, slim to
    GTDB_KEPT_COLUMNS, and stamp the version-info file.
    """

    print(color_text("\n    Downloading and parsing metadata tables from GTDB (only needs to be done once)...\n", "yellow"))

    arc_path = os.path.join(location, "ar53_metadata.tsv.gz")
    bac_path = os.path.join(location, "bac120_metadata.tsv.gz")
    metadata_path = os.path.join(location, GTDB_TABLE_FILENAME)
    version_info_path = os.path.join(location, GTDB_VERSION_FILENAME)

    # fail fast (instead of hanging) if the server can't be reached
    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)

    try:
        # getting archaea
        arc_link = f"{GTDB_BASE_URL}/ar53_metadata.tsv.gz"
        download_with_tqdm(arc_link, "        GTDB archaeal data", arc_path)
        arc_tab = pd.read_csv(arc_path, sep="\t", compression="gzip", on_bad_lines='skip', header=0, low_memory=False)
        arc_tab.rename(columns={arc_tab.columns[0]: "accession"}, inplace=True)
        arc_tab.dropna(inplace=True, how="all")

        # getting bacteria
        bac_link = f"{GTDB_BASE_URL}/bac120_metadata.tsv.gz"
        download_with_tqdm(bac_link, "        GTDB bacterial data", bac_path)
        bac_tab = pd.read_csv(bac_path, sep="\t", compression="gzip", on_bad_lines='skip', header=0, low_memory=False)
        bac_tab.rename(columns={bac_tab.columns[0]: "accession"}, inplace=True)
        bac_tab.dropna(inplace=True, how="all")
    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as err:
        report_gtdb_unreachable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)
        for p in (arc_path, bac_path):
            if os.path.exists(p):
                os.remove(p)

    # combining
    gtdb_tab = pd.concat([arc_tab, bac_tab])

    print("")

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

    # writing out, slimmed to only the columns GToTree uses (intersect with
    # present columns so a release missing one of them doesn't error)
    keep = [c for c in GTDB_KEPT_COLUMNS if c in gtdb_tab.columns]
    gtdb_tab = gtdb_tab[keep]
    gtdb_tab.to_csv(metadata_path, index=False, sep="\t")

    try:
        socket.setdefaulttimeout(30)
        urllib.request.urlretrieve(f"{GTDB_BASE_URL}/VERSION.txt", version_info_path)
    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as err:
        report_gtdb_unreachable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)


if __name__ == "__main__":
    main()
