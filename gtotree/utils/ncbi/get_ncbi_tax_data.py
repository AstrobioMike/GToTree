#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import tarfile
from datetime import datetime, timezone
from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm


TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

NAMES_FILENAME = "names.dmp"
NODES_FILENAME = "nodes.dmp"

DATE_FILENAME = "date-retrieved.txt"


################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Setup the NCBI taxonomy data",
        epilog="Example usage: gtt-get-ncbi-tax-data"
    )

    parser.add_argument("-f", "--force-update",
                        help="Re-download the NCBI taxonomy data even if it is "
                             "already present",
                        action="store_true")

    args = parser.parse_args()

    get_ncbi_tax_data(force_update=args.force_update)

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
    """
    True only if both names.dmp and nodes.dmp are present and non-empty. If
    either is missing/empty, any partial state is cleaned up and we return False
    so a fresh copy is pulled.
    """
    names_path = os.path.join(location, NAMES_FILENAME)
    nodes_path = os.path.join(location, NODES_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not (is_nonempty_file(names_path) and is_nonempty_file(nodes_path)):
        for p in (names_path, nodes_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False
    return True


def report_ncbi_tax_dl_failure(e):
    report_message(f"Downloading the NCBI taxonomy data failed with the following error:\n{e}", "red")
    report_message("You can check whether you can access this URL:")
    report_message(f"    {TAXDUMP_URL}")
    report_message("If this keeps failing, it may be a transient network issue - trying again "
                   "later often resolves it.")
    report_early_exit(None, copy_log=False)


def download_ncbi_tax_data(location):
    """
    Download the taxdump tarball and extract it atomically into `location`:
    download to a temp path, extract into a temp staging dir, validate the two
    needed files landed, then move the extracted contents into place. Avoids
    leaving a half-populated data dir if the download or extraction is
    interrupted.
    """
    taxdump_path = os.path.join(location, "taxdump.tar.gz")
    staging_dir = os.path.join(location, ".taxdump-staging")

    # clean any leftover staging from a prior interrupted run
    if os.path.isdir(staging_dir):
        shutil.rmtree(staging_dir)
    os.makedirs(staging_dir, exist_ok=True)

    try:
        download_with_tqdm(TAXDUMP_URL, "        NCBI Taxonomy data", taxdump_path)
    except Exception as e:
        _safe_rmtree(staging_dir)
        _safe_remove(taxdump_path)
        report_ncbi_tax_dl_failure(e)

    try:
        with tarfile.open(taxdump_path) as tarball:
            tarball.extractall(staging_dir)
    except Exception as e:
        _safe_rmtree(staging_dir)
        _safe_remove(taxdump_path)
        report_ncbi_tax_dl_failure(e)

    os.remove(taxdump_path)
    print()

    # validate the two files TaxonKit needs actually landed in staging
    staged_names = os.path.join(staging_dir, NAMES_FILENAME)
    staged_nodes = os.path.join(staging_dir, NODES_FILENAME)
    if not (os.path.isfile(staged_names) and os.path.isfile(staged_nodes)):
        _safe_rmtree(staging_dir)
        report_ncbi_tax_dl_failure(
            f"the downloaded archive did not contain the expected "
            f"{NAMES_FILENAME} and {NODES_FILENAME}")

    # move staged pieces into place (clearing any stale copies first)
    _clear_partial_state(location)
    for entry in os.listdir(staging_dir):
        shutil.move(os.path.join(staging_dir, entry), os.path.join(location, entry))
    shutil.rmtree(staging_dir)

    # stamp when this was pulled (the taxdump has no version string of its own)
    with open(os.path.join(location, DATE_FILENAME), "w") as f:
        f.write(datetime.now(timezone.utc).strftime("%Y,%m,%d") + "\n")


def _clear_partial_state(location):
    """Remove the checked-for files (and the date stamp) if they exist."""
    for fname in (NAMES_FILENAME, NODES_FILENAME, DATE_FILENAME):
        fpath = os.path.join(location, fname)
        if os.path.exists(fpath):
            os.remove(fpath)


def _safe_remove(path):
    if os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass


def _safe_rmtree(path):
    if os.path.isdir(path):
        try:
            shutil.rmtree(path)
        except OSError:
            pass


def get_ncbi_tax_data(force_update=False):
    ncbi_data_dir = check_tax_location_var_is_set()

    if force_update:
        print(color_text("\n    Re-downloading NCBI taxonomy data (force-update requested)...\n", "yellow"))
        _clear_partial_state(ncbi_data_dir)
        download_ncbi_tax_data(ncbi_data_dir)
        return

    if check_if_data_present(ncbi_data_dir):
        return

    print(color_text("\n    Downloading required NCBI taxonomy data (only needs to be done once)...\n", "yellow"))
    download_ncbi_tax_data(ncbi_data_dir)


################################################################################


if __name__ == "__main__":
    main()
