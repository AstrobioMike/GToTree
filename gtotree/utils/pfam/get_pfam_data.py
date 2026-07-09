#!/usr/bin/env python
"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
to download and setup the Pfam (https://www.ebi.ac.uk/interpro/) data files for use.

Default behavior: download the pre-built Pfam data packaged as a single .tar.gz
and rebuilt monthly by GToTree's refresh-pfam-data GitHub Action (if needed), hosted at a
fixed release asset.

`-f/--force-update` re-pulls that same hosted tarball.
"""

import sys
import os
import argparse
import shutil
import gzip
import tarfile
from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm


# pre-built Pfam data, packaged as a single .tar.gz and rebuilt monthly by
# GToTree's refresh-pfam-data GitHub Action, then re-hosted at this fixed asset
# URL. The archive holds exactly this at its root:
#   Pfam-A.hmm.gz        (the master HMM file, compressed)
#   pfamA.txt.gz         (the info table, compressed)
#   pfam-version.txt     (the resolved version string, e.g. "38.2")
#   date-retrieved.txt   (build date, YYYY,MM,DD)
PFAM_TARBALL_URL = "https://github.com/AstrobioMike/GToTree/releases/download/pfam-data-latest/pfam-data.tar.gz"

# files shipped in the asset
HMM_GZ_FILENAME = "Pfam-A.hmm.gz"
INFO_GZ_FILENAME = "pfamA.txt.gz"
VERSION_FILENAME = "pfam-version.txt"
DATE_FILENAME = "date-retrieved.txt"

# files produced by decompressing the two .gz files, ready for immediate use
HMM_FILENAME = "Pfam-A.hmm"
INFO_FILENAME = "pfamA.txt"

# (gz -> decompressed) pairs to expand after extraction
_DECOMPRESS_PAIRS = ((HMM_GZ_FILENAME, HMM_FILENAME),
                     (INFO_GZ_FILENAME, INFO_FILENAME))


################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Setup the Pfam data",
        epilog="Example usage: gtt-get-pfam-data"
    )

    parser.add_argument("-f", "--force-update",
                        help="Re-download the Pfam data even if it is already present",
                        action="store_true")

    args = parser.parse_args()

    get_pfam_data(force_update=args.force_update)

################################################################################


def check_location_var_is_set():
    """
    Ensure that the environment variable 'Pfam_data_dir' is set.
    Returns the directory specified by the environment variable.
    """
    try:
        pfam_data_dir = os.environ['Pfam_data_dir']
    except KeyError:
        wprint(color_text("The environment variable 'Pfam_data_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return pfam_data_dir


def check_if_data_present(location):
    """
    True only if all expected pieces are present, including the DECOMPRESSED
    files (that's the ready-to-use state). If anything is missing, remove
    whatever partial state exists so the subsequent download starts clean.
    """
    required = (HMM_FILENAME, INFO_FILENAME, VERSION_FILENAME)
    if all(os.path.isfile(os.path.join(location, f)) for f in required):
        return True

    _clear_partial_state(location)
    return False


def _clear_partial_state(location):
    """Remove any of the expected files (compressed or decompressed) that exist."""
    for fname in (HMM_GZ_FILENAME, INFO_GZ_FILENAME, HMM_FILENAME, INFO_FILENAME,
                  VERSION_FILENAME, DATE_FILENAME):
        fpath = os.path.join(location, fname)
        if os.path.exists(fpath):
            os.remove(fpath)


def get_stored_pfam_version(location):
    """Return the version string recorded in pfam-version.txt, or None."""
    version_path = os.path.join(location, VERSION_FILENAME)
    try:
        with open(version_path) as f:
            return f.read().strip()
    except (FileNotFoundError, OSError):
        return None


def report_pfam_dl_failure(e):
    report_message(f"Downloading the Pfam data failed with the following error:\n{e}", "red")
    report_message("GToTree pulls this data from a GitHub-hosted release asset over HTTPS.")
    report_message("You can check whether you can access this URL:")
    report_message(f"    {PFAM_TARBALL_URL}")
    report_message("The problem may be a transient network issue, and trying again "
                   "later often resolves it. If you're using this in a regular GToTree run, "
                   "you can also drop the additional target Pfams (being provided "
                   "to the `-p` flag) and run GToTree without them to get around this.")
    report_early_exit(None, copy_log=False)


def _gunzip(src_path, dest_path):
    """Decompress a .gz file to dest_path."""
    with gzip.open(src_path, "rb") as f_in, open(dest_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def download_pfam_data(location):
    """
    Pull the hosted tarball and set it up atomically into `location`:
    download to a temp path, extract into a temp staging dir, decompress the two
    .gz files there, validate, then move everything into place. Nothing lands in
    the real data dir until the decompressed files exist, so an interrupted
    download or gunzip never leaves a half-populated dir.
    """
    tarball_path = os.path.join(location, "pfam-data.tar.gz")
    staging_dir = os.path.join(location, ".pfam-staging")

    if os.path.isdir(staging_dir):
        shutil.rmtree(staging_dir)
    os.makedirs(staging_dir, exist_ok=True)

    try:
        download_with_tqdm(PFAM_TARBALL_URL, "        Pfam data", tarball_path, speed_gate=True)
    except Exception as e:
        _safe_rmtree(staging_dir)
        _safe_remove(tarball_path)
        report_pfam_dl_failure(e)

    try:
        with tarfile.open(tarball_path) as tarball:
            tarball.extractall(staging_dir)
    except Exception as e:
        _safe_rmtree(staging_dir)
        _safe_remove(tarball_path)
        report_pfam_dl_failure(e)

    os.remove(tarball_path)

    # confirm the shipped .gz files and version stamp landed in staging
    for fname in (HMM_GZ_FILENAME, INFO_GZ_FILENAME, VERSION_FILENAME):
        if not os.path.isfile(os.path.join(staging_dir, fname)):
            _safe_rmtree(staging_dir)
            report_pfam_dl_failure(
                f"the downloaded archive did not contain the expected {fname}")

    # decompress the two .gz files in staging, ready for immediate use
    print(color_text("    Decompressing...\n", "yellow"))
    for gz_name, out_name in _DECOMPRESS_PAIRS:
        gz_path = os.path.join(staging_dir, gz_name)
        out_path = os.path.join(staging_dir, out_name)
        try:
            _gunzip(gz_path, out_path)
        except Exception as e:
            _safe_rmtree(staging_dir)
            report_pfam_dl_failure(f"failed to decompress {gz_name}: {e}")

    # validate the decompressed files are present and non-empty
    for out_name in (HMM_FILENAME, INFO_FILENAME):
        out_path = os.path.join(staging_dir, out_name)
        if not (os.path.isfile(out_path) and os.path.getsize(out_path) > 0):
            _safe_rmtree(staging_dir)
            report_pfam_dl_failure(
                f"decompressed file {out_name} is missing or empty")

    # move staged pieces into place (clearing any stale copies first), dropping
    # the compressed originals -- only the decompressed files (ready for
    # immediate use), the version stamp, and the date stamp are kept.
    _clear_partial_state(location)
    for entry in os.listdir(staging_dir):
        if entry in (HMM_GZ_FILENAME, INFO_GZ_FILENAME):
            continue
        shutil.move(os.path.join(staging_dir, entry), os.path.join(location, entry))
    shutil.rmtree(staging_dir)
    print()


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


def get_pfam_data(force_update=False):
    pfam_data_dir = check_location_var_is_set()

    if force_update:
        print(color_text("\n    Re-downloading Pfam data (force-update requested)...\n", "yellow"))
        _clear_partial_state(pfam_data_dir)
        download_pfam_data(pfam_data_dir)
        # _report_version(pfam_data_dir)
        return

    if check_if_data_present(pfam_data_dir):
        return

    print(color_text("\n    Downloading required Pfam data (only needs to be done once)...\n", "yellow"))
    download_pfam_data(pfam_data_dir)
    # _report_version(pfam_data_dir)


def _report_version(location):
    version = get_stored_pfam_version(location)
    if version:
        print(color_text(f"\n    Pfam version: {version}\n", "yellow"))


if __name__ == "__main__":
    main()
