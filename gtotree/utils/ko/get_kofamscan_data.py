#!/usr/bin/env python
"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
to download and setup the KOFamScan (https://github.com/takaram/kofam_scan) data files for use.

Default behavior: download the pre-built KOFamScan data (README, ko_list,
profiles/, date-retrieved.txt) packaged as a single .tar.gz and rebuilt monthly
by GToTree's refresh-kofamscan-data GitHub Action, re-hosted at a fixed release
asset. This is served over HTTPS, so it sidesteps the ftp.genome.jp FTP access
problems that block many institutional / cluster networks.

`-f/--force-update` re-pulls that same hosted tarball.
"""

import sys
import os
import argparse
import shutil
import tarfile
from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm


# pre-built KOFamScan data, packaged as a single .tar.gz and rebuilt monthly by
# GToTree's refresh-kofamscan-data GitHub Action, then re-hosted at this fixed
# asset URL. The archive holds exactly this at its root:
#   README
#   ko_list
#   profiles.tar.gz
#   date-retrieved.txt   (build date, YYYY,MM,DD)
# Served over HTTPS to avoid the ftp.genome.jp FTP access problems.
KOFAMSCAN_TARBALL_URL = "https://github.com/AstrobioMike/GToTree/releases/download/kofamscan-data-latest/kofamscan-data.tar.gz"

# files/dirs expected at the data dir root after a successful setup
README_FILENAME = "README"
KO_LIST_FILENAME = "ko_list"
PROFILES_DIRNAME = "profiles"
PROFILES_TARBALL = "profiles.tar.gz"
DATE_FILENAME = "date-retrieved.txt"


################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Setup the KOFamScan (https://github.com/takaram/kofam_scan) data files for use.",
        epilog="Example usage: gtt-get-kofamscan-data"
    )

    parser.add_argument("-f", "--force-update",
                        help="Re-download the KOFamScan data even if it is already "
                             "present",
                        action="store_true")

    args = parser.parse_args()

    get_kofamscan_data(force_update=args.force_update)

################################################################################


def check_location_var_is_set():
    """
    Ensure that the environment variable 'KO_data_dir' is set.
    Returns the directory specified by the environment variable.
    """
    try:
        ko_data_dir = os.environ['KO_data_dir']
    except KeyError:
        wprint(color_text("The environment variable 'KO_data_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return ko_data_dir


def check_if_data_present(location):
    """
    True only if all expected pieces are present. If anything is missing, remove
    whatever partial state exists so the subsequent download starts clean.
    """
    readme_path = os.path.join(location, README_FILENAME)
    ko_list_path = os.path.join(location, KO_LIST_FILENAME)
    profiles_dir_path = os.path.join(location, PROFILES_DIRNAME)

    if (os.path.isfile(readme_path) and
            os.path.isfile(ko_list_path) and
            os.path.isdir(profiles_dir_path)):
        return True

    # incomplete - clear partial state
    _clear_partial_state(location)
    return False


def _clear_partial_state(location):
    """Remove any of the expected files/dirs that happen to exist."""
    for fname in (README_FILENAME, KO_LIST_FILENAME, DATE_FILENAME, PROFILES_TARBALL):
        fpath = os.path.join(location, fname)
        if os.path.exists(fpath):
            os.remove(fpath)
    profiles_dir_path = os.path.join(location, PROFILES_DIRNAME)
    if os.path.isdir(profiles_dir_path):
        shutil.rmtree(profiles_dir_path)


def report_kofam_dl_failure(e):
    report_message(f"Downloading the KOFamScan data failed with the following error:\n{e}", "red")
    report_message("GToTree pulls this data from a GitHub-hosted release asset over HTTPS.")
    report_message("You can check whether you can access this URL:")
    report_message(f"    {KOFAMSCAN_TARBALL_URL}")
    report_message("If this keeps failing, it may be a transient network issue - trying again "
                   "later often resolves it. You can also drop the additional target KOs (being "
                   "provided to the `-K` flag) and run GToTree without them.")
    report_early_exit(None, copy_log=False)


def download_kofamscan_data(location):
    """
    Pull the hosted tarball and extract it atomically into `location`:
    download to a temp path, extract the outer archive into a temp staging dir,
    extract the inner profiles.tar.gz, validate, then move the pieces into
    place. Avoids leaving a half-populated data dir if anything is interrupted.
    """
    tarball_path = os.path.join(location, "kofamscan-data.tar.gz")
    staging_dir = os.path.join(location, ".kofamscan-staging")

    # clean any leftover staging from a prior interrupted run
    if os.path.isdir(staging_dir):
        shutil.rmtree(staging_dir)
    os.makedirs(staging_dir, exist_ok=True)

    try:
        download_with_tqdm(KOFAMSCAN_TARBALL_URL, "        KOFamScan data", tarball_path)
    except Exception as e:
        _safe_rmtree(staging_dir)
        _safe_remove(tarball_path)
        report_kofam_dl_failure(e)

    # extract the OUTER archive into staging
    try:
        with tarfile.open(tarball_path) as tarball:
            tarball.extractall(staging_dir)
    except Exception as e:
        _safe_rmtree(staging_dir)
        _safe_remove(tarball_path)
        report_kofam_dl_failure(e)

    os.remove(tarball_path)

    # extract the INNER profiles.tar.gz (shipped compressed) into staging,
    # producing profiles/. Remove the inner tarball once extracted so it
    # doesn't get moved into the final data dir.
    inner_profiles_tar = os.path.join(staging_dir, PROFILES_TARBALL)
    if not os.path.isfile(inner_profiles_tar):
        _safe_rmtree(staging_dir)
        report_kofam_dl_failure(
            f"the downloaded archive did not contain {PROFILES_TARBALL}")
    try:
        with tarfile.open(inner_profiles_tar) as profiles_tar:
            profiles_tar.extractall(staging_dir)
    except Exception as e:
        _safe_rmtree(staging_dir)
        report_kofam_dl_failure(f"failed to extract {PROFILES_TARBALL}: {e}")
    os.remove(inner_profiles_tar)

    # validate expected contents landed in staging
    staged_readme = os.path.join(staging_dir, README_FILENAME)
    staged_ko_list = os.path.join(staging_dir, KO_LIST_FILENAME)
    staged_profiles = os.path.join(staging_dir, PROFILES_DIRNAME)
    if not (os.path.isfile(staged_readme) and
            os.path.isfile(staged_ko_list) and
            os.path.isdir(staged_profiles)):
        _safe_rmtree(staging_dir)
        report_kofam_dl_failure(
            "the downloaded archive did not contain the expected "
            f"{README_FILENAME}, {KO_LIST_FILENAME}, and {PROFILES_DIRNAME}/ "
            "(after extracting profiles)")

    # move staged pieces into place (clearing any stale copies first)
    _clear_partial_state(location)
    for entry in os.listdir(staging_dir):
        shutil.move(os.path.join(staging_dir, entry), os.path.join(location, entry))
    shutil.rmtree(staging_dir)


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


def get_kofamscan_data(force_update=False):
    ko_data_dir = check_location_var_is_set()

    if force_update:
        print(color_text("    Re-downloading KO data (force-update requested)...\n", "yellow"))
        _clear_partial_state(ko_data_dir)
        download_kofamscan_data(ko_data_dir)
        return

    if check_if_data_present(ko_data_dir):
        return

    print(color_text("    Downloading required KO data (only needs to be done once, "
                     "unless newer HMMs become available)...\n", "yellow"))
    download_kofamscan_data(ko_data_dir)


if __name__ == "__main__":
    main()
