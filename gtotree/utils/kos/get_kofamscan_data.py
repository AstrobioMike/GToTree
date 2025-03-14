#!/usr/bin/env python
"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
to download and setup the KOFamScan (https://github.com/takaram/kofam_scan) data files for use.

For examples, please visit the GToTree wiki here: https://github.com/AstrobioMike/GToTree/wiki/example-usage
"""

import sys
import os
import argparse
import shutil
import filecmp
import tarfile
import gzip
from gtotree.utils.messaging import wprint, color_text, report_message, report_early_exit
from gtotree.utils.general import download_with_tqdm



################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Setup the KOFamScan (https://github.com/takaram/kofam_scan) data files for use.",
        epilog="Example usage: gtt-get-kofamscan-data"
    )
    get_kofamscan_data()

################################################################################


def check_location_var_is_set():
    """
    Ensure that the environment variable 'KO_data_dir' is set.
    Returns the directory specified by the environment variable.
    """
    try:
        ko_data_dir = os.environ['KO_data_dir']
    except KeyError:
        wprint(color_text("The environment variable 'KO_data_dir' does not seem to be set :(", "yellow"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return ko_data_dir


def check_stored_data_up_to_date(location):
    """ checks if the stored kofamscan data is the latest """

    latest_readme_path = os.path.join(location, "README-latest")
    remote_url = "ftp://ftp.genome.jp/pub/db/kofam/README"
    try:
        download_with_tqdm(remote_url, latest_readme_path, "        Latest KOFam README")
    except Exception as e:
        report_KOFam_dl_failure(e)

    local_readme_path = os.path.join(location, "README")
    if os.path.isfile(local_readme_path) and filecmp.cmp(latest_readme_path, local_readme_path):
        os.remove(latest_readme_path)
        return True
    else:
        os.remove(latest_readme_path)
        wprint(color_text("A newer version of the KOFamScan data is available, updating...", "yellow"))
        return False


def report_KOFam_dl_failure(e):
    report_message(f"Downloading KOFam data failed with the following error:\n{e}", "red")
    report_message(f"Unfortunately, the KOFam data can only be accessed via ftp, which may be blocked on your system.")
    report_message(f"You can check and see if you can access this URL: ftp://ftp.genome.jp/pub/db/kofam/")
    report_message(f"If accessing that is a problem, you can drop the additional target KOs (being provided to the `-K` flag) and sadly run GToTree without them.")
    report_early_exit()


def check_if_data_present(location):

    readme_path = os.path.join(location, "README")
    ko_list_path = os.path.join(location, "ko_list")
    hmms_dir_path = os.path.join(location, "profiles")

    # If any file/directory is missing, remove what exists and signal that update is needed.
    if (not os.path.isfile(readme_path) or
        not os.path.isfile(ko_list_path) or
        not os.path.isdir(hmms_dir_path)):
        if os.path.exists(readme_path):
            os.remove(readme_path)
        if os.path.exists(ko_list_path):
            os.remove(ko_list_path)
        if os.path.isdir(hmms_dir_path):
            shutil.rmtree(hmms_dir_path)
        return False
    else:
        if check_stored_data_up_to_date(location):
            return True
        else:
            if os.path.exists(readme_path):
                os.remove(readme_path)
            if os.path.exists(ko_list_path):
                os.remove(ko_list_path)
            if os.path.isdir(hmms_dir_path):
                shutil.rmtree(hmms_dir_path)
            return False


def download_kofamscan_data(location):

    readme_url = "ftp://ftp.genome.jp/pub/db/kofam/README"
    ko_list_gz_url = "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz"
    hmms_tar_url = "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz"

    readme_path = os.path.join(location, "README")
    ko_list_gz_path = os.path.join(location, "ko_list.gz")
    ko_list_path = os.path.join(location, "ko_list")
    hmms_tar_path = os.path.join(location, "profiles.tar.gz")

    try:
        download_with_tqdm(readme_url, readme_path, "        KOFam README")
    except Exception as e:
        report_KOFam_dl_failure(e)
    try:
        download_with_tqdm(ko_list_gz_url, ko_list_gz_path, "        KOFam KO list")
        download_with_tqdm(hmms_tar_url, hmms_tar_path, "        KOFam HMMs")
    except Exception as e:
        report_KOFam_dl_failure(e)

    with gzip.open(ko_list_gz_path, 'rb') as f_in:
        with open(ko_list_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(ko_list_gz_path)

    with tarfile.open(hmms_tar_path) as tarball:
        tarball.extractall(location)
    os.remove(hmms_tar_path)


def get_kofamscan_data():
    ko_data_dir = check_location_var_is_set()
    if check_if_data_present(ko_data_dir):
        return
    else:
        print(color_text("    Downloading required KO data (only needs to be done once, unless newer HMMs become available)...\n", "yellow"))
        download_kofamscan_data(ko_data_dir)


if __name__ == "__main__":
    main()
