#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
to download and set up the slim NCBI assembly-info table if it is not present.

Default behavior: download the pre-built, pre-combined slim table (packaged as a
single .tar.gz and rebuilt weekly by bit's refresh GitHub Action, re-hosted at a
fixed release asset). This is fast and is what almost everyone gets.

`-f/--force-update` re-pulls that same slim tarball (it does NOT rebuild directly
from NCBI). A full rebuild from NCBI's GenBank + RefSeq summaries only happens as
a fallback, if the hosted tarball can't be retrieved or extracted.

There is no longer a 4-week freshness check: if the table is present it is left
alone unless `-f` is given.

Ex. usage: get_ncbi_assembly_tables.py
"""

import sys
import os
import socket
import shutil
import tarfile
import time
import urllib
import urllib.error
import argparse
from datetime import date

from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm
from gtotree.utils.ncbi.slim_ncbi_assembly_summary import (
    build_slim_assembly_summary,
    write_date_retrieved,
    GENBANK_URL,
    REFSEQ_URL,
    TABLE_FILENAME,
    DATE_FILENAME,
)


# pre-slimmed, pre-combined NCBI assembly-info table + date-retrieved.txt,
# packaged as a single .tar.gz and rebuilt weekly by bit's
# refresh-ncbi-assembly-info GitHub Action, then re-hosted at this fixed asset
# URL. The archive holds exactly two files at its root:
#   ncbi-assembly-info.tsv   (slim, clean-header; see slim_ncbi_assembly_summary)
#   date-retrieved.txt       (build date, YYYY,MM,DD)
# GToTree depends on the bit-hosted asset (identical layout to what GToTree needs).
NCBI_ASSEMBLY_TARBALL_URL = "https://github.com/AstrobioMike/bit/releases/download/ncbi-assembly-info-latest/ncbi-assembly-info.tar.gz"


################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="This is a helper program to download and set up the slim "
                    "NCBI assembly-info table if it is not present.",
        epilog="Ex. usage: get_ncbi_assembly_tables.py")

    parser.add_argument("-f", "--force-update",
                        help="Re-download the slim assembly-info table even if it "
                             "is already present (pulls the hosted slim tarball "
                             "again).",
                        action="store_true")

    args = parser.parse_args()

    get_ncbi_assembly_data(force_update=args.force_update)

################################################################################


def get_ncbi_assembly_summary_tab():
    """
    lazily resolve the path to the slim assembly-info table. Resolved on demand
    (not at import time) so importing this module doesn't require the
    NCBI_assembly_data_dir env variable to be set.
    """
    return os.path.join(check_ncbi_assembly_info_location_var_is_set(),
                        TABLE_FILENAME)


def check_ncbi_assembly_info_location_var_is_set():

    # making sure there is a NCBI_assembly_data_dir env variable
    try:
        ncbi_assembly_data_dir = os.environ['NCBI_assembly_data_dir']
    except KeyError:
        wprint(color_text("The environment variable 'NCBI_assembly_data_dir' does not seem to be set :(", "yellow"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        print("")
        sys.exit(0)

    return ncbi_assembly_data_dir


def check_if_data_present(location):
    """
    True if both the slim table and date-retrieved.txt are present and non-empty.
    If either is missing/empty, any stray copy is cleaned up and we return False
    so a fresh copy is pulled.
    """
    table_path = os.path.join(str(location), TABLE_FILENAME)
    date_retrieved_path = os.path.join(str(location), DATE_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(table_path) or not is_nonempty_file(date_retrieved_path):
        for p in (table_path, date_retrieved_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False

    return True


def _download_with_retries(url, label, dest, attempts=4, retry_wait=3):
    """
    download `url` to `dest`, retrying up to `attempts` times on transient
    failures (timeouts, connection resets, transient errors), with a short wait
    between tries. A 404 is raised immediately (not retried). Raises the last
    error if all attempts fail.
    """
    last_err = None
    for attempt in range(1, attempts + 1):
        try:
            download_with_tqdm(url, label, dest)
            return
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise
            last_err = err
        except (urllib.error.URLError, socket.timeout, TimeoutError,
                ConnectionError, OSError) as err:
            last_err = err

        if attempt < attempts:
            wprint(color_text(
                f"    download failed (attempt {attempt}/{attempts}); retrying...",
                "yellow"))
            time.sleep(retry_wait)

    raise last_err


def get_slim_ncbi_assembly_data(location):
    """
    default path: download the pre-built slim NCBI assembly-info tarball and
    extract its two files (the slim table and date-retrieved.txt) into
    `location`. Falls back to rebuilding directly from NCBI
    (download_ncbi_assembly_summary_data) if no tarball URL is configured or the
    download/extract fails, so the user always ends up with a usable table.
    """
    table_path = os.path.join(location, TABLE_FILENAME)
    date_path = os.path.join(location, DATE_FILENAME)

    if not NCBI_ASSEMBLY_TARBALL_URL:
        download_ncbi_assembly_summary_data(location)
        return

    print(color_text("\n    Downloading the prepared NCBI assembly-info table (only needs to be done once)...\n", "yellow"))

    tarball_path = os.path.join(location, "ncbi-assembly-info.tar.gz")
    expected = {TABLE_FILENAME, DATE_FILENAME}

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        _download_with_retries(
            NCBI_ASSEMBLY_TARBALL_URL, "        NCBI prepared data", tarball_path)

        # a truncated/corrupt download trips here (full-stream read via
        # getmembers), triggering the fallback below rather than writing a
        # corrupt table.
        with tarfile.open(tarball_path, "r:gz") as tar:
            members = tar.getmembers()
            # take only the two expected files, matched by basename, guarding
            # against path traversal / nested dirs by extracting to flat names.
            wanted = {}
            for m in members:
                base = os.path.basename(m.name)
                if base in expected and m.isfile():
                    wanted[base] = m
            missing = expected - set(wanted)
            if missing:
                raise ValueError(
                    f"prepared NCBI archive is missing expected file(s): {', '.join(sorted(missing))}")

            for base, m in wanted.items():
                src = tar.extractfile(m)
                if src is None:
                    raise ValueError(f"could not read '{base}' from prepared NCBI archive")
                dest = os.path.join(location, base)
                with src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            tarfile.TarError, ValueError, OSError) as err:
        # clean up partial artifacts, then fall back to rebuilding from NCBI
        for p in (tarball_path, table_path, date_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        print("")
        wprint(color_text("  Couldn't get the prepared NCBI assembly-info table; "
                          "rebuilding directly from NCBI instead.", "yellow"))
        report_message(f"Underlying issue: {err}", ii="    ",
                       si="    ")
        print("")
        download_ncbi_assembly_summary_data(location)
        return
    finally:
        socket.setdefaulttimeout(default_timeout)
        if os.path.exists(tarball_path):
            try:
                os.remove(tarball_path)
            except OSError:
                pass
    print("")


def download_ncbi_assembly_summary_data(location):
    """
    rebuild fallback (used only if the hosted slim tarball can't be retrieved):
    download the GenBank + RefSeq assembly summaries directly from NCBI, combine
    and slim them to the columns GToTree uses (identical layout to the hosted
    asset), and stamp date-retrieved.txt with today's date.
    """
    genbank_temp = os.path.join(str(location), "assembly_summary_genbank.txt")
    refseq_temp = os.path.join(str(location), "assembly_summary_refseq.txt")
    table_path = os.path.join(str(location), TABLE_FILENAME)
    date_path = os.path.join(str(location), DATE_FILENAME)

    print(color_text("\n    Downloading NCBI assembly summaries (only needs to be done once)...\n", "yellow"))

    try:
        download_with_tqdm(GENBANK_URL, "        Genbank assemblies summary", genbank_temp)
        download_with_tqdm(REFSEQ_URL, "        RefSeq assemblies summary", refseq_temp)
        print("")
    except Exception as e:
        report_message(f"Downloading the NCBI assembly summary tables failed with the following error:\n{e}", "red")
        report_early_exit(None, copy_log=False)
        return

    # combine + slim to the columns GToTree uses (clean-header layout)
    try:
        build_slim_assembly_summary(genbank_temp, refseq_temp, table_path)
    except Exception as e:
        report_message(f"Combining/slimming the NCBI assembly summary tables failed with the following error:\n{e}", "red")
        report_early_exit(None, copy_log=False)
        return
    finally:
        for p in (genbank_temp, refseq_temp):
            if os.path.exists(p):
                os.remove(p)

    # storing date retrieved (YYYY,MM,DD), matching the asset's build stamp
    write_date_retrieved(date_path, when=date.today())


def get_ncbi_assembly_data(force_update=False):
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    data_present = check_if_data_present(ncbi_dir)

    if data_present and not force_update:
        return

    get_slim_ncbi_assembly_data(ncbi_dir)

################################################################################

if __name__ == "__main__":
    main()
