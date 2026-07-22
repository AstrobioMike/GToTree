#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
for setting up the NCBI assembly-info table if it is not present.

GToTree uses bit's hosted, pre-built NCBI Parquet asset (ncbi-data.parquet),
which is a combined, slimmed GenBank + RefSeq assembly summary with taxonomy resolved
"""

import sys
import os
import socket
import urllib
import urllib.error

from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm


_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/ncbi-assembly-info-latest"

PARQUET_FILENAME = "ncbi-data.parquet"
DATE_FILENAME = "date-retrieved.txt"

NCBI_DATA_URL = f"{_RELEASE_BASE}/{PARQUET_FILENAME}"
NCBI_DATE_URL = f"{_RELEASE_BASE}/{DATE_FILENAME}"


def get_ncbi_assembly_data(force_update=False):
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    data_present = check_if_data_present(ncbi_dir)

    if data_present and not force_update:
        return

    get_slim_ncbi_assembly_data(ncbi_dir)


def check_ncbi_assembly_info_location_var_is_set():

    # making sure there is a NCBI_ASSEMBLY_DATA_DIR env variable
    try:
        ncbi_assembly_data_dir = os.environ['NCBI_ASSEMBLY_DATA_DIR']
    except KeyError:
        wprint(color_text("The environment variable 'NCBI_ASSEMBLY_DATA_DIR' does not seem to be set :(", "yellow"))
        wprint("This shouldn't happen, check on things with `gtt data locations check`.")
        print("")
        sys.exit(0)

    return ncbi_assembly_data_dir


def ncbi_data_table_path(location=None):
    """Path to the local NCBI Parquet asset (resolving the location if not given)."""
    if location is None:
        location = check_ncbi_assembly_info_location_var_is_set()
    return os.path.join(str(location), PARQUET_FILENAME)


def get_ncbi_assembly_summary_tab():
    """
    lazily resolve the path to the NCBI Parquet assembly-info table. Resolved on demand
    (not at import time) so importing this module doesn't require the
    NCBI_ASSEMBLY_DATA_DIR env variable to be set.
    """
    return ncbi_data_table_path()


def check_if_data_present(location):
    """
    True if both the Parquet table and date-retrieved.txt are present and non-empty.
    If either is missing/empty, any stray copy is cleaned up and we return False so a
    fresh copy is pulled.
    """
    table_path = os.path.join(str(location), PARQUET_FILENAME)
    date_retrieved_path = os.path.join(str(location), DATE_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(table_path) or not is_nonempty_file(date_retrieved_path):
        for p in (table_path, date_retrieved_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False

    return True


def report_ncbi_unavailable(err):

    print("")
    wprint(color_text("  Couldn't download the prepared NCBI assembly-info table :(", "yellow"))
    report_message(f"Underlying issue: {err}", color=None, ii="    ", si="    ")
    print("")
    report_message("This is usually a transient network problem, and trying again in a few "
                   "minutes often works. If it persists, the table can be fetched manually "
                   "from:", color=None, ii="    ", si="    ")
    print(f"        {color_text(NCBI_DATA_URL)}")
    print(f"        {color_text(NCBI_DATE_URL)}")
    report_message(f"and placed (as '{PARQUET_FILENAME}' and '{DATE_FILENAME}') in the "
                   "directory shown by `gtt data locations check`.", color=None,
                   ii="    ", si="    ")
    print("")
    report_early_exit(None, copy_log=False)


def get_slim_ncbi_assembly_data(location):
    """
    Download bit's prepared NCBI Parquet asset and its date-retrieved file into
    `location`. The Parquet footer is verified before we trust the table, and the date
    file is written atomically. On any network/integrity failure the partial artifacts
    are cleaned up and we exit with a helpful message -- there is no NCBI-rebuild
    fallback, since the hosted asset is the prepared table.
    """
    table_path = os.path.join(location, PARQUET_FILENAME)
    date_path = os.path.join(location, DATE_FILENAME)

    print(color_text("\n    Downloading the prepared NCBI assembly-info table (only needs to be done once)...\n", "yellow"))

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        download_with_tqdm(NCBI_DATA_URL, "        NCBI prepared data", table_path,
                           speed_gate=True)

        # confirm the file is a readable Parquet before we trust it (catches a
        # truncated/corrupt download without reading the whole table)
        _verify_parquet(table_path)

        _download_date_file(date_path)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            ValueError, OSError) as err:
        for p in (table_path, date_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        report_ncbi_unavailable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)

    print("")


def _verify_parquet(path):
    import pyarrow.parquet as pq # type: ignore
    md = pq.ParquetFile(path).metadata
    if md.num_columns == 0 or md.num_row_groups == 0:
        raise ValueError("downloaded NCBI table has no data (truncated download?)")


def _download_date_file(date_path):
    tmp = date_path + ".part"
    try:
        download_with_tqdm(NCBI_DATE_URL, "        date stamp", tmp, leave=False)
        _validate_date_file(tmp)
        os.replace(tmp, date_path)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


def _validate_date_file(path):
    with open(path) as fh:
        first = fh.readline().strip()
    parts = first.split(",")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise ValueError(f"date-retrieved.txt is not a 'YYYY,MM,DD' stamp: {first!r}")


def read_date_retrieved(location):
    """
    Read date-retrieved.txt (a 'YYYY,MM,DD' stamp) from `location` and return it
    formatted like 'Jan 05, 2026'. Returns the raw string if it can't be parsed.
    """
    import datetime

    with open(os.path.join(location, DATE_FILENAME)) as fh:
        stamp = fh.readline().strip()
    try:
        y, m, d = (int(p) for p in stamp.split(","))
        return datetime.date(y, m, d).strftime("%b %d, %Y")
    except (ValueError, TypeError):
        return stamp

