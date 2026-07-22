#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
for setting up reference files for the glorious Genome Taxonomy Database
(gtdb.ecogenomic.org/).

It ensures the prepared GTDB metadata table is present, downloading bit's hosted
Parquet asset (gtdb-data.parquet) if it isn't. The Parquet asset is already slimmed
to the columns GToTree uses and already has the 7 taxonomic ranks split into their
own columns
"""

import sys
import os
import socket
import urllib
import urllib.error

from gtotree.utils.messaging import (wprint, color_text, report_message,
                                     report_early_exit)
from gtotree.utils.general import download_with_tqdm


# bit hosts the prepared GTDB assets on a rolling "-latest" GitHub release. GToTree
# consumes them directly; it does not build its own. The Parquet file is the slimmed,
# rank-split metadata table; VERSION.txt carries the GTDB release + date lines.
_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/gtdb-metadata-latest"

PARQUET_FILENAME = "gtdb-data.parquet"
VERSION_FILENAME = "VERSION.txt"

GTDB_DATA_URL = f"{_RELEASE_BASE}/{PARQUET_FILENAME}"
GTDB_VERSION_URL = f"{_RELEASE_BASE}/{VERSION_FILENAME}"


################################################################################


def get_gtdb_data(force_update=False):
    """
    Ensure the GTDB Parquet table is present locally, and return the GTDB data dir.
    """
    gtdb_dir = check_gtdb_location_var_is_set()
    data_present = check_if_gtdb_data_present(gtdb_dir)

    if data_present and not force_update:
        return gtdb_dir

    get_slim_gtdb_tab(gtdb_dir)
    return gtdb_dir


def check_gtdb_location_var_is_set():
    try:
        gtdb_data_dir = os.environ['GTDB_DIR']
    except KeyError:
        wprint(color_text("The environment variable 'GTDB_DIR' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt data locations check`.")
        sys.exit(1)
    return gtdb_data_dir


def gtdb_data_table_path(location=None):
    """Path to the local GTDB Parquet asset (resolving the location if not given)."""
    if location is None:
        location = check_gtdb_location_var_is_set()
    return os.path.join(str(location), PARQUET_FILENAME)


def check_if_gtdb_data_present(location):
    """
    True if both the Parquet table and version-info file are present and non-empty.
    If either is missing/empty, any stray copy is cleaned up and we return False so a
    fresh copy is pulled.
    """
    table_path = os.path.join(str(location), PARQUET_FILENAME)
    version_info_path = os.path.join(str(location), VERSION_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(table_path) or not is_nonempty_file(version_info_path):
        for p in (table_path, version_info_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False

    return True


def report_gtdb_unreachable(err):

    print("")
    wprint(color_text("  Couldn't download the prepared GTDB table :(", "red"))
    report_message(f"Underlying issue: {err}", color=None, ii="    ", si="    ")
    print("")
    report_message("This is usually a transient network problem, and trying again in a few "
                   "minutes often works. If it persists, the table can be fetched manually "
                   "from:", color=None, ii="    ", si="    ")
    print(f"        {color_text(GTDB_DATA_URL)}")
    print(f"        {color_text(GTDB_VERSION_URL)}")
    report_message(f"and placed (as '{PARQUET_FILENAME}' and '{VERSION_FILENAME}') in the "
                   "directory shown by `gtt data locations check`.", color=None,
                   ii="    ", si="    ")
    print("")
    report_early_exit(None, copy_log=False)


def get_slim_gtdb_tab(location):
    """
    Download bit's prepared GTDB Parquet asset and its version-info file into
    `location`. The Parquet footer is verified before we trust the table, and the
    version file is written atomically. On any network/integrity failure the partial
    artifacts are cleaned up and we exit with a helpful message
    """
    table_path = os.path.join(location, PARQUET_FILENAME)
    version_path = os.path.join(location, VERSION_FILENAME)

    print(color_text("\n    Downloading the prepared GTDB table (only needs to be done once)...\n", "yellow"))

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        download_with_tqdm(GTDB_DATA_URL, "        GTDB prepared data", table_path,
                           speed_gate=True)

        # confirm the file is a readable Parquet before we trust it (catches a
        # truncated/corrupt download without reading the whole table)
        _verify_parquet(table_path)

        _download_version_file(version_path)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            ValueError, OSError) as err:
        for p in (table_path, version_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        report_gtdb_unreachable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)

    print("")


def _verify_parquet(path):
    """
    Cheap integrity check: open the Parquet footer and confirm the file has a schema
    and at least one row group. Reads only the footer, not the whole table.
    """
    import pyarrow.parquet as pq # type: ignore
    md = pq.ParquetFile(path).metadata
    if md.num_columns == 0 or md.num_row_groups == 0:
        raise ValueError("downloaded GTDB table has no data (truncated download?)")


def _download_version_file(version_path):
    """Download VERSION.txt atomically: write to .part, validate, then os.replace()."""
    tmp = version_path + ".part"
    try:
        download_with_tqdm(GTDB_VERSION_URL, "        version info", tmp, leave=False)
        _validate_version_file(tmp)
        os.replace(tmp, version_path)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


def _validate_version_file(path):
    lines = [ln.strip() for ln in open(path) if ln.strip()]
    if len(lines) < 2:
        raise ValueError("GTDB version file doesn't have the expected version + date lines")


def report_gtdb_version_info(location):
    """Return (version, release_date) from the local VERSION.txt (first two lines)."""
    version_info = []
    with open(os.path.join(location, VERSION_FILENAME)) as version_info_file:
        for line in version_info_file:
            line = line.strip()
            if line != "":
                version_info.append(line)
    return version_info[0], version_info[1]

