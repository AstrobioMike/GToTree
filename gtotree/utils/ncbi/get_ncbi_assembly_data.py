#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
for setting up the NCBI assembly-info table if it is not present.

GToTree uses bit's hosted, pre-built NCBI Parquet asset (ncbi-data.parquet),
which is a combined, slimmed GenBank + RefSeq assembly summary with taxonomy resolved

The download/verify/cleanup machinery is shared with the GTDB table, see
gtotree/utils/misc/hosted_parquet_asset.py. What lives here is the NCBI-specific
configuration and the names the rest of the codebase imports.
"""

import os

from gtotree.utils.misc.hosted_parquet_asset import (HostedParquetAsset,
                                                     validate_date_stamp)


_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/ncbi-assembly-info-latest"

PARQUET_FILENAME = "ncbi-data.parquet"
DATE_FILENAME = "date-retrieved.txt"

NCBI_ASSET = HostedParquetAsset(
    env_var="NCBI_ASSEMBLY_DATA_DIR",
    release_base=_RELEASE_BASE,
    parquet_filename=PARQUET_FILENAME,
    sidecar_filename=DATE_FILENAME,
    sidecar_validator=validate_date_stamp,
    display_name="NCBI assembly-info table",
    download_label="NCBI prepared data",
    sidecar_label="date stamp",
)

NCBI_DATA_URL = NCBI_ASSET.data_url
NCBI_DATE_URL = NCBI_ASSET.sidecar_url


def get_ncbi_assembly_data(force_update=False):
    """Ensure the NCBI Parquet table is present locally."""
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    data_present = check_if_data_present(ncbi_dir)

    if data_present and not force_update:
        return

    get_slim_ncbi_assembly_data(ncbi_dir)


def check_ncbi_assembly_info_location_var_is_set():
    return NCBI_ASSET.location()


def ncbi_data_table_path(location=None):
    """Path to the local NCBI Parquet asset (resolving the location if not given)."""
    return NCBI_ASSET.table_path(location)


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
    return NCBI_ASSET.is_present(location)


def get_slim_ncbi_assembly_data(location):
    """
    Download bit's prepared NCBI Parquet asset and its date-retrieved file into
    `location`. The Parquet footer is verified before we trust the table, and the date
    file is written atomically. On any network/integrity failure the partial artifacts
    are cleaned up and we exit with a helpful message -- there is no NCBI-rebuild
    fallback, since the hosted asset is the prepared table.
    """
    NCBI_ASSET.download(location)


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
