#!/usr/bin/env python

"""
This is a helper program of GToTree (https://github.com/AstrobioMike/GToTree/wiki)
for setting up reference files for the glorious Genome Taxonomy Database
(gtdb.ecogenomic.org/).

It ensures the prepared GTDB metadata table is present, downloading bit's hosted
Parquet asset (gtdb-data.parquet) if it isn't. The Parquet asset is already slimmed
to the columns GToTree uses and already has the 7 taxonomic ranks split into their
own columns

The download/verify/cleanup machinery is shared with the NCBI table, see
gtotree/utils/misc/hosted_parquet_asset.py. What lives here is the GTDB-specific
configuration and the names the rest of the codebase imports.
"""

import os

from gtotree.utils.misc.hosted_parquet_asset import (HostedParquetAsset,
                                                     validate_version_lines)


# bit hosts the prepared GTDB assets on a rolling "-latest" GitHub release. GToTree
# consumes them directly; it does not build its own. The Parquet file is the slimmed,
# rank-split metadata table; VERSION.txt carries the GTDB release + date lines.
_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/gtdb-metadata-latest"

PARQUET_FILENAME = "gtdb-data.parquet"
VERSION_FILENAME = "VERSION.txt"

GTDB_ASSET = HostedParquetAsset(
    env_var="GTDB_DIR",
    release_base=_RELEASE_BASE,
    parquet_filename=PARQUET_FILENAME,
    sidecar_filename=VERSION_FILENAME,
    sidecar_validator=validate_version_lines,
    display_name="GTDB table",
    download_label="GTDB prepared data",
    sidecar_label="version info",
)

GTDB_DATA_URL = GTDB_ASSET.data_url
GTDB_VERSION_URL = GTDB_ASSET.sidecar_url


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
    return GTDB_ASSET.location()


def gtdb_data_table_path(location=None):
    """Path to the local GTDB Parquet asset (resolving the location if not given)."""
    return GTDB_ASSET.table_path(location)


def check_if_gtdb_data_present(location):
    """
    True if both the Parquet table and version-info file are present and non-empty.
    If either is missing/empty, any stray copy is cleaned up and we return False so a
    fresh copy is pulled.
    """
    return GTDB_ASSET.is_present(location)


def get_slim_gtdb_tab(location):
    """
    Download bit's prepared GTDB Parquet asset and its version-info file into
    `location`. The Parquet footer is verified before we trust the table, and the
    version file is written atomically. On any network/integrity failure the partial
    artifacts are cleaned up and we exit with a helpful message
    """
    GTDB_ASSET.download(location)


def report_gtdb_version_info(location):
    """Return (version, release_date) from the local VERSION.txt (first two lines)."""
    version_info = []
    with open(os.path.join(location, VERSION_FILENAME)) as version_info_file:
        for line in version_info_file:
            line = line.strip()
            if line != "":
                version_info.append(line.replace("Released ", ""))
    return version_info[0], version_info[1]
