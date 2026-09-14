#!/usr/bin/env python

"""
The one implementation behind GToTree's hosted Parquet reference assets.

Everything source-specific lives in a HostedParquetAsset instance:

    env_var           the variable naming the local directory
    release_base      the bit release the two files come from
    parquet_filename  the table
    sidecar_filename  the small stamp alongside it, whatever shape it takes
    sidecar_validator raises ValueError if the downloaded sidecar is malformed
    display_name      how the asset is named to the user in messages

`gtotree/utils/ncbi/get_ncbi_assembly_data.py` and `gtotree/utils/gtdb/get_gtdb_data.py`
are thin wrappers over an instance each. They keep their own public function names
because a lot of call sites and tests reach for those directly.

This mirrors bit/modules/hosted_parquet_asset.py
"""

import os
import socket
import sys
import urllib
import urllib.error

from gtotree.utils.misc.messaging import (wprint, color_text, report_message,
                                          report_early_exit)
from gtotree.utils.misc.general import download_with_tqdm


#: how long to let a stalled connection hang before giving up, in seconds
DOWNLOAD_TIMEOUT = 30


def _is_nonempty_file(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def validate_date_stamp(path):
    """A 'YYYY,MM,DD' single-line stamp, as published beside the NCBI table."""
    with open(path) as fh:
        first = fh.readline().strip()
    parts = first.split(",")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise ValueError(f"date-retrieved.txt is not a 'YYYY,MM,DD' stamp: {first!r}")


def validate_version_lines(path):
    """A version line and a release-date line, as published beside the GTDB table."""
    with open(path) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    if len(lines) < 2:
        raise ValueError("GTDB version file doesn't have the expected version + date lines")


class HostedParquetAsset:
    """
    One hosted Parquet table plus its sidecar stamp, and everything done to it.

    The methods map onto what the per-source wrappers expose: location() resolves the
    env var, table_path() names the local table, is_present() decides whether a fetch
    is needed, and download() fetches both files or exits with something readable.
    """

    def __init__(self, *, env_var, release_base, parquet_filename, sidecar_filename,
                 sidecar_validator, display_name, download_label, sidecar_label):
        self.env_var = env_var
        self.release_base = release_base
        self.parquet_filename = parquet_filename
        self.sidecar_filename = sidecar_filename
        self.sidecar_validator = sidecar_validator
        self.display_name = display_name
        self.download_label = download_label
        self.sidecar_label = sidecar_label

    @property
    def data_url(self):
        return f"{self.release_base}/{self.parquet_filename}"

    @property
    def sidecar_url(self):
        return f"{self.release_base}/{self.sidecar_filename}"

    # -- locating it -------------------------------------------------------

    def location(self):
        """
        The local directory for this asset, from its environment variable.

        The variable is set by the conda activate script, so an unset one means the
        install is in a state the user can't be expected to reason about.
        """
        try:
            return os.environ[self.env_var]
        except KeyError:
            wprint(color_text(f"The environment variable '{self.env_var}' does not "
                              "seem to be set :(", "yellow"))
            wprint("This shouldn't happen, check on things with "
                   "`gtt data locations check`.")
            print("")
            sys.exit(1)

    def table_path(self, location=None):
        if location is None:
            location = self.location()
        return os.path.join(str(location), self.parquet_filename)

    def sidecar_path(self, location=None):
        if location is None:
            location = self.location()
        return os.path.join(str(location), self.sidecar_filename)

    # -- deciding whether to fetch -----------------------------------------

    def is_present(self, location):
        """
        True when both files are there and non-empty.

        A half-present pair is actively harmful: the table alone will be read and
        trusted. So anything left over from an interrupted fetch is cleared out here
        and False returned, which sends us down the download path.
        """
        table_path = self.table_path(location)
        sidecar_path = self.sidecar_path(location)

        if not _is_nonempty_file(table_path) or not _is_nonempty_file(sidecar_path):
            for p in (table_path, sidecar_path):
                if os.path.exists(p) and os.path.isfile(p):
                    os.remove(p)
            return False

        return True

    # -- fetching it -------------------------------------------------------

    def download(self, location):
        """
        Pull both files into `location`.

        The Parquet footer is checked and the sidecar validated before either is
        trusted, and both are removed on any failure. There is no rebuild fallback,
        so a half-written pair here would be a table that silently reads as empty on
        the next run.
        """
        table_path = self.table_path(location)
        sidecar_path = self.sidecar_path(location)

        print(color_text(f"\n    Downloading the prepared {self.display_name} "
                         "(only needs to be done once)...\n", "yellow"))

        default_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(DOWNLOAD_TIMEOUT)
        try:
            download_with_tqdm(self.data_url, f"        {self.download_label}",
                               table_path, speed_gate=True)

            # catches a truncated/corrupt download without reading the whole table
            self.verify_parquet(table_path)

            self.download_sidecar(sidecar_path)

        except (urllib.error.URLError, TimeoutError, ConnectionError, ValueError,
                OSError) as err:
            for p in (table_path, sidecar_path):
                if os.path.exists(p):
                    try:
                        os.remove(p)
                    except OSError:
                        pass
            self.report_unavailable(err)
        finally:
            socket.setdefaulttimeout(default_timeout)

        print("")

    def verify_parquet(self, path):
        """
        Cheap integrity check: open the Parquet footer and confirm the file has a
        schema and at least one row group. Reads only the footer, not the table.
        """
        import pyarrow.parquet as pq  # type: ignore
        md = pq.ParquetFile(path).metadata
        if md.num_columns == 0 or md.num_row_groups == 0:
            raise ValueError(f"downloaded {self.display_name} has no data "
                             "(truncated download?)")

    def download_sidecar(self, sidecar_path):
        """Atomic: write to .part, validate, then os.replace() into place."""
        tmp = sidecar_path + ".part"
        try:
            download_with_tqdm(self.sidecar_url, f"        {self.sidecar_label}", tmp,
                               leave=False)
            self.sidecar_validator(tmp)
            os.replace(tmp, sidecar_path)
        finally:
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass

    def report_unavailable(self, err):
        """Explain the failure, name the manual route, and exit."""
        print("")
        wprint(color_text(f"  Couldn't download the prepared {self.display_name} :(",
                          "yellow"))
        report_message(f"Underlying issue: {err}", color=None, ii="    ", si="    ")
        print("")
        report_message("This is usually a transient network problem, and trying again "
                       "in a few minutes often works. If it persists, the table can be "
                       "fetched manually from:", color=None, ii="    ", si="    ")
        print(f"        {color_text(self.data_url)}")
        print(f"        {color_text(self.sidecar_url)}")
        report_message(f"and placed (as '{self.parquet_filename}' and "
                       f"'{self.sidecar_filename}') in the directory shown by "
                       "`gtt data locations check`.", color=None, ii="    ", si="    ")
        print("")
        report_early_exit(None, copy_log=False)
