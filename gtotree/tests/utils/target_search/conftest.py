"""
Fixtures for the `gtt search-pfams` / `gtt search-kos` tests.

The Pfam fixtures reuse the mock profile set the gen-scg-hmms tests already run
against, staged into a directory laid out like a real managed `Pfam_data_dir` so the
production `get_additional_pfam_targets` runs unmodified against it.
"""

import os
import shutil
import dataclasses

import pytest  # type: ignore

from gtotree.tests.paths import MOCK_PFAM_HMM, MOCK_PFAM_INFO
from gtotree.utils.pfam.get_pfam_data import (HMM_FILENAME, INFO_FILENAME,
                                              VERSION_FILENAME)
from gtotree.utils.target_search.target_search_spec import get_spec


# each mock profile's exact matching peptide, from the gen-scg-hmms fixtures
MOTIFS = {
    "PF90001": "MKVLAAAL",
    "PF90002": "MARTKQTA",
    "PF90003": "MSDKIIHL",
    "PF90004": "MAHHWWGS",
}

MOCK_PFAM_VERSION = "mock-1.0"


@pytest.fixture
def mock_pfam_data_dir(tmp_path, monkeypatch):
    """
    A directory shaped like a real Pfam_data_dir, with `Pfam_data_dir` pointed at it.

    monkeypatch.setenv rather than a plain assignment so the real variable (if the
    developer running the tests has one) is restored afterwards.
    """
    pfam_dir = tmp_path / "pfam-data"
    pfam_dir.mkdir()
    shutil.copy(MOCK_PFAM_HMM, pfam_dir / HMM_FILENAME)
    shutil.copy(MOCK_PFAM_INFO, pfam_dir / INFO_FILENAME)
    (pfam_dir / VERSION_FILENAME).write_text(MOCK_PFAM_VERSION + "\n")

    monkeypatch.setenv("Pfam_data_dir", str(pfam_dir))
    return pfam_dir


@pytest.fixture
def pfam_spec(mock_pfam_data_dir):
    """
    The real Pfam spec with only the managed-dataset download stubbed out.

    Everything else -- target collection, the search worker, the counts writer, the
    hit-combining -- is the production code path. `dataclasses.replace` because the
    spec is frozen, so this can't accidentally leak a mutation into another test.
    """
    return dataclasses.replace(
        get_spec("pfam"),
        ensure_data=lambda **kwargs: str(mock_pfam_data_dir),
    )


@pytest.fixture
def write_genome(tmp_path):
    """
    Returns a callable making a one-genome amino-acid fasta out of Pfam IDs.

    Repeating an ID gives that genome a second copy of the profile, which is how the
    multi-copy cases below are built.
    """
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()

    def _write(name, pfam_ids):
        path = genome_dir / f"{name}.faa"
        with open(path, "w") as f:
            for i, pfam_id in enumerate(pfam_ids, 1):
                f.write(f">p{i}\n{MOTIFS[pfam_id]}\n")
        return path

    return _write


@pytest.fixture
def listing(tmp_path):
    """Returns a callable writing a single-column listing file and returning its path."""

    def _listing(name, entries):
        path = tmp_path / name
        path.write_text("\n".join(str(e) for e in entries) + "\n")
        return str(path)

    return _listing


@pytest.fixture
def in_tmp_cwd(tmp_path, monkeypatch):
    """Run inside tmp_path so relative output directories don't litter the repo."""
    monkeypatch.chdir(tmp_path)
    return tmp_path
