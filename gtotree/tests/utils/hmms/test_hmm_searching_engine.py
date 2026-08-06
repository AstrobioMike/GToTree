"""Unit tests for gtotree/utils/hmms/hmm_searching_engine.py."""

import os
import shutil
import subprocess

import pyhmmer.errors  # type: ignore
import pytest  # type: ignore

from gtotree.tests.paths import DATA_DIR, MOCK_HMM
from gtotree.utils.hmms.hmm_searching_engine import (
    press_profiles,
    profiles_missing_gathering_cutoffs,
    search_one_genome,
)

GENOME = DATA_DIR / "mock-1.faa"


def _rows(path):
    """Data rows of a tblout, as field lists (comments and blanks dropped)."""
    return [line.split()[:6] for line in open(path)
            if line.strip() and not line.startswith("#")]


class TestGatheringCutoffs:

    def test_packaged_hmms_have_gathering_thresholds(self):
        assert profiles_missing_gathering_cutoffs(str(MOCK_HMM)) == []

    @pytest.mark.skipif(shutil.which("hmmbuild") is None, reason="hmmbuild not on PATH")
    def test_profile_without_ga_is_reported(self, tmp_path):
        """Preflight uses this to reject a GA-less HMM once, before any work starts."""
        hmm = tmp_path / "noga.hmm"
        subprocess.run(["hmmbuild", "--amino", "-n", "NOGA", str(hmm), str(GENOME)],
                       capture_output=True, check=True)

        assert profiles_missing_gathering_cutoffs(str(hmm)) == ["NOGA"]

    @pytest.mark.skipif(shutil.which("hmmbuild") is None, reason="hmmbuild not on PATH")
    def test_searching_without_ga_raises_rather_than_falling_back(self, tmp_path):
        """Falling back to an E-value cutoff would change which genes are called."""
        hmm = tmp_path / "noga.hmm"
        subprocess.run(["hmmbuild", "--amino", "-n", "NOGA", str(hmm), str(GENOME)],
                       capture_output=True, check=True)

        with pytest.raises(pyhmmer.errors.MissingCutoffs):
            search_one_genome(str(hmm), str(GENOME), str(tmp_path / "out.tbl"))


class TestPressProfiles:

    def test_pressing_produces_a_usable_base(self, tmp_path):
        base = press_profiles(str(MOCK_HMM), str(tmp_path), "pressed")
        assert base is not None
        # hmmpress emits a set of sidecar files alongside the base path
        produced = [f for f in os.listdir(tmp_path) if f.startswith("pressed")]
        assert produced, "hmmpress produced no files"

    def test_pressing_a_bad_path_returns_none_rather_than_raising(self, tmp_path):
        """Pressing is an optimization: failing it degrades to reading the plain HMM."""
        assert press_profiles(str(tmp_path / "nope.hmm"), str(tmp_path), "x") is None

    def test_pressed_and_unpressed_searches_agree(self, tmp_path):
        """
        Both paths are live -- `_open_profiles` falls back to the plain HMM when no
        pressed set is available -- so they have to keep producing the same hits.
        """
        plain = tmp_path / "plain.tbl"
        pressed_out = tmp_path / "pressed.tbl"

        search_one_genome(str(MOCK_HMM), str(GENOME), str(plain))
        base = press_profiles(str(MOCK_HMM), str(tmp_path), "pressed")
        search_one_genome(str(MOCK_HMM), str(GENOME), str(pressed_out),
                          pressed_base=base)

        assert _rows(str(plain)) == _rows(str(pressed_out))

    def test_falls_back_when_pressed_base_is_missing(self, tmp_path):
        """A bogus pressed base must fall through to the plain HMM, still returning hits."""
        out = tmp_path / "out.tbl"
        search_one_genome(str(MOCK_HMM), str(GENOME), str(out),
                          pressed_base=str(tmp_path / "does-not-exist"))
        assert _rows(str(out)), "fallback produced no hits"


class TestSearchOutput:

    def test_writes_a_tblout_and_returns_hits(self, tmp_path):
        out = tmp_path / "out.tbl"
        results = search_one_genome(str(MOCK_HMM), str(GENOME), str(out))

        assert out.exists()
        assert results, "no TopHits returned"
        assert _rows(str(out)), "tblout had no data rows"

    def test_output_is_tblout_format(self, tmp_path):
        """
        Comment lines, then one row per hit whose leading fields are target name,
        target accession, query name, query accession, E-value, bit score. This is what
        `read_hmmer_results` parses and what users read in the per-genome Pfam results.
        """
        out = tmp_path / "out.tbl"
        search_one_genome(str(MOCK_HMM), str(GENOME), str(out))

        text = out.read_text()
        assert text.startswith("#"), "no comment header"

        rows = _rows(str(out))
        assert len(rows) == 1, "expected exactly one hit from the mock profile"
        target, target_acc, query, query_acc, evalue, score = rows[0]
        assert query == "MOCK"             # the profile name
        assert query_acc == "Mike.1"       # the profile accession
        assert target_acc == "-"           # sequences carry no accession
        float(evalue), float(score)        # numeric, in that order

    def test_genome_with_no_hits_yields_an_empty_table(self, tmp_path):
        """A genome with no hits still produces a parseable, header-only file."""
        empty = tmp_path / "nothing.faa"
        empty.write_text(">gene1\nWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n")

        out = tmp_path / "out.tbl"
        search_one_genome(str(MOCK_HMM), str(empty), str(out))

        assert out.exists()
        assert _rows(str(out)) == []
