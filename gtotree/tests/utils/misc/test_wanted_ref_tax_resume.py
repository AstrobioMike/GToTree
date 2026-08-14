"""
A resume reuses the `-w` accessions the previous run selected instead of re-resolving.

Re-resolving isn't only wasted work (a full taxonomy-table read, the dereplication, and
the NCBI liveness screen, on every resume). It's a hole in the resume contract: the
fingerprint covers `-w`, `--target-rank`, `--derep-rank` and `--source`, but not the
version of the GTDB/NCBI asset they resolve against. Update that asset between two
invocations and re-resolution returns a different accession set, which gets folded into
`ncbi_accs` after the fingerprint check has already passed. `process_genomes` picks up
work by per-genome flags rather than a stage guard, so the new genomes would quietly be
downloaded and searched into a run whose finished stages were counted against a
smaller set.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import RunData, GenomeData
from gtotree.utils.misc import preflight_checks
from gtotree.utils.misc.preflight_checks import (select_wanted_ref_tax,
                                            wanted_ref_tax_already_resolved,
                                            merge_wanted_ref_tax)


class _Args:
    def __init__(self, **kw):
        self.wanted_ref_tax = None
        self.source = "GTDB"
        self.target_rank = None
        self.derep_rank = "auto"
        self.__dict__.update(kw)


def _previous_run_data(accs=("GCF_1", "GCF_2"), user_accs=()):
    rd = RunData()
    for acc in user_accs:
        rd.ncbi_accs.append(GenomeData.from_acc(acc))
    for acc in accs:
        gd = GenomeData.from_acc(acc)
        gd.from_wanted_ref_tax = True
        rd.ncbi_accs.append(gd)
    rd.update_all_input_genomes()
    return rd


@pytest.fixture
def _never_resolves(monkeypatch):
    """Make any actual resolution attempt a loud failure rather than a slow one."""
    def _boom(*a, **kw):
        raise AssertionError("resolved `-w` when it should have reused the previous run")
    monkeypatch.setattr(preflight_checks, "resolve_wanted_ref_tax_accessions", _boom)


class TestReuseDetection:

    def test_a_fresh_run_has_nothing_to_reuse(self):
        assert not wanted_ref_tax_already_resolved(None)

    def test_a_previous_run_with_selected_accessions_can_be_reused(self):
        assert wanted_ref_tax_already_resolved(_previous_run_data())

    def test_a_previous_run_with_only_user_accessions_cannot(self):
        """
        Guarded on the accessions rather than on `-R` alone: if the run being resumed
        never got as far as selecting any, there's nothing to reuse and we have to
        resolve after all.
        """
        rd = _previous_run_data(accs=(), user_accs=("GCF_9",))
        assert not wanted_ref_tax_already_resolved(rd)


class TestSelectWantedRefTax:

    def test_resuming_skips_resolution_entirely(self, _never_resolves):
        args = _Args(wanted_ref_tax="Alteromonas")
        assert select_wanted_ref_tax(args, _previous_run_data()) == []

    def test_no_w_flag_skips_resolution_regardless(self, _never_resolves):
        assert select_wanted_ref_tax(_Args(), None) == []

    def test_a_fresh_run_still_resolves(self, monkeypatch):
        called = {}

        class _Selection:
            accessions = ["GCF_1"]
            warnings = []

        def _fake(*a, **kw):
            called["yes"] = True
            return ["GCF_1"], _Selection()

        monkeypatch.setattr(preflight_checks, "resolve_wanted_ref_tax_accessions", _fake)

        args = _Args(wanted_ref_tax="Alteromonas")
        assert len(select_wanted_ref_tax(args, None)) == 1
        assert called


class TestReusedAccessionsSurviveTheSkip:

    def test_the_previous_run_s_accessions_are_still_there(self):
        """
        This is what makes the skip invisible: the RUN INFO block reads the count off
        run_data, not off the selection, so the `-w` line prints the same either way.
        """
        rd = _previous_run_data(accs=("GCF_1", "GCF_2", "GCF_3"))
        assert len(rd.get_wanted_ref_tax_accs()) == 3

        merge_wanted_ref_tax(rd, None)   # what preflight does with a skipped selection

        assert len(rd.get_wanted_ref_tax_accs()) == 3

    def test_reference_and_user_accessions_stay_distinguishable(self):
        rd = _previous_run_data(accs=("GCF_1",), user_accs=("GCF_9",))
        merge_wanted_ref_tax(rd, None)
        assert rd.get_wanted_ref_tax_accs() == ["GCF_1"]
        assert rd.get_user_provided_ncbi_accs() == ["GCF_9"]
