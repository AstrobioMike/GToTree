"""
Settling `-H` in preflight: explicit, auto-selected from `-w`, or carried over.

These cover the seam rather than the picking itself (that's
test_scg_set_autopick.py): which of the three routes a given invocation takes, and
that `-H` always comes out canonicalized whichever one it was.
"""

import argparse
import pytest  # type: ignore

from gtotree.utils.misc import preflight_checks as P
from gtotree.utils.hmms.scg_hmms_setup import AutoPickedSCGSet


def make_args(**overrides):
    defaults = dict(hmm=None, hmm_auto_selected=None, source="gtdb",
                    wanted_ref_tax=None, ncbi_accessions=None, genbank_files=None,
                    fasta_files=None, amino_acid_files=None)
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


class FakeRunData:
    def __init__(self, hmm):
        self.fingerprint = {"hmm": hmm}


@pytest.fixture(autouse=True)
def packaged_sets(monkeypatch):
    """
    Stand in for hmm-sources-and-info.tsv, which resolve_hmm reads through
    canonicalize_hmm_arg, so these don't need GToTree_HMM_dir set.
    """
    monkeypatch.setattr(
        "gtotree.utils.hmms.scg_hmms_setup.packaged_scg_set_names",
        lambda: ["Bacteria", "Archaea", "Nitrospirota", "Universal-Hug-et-al"])


@pytest.fixture
def picks(monkeypatch):
    """Record what auto-selection was asked for, and dictate what it answers."""
    calls = []

    def install(name, reason="of a reason"):
        def fake(source, selection):
            calls.append((source, selection))
            return AutoPickedSCGSet(name, reason=reason)
        monkeypatch.setattr(P, "autopick_scg_set", fake)
        return calls

    return install


################################################################################
# which route -H takes
################################################################################

def test_an_explicit_hmm_is_never_second_guessed(picks):
    calls = picks("Bacteria")
    args = make_args(hmm="Universal-Hug-et-al", wanted_ref_tax="Nitrospirota")

    args = P.resolve_hmm(args, selection=object())

    assert args.hmm == "Universal-Hug-et-al"
    assert args.hmm_auto_selected is None
    assert calls == []


def test_a_missing_hmm_is_auto_selected_from_the_taxon(picks, capsys):
    picks("Nitrospirota", reason="'Nitrospirota' has a pre-built set of its own")
    args = make_args(wanted_ref_tax="Nitrospirota")
    selection = object()

    args = P.resolve_hmm(args, selection=selection)

    assert args.hmm == "Nitrospirota"
    assert args.hmm_auto_selected == "'Nitrospirota' has a pre-built set of its own"
    # and it says so, since the user didn't ask for this set by name
    assert "auto-selected" in capsys.readouterr().out


def test_the_selection_and_source_are_what_gets_handed_to_the_picker(picks):
    calls = picks("Archaea")
    selection = object()

    P.resolve_hmm(make_args(source="ncbi", wanted_ref_tax="Archaea"), selection=selection)

    assert calls == [("ncbi", selection)]


################################################################################
# resuming
################################################################################

def test_a_resume_carries_the_original_runs_set_forward(picks):
    """
    Re-deriving would be wasted work at best and wrong at worst -- the reference assets
    can be updated between runs, and a run should finish with the set it started with
    rather than tripping its own fingerprint check.
    """
    calls = picks("Bacteria")
    args = make_args(wanted_ref_tax="Nitrospirota")

    args = P.resolve_hmm(args, selection=object(),
                         previous_run_data=FakeRunData("Nitrospirota"))

    assert args.hmm == "Nitrospirota"
    assert "resumed" in args.hmm_auto_selected
    assert calls == []


def test_a_resume_with_nothing_recorded_falls_back_to_picking(picks):
    calls = picks("Bacteria")

    args = P.resolve_hmm(make_args(wanted_ref_tax="Bacteria"), selection=object(),
                         previous_run_data=FakeRunData(None))

    assert args.hmm == "Bacteria"
    assert len(calls) == 1


################################################################################
# canonicalization applies either way
################################################################################

@pytest.mark.parametrize("given, expected", [
    ("universal", "Universal-Hug-et-al"),
    ("bacteria", "Bacteria"),
    ("Bacteria.HMM", "Bacteria"),
])
def test_an_explicit_hmm_still_gets_canonicalized(given, expected):
    args = P.resolve_hmm(make_args(hmm=given))

    assert args.hmm == expected


def test_an_auto_selected_name_goes_through_the_same_canonicalization(picks):
    """
    The fingerprint, the reporting and the filename under GToTree_HMM_dir all key off
    this string, so an auto-selected set has to spell it the same way an explicit one
    would.
    """
    picks("universal")

    args = P.resolve_hmm(make_args(wanted_ref_tax="Ascomycota"), selection=object())

    assert args.hmm == "Universal-Hug-et-al"


################################################################################
# when -H is required at all
################################################################################

def test_no_hmm_and_no_wanted_ref_tax_is_an_error(capsys):
    with pytest.raises(SystemExit):
        P.check_for_minimum_args(make_args(fasta_files="genomes.txt"))

    assert "`-w`/`--wanted-ref-tax`" in capsys.readouterr().out


def test_no_hmm_is_fine_when_a_taxon_was_given():
    P.check_for_minimum_args(make_args(wanted_ref_tax="Nitrospirota"))


def test_no_input_genomes_at_all_is_still_an_error():
    with pytest.raises(SystemExit):
        P.check_for_minimum_args(make_args(hmm="Bacteria"))
