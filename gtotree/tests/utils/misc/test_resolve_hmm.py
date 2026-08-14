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
        def fake(source, selections):
            calls.append((source, selections))
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

    args = P.resolve_hmm(args, selections=object())

    assert args.hmm == "Universal-Hug-et-al"
    assert args.hmm_auto_selected is None
    assert calls == []


def test_a_missing_hmm_is_auto_selected_from_the_taxon(picks, capsys):
    picks("Nitrospirota", reason="")
    args = make_args(wanted_ref_tax="Nitrospirota")
    selection = object()

    args = P.resolve_hmm(args, selections=selection)

    assert args.hmm == "Nitrospirota"


def test_the_selection_and_source_are_what_gets_handed_to_the_picker(picks):
    calls = picks("Archaea")
    selection = object()

    P.resolve_hmm(make_args(source="ncbi", wanted_ref_tax="Archaea"), selections=selection)

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

    args = P.resolve_hmm(args, selections=object(),
                         previous_run_data=FakeRunData("Nitrospirota"))

    assert args.hmm == "Nitrospirota"
    assert "resumed" in args.hmm_auto_selected
    assert calls == []


def test_a_resume_with_nothing_recorded_falls_back_to_picking(picks):
    calls = picks("Bacteria")

    args = P.resolve_hmm(make_args(wanted_ref_tax="Bacteria"), selections=object(),
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

    args = P.resolve_hmm(make_args(wanted_ref_tax="Ascomycota"), selections=object())

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


################################################################################
# viral targets never reach the picker
################################################################################

class FakeSelection:
    def __init__(self, canonical, rows):
        self.canonical = canonical
        self.resolved_rank = "class"
        self.rows = rows
        self.accessions = []


def test_a_viral_target_exits_before_anything_is_picked(picks, capsys):
    """
    The refusal has to come ahead of BOTH routes through resolve_hmm, so it lands the
    same whether -H was left off or pointed at a pre-built set.
    """
    calls = picks("Bacteria")
    args = make_args(wanted_ref_tax="Caudoviricetes")
    selection = FakeSelection("Caudoviricetes", [{"domain": "Viruses"}])

    with pytest.raises(SystemExit):
        P.resolve_hmm(args, selections=[selection])

    assert calls == []
    assert "Caudoviricetes" in capsys.readouterr().out


def test_a_viral_target_with_an_explicit_pre_built_set_exits_too(capsys):
    args = make_args(hmm="Universal-Hug-et-al", wanted_ref_tax="Caudoviricetes")
    selection = FakeSelection("Caudoviricetes", [{"domain": "Viruses"}])

    with pytest.raises(SystemExit):
        P.resolve_hmm(args, selections=[selection])

    assert "Universal-Hug-et-al" in capsys.readouterr().out


def test_a_viral_target_with_the_users_own_hmm_file_carries_on(tmp_path):
    hmm_file = tmp_path / "viral-SCGs.hmm"
    hmm_file.write_text("HMMER3/f\n")

    args = make_args(hmm=str(hmm_file), wanted_ref_tax="Caudoviricetes")
    selection = FakeSelection("Caudoviricetes", [{"domain": "Viruses"}])

    args = P.resolve_hmm(args, selections=[selection])

    assert args.hmm == str(hmm_file)
