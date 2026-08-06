import argparse
import os
import pytest  # type: ignore
from gtotree.utils.hmms import scg_hmms_setup as S


PACKAGED = ["Actinobacteria", "Bacteria", "Bacteria_and_Archaea",
            "Gammaproteobacteria", "Universal-Hug-et-al"]


@pytest.fixture(autouse=True)
def packaged_sets(monkeypatch):
    """Stand in for hmm-sources-and-info.tsv so tests don't need GToTree_HMM_dir set."""
    monkeypatch.setattr(S, "packaged_scg_set_names", lambda: list(PACKAGED))


@pytest.mark.parametrize("given, expected", [
    ("Bacteria", "Bacteria"),
    ("bacteria", "Bacteria"),
    ("BACTERIA", "Bacteria"),
    ("BaCtErIa", "Bacteria"),
    ("bacteria.hmm", "Bacteria"),
    ("Bacteria.HMM", "Bacteria"),
    ("bacteria_and_archaea", "Bacteria_and_Archaea"),
    ("GAMMAPROTEOBACTERIA", "Gammaproteobacteria"),
    ("universal-hug-et-al", "Universal-Hug-et-al"),
])
def test_packaged_names_resolve_regardless_of_case(given, expected):
    assert S.canonicalize_hmm_arg(given) == expected


@pytest.mark.parametrize("given", ["Universal", "universal", "UNIVERSAL"])
def test_the_universal_alias_still_works_in_any_case(given):
    assert S.canonicalize_hmm_arg(given) == "Universal-Hug-et-al"


@pytest.mark.parametrize("given", ["Bacterai", "not-a-set", "/some/path/mine.hmm"])
def test_unrecognized_names_pass_through_untouched(given):
    """Typos and paths must reach the existing not-found reporting unchanged."""
    assert S.canonicalize_hmm_arg(given) == given


def test_resolution_survives_an_unreadable_summary_table(monkeypatch):
    # the autouse fixture stubs packaged_scg_set_names; this test needs the real one
    monkeypatch.undo()

    def boom():
        raise OSError("no table here")
    monkeypatch.setattr(S, "read_in_hmm_summary_table", boom)
    assert S.packaged_scg_set_names() == []
    assert S.canonicalize_hmm_arg("bacteria") == "bacteria"


# ---------------------------------------------------------------------------
# user-supplied files stay case-sensitive
# ---------------------------------------------------------------------------

def test_a_user_file_is_found_and_left_alone(tmp_path, monkeypatch):
    mine = tmp_path / "MyOwnSet.hmm"
    mine.write_text("HMMER3/f\n")
    monkeypatch.chdir(tmp_path)

    assert S.find_local_hmm_file("MyOwnSet") == str(mine)
    assert S.find_local_hmm_file("MyOwnSet.hmm") == str(mine)


def test_a_user_file_is_not_case_folded(tmp_path, monkeypatch):
    """
    Filesystems are case-sensitive nearly everywhere GToTree runs, so a mistyped path
    must fail rather than silently resolve to a different file.
    """
    if os.path.exists(str(tmp_path).swapcase()):
        pytest.skip("case-insensitive filesystem")

    (tmp_path / "MyOwnSet.hmm").write_text("HMMER3/f\n")
    monkeypatch.chdir(tmp_path)

    assert S.find_local_hmm_file("myownset") is None


def test_a_local_file_wins_over_a_packaged_name(tmp_path, monkeypatch):
    """A Bacteria.hmm in the working directory is used ahead of the packaged set."""
    local = tmp_path / "Bacteria.hmm"
    local.write_text("HMMER3/f\n")
    monkeypatch.chdir(tmp_path)

    args = argparse.Namespace(hmm="Bacteria")
    run_data = argparse.Namespace(hmm_path=None, SCG_targets=None)
    monkeypatch.setattr(S, "check_gathering_cutoffs", lambda *a, **k: None)
    monkeypatch.setattr(S, "get_SCG_hmm_targets", lambda p: ["PF00001"])
    monkeypatch.setattr(S, "get_hmm_path",
                        lambda *a: pytest.fail("packaged lookup should not be reached"))

    S.check_hmm_file(args, run_data)
    assert run_data.hmm_path == str(local)


def test_check_hmm_file_rewrites_args_hmm_to_official_casing(monkeypatch, tmp_path):
    """
    preflight_checks keys `tools_used.universal_SCGs_used` off args.hmm, and messaging
    prints it, so the canonical name has to land back on args.
    """
    monkeypatch.chdir(tmp_path)
    seen = {}

    def fake_get_hmm_path(hmm_file, hmm_arg):
        seen["file"], seen["arg"] = hmm_file, hmm_arg
        return "/fake/Bacteria.hmm"

    monkeypatch.setattr(S, "get_hmm_path", fake_get_hmm_path)
    monkeypatch.setattr(S, "check_gathering_cutoffs", lambda *a, **k: None)
    monkeypatch.setattr(S, "get_SCG_hmm_targets", lambda p: ["PF00001"])

    args = argparse.Namespace(hmm="bACTERIA")
    S.check_hmm_file(args, argparse.Namespace(hmm_path=None, SCG_targets=None))

    assert args.hmm == "Bacteria"
    assert seen["file"] == "Bacteria.hmm", "a lowercase name would download a second copy"
