"""
Unit tests for the multi-target-type logic in the `search-annotations` driver.

The end-to-end behavior (running a search, the per-type layout, the summaries) is
covered against the real code path in test_target_search_integration.py, which drives
the same driver with a one-spec list. What can't be exercised there is the logic that
only matters when *both* target types are in play -- and which can't run a real search
because the KO worker needs the external `exec_annotation` binary. Those pieces are
pure functions over args/specs, so they're tested directly here:

  * requested-spec selection (which of `-p` / `-K` a run actually acts on),
  * combined-argument validation (at least one target type, genomes required),
  * SearchPlan construction for each flag combination, and
  * the independence of the two target types in the resume fingerprint.
"""

import pytest  # type: ignore

from gtotree.utils.target_search import target_search_cli as driver
from gtotree.utils.target_search.target_search_spec import get_spec
from gtotree.utils.target_search.target_search_setup import (TargetSearchError,
                                                             make_args,
                                                             check_args_multi,
                                                             requested_specs)


################################################################################
# requested-spec selection
################################################################################

def test_requested_specs_picks_only_flagged_types():
    specs = [get_spec("pfam"), get_spec("ko")]

    both = make_args(target_pfams_file="p.txt", target_kos_file="k.txt")
    assert [s.subcommand for s in requested_specs(both, specs)] == \
        ["search-pfams", "search-kos"]

    pfam_only = make_args(target_pfams_file="p.txt")
    assert [s.subcommand for s in requested_specs(pfam_only, specs)] == \
        ["search-pfams"]

    ko_only = make_args(target_kos_file="k.txt")
    assert [s.subcommand for s in requested_specs(ko_only, specs)] == \
        ["search-kos"]


################################################################################
# combined-argument validation
################################################################################

def test_check_args_multi_requires_at_least_one_target_type():
    specs = [get_spec("pfam"), get_spec("ko")]

    args = make_args(amino_acid_files="aa.txt")  # genomes but no targets
    with pytest.raises(TargetSearchError) as excinfo:
        check_args_multi(args, specs)
    assert "-p" in str(excinfo.value) and "-K" in str(excinfo.value)


def test_check_args_multi_requires_genomes():
    specs = [get_spec("pfam"), get_spec("ko")]

    args = make_args(target_pfams_file="p.txt")  # targets but no genomes
    with pytest.raises(TargetSearchError) as excinfo:
        check_args_multi(args, specs)
    assert "input genomes" in str(excinfo.value).lower()


def test_check_args_multi_accepts_one_or_both_types():
    specs = [get_spec("pfam"), get_spec("ko")]

    for overrides in (
        {"target_pfams_file": "p.txt"},
        {"target_kos_file": "k.txt"},
        {"target_pfams_file": "p.txt", "target_kos_file": "k.txt"},
    ):
        args = make_args(amino_acid_files="aa.txt", **overrides)
        assert check_args_multi(args, specs) is args  # should not raise


################################################################################
# combined SearchPlan
################################################################################

def test_combined_plan_turns_on_only_requested_searches():
    pfam = get_spec("pfam")
    ko = get_spec("ko")

    both = driver._combined_plan(make_args(), [pfam, ko])
    assert both.do_pfam and both.do_ko and not both.do_scg

    pfam_only = driver._combined_plan(make_args(), [pfam])
    assert pfam_only.do_pfam and not pfam_only.do_ko

    ko_only = driver._combined_plan(make_args(), [ko])
    assert ko_only.do_ko and not ko_only.do_pfam


def test_combined_plan_forwards_keep_working_dir():
    plan = driver._combined_plan(make_args(keep_working_dir=True), [get_spec("pfam")])
    assert plan.keep_genome_files is True


################################################################################
# fingerprint independence of the two target types
################################################################################

class _FakeRunData:
    """Minimal stand-in exposing what build_fingerprint reads off RunData."""

    genbank_files = ()
    fasta_files = ()
    amino_acid_files = ()

    def get_input_ncbi_accs(self):
        return []


def _fingerprint(args, tmp_path, pfam_ids=None, ko_ids=None):
    # build_fingerprint hashes the target files by contents, so real files are needed
    if pfam_ids is not None:
        p = tmp_path / "p.txt"
        p.write_text("\n".join(pfam_ids) + "\n")
        args.target_pfams_file = str(p)
    if ko_ids is not None:
        k = tmp_path / "k.txt"
        k.write_text("\n".join(ko_ids) + "\n")
        args.target_kos_file = str(k)

    specs = [get_spec("pfam"), get_spec("ko")]
    return driver.build_fingerprint(_FakeRunData(), args, specs)


def test_editing_pfam_targets_changes_only_the_pfam_hash(tmp_path):
    args1 = make_args(amino_acid_files="aa.txt")
    fp1 = _fingerprint(args1, tmp_path, pfam_ids=["PF00001"], ko_ids=["K00001"])

    args2 = make_args(amino_acid_files="aa.txt")
    fp2 = _fingerprint(args2, tmp_path, pfam_ids=["PF00002"], ko_ids=["K00001"])

    assert fp1["pfam_targets_sha256"] != fp2["pfam_targets_sha256"]
    assert fp1["ko_targets_sha256"] == fp2["ko_targets_sha256"]


def test_adding_a_target_type_changes_the_fingerprint(tmp_path):
    args1 = make_args(amino_acid_files="aa.txt")
    fp_pfam_only = _fingerprint(args1, tmp_path, pfam_ids=["PF00001"])

    args2 = make_args(amino_acid_files="aa.txt")
    fp_both = _fingerprint(args2, tmp_path, pfam_ids=["PF00001"], ko_ids=["K00001"])

    # adding -K is a real change to what the run produces
    assert fp_pfam_only["target_types"] != fp_both["target_types"]
    assert fp_pfam_only["ko_targets_sha256"] != fp_both["ko_targets_sha256"]
