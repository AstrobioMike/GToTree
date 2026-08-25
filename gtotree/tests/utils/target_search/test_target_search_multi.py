"""
Tests for the combined `gtt search-annotations` command.

Two layers:

  * End-to-end runs through `target_search_multi.run_search` with only the managed-Pfam
    *download* stubbed out -- target collection, the fused preprocess-and-search stage,
    the counts matrix, and hit-seq combining are all production code. These use the
    Pfam path only, because the KO search worker needs the external `exec_annotation`
    binary (the single-command integration suite is Pfam-only for the same reason). A
    Pfam-only run through the *combined* driver still exercises everything specific to
    this command: multi-spec argument validation, per-type subdirectory layout, the
    combined SearchPlan, the run-level summary, and the combined finish banner.

  * Unit tests for the multi-specific logic that doesn't need a search to run:
    requested-spec selection, plan construction for each flag combination, and the
    independence of the two target types in the resume fingerprint.

The interesting failure mode for this feature is wiring -- every per-genome helper is
the main driver's and already tested -- so an interesting bug is a results directory
pointed at the wrong place, a search that silently never runs, or a fingerprint that
conflates the two target types.
"""

import os
import pytest  # type: ignore

from gtotree.utils.target_search import target_search_multi as multi
from gtotree.utils.target_search.target_search_spec import get_spec
from gtotree.utils.target_search.target_search_setup import (TargetSearchError,
                                                             make_args,
                                                             WORKING_DIR_NAME)
from gtotree.utils.target_search.target_search_outputs import SUMMARY_FILENAME


def _read_tsv(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    header, *body = rows
    return header, [dict(zip(header, row, strict=True)) for row in body]


@pytest.fixture
def run_annotation_search(pfam_spec, tmp_path, listing, write_genome, capsys,
                          monkeypatch):
    """
    Run the combined driver over a {genome: [pfam ids]} map, Pfam side only.

    `pfam_spec` (from the shared conftest) is the production Pfam spec with only the
    dataset download stubbed. It's substituted for the real Pfam spec everywhere the
    combined driver reaches for it, so the run uses the mock Pfam_data_dir.
    """
    monkeypatch.setattr(multi, "_all_specs", lambda: [pfam_spec])

    def fake_get_spec(name):
        return pfam_spec if name == "pfam" else get_spec(name)

    monkeypatch.setattr(multi, "get_spec", fake_get_spec)

    def _run(genomes, pfam_targets, **arg_overrides):
        paths = [write_genome(name, ids) for name, ids in genomes.items()]
        args = make_args(
            amino_acid_files=listing("aa.txt", paths),
            target_pfams_file=listing("pfam-targets.txt", pfam_targets),
            output_dir=str(tmp_path / "out"),
            num_jobs=2,
            **arg_overrides,
        )
        run_data = multi.run_search(args, specs=[pfam_spec])
        capsys.readouterr()
        return run_data, str(tmp_path / "out")

    return _run


################################################################################
# end-to-end: single target type through the combined command
################################################################################

@pytest.fixture
def basic_run(run_annotation_search):
    return run_annotation_search(
        genomes={
            "g1": ["PF90001", "PF90002"],
            "g2": ["PF90001", "PF90002"],
            "g3": ["PF90001", "PF90002", "PF90001"],
        },
        pfam_targets=["PF90001", "PF90002", "PF99999"],
    )


def test_pfam_results_land_in_the_pfam_subdirectory(basic_run, pfam_spec):
    """The whole point of the layout: `<out>/pfam/`, not the flat top level."""
    _run_data, out_dir = basic_run

    pfam_dir = os.path.join(out_dir, "pfam")
    assert os.path.isdir(pfam_dir)
    assert os.path.isfile(os.path.join(pfam_dir, pfam_spec.counts_filename))
    # the flat single-command location must NOT be used by the combined command
    assert not os.path.isfile(os.path.join(out_dir, pfam_spec.counts_filename))


def test_hit_counts_matrix_matches_the_planted_copies(basic_run, pfam_spec):
    _run_data, out_dir = basic_run
    counts_path = os.path.join(out_dir, "pfam", pfam_spec.counts_filename)
    header, rows = _read_tsv(counts_path)

    assert header == ["genome_id", "total_gene_count", "PF90001", "PF90002"]
    by_genome = {row["genome_id"]: row for row in rows}

    assert by_genome["g1"]["PF90001"] == "1"
    assert by_genome["g2"]["PF90001"] == "1"
    assert by_genome["g3"]["PF90001"] == "2"
    assert by_genome["g3"]["PF90002"] == "1"


def test_run_level_summary_is_written_at_the_top(basic_run):
    """One genomes-summary for the whole run, not one per target type."""
    _run_data, out_dir = basic_run

    top_summary = os.path.join(out_dir, SUMMARY_FILENAME)
    assert os.path.isfile(top_summary)

    header, rows = _read_tsv(top_summary)
    assert {row["genome_id"] for row in rows} == {"g1", "g2", "g3"}
    assert all(row["search_completed"] == "Yes" for row in rows)


def test_hit_seqs_written_under_the_pfam_subdir(basic_run, pfam_spec):
    _run_data, out_dir = basic_run
    hit_seqs = os.path.join(out_dir, "pfam", pfam_spec.hit_seqs_subdir)
    assert os.path.isdir(hit_seqs)
    # PF90001 and PF90002 both hit; PF99999 didn't exist, so no file for it
    files = os.listdir(hit_seqs)
    assert any("PF90001" in f for f in files)
    assert any("PF90002" in f for f in files)
    assert not any("PF99999" in f for f in files)


def test_working_dir_is_cleaned_up(basic_run):
    _run_data, out_dir = basic_run
    assert not os.path.isdir(os.path.join(out_dir, WORKING_DIR_NAME))


################################################################################
# end-to-end: resume through the combined command
################################################################################

def test_resume_reuses_finished_genomes(run_annotation_search, tmp_path):
    genomes = {"g1": ["PF90001"], "g2": ["PF90002"]}
    targets = ["PF90001", "PF90002"]

    # keep_working_dir so the resume state and run-data survive to the second run,
    # matching the single-command resume tests
    run_annotation_search(genomes, targets, keep_working_dir=True)

    run_data, out_dir = run_annotation_search(genomes, targets, keep_working_dir=True,
                                              resume=True)

    # both genomes finished the first time, so the resumed run re-derives the same
    # completed set and produces the same counts matrix
    counts_path = os.path.join(out_dir, "pfam",
                               get_spec("pfam").counts_filename)
    assert os.path.isfile(counts_path)
    _header, rows = _read_tsv(counts_path)
    assert {row["genome_id"] for row in rows} == {"g1", "g2"}


def test_resume_refuses_when_pfam_targets_change(run_annotation_search):
    genomes = {"g1": ["PF90001", "PF90002"]}

    run_annotation_search(genomes, ["PF90001"], keep_working_dir=True)

    with pytest.raises(TargetSearchError) as excinfo:
        run_annotation_search(genomes, ["PF90002"], keep_working_dir=True, resume=True)

    # the refusal should name the changed thing
    assert "target" in str(excinfo.value).lower()


################################################################################
# unit: requested-spec selection
################################################################################

def test_requested_specs_picks_only_flagged_types():
    from gtotree.utils.target_search.target_search_setup import requested_specs

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


def test_check_args_multi_requires_at_least_one_target_type():
    from gtotree.utils.target_search.target_search_setup import check_args_multi

    specs = [get_spec("pfam"), get_spec("ko")]

    # genomes but no targets at all
    args = make_args(amino_acid_files="aa.txt")
    with pytest.raises(TargetSearchError) as excinfo:
        check_args_multi(args, specs)
    assert "-p" in str(excinfo.value) and "-K" in str(excinfo.value)


def test_check_args_multi_requires_genomes():
    from gtotree.utils.target_search.target_search_setup import check_args_multi

    specs = [get_spec("pfam"), get_spec("ko")]

    args = make_args(target_pfams_file="p.txt")
    with pytest.raises(TargetSearchError) as excinfo:
        check_args_multi(args, specs)
    assert "input genomes" in str(excinfo.value).lower()


def test_check_args_multi_accepts_one_or_both_types():
    from gtotree.utils.target_search.target_search_setup import check_args_multi

    specs = [get_spec("pfam"), get_spec("ko")]

    for overrides in (
        {"target_pfams_file": "p.txt"},
        {"target_kos_file": "k.txt"},
        {"target_pfams_file": "p.txt", "target_kos_file": "k.txt"},
    ):
        args = make_args(amino_acid_files="aa.txt", **overrides)
        # should not raise
        assert check_args_multi(args, specs) is args


################################################################################
# unit: combined SearchPlan
################################################################################

def test_combined_plan_turns_on_only_requested_searches():
    pfam = get_spec("pfam")
    ko = get_spec("ko")

    both = multi._combined_plan(make_args(), [pfam, ko])
    assert both.do_pfam and both.do_ko and not both.do_scg

    pfam_only = multi._combined_plan(make_args(), [pfam])
    assert pfam_only.do_pfam and not pfam_only.do_ko

    ko_only = multi._combined_plan(make_args(), [ko])
    assert ko_only.do_ko and not ko_only.do_pfam


def test_combined_plan_forwards_keep_working_dir():
    pfam = get_spec("pfam")
    plan = multi._combined_plan(make_args(keep_working_dir=True), [pfam])
    assert plan.keep_genome_files is True


################################################################################
# unit: fingerprint independence of the two target types
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
    return multi.build_fingerprint(_FakeRunData(), args, specs)


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
