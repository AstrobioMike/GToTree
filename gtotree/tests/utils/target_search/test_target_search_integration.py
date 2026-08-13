"""
End-to-end tests for `gtt search-pfams`.

These drive `run_search` with only the managed-Pfam *download* stubbed out; target
collection, the fused preprocess-and-search stage, the counts matrix, and the hit-seq
combining are all the production code paths -- the same functions a full GToTree run
with `-p` calls.

The point of running the whole thing rather than each piece is that the interesting
failure mode for this feature is wiring: every individual helper works (they're the
main driver's), and what could break is a directory attribute pointed somewhere wrong
or a search that silently never runs.
"""

import os
import time
import pytest  # type: ignore
from gtotree.utils.target_search import target_search_cli
from gtotree.utils.target_search import target_search_stages as stages
from gtotree.utils.target_search import target_search_outputs as outputs
from gtotree.utils.target_search.target_search_setup import (TargetSearchError,
                                                             make_args,
                                                             WORKING_DIR_NAME)
from gtotree.utils.misc.messaging import REMOVED_GENOMES_FILENAME
from gtotree.utils.misc.stages import GenomeRemovalStage
from gtotree.tests.utils.target_search.conftest import MOCK_PFAM_VERSION


def _read_tsv(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    header, *body = rows
    return header, [dict(zip(header, row, strict=True)) for row in body]


@pytest.fixture
def run_pfam_search(pfam_spec, tmp_path, listing, write_genome, capsys):
    """
    Returns a callable that runs a full search over the given {genome: [pfam ids]} map.
    """

    def _run(genomes, targets, **arg_overrides):
        paths = [write_genome(name, ids) for name, ids in genomes.items()]
        args = make_args(
            amino_acid_files=listing("aa.txt", paths),
            target_pfams_file=listing("targets.txt", targets),
            output_dir=str(tmp_path / "out"),
            num_jobs=2,
            **arg_overrides,
        )
        run_data = target_search_cli.run_search(args, pfam_spec)
        capsys.readouterr()
        return run_data, str(tmp_path / "out")

    return _run


################################################################################
# the happy path
################################################################################

@pytest.fixture
def basic_run(run_pfam_search):
    """
    g1 and g2 each carry one copy of PF90001 and PF90002; g3 carries two copies of
    PF90001 and one of PF90002. PF99999 doesn't exist, so it exercises the
    partially-failed-targets path in the same run.
    """
    return run_pfam_search(
        genomes={
            "g1": ["PF90001", "PF90002"],
            "g2": ["PF90001", "PF90002"],
            "g3": ["PF90001", "PF90002", "PF90001"],
        },
        targets=["PF90001", "PF90002", "PF99999"],
    )


def test_hit_counts_matrix_matches_the_planted_copies(basic_run, pfam_spec):
    _run_data, out_dir = basic_run
    header, rows = _read_tsv(os.path.join(out_dir, pfam_spec.counts_filename))

    assert header == ["genome_id", "total_gene_count", "PF90001", "PF90002"]
    by_genome = {row["genome_id"]: row for row in rows}

    assert by_genome["g1"]["PF90001"] == "1"
    assert by_genome["g2"]["PF90001"] == "1"
    # the planted second copy
    assert by_genome["g3"]["PF90001"] == "2"
    assert by_genome["g3"]["PF90002"] == "1"
    assert by_genome["g3"]["total_gene_count"] == "3"


def test_hit_seqs_are_combined_per_target(basic_run, pfam_spec):
    _run_data, out_dir = basic_run
    hit_dir = os.path.join(out_dir, pfam_spec.hit_seqs_subdir)

    assert sorted(os.listdir(hit_dir)) == ["PF90001.faa", "PF90002.faa"]

    headers = [line.strip() for line in open(os.path.join(hit_dir, "PF90001.faa"))
               if line.startswith(">")]
    # four hits total: one each from g1 and g2, two from g3
    assert len(headers) == 4
    # headers are genome-prefixed by the shared processing, so hits stay traceable
    assert sum(h.startswith(">g3_") for h in headers) == 2


def test_per_genome_result_files_are_written(basic_run):
    _run_data, out_dir = basic_run
    for genome in ("g1", "g2", "g3"):
        path = os.path.join(out_dir, "individual-genome-results", genome,
                            "pfam-hmmsearch.txt")
        assert os.path.isfile(path), f"no per-genome tblout for {genome}"


def test_failed_targets_are_reported(basic_run, pfam_spec):
    _run_data, out_dir = basic_run
    path = os.path.join(out_dir, pfam_spec.failed_targets_filename)
    assert open(path).read().split() == ["PF99999"]


def test_requested_and_pulled_table_records_every_request(basic_run):
    _run_data, out_dir = basic_run
    header, rows = _read_tsv(os.path.join(out_dir, "info",
                                          "requested-and-pulled-pfams.tsv"))
    assert header == ["specified_pfam", "pulled_pfam"]
    pulled = {row["specified_pfam"]: row["pulled_pfam"] for row in rows}
    # the versioned accession as it appears in the master HMM
    assert pulled["PF90001"] == "PF90001.3"
    assert pulled["PF99999"] == "NA"


def test_summary_table_covers_every_input_genome(basic_run):
    _run_data, out_dir = basic_run
    header, rows = _read_tsv(os.path.join(out_dir, "genomes-summary-info.tsv"))

    assert "search_completed" in header
    assert [row["genome_id"] for row in rows] == ["g1", "g2", "g3"]
    assert all(row["search_completed"] == "Yes" for row in rows)
    assert all(row["source"] == "amino-acid-fasta" for row in rows)
    assert rows[2]["num_genes"] == "3"


def test_no_nested_search_results_directory(basic_run):
    """The flattened layout: no `pfam-search-results/` wrapper."""
    _run_data, out_dir = basic_run
    assert not os.path.exists(os.path.join(out_dir, "pfam-search-results"))


def test_working_dir_is_cleaned_up_on_success(basic_run):
    _run_data, out_dir = basic_run
    assert not os.path.exists(os.path.join(out_dir, WORKING_DIR_NAME))


def test_no_scg_or_tree_artifacts_are_produced(basic_run):
    """
    `do_scg=False` has to actually stop the SCG search from running -- if it leaked
    through it would fail on the absent `run_data.hmm_path` rather than no-op.
    """
    _run_data, out_dir = basic_run
    names = os.listdir(out_dir)
    assert not any("SCG" in n or "aligned" in n or n.endswith(".tre")
                   for n in names), names


def test_genomes_are_marked_searched(basic_run, pfam_spec):
    run_data, _out_dir = basic_run
    for gd in run_data.all_input_genomes:
        assert getattr(gd, pfam_spec.search_done_flag) is True
        assert gd.hmm_search_done is False, \
            "the SCG search ran even though do_scg was False"


################################################################################
# failure paths
################################################################################

def test_no_findable_targets_is_a_friendly_error(run_pfam_search):
    with pytest.raises(TargetSearchError, match="none of the requested pfams"):
        try:
            run_pfam_search(genomes={"g1": ["PF90001"]}, targets=["PF99998", "PF99999"])
        except TargetSearchError as e:
            raise TargetSearchError(str(e).lower()) from e


def test_a_genome_with_no_hits_still_appears_everywhere(run_pfam_search, pfam_spec):
    """
    A zero-hit genome is a result, not an absence -- it needs a row in the counts
    matrix and the summary, or the output silently under-reports the input set.
    """
    _run_data, out_dir = run_pfam_search(
        genomes={"hit": ["PF90001"], "nohit": ["PF90004"]},
        targets=["PF90001"],
    )

    _header, rows = _read_tsv(os.path.join(out_dir, pfam_spec.counts_filename))
    by_genome = {row["genome_id"]: row for row in rows}
    assert by_genome["nohit"]["PF90001"] == "0"

    _header, summary = _read_tsv(os.path.join(out_dir, "genomes-summary-info.tsv"))
    assert {row["genome_id"] for row in summary} == {"hit", "nohit"}


def test_single_genome_run_works(run_pfam_search):
    """No four-genome floor."""
    run_data, out_dir = run_pfam_search(genomes={"solo": ["PF90001"]},
                                        targets=["PF90001"])
    assert len(run_data.all_input_genomes) == 1
    assert os.path.isfile(os.path.join(out_dir, "pfam-hit-seqs", "PF90001.faa"))


################################################################################
# a search that fails
################################################################################

@pytest.fixture
def fail_search_for(monkeypatch):
    """
    Returns a callable making the Pfam search fail for the named genomes.

    Patched at `processing_genomes`, which is where the fused worker looks the search
    function up, and returning the same shape the real worker returns when it catches
    something -- so this exercises the production failure path rather than a synthetic
    one bolted on beside it.
    """
    from gtotree.main_stages import processing_genomes

    real_worker = processing_genomes._pfam_search_worker

    def _install(*genome_ids):
        wanted = set(genome_ids)

        def worker(genome, run_data, **kwargs):
            if genome.id in wanted:
                return {"pfam_search_failed": True, "error": "RuntimeError: planted"}
            return real_worker(genome, run_data, **kwargs)

        monkeypatch.setattr(processing_genomes, "_pfam_search_worker", worker)

    return _install


def test_a_failed_search_drops_the_genome_from_the_run(run_pfam_search, fail_search_for):
    """
    In a full GToTree run a failed Pfam search leaves the genome in place -- it still
    has its SCG hits and still belongs in the tree. Here the search is the whole run,
    so a genome that failed it has nothing left to contribute and has to leave.
    """
    fail_search_for("broken")
    run_data, _out_dir = run_pfam_search(
        genomes={"fine": ["PF90001"], "broken": ["PF90001"]},
        targets=["PF90001"],
    )

    by_id = {gd.id: gd for gd in run_data.all_input_genomes}

    assert by_id["broken"].removed
    assert by_id["broken"].removed_at == GenomeRemovalStage.HMM_SEARCH
    assert "search failed" in by_id["broken"].reason_removed
    assert not by_id["fine"].removed


def test_a_failed_search_is_recorded_in_the_removals_table(run_pfam_search,
                                                           fail_search_for):
    fail_search_for("broken")
    _run_data, out_dir = run_pfam_search(
        genomes={"fine": ["PF90001"], "broken": ["PF90001"]},
        targets=["PF90001"],
    )

    header, rows = _read_tsv(os.path.join(out_dir, REMOVED_GENOMES_FILENAME))
    assert header == ["genome_id", "input", "source", "stage_removed",
                      "reason_removed"]
    assert [row["genome_id"] for row in rows] == ["broken"]
    assert rows[0]["stage_removed"] == GenomeRemovalStage.HMM_SEARCH


def test_a_failed_search_stays_out_of_the_counts_matrix(run_pfam_search,
                                                        fail_search_for, pfam_spec):
    """
    The regression this guards: `mark_pfam_search_failed` sets the *done* flag too, and
    the counts writer selects on `search_done and not removed` -- so before the genome
    was removed it landed in the matrix as a row of zeros, indistinguishable from a
    genome that really was searched and hit nothing.
    """
    fail_search_for("broken")
    _run_data, out_dir = run_pfam_search(
        genomes={"fine": ["PF90001"], "broken": ["PF90001"]},
        targets=["PF90001"],
    )

    _header, rows = _read_tsv(os.path.join(out_dir, pfam_spec.counts_filename))
    assert [row["genome_id"] for row in rows] == ["fine"]

    # but it still gets a summary row, saying plainly that it didn't finish
    _header, summary = _read_tsv(os.path.join(out_dir, "genomes-summary-info.tsv"))
    by_id = {row["genome_id"]: row for row in summary}
    assert by_id["broken"]["search_completed"] == "No"
    assert by_id["fine"]["search_completed"] == "Yes"


def test_a_failed_search_is_not_counted_as_searched(run_pfam_search, fail_search_for,
                                                    pfam_spec):
    fail_search_for("broken")
    run_data, _out_dir = run_pfam_search(
        genomes={"fine": ["PF90001"], "broken": ["PF90001"]},
        targets=["PF90001"],
    )

    assert stages.count_searched_genomes(run_data, pfam_spec) == 1

    searched, removed, failed = outputs.summarize_counts(run_data, pfam_spec)
    assert (searched, removed, failed) == (1, 1, 0)


def test_every_genome_failing_its_search_is_a_friendly_error(run_pfam_search,
                                                             fail_search_for):
    fail_search_for("a", "b")
    with pytest.raises(TargetSearchError, match="completed search"):
        run_pfam_search(genomes={"a": ["PF90001"], "b": ["PF90001"]},
                        targets=["PF90001"])


################################################################################
# which removals the search phase reports
################################################################################

def test_search_phase_stages_exclude_the_phase_one_lookup():
    """
    The NCBI lookup happens in phase 1 and is reported there. Counting it again in
    phase 3 was the bug: total-input minus searched swept in genomes the search phase
    never touched and blamed them on it.
    """
    assert GenomeRemovalStage.NCBI_LOOKUP not in stages.SEARCH_PHASE_REMOVAL_STAGES
    assert GenomeRemovalStage.NCBI_DOWNLOAD in stages.SEARCH_PHASE_REMOVAL_STAGES
    assert GenomeRemovalStage.HMM_SEARCH in stages.SEARCH_PHASE_REMOVAL_STAGES


def test_only_this_phases_removals_are_counted():
    """
    A genome lost at the lookup and one lost at the search sit in the same removals
    table; only the second is the search phase's to report.
    """
    from gtotree.utils.misc.general import RunData, GenomeData

    lost_at_lookup = GenomeData.from_acc("GCF_000000001.1")
    lost_at_lookup.mark_removed("not found at NCBI",
                                GenomeRemovalStage.NCBI_LOOKUP)

    lost_at_search = GenomeData.from_acc("GCF_000000002.1")
    lost_at_search.mark_removed("Pfam search failed",
                                GenomeRemovalStage.HMM_SEARCH)

    survivor = GenomeData.from_acc("GCF_000000003.1")

    run_data = RunData()
    run_data.ncbi_accs = [lost_at_lookup, lost_at_search, survivor]
    run_data.update_all_input_genomes()

    reported = run_data.genomes_removed_at(*stages.SEARCH_PHASE_REMOVAL_STAGES)
    assert [gd.id for gd in reported] == ["GCF_000000002.1"]


################################################################################
# resume
################################################################################

def test_rewriting_identical_genome_contents_leaves_the_file_alone(write_genome):
    """
    Guards the resume tests below. `hash_local_genomes` fingerprints local inputs by
    size and `int(st_mtime)`, so re-saving a byte-identical genome between two runs
    changes the fingerprint whenever the two writes straddle a whole-second boundary --
    which made `test_resume_reuses_finished_genomes` fail intermittently. The fixture
    now skips the write when nothing changed; if that's ever undone, this fails
    deterministically instead of the resume test failing one run in twelve.
    """
    path = write_genome("stable", ["PF90001"])
    before = path.stat().st_mtime_ns

    time.sleep(0.01)
    assert write_genome("stable", ["PF90001"]) == path
    assert path.stat().st_mtime_ns == before

    # a genuinely different genome is still rewritten, so tests can model an edit
    time.sleep(0.01)
    write_genome("stable", ["PF90002"])
    assert path.stat().st_mtime_ns != before


def test_resume_reuses_finished_genomes(run_pfam_search, tmp_path, pfam_spec, capsys):
    genomes = {"g1": ["PF90001"], "g2": ["PF90001"]}
    targets = ["PF90001"]

    run_pfam_search(genomes=genomes, targets=targets, keep_working_dir=True)

    run_data, out_dir = run_pfam_search(genomes=genomes, targets=targets,
                                        keep_working_dir=True, resume=True)

    # the results are still complete after a run that did no new work
    _header, rows = _read_tsv(os.path.join(out_dir, pfam_spec.counts_filename))
    assert len(rows) == 2
    assert all(gd.pfam_search_done for gd in run_data.all_input_genomes)


def test_resume_refuses_a_changed_target_list(run_pfam_search):
    run_pfam_search(genomes={"g1": ["PF90001"]}, targets=["PF90001"],
                    keep_working_dir=True)

    with pytest.raises(TargetSearchError, match="list of search targets changed"):
        run_pfam_search(genomes={"g1": ["PF90001"]},
                        targets=["PF90001", "PF90002"],
                        keep_working_dir=True, resume=True)


def test_resume_refuses_a_changed_genome_set(run_pfam_search):
    run_pfam_search(genomes={"g1": ["PF90001"]}, targets=["PF90001"],
                    keep_working_dir=True)

    with pytest.raises(TargetSearchError, match="local genome files"):
        run_pfam_search(genomes={"g1": ["PF90001"], "g2": ["PF90001"]},
                        targets=["PF90001"], keep_working_dir=True, resume=True)


def test_resume_without_a_previous_run_starts_fresh(run_pfam_search, pfam_spec):
    _run_data, out_dir = run_pfam_search(genomes={"g1": ["PF90001"]},
                                         targets=["PF90001"], resume=True)
    assert os.path.isfile(os.path.join(out_dir, pfam_spec.counts_filename))


def test_pfam_version_is_recorded_in_the_run_state(run_pfam_search, tmp_path):
    import json
    run_pfam_search(genomes={"g1": ["PF90001"]}, targets=["PF90001"],
                    keep_working_dir=True)
    state_path = tmp_path / "out" / WORKING_DIR_NAME / "run-state.json"
    state = json.loads(state_path.read_text())
    assert state["fingerprint"]["data_version"] == MOCK_PFAM_VERSION
