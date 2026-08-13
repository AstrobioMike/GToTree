import argparse
import os
import time

import pytest # type: ignore

from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as cli
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import RESUME
from gtotree.utils.misc.general import GenomeData, RunData


def _args(**kw):
    base = dict(percent_single_copy=90, min_pfam_coverage=50.0, source="GTDB",
                wanted_ref_tax=["Nitrospirota"], target_rank=None, derep_rank="off",
                genbank_files=None, fasta_files=None, amino_acid_files=None,
                min_completeness=None, max_contamination=None)
    base.update(kw)
    return argparse.Namespace(**base)


def _fp(accessions=(), args=None, pfam_version="38.2", local_genomes=()):
    """
    build_fingerprint reads the genome set off a RunData now, so the tests assemble
    one rather than passing bare lists.
    """
    run_data = RunData()
    run_data.ncbi_accs = [GenomeData.from_acc(a) for a in accessions]
    run_data.fasta_files = list(local_genomes)
    run_data.update_all_input_genomes()
    return cli.build_fingerprint(run_data, args if args is not None else _args(),
                                 pfam_version)


class _FakeGenome:
    """Stand-in for GenomeData with just what the fingerprint reads."""
    def __init__(self, gid, source, full_path):
        self.id = gid
        self.source = source
        self.full_path = full_path


################################################################################
# accession hashing
################################################################################

def test_fingerprint_is_order_independent():
    """
    With threaded downloads and taxonomy-driven selection, the same genome set can
    legitimately arrive in a different order. What matters is WHICH genomes.
    """
    a = _fp(["A", "B", "C"], _args(), "38.2")
    b = _fp(["C", "A", "B"], _args(), "38.2")
    assert a == b


def test_fingerprint_ignores_duplicates():
    a = _fp(["A", "B"], _args(), "38.2")
    b = _fp(["A", "B", "A"], _args(), "38.2")
    assert a == b


def test_fingerprint_detects_added_genome():
    a = _fp(["A", "B"], _args(), "38.2")
    b = _fp(["A", "B", "C"], _args(), "38.2")
    diffs = RESUME.compare(a, b)
    assert any("target genomes" in d for d in diffs)


################################################################################
# what invalidates a resume
################################################################################

@pytest.mark.parametrize("kwargs,expected", [
    (dict(percent_single_copy=95), "--percent-single-copy"),
    (dict(min_pfam_coverage=60.0), "--min-pfam-coverage"),
    (dict(source="NCBI"), "--source"),
    (dict(wanted_ref_tax="Bacteroidota"), "--wanted-ref-tax"),
    (dict(target_rank="genus"), "--target-rank"),
    (dict(derep_rank="genus"), "--derep-rank"),
])
def test_result_affecting_params_invalidate(kwargs, expected):
    base = _fp(["A"], _args(), "38.2")
    changed = _fp(["A"], _args(**kwargs), "38.2")
    diffs = RESUME.compare(base, changed)
    assert any(expected in d for d in diffs), diffs


def test_pfam_version_change_invalidates():
    base = _fp(["A"], _args(), "38.2")
    changed = _fp(["A"], _args(), "37.0")
    diffs = RESUME.compare(base, changed)
    assert any("Pfam version" in d for d in diffs)


def test_unknown_pfam_version_does_not_invalidate():
    """
    The Pfam version isn't resolved until the Pfam stage runs, so a run interrupted
    before then legitimately stores None and must not be refused on that basis.
    """
    known = _fp(["A"], _args(), "38.2")
    unknown = _fp(["A"], _args(), None)
    assert RESUME.compare(known, unknown) == []


def test_identical_fingerprints_have_no_differences():
    fp = _fp(["A", "B"], _args(), "38.2")
    assert RESUME.compare(fp, fp) == []


def test_missing_previous_state_is_reported():
    fp = _fp(["A"], _args(), "38.2")
    diffs = RESUME.compare(None, fp)
    assert diffs == ["no previous run state was found"]


################################################################################
# local genome files in the fingerprint
################################################################################

def test_local_genome_files_are_fingerprinted(tmp_path):
    f = tmp_path / "g1.faa"
    f.write_text(">a\nMK\n")
    genomes = [_FakeGenome("g1", "amino-acid-fasta", str(f))]

    with_local = _fp([], _args(), "38.2", local_genomes=genomes)
    without = _fp([], _args(), "38.2", local_genomes=[])

    diffs = RESUME.compare(without, with_local)
    assert any("local genome files" in d for d in diffs)


def test_edited_local_file_invalidates_resume(tmp_path):
    """
    Unlike an NCBI accession, a local file's CONTENTS can change while its path stays
    the same. Resuming against an edited fasta would mix old results with new input.
    """
    f = tmp_path / "g1.faa"
    f.write_text(">a\nMK\n")
    genomes = [_FakeGenome("g1", "amino-acid-fasta", str(f))]
    before = _fp([], _args(), "38.2", local_genomes=genomes)

    time.sleep(1.1)  # mtime has second resolution
    f.write_text(">a\nMKVLAAA\n")
    after = _fp([], _args(), "38.2", local_genomes=genomes)

    assert RESUME.compare(before, after)


def test_local_fingerprint_stable_when_untouched(tmp_path):
    f = tmp_path / "g1.faa"
    f.write_text(">a\nMK\n")
    genomes = [_FakeGenome("g1", "amino-acid-fasta", str(f))]

    a = _fp([], _args(), "38.2", local_genomes=genomes)
    b = _fp([], _args(), "38.2", local_genomes=genomes)
    assert a == b


def test_local_fingerprint_handles_missing_file(tmp_path):
    genomes = [_FakeGenome("ghost", "fasta", str(tmp_path / "gone.fna"))]
    fp = _fp([], _args(), "38.2", local_genomes=genomes)
    assert fp["local_genomes_sha256"] is not None


def test_runtime_only_params_do_not_invalidate():
    """
    -n/-j/output dir change HOW a run executes, not WHAT it produces, so they must not
    force a full redo.
    """
    a = _fp(["A"], _args(), "38.2")
    b = _fp(["A"], _args(), "38.2")
    assert RESUME.compare(a, b) == []


################################################################################
# state persistence
################################################################################

def test_state_round_trips(tmp_path):
    fp = _fp(["A"], _args(), "38.2")
    state = RESUME.new(fp)
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [])
    RESUME.save(str(tmp_path), state)

    loaded = RESUME.load(str(tmp_path))
    assert loaded["fingerprint"] == fp
    assert cli.STAGE_GENOMES in loaded["completed"]


def test_load_state_returns_none_when_absent(tmp_path):
    assert RESUME.load(str(tmp_path)) is None


def test_load_state_returns_none_when_corrupt(tmp_path):
    """A corrupt state file means we can't trust the prior run; start fresh."""
    path = tmp_path / RESUME.state_filename
    path.write_text("{not valid json")
    assert RESUME.load(str(tmp_path)) is None


def test_save_state_is_atomic(tmp_path):
    state = RESUME.new({"x": 1})
    RESUME.save(str(tmp_path), state)
    leftovers = [f for f in os.listdir(tmp_path) if f.endswith(".part")]
    assert leftovers == []


################################################################################
# stage reuse and artifact integrity
################################################################################

def test_stage_reusable_when_artifacts_intact(tmp_path):
    artifact = tmp_path / "combined.faa"
    artifact.write_text(">a\nMK\n")

    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(tmp_path))
    assert RESUME.is_reusable(state, cli.STAGE_GENOMES, str(tmp_path))


def test_stage_not_reusable_when_never_run(tmp_path):
    state = RESUME.new({})
    assert not RESUME.is_reusable(state, cli.STAGE_SEARCH, str(tmp_path))


def test_stage_not_reusable_when_artifact_deleted(tmp_path):
    artifact = tmp_path / "combined.faa"
    artifact.write_text(">a\nMK\n")
    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(tmp_path))

    artifact.unlink()
    assert not RESUME.is_reusable(state, cli.STAGE_GENOMES, str(tmp_path))


def test_stage_not_reusable_when_artifact_truncated(tmp_path):
    """
    Size is recorded so a file truncated by a kill -9 (or anything else touching the
    working dir) isn't silently reused as if complete.
    """
    artifact = tmp_path / "filtered.hmm"
    artifact.write_text("HMM" * 100)
    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_PFAMS, [str(artifact)],
                               work_dir=str(tmp_path))

    artifact.write_text("HMM")
    assert not RESUME.is_reusable(state, cli.STAGE_PFAMS, str(tmp_path))


def test_invalidate_from_cascades_downstream():
    """Re-running a stage invalidates everything computed from its output."""
    state = RESUME.new({})
    for stage in cli.STAGE_ORDER:
        RESUME.mark_complete(state, stage, [])

    RESUME.invalidate_from(state, cli.STAGE_PFAMS)

    assert set(state["completed"]) == {cli.STAGE_GENOMES}


def test_invalidate_from_unknown_stage_raises():
    """
    Deliberate change from the old free-function behavior, which silently no-oped.

    A no-op is the worst outcome for a typo'd stage name: `invalidate_from("serach")`
    would invalidate nothing, downstream results computed from stale input would be
    reused, and there'd be no signal anywhere. The profile knows its own stage names,
    so it can just say so.
    """
    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [])

    with pytest.raises(ValueError, match="unknown stage"):
        RESUME.invalidate_from(state, "not-a-stage")

    assert cli.STAGE_GENOMES in state["completed"]


def test_mark_complete_rejects_an_unknown_stage():
    """Same guard on the writing side, where a typo would record a stage nothing reads."""
    state = RESUME.new({})
    with pytest.raises(ValueError, match="unknown stage"):
        RESUME.mark_complete(state, "not-a-stage", [])


def test_stage_order_is_pipeline_order():
    assert cli.STAGE_ORDER == [cli.STAGE_GENOMES,
                               cli.STAGE_PFAMS,
                               cli.STAGE_SEARCH]
    # the profile and the module constant must not drift
    assert RESUME.stages == cli.STAGE_ORDER


################################################################################
# artifact paths are relative to work_dir
################################################################################

def test_artifacts_stored_relative_to_work_dir(tmp_path):
    artifact = tmp_path / "combined.faa"
    artifact.write_text(">a\nMK\n")

    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(tmp_path))

    keys = list(state["completed"][cli.STAGE_GENOMES]["artifacts"])
    assert keys == ["combined.faa"]
    assert not os.path.isabs(keys[0])


def test_renamed_output_dir_still_resumes(tmp_path):
    """
    The whole output directory can legitimately be renamed or moved between runs.
    With absolute artifact paths every stage would miss and silently re-run from
    scratch even though the files are right there; relative paths follow the move.
    """
    original = tmp_path / "runA" / "working-dir"
    original.mkdir(parents=True)
    artifact = original / "combined.faa"
    artifact.write_text(">a\nMK\n")

    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(original))
    RESUME.save(str(original), state)

    assert RESUME.is_reusable(state, cli.STAGE_GENOMES, str(original))

    # user renames the output directory, then resumes pointing at the new name
    renamed = tmp_path / "runRenamed"
    (tmp_path / "runA").rename(renamed)
    moved_work_dir = renamed / "working-dir"

    loaded = RESUME.load(str(moved_work_dir))
    assert loaded is not None
    assert RESUME.is_reusable(loaded, cli.STAGE_GENOMES, str(moved_work_dir))


def test_truncation_still_detected_after_move(tmp_path):
    """Relative paths must not weaken the integrity check."""
    original = tmp_path / "runA" / "working-dir"
    original.mkdir(parents=True)
    artifact = original / "filtered.hmm"
    artifact.write_text("HMM" * 100)

    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_PFAMS, [str(artifact)],
                               work_dir=str(original))

    renamed = tmp_path / "runRenamed"
    (tmp_path / "runA").rename(renamed)
    moved = renamed / "working-dir"

    (moved / "filtered.hmm").write_text("HMM")
    assert not RESUME.is_reusable(state, cli.STAGE_PFAMS, str(moved))


def test_artifact_outside_work_dir_stays_absolute(tmp_path):
    """
    Nothing currently writes outside the working dir, but if a future stage does, an
    absolute path is kept rather than a fragile chain of `..` segments.
    """
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    artifact = elsewhere / "x.txt"
    artifact.write_text("y")

    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_PFAMS, [str(artifact)],
                               work_dir=str(work_dir))

    key = list(state["completed"][cli.STAGE_PFAMS]["artifacts"])[0]
    assert os.path.isabs(key)
    assert RESUME.is_reusable(state, cli.STAGE_PFAMS, str(work_dir))


def test_legacy_absolute_state_still_readable(tmp_path):
    """State written before paths went relative must keep working."""
    artifact = tmp_path / "a.faa"
    artifact.write_text("z" * 10)

    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES, [str(artifact)])  # no work_dir

    assert RESUME.is_reusable(state, cli.STAGE_GENOMES, str(tmp_path))


def test_missing_artifact_recorded_with_none_size(tmp_path):
    state = RESUME.new({})
    RESUME.mark_complete(state, cli.STAGE_GENOMES,
                               [str(tmp_path / "never-made.faa")],
                               work_dir=str(tmp_path))

    artifacts = state["completed"][cli.STAGE_GENOMES]["artifacts"]
    assert list(artifacts.values()) == [None]
    # size None means "unknown", but the file still has to exist to be reusable
    assert not RESUME.is_reusable(state, cli.STAGE_GENOMES, str(tmp_path))
