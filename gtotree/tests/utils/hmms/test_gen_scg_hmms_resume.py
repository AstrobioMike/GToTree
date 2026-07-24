import argparse
import json
import os
import time

import pytest # type: ignore

from gtotree.utils.hmms import gen_scg_hmms_resume as resume


def _args(**kw):
    base = dict(percent_single_copy=90, min_pfam_coverage=50.0, source="GTDB",
                wanted_ref_tax="Nitrospirota", target_rank=None, derep_rank="off",
                genbank_files=None, fasta_files=None, amino_acid_files=None)
    base.update(kw)
    return argparse.Namespace(**base)


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
    a = resume.build_fingerprint(["A", "B", "C"], _args(), "38.2")
    b = resume.build_fingerprint(["C", "A", "B"], _args(), "38.2")
    assert a == b


def test_fingerprint_ignores_duplicates():
    a = resume.build_fingerprint(["A", "B"], _args(), "38.2")
    b = resume.build_fingerprint(["A", "B", "A"], _args(), "38.2")
    assert a == b


def test_fingerprint_detects_added_genome():
    a = resume.build_fingerprint(["A", "B"], _args(), "38.2")
    b = resume.build_fingerprint(["A", "B", "C"], _args(), "38.2")
    diffs = resume.compare_fingerprints(a, b)
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
    base = resume.build_fingerprint(["A"], _args(), "38.2")
    changed = resume.build_fingerprint(["A"], _args(**kwargs), "38.2")
    diffs = resume.compare_fingerprints(base, changed)
    assert any(expected in d for d in diffs), diffs


def test_pfam_version_change_invalidates():
    base = resume.build_fingerprint(["A"], _args(), "38.2")
    changed = resume.build_fingerprint(["A"], _args(), "37.0")
    diffs = resume.compare_fingerprints(base, changed)
    assert any("Pfam version" in d for d in diffs)


def test_unknown_pfam_version_does_not_invalidate():
    """
    The Pfam version isn't resolved until the Pfam stage runs, so a run interrupted
    before then legitimately stores None and must not be refused on that basis.
    """
    known = resume.build_fingerprint(["A"], _args(), "38.2")
    unknown = resume.build_fingerprint(["A"], _args(), None)
    assert resume.compare_fingerprints(known, unknown) == []


def test_identical_fingerprints_have_no_differences():
    fp = resume.build_fingerprint(["A", "B"], _args(), "38.2")
    assert resume.compare_fingerprints(fp, fp) == []


def test_missing_previous_state_is_reported():
    fp = resume.build_fingerprint(["A"], _args(), "38.2")
    diffs = resume.compare_fingerprints(None, fp)
    assert diffs == ["no previous run state was found"]


################################################################################
# local genome files in the fingerprint
################################################################################

def test_local_genome_files_are_fingerprinted(tmp_path):
    f = tmp_path / "g1.faa"
    f.write_text(">a\nMK\n")
    genomes = [_FakeGenome("g1", "amino-acid", str(f))]

    with_local = resume.build_fingerprint([], _args(), "38.2", local_genomes=genomes)
    without = resume.build_fingerprint([], _args(), "38.2", local_genomes=[])

    diffs = resume.compare_fingerprints(without, with_local)
    assert any("local genome files" in d for d in diffs)


def test_edited_local_file_invalidates_resume(tmp_path):
    """
    Unlike an NCBI accession, a local file's CONTENTS can change while its path stays
    the same. Resuming against an edited fasta would mix old results with new input.
    """
    f = tmp_path / "g1.faa"
    f.write_text(">a\nMK\n")
    genomes = [_FakeGenome("g1", "amino-acid", str(f))]
    before = resume.build_fingerprint([], _args(), "38.2", local_genomes=genomes)

    time.sleep(1.1)  # mtime has second resolution
    f.write_text(">a\nMKVLAAA\n")
    after = resume.build_fingerprint([], _args(), "38.2", local_genomes=genomes)

    assert resume.compare_fingerprints(before, after)


def test_local_fingerprint_stable_when_untouched(tmp_path):
    f = tmp_path / "g1.faa"
    f.write_text(">a\nMK\n")
    genomes = [_FakeGenome("g1", "amino-acid", str(f))]

    a = resume.build_fingerprint([], _args(), "38.2", local_genomes=genomes)
    b = resume.build_fingerprint([], _args(), "38.2", local_genomes=genomes)
    assert a == b


def test_local_fingerprint_handles_missing_file(tmp_path):
    genomes = [_FakeGenome("ghost", "fasta", str(tmp_path / "gone.fna"))]
    fp = resume.build_fingerprint([], _args(), "38.2", local_genomes=genomes)
    assert fp["local_genomes_sha256"] is not None


def test_runtime_only_params_do_not_invalidate():
    """
    -n/-j/output dir change HOW a run executes, not WHAT it produces, so they must not
    force a full redo.
    """
    a = resume.build_fingerprint(["A"], _args(), "38.2")
    b = resume.build_fingerprint(["A"], _args(), "38.2")
    assert resume.compare_fingerprints(a, b) == []


################################################################################
# state persistence
################################################################################

def test_state_round_trips(tmp_path):
    fp = resume.build_fingerprint(["A"], _args(), "38.2")
    state = resume.new_state(fp)
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [])
    resume.save_state(str(tmp_path), state)

    loaded = resume.load_state(str(tmp_path))
    assert loaded["fingerprint"] == fp
    assert resume.STAGE_GENOMES in loaded["completed"]


def test_load_state_returns_none_when_absent(tmp_path):
    assert resume.load_state(str(tmp_path)) is None


def test_load_state_returns_none_when_corrupt(tmp_path):
    """A corrupt state file means we can't trust the prior run; start fresh."""
    path = tmp_path / resume.STATE_FILENAME
    path.write_text("{not valid json")
    assert resume.load_state(str(tmp_path)) is None


def test_save_state_is_atomic(tmp_path):
    state = resume.new_state({"x": 1})
    resume.save_state(str(tmp_path), state)
    leftovers = [f for f in os.listdir(tmp_path) if f.endswith(".part")]
    assert leftovers == []


################################################################################
# stage reuse and artifact integrity
################################################################################

def test_stage_reusable_when_artifacts_intact(tmp_path):
    artifact = tmp_path / "combined.faa"
    artifact.write_text(">a\nMK\n")

    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(tmp_path))
    assert resume.stage_is_reusable(state, resume.STAGE_GENOMES, str(tmp_path))


def test_stage_not_reusable_when_never_run(tmp_path):
    state = resume.new_state({})
    assert not resume.stage_is_reusable(state, resume.STAGE_SEARCH, str(tmp_path))


def test_stage_not_reusable_when_artifact_deleted(tmp_path):
    artifact = tmp_path / "combined.faa"
    artifact.write_text(">a\nMK\n")
    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(tmp_path))

    artifact.unlink()
    assert not resume.stage_is_reusable(state, resume.STAGE_GENOMES, str(tmp_path))


def test_stage_not_reusable_when_artifact_truncated(tmp_path):
    """
    Size is recorded so a file truncated by a kill -9 (or anything else touching the
    working dir) isn't silently reused as if complete.
    """
    artifact = tmp_path / "filtered.hmm"
    artifact.write_text("HMM" * 100)
    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_PFAMS, [str(artifact)],
                               work_dir=str(tmp_path))

    artifact.write_text("HMM")
    assert not resume.stage_is_reusable(state, resume.STAGE_PFAMS, str(tmp_path))


def test_invalidate_from_cascades_downstream():
    """Re-running a stage invalidates everything computed from its output."""
    state = resume.new_state({})
    for stage in resume.STAGE_ORDER:
        resume.mark_stage_complete(state, stage, [])

    resume.invalidate_from(state, resume.STAGE_PFAMS)

    assert set(state["completed"]) == {resume.STAGE_GENOMES}


def test_invalidate_from_unknown_stage_is_noop():
    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [])
    resume.invalidate_from(state, "not-a-stage")
    assert resume.STAGE_GENOMES in state["completed"]


def test_stage_order_is_pipeline_order():
    assert resume.STAGE_ORDER == [resume.STAGE_GENOMES,
                                  resume.STAGE_PFAMS,
                                  resume.STAGE_SEARCH]


################################################################################
# artifact paths are relative to work_dir
################################################################################

def test_artifacts_stored_relative_to_work_dir(tmp_path):
    artifact = tmp_path / "combined.faa"
    artifact.write_text(">a\nMK\n")

    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(tmp_path))

    keys = list(state["completed"][resume.STAGE_GENOMES]["artifacts"])
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

    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [str(artifact)],
                               work_dir=str(original))
    resume.save_state(str(original), state)

    assert resume.stage_is_reusable(state, resume.STAGE_GENOMES, str(original))

    # user renames the output directory, then resumes pointing at the new name
    renamed = tmp_path / "runRenamed"
    (tmp_path / "runA").rename(renamed)
    moved_work_dir = renamed / "working-dir"

    loaded = resume.load_state(str(moved_work_dir))
    assert loaded is not None
    assert resume.stage_is_reusable(loaded, resume.STAGE_GENOMES, str(moved_work_dir))


def test_truncation_still_detected_after_move(tmp_path):
    """Relative paths must not weaken the integrity check."""
    original = tmp_path / "runA" / "working-dir"
    original.mkdir(parents=True)
    artifact = original / "filtered.hmm"
    artifact.write_text("HMM" * 100)

    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_PFAMS, [str(artifact)],
                               work_dir=str(original))

    renamed = tmp_path / "runRenamed"
    (tmp_path / "runA").rename(renamed)
    moved = renamed / "working-dir"

    (moved / "filtered.hmm").write_text("HMM")
    assert not resume.stage_is_reusable(state, resume.STAGE_PFAMS, str(moved))


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

    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_PFAMS, [str(artifact)],
                               work_dir=str(work_dir))

    key = list(state["completed"][resume.STAGE_PFAMS]["artifacts"])[0]
    assert os.path.isabs(key)
    assert resume.stage_is_reusable(state, resume.STAGE_PFAMS, str(work_dir))


def test_legacy_absolute_state_still_readable(tmp_path):
    """State written before paths went relative must keep working."""
    artifact = tmp_path / "a.faa"
    artifact.write_text("z" * 10)

    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES, [str(artifact)])  # no work_dir

    assert resume.stage_is_reusable(state, resume.STAGE_GENOMES, str(tmp_path))


def test_missing_artifact_recorded_with_none_size(tmp_path):
    state = resume.new_state({})
    resume.mark_stage_complete(state, resume.STAGE_GENOMES,
                               [str(tmp_path / "never-made.faa")],
                               work_dir=str(tmp_path))

    artifacts = state["completed"][resume.STAGE_GENOMES]["artifacts"]
    assert list(artifacts.values()) == [None]
    # size None means "unknown", but the file still has to exist to be reusable
    assert not resume.stage_is_reusable(state, resume.STAGE_GENOMES, str(tmp_path))
