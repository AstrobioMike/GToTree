import os
from datetime import datetime
import pytest  # type: ignore
from gtotree.main_stages import concatenating_SCG_sets
from gtotree.main_stages.concatenating_SCG_sets import (concatenate_SCG_sets,
                                                        _concatenation_can_be_skipped)
from gtotree.utils.misc.general import GenomeData, RunData, SCGset, read_run_data
from gtotree.utils.misc.stages import PipelineStage


SCG_IDS = ("A", "B")
GENOME_IDS = ("G0", "G1")


def _run_data(tmp_path):
    out_dir = tmp_path / "output"
    run_files_dir = out_dir / "run-files"
    scg_dir = run_files_dir / "SCG-seqs"
    for d in (out_dir, run_files_dir, scg_dir):
        d.mkdir(parents=True, exist_ok=True)

    rd = RunData()
    rd.start_time = datetime.now()
    rd.output_dir = str(out_dir)
    rd.run_files_dir = str(run_files_dir)
    rd.found_SCG_seqs_dir = str(scg_dir)
    rd.run_data_path = str(run_files_dir / "run-data.json")
    rd.general_ext = ".faa"

    rd.ncbi_accs = [GenomeData.from_acc(g) for g in GENOME_IDS]
    rd.update_all_input_genomes()

    rd.SCG_targets = []
    for scg_id in SCG_IDS:
        scg = SCGset.from_id(scg_id)
        scg.aligned = True
        scg.trimmed = True
        scg.ready_for_cat = True
        rd.SCG_targets.append(scg)
        # 6 aligned columns per SCG, distinct per SCG so the join order is checkable
        (scg_dir / f"{scg_id}-final.faa").write_text(
            "".join(f">{g}\n{scg_id * 6}\n" for g in GENOME_IDS))

    return rd


def _alignment_path(rd):
    return os.path.join(rd.output_dir, "aligned-SCGs.faa")


class TestConcatenation:

    def test_writes_the_alignment_and_both_partitions_files(self, tmp_path):
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)

        assert rd.concatenated_alignment_path == _alignment_path(rd)
        assert os.path.getsize(rd.concatenated_alignment_path) > 0
        assert os.path.getsize(rd.run_files_dir + "/partitions.txt") > 0
        assert os.path.getsize(rd.run_files_dir + "/partitions.nex") > 0

    def test_records_the_alignment_length(self, tmp_path):
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)
        # 6 columns per SCG + a 5-residue spacer between the two
        assert rd.final_alignment_length == 6 + 5 + 6

    def test_persists_its_results_to_run_data(self, tmp_path):
        """
        The whole point of the write: a resume must be able to see that concatenation
        happened and where the alignment landed.
        """
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)

        back = read_run_data(rd.run_data_path)
        assert back is not None
        assert back.stage_is_complete(PipelineStage.CONCATENATE_SCG_SETS)
        assert back.concatenated_alignment_path == rd.concatenated_alignment_path
        assert back.final_alignment_length == rd.final_alignment_length

    def test_leaves_no_part_files_behind(self, tmp_path):
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)
        leftovers = [p for p in tmp_path.rglob("*.part")]
        assert leftovers == []


# ---------------------------------------------------------------------------
# resume
# ---------------------------------------------------------------------------

class TestResumeSkip:

    def test_a_finished_stage_is_skipped(self, tmp_path, monkeypatch):
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)

        def explode(*a, **kw):
            raise AssertionError("concatenation re-ran on a completed stage")

        monkeypatch.setattr(concatenating_SCG_sets, "concatenate_alignments", explode)
        concatenate_SCG_sets(rd)

    def test_a_fresh_run_is_not_skipped(self, tmp_path):
        rd = _run_data(tmp_path)
        assert _concatenation_can_be_skipped(rd) is False

    def test_a_missing_alignment_forces_a_redo(self, tmp_path):
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)
        os.remove(rd.concatenated_alignment_path)

        assert _concatenation_can_be_skipped(rd) is False

        rd = concatenate_SCG_sets(rd)
        assert os.path.getsize(rd.concatenated_alignment_path) > 0

    def test_an_empty_alignment_forces_a_redo_and_is_cleared(self, tmp_path):
        """A zero-byte alignment is a killed run, not a finished one."""
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)
        open(rd.concatenated_alignment_path, "w").close()

        assert _concatenation_can_be_skipped(rd) is False
        assert not os.path.exists(rd.concatenated_alignment_path)

    def test_skipped_when_the_alignment_has_been_moved_into_run_files(self, tmp_path):
        """
        swap_labels_in_alignment moves the concatenated alignment into run-files/
        and puts the relabeled copy at the original path. Re-concatenating then
        would leave a stray un-relabeled alignment in the output dir.
        """
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)

        moved_to = os.path.join(rd.run_files_dir,
                                os.path.basename(rd.concatenated_alignment_path))
        os.replace(rd.concatenated_alignment_path, moved_to)
        rd.mark_stage_complete(PipelineStage.UPDATE_HEADERS)

        assert _concatenation_can_be_skipped(rd) is True

        concatenate_SCG_sets(rd)
        assert not os.path.exists(rd.concatenated_alignment_path)

    def test_not_skipped_when_the_alignment_is_gone_from_both_places(self, tmp_path):
        rd = _run_data(tmp_path)
        rd = concatenate_SCG_sets(rd)
        os.remove(rd.concatenated_alignment_path)
        rd.mark_stage_complete(PipelineStage.UPDATE_HEADERS)

        assert _concatenation_can_be_skipped(rd) is False

    def test_an_unset_alignment_path_is_not_mistaken_for_a_directory(self, tmp_path):
        """
        Guard against os.path.join(run_files_dir, "") resolving to the directory
        itself, which getsize() would happily report as non-empty.
        """
        rd = _run_data(tmp_path)
        rd.mark_stage_complete(PipelineStage.CONCATENATE_SCG_SETS)
        rd.concatenated_alignment_path = ""

        assert _concatenation_can_be_skipped(rd) is False

    def test_headers_updated_alone_does_not_skip_an_unconcatenated_run(self, tmp_path):
        """The flag gates everything; file presence alone is never enough."""
        rd = _run_data(tmp_path)
        rd.mark_stage_complete(PipelineStage.UPDATE_HEADERS)
        assert _concatenation_can_be_skipped(rd) is False


# ---------------------------------------------------------------------------
# atomicity
# ---------------------------------------------------------------------------

class TestAtomicity:

    @pytest.mark.parametrize("target", ["alignment", "partitions"])
    def test_an_interrupted_write_leaves_no_usable_output(self, tmp_path, monkeypatch,
                                                          target):
        """
        The resume skip trusts these files, so a killed run must not leave a truncated
        one that looks complete.
        """
        rd = _run_data(tmp_path)

        if target == "partitions":
            monkeypatch.setattr(concatenating_SCG_sets, "gen_partitions_file",
                                _exploding)
        else:
            monkeypatch.setattr(concatenating_SCG_sets, "concatenate_alignments",
                                _exploding)

        with pytest.raises(RuntimeError):
            concatenate_SCG_sets(rd)

        monkeypatch.undo()

        assert not rd.stage_is_complete(PipelineStage.CONCATENATE_SCG_SETS)
        assert list(tmp_path.rglob("*.part")) == []
        assert _concatenation_can_be_skipped(rd) is False


def _exploding(*args, **kwargs):
    raise RuntimeError("killed mid-stage")
