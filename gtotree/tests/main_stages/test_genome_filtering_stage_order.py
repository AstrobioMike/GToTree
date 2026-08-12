"""
`filter_genomes` marks its stage complete only once its own bookkeeping has landed.

It used to mark the stage complete and write run-data.json *before* calling
`capture_removed_genomes` and the GENOME_FILTER `check_target_SCGs_have_seqs`. A crash in
that window left a resume seeing the stage as done and skipping the whole block -- so the
zero-byte '-genome-filtered' files it never got to clear went straight into alignment,
where `add_needed_gap_seqs` raises on them, the worker's blanket `except BaseException`
swallows it, and the set gets recorded as `alignment` / "alignment failed". Wrong stage,
wrong reason, and no sign of the real cause anywhere.
"""

import pytest  # type: ignore

from gtotree.main_stages import filtering_genomes
from gtotree.main_stages.filtering_genomes import filter_genomes
from gtotree.utils.misc.general import GenomeData, RunData, SCGset
from gtotree.utils.misc.stages import PipelineStage, SCGRemovalStage


class _Args:
    best_hit_mode = False
    genome_hits_cutoff = 0.5
    gene_representation_cutoff = 0.5
    num_jobs = 1


def _run_data(tmp_path, scg_hits):
    """
    `scg_hits` maps SCG id -> the genome ids with a surviving gene-filtered hit.
    """
    rd = RunData()
    rd.run_files_dir = str(tmp_path)
    rd.run_files_dir_rel = str(tmp_path)
    rd.found_SCG_seqs_dir = str(tmp_path)
    rd.run_data_path = str(tmp_path / "run-data.json")
    rd.general_ext = ".faa"
    rd.gene_representation_cutoff = 0.5
    rd.SCG_targets = [SCGset.from_id(s) for s in scg_hits]

    all_genomes = sorted({g for ids in scg_hits.values() for g in ids})
    for gid in all_genomes:
        gd = GenomeData.from_acc(gid)
        gd.processing_done = True
        gd.hmm_search_done = True
        gd.num_SCG_hits_after_filtering = sum(1 for ids in scg_hits.values() if gid in ids)
        rd.ncbi_accs.append(gd)
    rd.update_all_input_genomes()

    for scg_id, ids in scg_hits.items():
        with open(tmp_path / f"{scg_id}-gene-filtered.faa", "w") as f:
            for gid in ids:
                f.write(f">{gid}\nMAAAA\n")

    return rd


@pytest.fixture(autouse=True)
def _quiet(monkeypatch):
    """Silence the reports; their content is specified in the reporting tests."""
    exits = []
    for name in ("report_processing_stage", "report_message",
                 "report_genome_filtering_update", "report_SCG_genome_filtering_update"):
        monkeypatch.setattr(filtering_genomes, name, lambda *a, **kw: None, raising=False)
    return exits


class TestBookkeepingLandsBeforeTheStageIsMarkedComplete:

    def test_a_set_that_lost_every_genome_is_recorded_at_the_genome_filter_stage(
            self, tmp_path):
        # G3/G4 hit only SCG_C, so they fall under `-G` and take SCG_C with them
        rd = _run_data(tmp_path, {
            "SCG_A": ["G1", "G2", "G3", "G4"],
            "SCG_B": ["G1", "G2", "G3", "G4"],
            "SCG_C": ["G5", "G6"],
        })
        # G5/G6 have 1 of 3 sets -> below the 50% cutoff -> removed
        filter_genomes(_Args(), rd)

        scg_c = next(s for s in rd.SCG_targets if s.id == "SCG_C")
        assert scg_c.removed_at == SCGRemovalStage.GENOME_FILTER
        assert "`-G`" in scg_c.reason_removed

    def test_the_stage_is_complete_only_after_those_removals_are_recorded(self, tmp_path):
        """
        The ordering itself: at the moment FILTER_GENOMES is marked, the removals a
        resume would otherwise skip must already be on the object being written out.
        """
        rd = _run_data(tmp_path, {
            "SCG_A": ["G1", "G2", "G3", "G4"],
            "SCG_B": ["G1", "G2", "G3", "G4"],
            "SCG_C": ["G5", "G6"],
        })

        seen = {}
        real_mark = RunData.mark_stage_complete

        def _spy(self, stage):
            if stage == PipelineStage.FILTER_GENOMES:
                seen["removed_at_mark"] = [s.id for s in self.SCG_targets if s.removed]
            return real_mark(self, stage)

        RunData.mark_stage_complete = _spy
        try:
            filter_genomes(_Args(), rd)
        finally:
            RunData.mark_stage_complete = real_mark

        assert seen["removed_at_mark"] == ["SCG_C"]

    def test_the_info_table_is_on_disk_by_the_time_the_stage_is_complete(self, tmp_path):
        rd = _run_data(tmp_path, {
            "SCG_A": ["G1", "G2", "G3", "G4"],
            "SCG_B": ["G1", "G2", "G3", "G4"],
            "SCG_C": ["G5", "G6"],
        })

        seen = {}
        real_mark = RunData.mark_stage_complete

        def _spy(self, stage):
            if stage == PipelineStage.FILTER_GENOMES:
                seen["exists"] = (tmp_path / "SCG-info.tsv").exists()
            return real_mark(self, stage)

        RunData.mark_stage_complete = _spy
        try:
            filter_genomes(_Args(), rd)
        finally:
            RunData.mark_stage_complete = real_mark

        assert seen["exists"]


class TestPostGenomeFilterCounts:

    def test_each_surviving_set_records_how_many_genomes_it_kept(self, tmp_path):
        rd = _run_data(tmp_path, {
            "SCG_A": ["G1", "G2", "G3", "G4"],
            "SCG_B": ["G1", "G2", "G3", "G4"],
            "SCG_C": ["G5", "G6"],
        })

        filter_genomes(_Args(), rd)

        by_id = {s.id: s for s in rd.SCG_targets}
        assert by_id["SCG_A"].num_genomes_with_hits_after_genome_filtering == 4
        assert by_id["SCG_C"].num_genomes_with_hits_after_genome_filtering == 0

    def test_the_no_removals_fast_path_still_records_a_count(self, tmp_path):
        """
        With nothing to remove the worker copies the file rather than filtering it, so
        the count comes from a different branch and has to mean the same thing.
        """
        rd = _run_data(tmp_path, {
            "SCG_A": ["G1", "G2", "G3", "G4"],
            "SCG_B": ["G1", "G2", "G3", "G4"],
        })

        filter_genomes(_Args(), rd)

        assert not any(g.removed for g in rd.all_input_genomes)
        assert all(s.num_genomes_with_hits_after_genome_filtering == 4
                   for s in rd.SCG_targets)
