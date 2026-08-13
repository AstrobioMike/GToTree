"""
Why an SCG-set left the run, as recorded in SCG-info.tsv.

`check_target_SCGs_have_seqs` is the single removal site for three different stages, and
it used to hard-code one reason for all of them: "no seqs found or no seqs remaining
after length-filtering". That was wrong at NO_HITS (nothing has been length-filtered
yet) and wrong at GENOME_FILTER (the cause there is `-G`, not `-c`), so `stage_removed`
was right in the table while `reason_removed` disagreed with it in two cases out of
three.

The multi-copy case was invisible on top of that. With `-B` off, `parse_hmmer_results`
discards a genome's hit when the same profile hits more than once, so a profile that
matched every genome several times produces an empty combined file and reads as "no
seqs found" -- with nothing to tell the user that `-B` would have kept it.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData, SCGset
from gtotree.utils.misc.seqs import NO_SEQS_REASONS, check_target_SCGs_have_seqs
from gtotree.utils.misc.stages import SCGRemovalStage, SCG_REMOVAL_STAGE_ORDER
from gtotree.utils.hmms.hmm_searching import no_hits_reason


def _run_data(tmp_path, scg_ids=("A", "B")):
    rd = RunData()
    rd.run_files_dir = str(tmp_path)
    rd.found_SCG_seqs_dir = str(tmp_path)
    rd.general_ext = ".faa"
    rd.SCG_targets = [SCGset.from_id(s) for s in scg_ids]
    rd.ncbi_accs = [GenomeData.from_acc(f"G{i}") for i in range(4)]
    rd.update_all_input_genomes()
    return rd


class TestStageSpecificReasons:

    @pytest.mark.parametrize("stage", [
        SCGRemovalStage.NO_HITS,
        SCGRemovalStage.GENE_FILTER,
        SCGRemovalStage.GENOME_FILTER,
    ])
    def test_each_stage_records_its_own_reason(self, tmp_path, stage):
        rd = _run_data(tmp_path, scg_ids=("A",))
        (tmp_path / "A.faa").write_text("")   # empty -> nothing usable

        check_target_SCGs_have_seqs(rd, ".faa", stage)

        assert rd.SCG_targets[0].removed_at == stage
        assert rd.SCG_targets[0].reason_removed == NO_SEQS_REASONS[stage]

    def test_the_three_reasons_are_actually_distinct(self):
        """
        The whole point of the change -- one shared string meant the stage column was
        the only thing carrying the distinction.
        """
        assert len(set(NO_SEQS_REASONS.values())) == 3

    def test_every_stage_this_check_runs_at_has_a_reason(self):
        """
        A missing key would be a KeyError at the exact moment a set is being removed,
        i.e. only on the failure path, where it's least likely to be noticed early.
        """
        for stage in NO_SEQS_REASONS:
            assert stage in SCG_REMOVAL_STAGE_ORDER

    def test_a_set_with_seqs_is_left_alone(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=("A",))
        (tmp_path / "A.faa").write_text(">G0\nMAAA\n")

        check_target_SCGs_have_seqs(rd, ".faa", SCGRemovalStage.NO_HITS)

        assert not rd.SCG_targets[0].removed

    def test_an_already_removed_set_keeps_its_original_reason(self, tmp_path):
        """
        A set dropped at NO_HITS has no files at any later extension, so without the
        skip every later check would overwrite the real reason with its own.
        """
        rd = _run_data(tmp_path, scg_ids=("A",))
        rd.SCG_targets[0].mark_removed("no hits in any genome", SCGRemovalStage.NO_HITS)

        check_target_SCGs_have_seqs(rd, "-genome-filtered.faa",
                                    SCGRemovalStage.GENOME_FILTER)

        assert rd.SCG_targets[0].removed_at == SCGRemovalStage.NO_HITS
        assert rd.SCG_targets[0].reason_removed == "no hits in any genome"

    def test_reason_fn_overrides_the_stage_default(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=("A",))
        (tmp_path / "A.faa").write_text("")

        check_target_SCGs_have_seqs(rd, ".faa", SCGRemovalStage.NO_HITS,
                                    reason_fn=lambda scg: f"custom for {scg.id}")

        assert rd.SCG_targets[0].reason_removed == "custom for A"


class TestNoHitsReason:

    def test_nothing_hit_it_at_all(self):
        scg = SCGset.from_id("A")
        scg.num_genomes_with_any_hit = 0
        assert no_hits_reason(scg, best_hit_mode=False) == "no hits in any genome"

    def test_hits_present_but_never_single_copy_points_at_best_hit_mode(self):
        scg = SCGset.from_id("A")
        scg.num_genomes_with_any_hit = 37

        reason = no_hits_reason(scg, best_hit_mode=False)

        assert "37 genomes" in reason
        assert "single copy" in reason
        assert "-B" in reason

    def test_best_hit_mode_blames_extraction_instead(self):
        """
        With `-B` on, a multi-copy hit is still taken, so hits-but-no-seqs can't be the
        multi-copy case -- pointing the user at `-B` there would be misleading.
        """
        scg = SCGset.from_id("A")
        scg.num_genomes_with_any_hit = 37

        reason = no_hits_reason(scg, best_hit_mode=True)

        assert "-B" not in reason
        assert "extracted" in reason

    def test_a_single_genome_is_not_pluralized(self):
        scg = SCGset.from_id("A")
        scg.num_genomes_with_any_hit = 1
        assert "1 genome," in no_hits_reason(scg, best_hit_mode=False)

    def test_an_unset_count_reads_as_no_hits(self):
        """
        None means the search-stage tally never ran for this set (a run started before
        the counts existed, resumed after). "No hits" is the safe reading -- it makes no
        claim the data doesn't support.
        """
        scg = SCGset.from_id("A")
        assert no_hits_reason(scg, best_hit_mode=False) == "no hits in any genome"
