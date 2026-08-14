"""
What a stage reports must depend only on that stage.

The bug these specify: `removed` was a bare boolean, so every report counted every
removal, and program order was the only thing keeping earlier stages from seeing later
ones. A resume loads all of run-data.json at once, so that ordering disappears -- an
NCBI-stage report started counting SCG-hit-filter removals, and the processing
summary announced that genomes "failed processing" when they'd been dropped several
stages later, then bailed out on that basis three stages before the fresh run did.

The specification is: a report reads only stage-scoped state, so its output is the same
whether the removals downstream of it have happened yet or not. That's what
`TestTranscriptEquivalence` pins down -- it renders each stage's report against a
"fresh" run (only that stage's removals recorded so far) and against a "resumed" one
(the finished run's full state) and requires them to match.
"""

import io
import contextlib

import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData, SCGset
from gtotree.utils.misc.stages import (GenomeRemovalStage, SCGRemovalStage,
                                  PipelineStage, UnknownStage)
from gtotree.utils.misc import messaging


def _genome(gid, source="ncbi-accession"):
    gd = GenomeData.from_acc(gid) if source == "ncbi-accession" \
        else GenomeData.from_path(f"/in/{gid}.fa", source)
    gd.id = gid
    gd.processing_done = True
    gd.hmm_search_done = True
    return gd


def _run_data(n_accs=5, n_scgs=16):
    rd = RunData()
    rd.run_files_dir_rel = "gtotree-output/run-files"
    rd.ncbi_accs = [_genome(f"GCF_{i}") for i in range(n_accs)]
    rd.SCG_targets = [SCGset.from_id(f"SCG{i}") for i in range(n_scgs)]
    rd.update_all_input_genomes()
    return rd


def _render(report_fn, run_data):
    """Capture a report's output, with the exit path neutralized."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        report_fn(run_data)
    return buf.getvalue()


def _searched(report_fn):
    """The summary as the main run calls it, with the SCG search folded in."""
    return lambda run_data: report_fn(run_data, searched=True)


@pytest.fixture(autouse=True)
def _no_early_exit(monkeypatch):
    """
    The `< 4 genomes remaining` guards call sys.exit via report_too_few_genomes.
    Record the call instead so a test can assert on whether a stage decided to bail.
    """
    calls = []
    monkeypatch.setattr(messaging, "report_too_few_genomes",
                        lambda run_data: calls.append(run_data))
    monkeypatch.setattr(messaging, "report_no_SCGs_remaining",
                        lambda run_data: calls.append(run_data))
    return calls


class TestStageScopedCounts:
    """The accessors themselves, independent of any message wording."""

    def test_a_later_removal_is_not_counted_by_an_earlier_stage(self):
        rd = _run_data()
        rd.ncbi_accs[0].mark_removed("too few unique SCG hits",
                                     GenomeRemovalStage.SCG_HIT_FILTER)

        assert rd.get_removed_ncbi_accs() == []
        assert rd.get_genomes_removed_during_processing() == []
        assert len(rd.genomes_removed_at(GenomeRemovalStage.SCG_HIT_FILTER)) == 1

    def test_alive_through_ignores_removals_from_later_stages(self):
        rd = _run_data()
        rd.ncbi_accs[0].mark_removed("acc download failed",
                                     GenomeRemovalStage.NCBI_DOWNLOAD)
        rd.ncbi_accs[1].mark_removed("HMM search failed",
                                     GenomeRemovalStage.HMM_SEARCH)
        rd.ncbi_accs[2].mark_removed("too few unique SCG hits",
                                     GenomeRemovalStage.SCG_HIT_FILTER)

        # as of the end of processing only the download failure had happened
        assert len(rd.genomes_alive_through(
            GenomeRemovalStage.AMINO_ACID_PREP)) == 4
        assert len(rd.genomes_alive_through(GenomeRemovalStage.HMM_SEARCH)) == 3
        assert len(rd.genomes_alive_through(GenomeRemovalStage.SCG_HIT_FILTER)) == 2

    def test_source_scoping_keeps_one_input_types_failures_out_of_anothers(self):
        rd = RunData()
        acc = _genome("GCF_1")
        acc.mark_removed("acc download failed", GenomeRemovalStage.NCBI_DOWNLOAD)
        gb = _genome("gb1", source="genbank-file")
        gb.mark_removed("genbank-file processing failed",
                        GenomeRemovalStage.GENBANK_PREP)
        rd.ncbi_accs, rd.genbank_files = [acc], [gb]
        rd.update_all_input_genomes()

        assert rd.get_removed_ncbi_accs() == ["GCF_1"]
        assert rd.get_failed_genbank_ids() == ["gb1"]

    def test_removed_is_derived_and_cannot_drift_from_the_stage(self):
        gd = _genome("GCF_1")
        assert gd.removed is False
        gd.mark_removed("HMM search failed", GenomeRemovalStage.HMM_SEARCH)
        assert gd.removed is True
        with pytest.raises(AttributeError):
            gd.removed = False

    def test_an_unknown_removal_stage_raises(self):
        gd = _genome("GCF_1")
        with pytest.raises(UnknownStage):
            gd.mark_removed("something", "made-up-stage")
        assert gd.removed is False


class TestTranscriptEquivalence:
    """
    Each stage's report must read the same on a resumed run as on a fresh one.

    "Fresh" here is the state that stage would have seen when the run first reached it;
    "resumed" is the finished run's state, which is what read_run_data hands back.
    """

    def _fresh_and_resumed(self, apply_up_to, apply_rest):
        fresh = _run_data()
        apply_up_to(fresh)

        resumed = _run_data()
        apply_up_to(resumed)
        apply_rest(resumed)

        return fresh, resumed

    def test_ncbi_update_is_unchanged_by_later_removals(self):
        fresh, resumed = self._fresh_and_resumed(
            lambda rd: None,
            lambda rd: [rd.ncbi_accs[i].mark_removed(
                "too few unique SCG hits", GenomeRemovalStage.SCG_HIT_FILTER)
                for i in (0, 1)])

        out = _render(messaging.report_ncbi_update, fresh)
        assert out == _render(messaging.report_ncbi_update, resumed)
        # and it's the all-succeeded wording, not the itemized-failure one
        assert "All 5 input accessions were successfully downloaded" in out
        assert "Of the input genomes provided as" not in out

    def test_processing_summary_is_unchanged_by_later_removals(self, _no_early_exit):
        fresh, resumed = self._fresh_and_resumed(
            lambda rd: None,
            lambda rd: [rd.ncbi_accs[i].mark_removed(
                "too few unique SCG hits", GenomeRemovalStage.SCG_HIT_FILTER)
                for i in (0, 1)])

        out = _render(messaging.report_genome_processing_update, fresh)
        assert out == _render(messaging.report_genome_processing_update, resumed)
        assert "All 5 input genomes were successfully processed" in out
        assert "failed processing" not in out

    def test_processing_summary_does_not_bail_out_on_later_removals(self, _no_early_exit):
        """
        The regression that ended the resumed run early: 5 inputs, 2 dropped much later
        for too few SCG hits, and this stage concluded only 3 remained and quit.
        """
        rd = _run_data()
        for i in (0, 1):
            rd.ncbi_accs[i].mark_removed("too few unique SCG hits",
                                         GenomeRemovalStage.SCG_HIT_FILTER)

        _render(messaging.report_genome_processing_update, rd)

        assert _no_early_exit == [], "the run should not have exited at processing"

    def test_processing_summary_still_reports_real_processing_failures(self, _no_early_exit):
        rd = _run_data()
        rd.ncbi_accs[0].mark_removed("acc download failed",
                                     GenomeRemovalStage.NCBI_DOWNLOAD)

        out = _render(messaging.report_genome_processing_update, rd)

        assert "1 failed processing" in out
        assert "4 of the input 5" in out

    def test_searched_summary_is_unchanged_by_later_removals(self, _no_early_exit):
        fresh, resumed = self._fresh_and_resumed(
            lambda rd: rd.ncbi_accs[4].mark_removed(
                "HMM search failed", GenomeRemovalStage.HMM_SEARCH),
            lambda rd: [rd.ncbi_accs[i].mark_removed(
                "too few unique SCG hits", GenomeRemovalStage.SCG_HIT_FILTER)
                for i in (0, 1)])

        out = _render(_searched(messaging.report_genome_processing_update), fresh)
        assert out == _render(_searched(messaging.report_genome_processing_update),
                              resumed)
        assert "1 failed the target-gene search" in out
        assert "4 of the input 5" in out

    def test_genome_filtering_update_reports_only_its_own_stage(self, _no_early_exit):
        rd = _run_data()
        rd.ncbi_accs[4].mark_removed("HMM search failed",
                                     GenomeRemovalStage.HMM_SEARCH)
        for i in (0, 1):
            rd.ncbi_accs[i].mark_removed("too few unique SCG hits",
                                         GenomeRemovalStage.SCG_HIT_FILTER)

        out = _render(messaging.report_genome_filtering_update, rd)

        assert "2 genome(s) removed due to having too few unique SCG hits" in out

    def test_scg_set_filtering_update_is_unchanged_by_alignment_failures(self, _no_early_exit):
        fresh, resumed = self._fresh_and_resumed(
            lambda rd: [rd.SCG_targets[i].mark_removed(
                "too few genomes with hits (1 < 3 required)",
                SCGRemovalStage.GENE_FILTER) for i in (0, 1, 2)],
            lambda rd: rd.SCG_targets[5].mark_removed(
                "alignment failed", SCGRemovalStage.ALIGNMENT))

        out = _render(messaging.report_SCG_set_filtering_update, fresh)
        assert out == _render(messaging.report_SCG_set_filtering_update, resumed)
        assert "3 had no hits or were filtered out" in out
        assert "remaining 13" in out


class TestPipelineStageCompletion:
    """The completed_stages map that replaced the loose booleans."""

    def test_marking_and_checking(self):
        rd = RunData()
        assert not rd.stage_is_complete(PipelineStage.ALIGN_SCG_SETS)
        rd.mark_stage_complete(PipelineStage.ALIGN_SCG_SETS)
        assert rd.stage_is_complete(PipelineStage.ALIGN_SCG_SETS)

    def test_invalidating_a_stage_drops_everything_downstream(self):
        rd = RunData()
        for stage in (PipelineStage.PROCESS_GENOMES, PipelineStage.FILTER_GENES,
                      PipelineStage.FILTER_GENOMES, PipelineStage.ALIGN_SCG_SETS):
            rd.mark_stage_complete(stage)

        rd.invalidate_stages_from(PipelineStage.FILTER_GENOMES)

        assert rd.stage_is_complete(PipelineStage.PROCESS_GENOMES)
        assert rd.stage_is_complete(PipelineStage.FILTER_GENES)
        assert not rd.stage_is_complete(PipelineStage.FILTER_GENOMES)
        assert not rd.stage_is_complete(PipelineStage.ALIGN_SCG_SETS)

    def test_an_unknown_pipeline_stage_raises(self):
        rd = RunData()
        with pytest.raises(UnknownStage):
            rd.mark_stage_complete("align-everything")


class TestSerializationRoundTrip:
    def test_stage_attribution_and_completion_survive_run_data_json(self, tmp_path):
        from gtotree.utils.misc.general import write_run_data, read_run_data

        rd = _run_data()
        rd.run_data_path = str(tmp_path / "run-data.json")
        rd.ncbi_accs[0].mark_removed("too few unique SCG hits",
                                     GenomeRemovalStage.SCG_HIT_FILTER)
        rd.SCG_targets[0].mark_removed("alignment failed", SCGRemovalStage.ALIGNMENT)
        rd.mark_stage_complete(PipelineStage.FILTER_GENOMES)
        write_run_data(rd)

        back = read_run_data(rd.run_data_path)

        assert back.ncbi_accs[0].removed_at == GenomeRemovalStage.SCG_HIT_FILTER
        assert back.ncbi_accs[0].removed is True
        assert back.SCG_targets[0].removed_at == SCGRemovalStage.ALIGNMENT
        assert back.stage_is_complete(PipelineStage.FILTER_GENOMES)
        assert back.get_removed_ncbi_accs() == []

    def test_a_run_data_json_from_an_older_format_still_loads(self, tmp_path):
        """
        So the fingerprint's state_version can refuse the resume with an explanation,
        rather than the load dying first and reporting the file as corrupt.
        """
        import json
        from gtotree.utils.misc.general import read_run_data

        path = tmp_path / "run-data.json"
        path.write_text(json.dumps({
            "ncbi_accs": [{"id": "GCF_1", "source": "ncbi-accession", "full_path": None,
                           "provided_path": None, "basename": "GCF_1",
                           "removed": True, "processing_failed": True}],
            "run_complete": True,
            "SCG_hits_filtered": True,
            "fingerprint": {"state_version": 1},
        }))

        back = read_run_data(str(path))

        assert back is not None
        assert back.fingerprint["state_version"] == 1
        assert not back.stage_is_complete(PipelineStage.FINALIZE)


class TestAlignmentFailureReporting:
    """
    Alignment-stage drops are error paths, not filtering: the only ways an SCG-set is
    removed there are muscle exiting non-zero, trimal exiting non-zero, or an
    unexpected error in the worker. They used to be reported nowhere -- the set just
    disappeared from the concatenated alignment and the partitions file.
    """

    def test_a_clean_alignment_stage_says_so(self, _no_early_exit):
        rd = _run_data(n_scgs=13)
        out = _render(messaging.report_SCG_alignment_update, rd)

        assert "All 13 SCG-sets were successfully aligned and trimmed" in out
        assert _no_early_exit == []

    def test_failures_are_reported_and_named_as_failures(self, _no_early_exit):
        rd = _run_data(n_scgs=13)
        rd.logs_dir_rel = "gtotree-output/logs"
        for i in (0, 1):
            rd.SCG_targets[i].mark_removed("alignment failed",
                                           SCGRemovalStage.ALIGNMENT)

        out = _render(messaging.report_SCG_alignment_update, rd)

        assert "2 SCG-sets failed to align or trim" in out
        assert "SCG-info.tsv" in out
        assert "gtotree-output/logs/failed-SCG-alignments/" in out
        assert "remaining 11 target gene(s)" in out
        # not fatal -- the run finishes and the final summary states the real counts
        assert _no_early_exit == []

    def test_earlier_stage_drops_are_not_counted_as_alignment_failures(self, _no_early_exit):
        rd = _run_data(n_scgs=16)
        for i in range(3):
            rd.SCG_targets[i].mark_removed("too few genomes with hits (1 < 3 required)",
                                           SCGRemovalStage.GENE_FILTER)

        out = _render(messaging.report_SCG_alignment_update, rd)

        assert "All 13 SCG-sets were successfully aligned and trimmed" in out

    def test_losing_every_set_still_stops_the_run(self, _no_early_exit):
        rd = _run_data(n_scgs=2)
        for scg in rd.SCG_targets:
            scg.mark_removed("alignment failed", SCGRemovalStage.ALIGNMENT)

        _render(messaging.report_SCG_alignment_update, rd)

        assert len(_no_early_exit) == 1, "nothing left to concatenate is still fatal"


class TestFailedAlignmentLogCapture:
    def test_logs_for_failed_sets_are_copied_out_of_the_temp_dir(self, tmp_path):
        """
        The align/trimal logs live under tmp_dir, which is deleted at the end of a run
        unless --keep-working-dir. Without copying them the report would point at a directory that
        no longer exists by the time anyone reads it.
        """
        from gtotree.main_stages.aligning_and_preparing_SCG_sets import (
            capture_failed_alignment_logs)

        scg_dir = tmp_path / "tmp" / "found-SCG-seqs"
        scg_dir.mkdir(parents=True)
        (scg_dir / "SCG0-align.log").write_text("muscle blew up\n")
        (scg_dir / "SCG0-trimmal.log").write_text("")
        (scg_dir / "SCG1-align.log").write_text("this one was fine\n")

        rd = _run_data(n_scgs=2)
        rd.found_SCG_seqs_dir = str(scg_dir)
        rd.logs_dir = str(tmp_path / "logs")
        rd.general_ext = ".faa"
        rd.SCG_targets[0].mark_removed("alignment failed", SCGRemovalStage.ALIGNMENT)

        capture_failed_alignment_logs(rd)

        kept = sorted(p.name for p in (tmp_path / "logs" / "failed-SCG-alignments").iterdir())
        assert kept == ["SCG0-align.log", "SCG0-trimmal.log"]

    def test_nothing_is_created_when_no_set_failed(self, tmp_path):
        from gtotree.main_stages.aligning_and_preparing_SCG_sets import (
            capture_failed_alignment_logs)

        rd = _run_data(n_scgs=2)
        rd.logs_dir = str(tmp_path / "logs")

        capture_failed_alignment_logs(rd)

        assert not (tmp_path / "logs" / "failed-SCG-alignments").exists()


class TestProcessingSummaryWording:
    """
    One count for the whole input-genome funnel.

    Two reports used to share this section, so a run that lost a genome during
    processing printed "5 of the input 6 ..." and then "All 5 genomes were
    successfully searched ..." right under it -- a second, smaller total that read
    like it was the number of genomes the user handed in.
    """

    def test_the_search_is_folded_into_the_one_overall_count(self, _no_early_exit):
        rd = _run_data(n_accs=6)
        rd.ncbi_accs[0].mark_removed("acc download failed",
                                     GenomeRemovalStage.NCBI_DOWNLOAD)

        out = _render(_searched(messaging.report_genome_processing_update), rd)

        assert "1 failed processing" in out
        assert "5 of the input 6 genomes were successfully processed and searched" in out
        assert "Moving forward with those :)" in out
        # no second total underneath
        assert "successfully searched for target genes" not in out

    def test_search_failures_get_their_own_line_and_the_file_they_are_in(
            self, _no_early_exit):
        rd = _run_data(n_accs=6)
        rd.ncbi_accs[0].mark_removed("acc download failed",
                                     GenomeRemovalStage.NCBI_DOWNLOAD)
        rd.ncbi_accs[1].mark_removed("HMM search failed",
                                     GenomeRemovalStage.HMM_SEARCH)

        out = _render(_searched(messaging.report_genome_processing_update), rd)

        assert "1 failed processing" in out
        assert "1 failed the target-gene search" in out
        assert "removed-genomes.tsv" in out
        assert f"stage_removed = {GenomeRemovalStage.HMM_SEARCH}" in out
        assert "4 of the input 6 genomes were successfully processed and searched" in out

    def test_a_clean_run_says_so_once(self, _no_early_exit):
        out = _render(_searched(messaging.report_genome_processing_update),
                      _run_data(n_accs=6))

        assert "All 6 input genomes were successfully processed and searched." in out
        assert "Of all the input genomes provided" not in out

    def test_without_a_search_the_wording_does_not_claim_one(self, _no_early_exit):
        """`gtt search-pfams` / `search-kos` run this stage with do_scg False."""
        rd = _run_data(n_accs=6)
        rd.ncbi_accs[0].mark_removed("acc download failed",
                                     GenomeRemovalStage.NCBI_DOWNLOAD)

        out = _render(messaging.report_genome_processing_update, rd)

        assert "5 of the input 6 genomes were successfully processed." in out
        assert "searched" not in out

    def test_it_bails_out_on_the_post_search_count(self, _no_early_exit):
        """Losing genomes to the search can drop the run under the threshold too."""
        rd = _run_data(n_accs=6)
        for i in (0, 1, 2):
            rd.ncbi_accs[i].mark_removed("HMM search failed",
                                         GenomeRemovalStage.HMM_SEARCH)

        _render(_searched(messaging.report_genome_processing_update), rd)

        assert len(_no_early_exit) == 1


class TestSCGGenomeFilteringReporting:
    """
    What `-G` does to the SCG-sets had no report at all: `report_genome_filtering_update`
    covers genomes only. Two things were invisible as a result -- sets that lost every
    genome with a hit and left the run, and sets that survived but whose representation
    fell below the `-r` the user asked for, since `-r` is evaluated once before `-G` and
    never re-checked afterwards.
    """

    def _eroded(self, n_genomes=100, n_scgs=4, kept=5, dropped_genomes=45):
        rd = _run_data(n_accs=n_genomes, n_scgs=n_scgs)
        rd.gene_representation_cutoff = 0.5
        for gd in rd.all_input_genomes[:dropped_genomes]:
            gd.mark_removed("too few unique SCG hits",
                            GenomeRemovalStage.SCG_HIT_FILTER)
        for scg in rd.SCG_targets:
            scg.num_genomes_after_length_filtering = 50
            scg.num_genomes_after_genome_filtering = 55
        rd.SCG_targets[0].num_genomes_after_genome_filtering = kept
        return rd

    def test_a_clean_genome_filtering_stage_says_nothing(self, _no_early_exit):
        rd = self._eroded(kept=55)
        assert _render(messaging.report_SCG_genome_filtering_update, rd) == ""
        assert _no_early_exit == []

    def test_sets_that_lost_every_genome_are_reported(self, _no_early_exit):
        rd = self._eroded(kept=55)
        rd.SCG_targets[0].mark_removed("no hits remaining after genome filtering (`-G`)",
                                       SCGRemovalStage.GENOME_FILTER)

        out = _render(messaging.report_SCG_genome_filtering_update, rd)

        assert "1 SCG-set lost every genome with a hit" in out
        assert "SCG-info.tsv" in out
        assert _no_early_exit == []

    def test_erosion_below_the_r_cutoff_is_flagged_without_removing_anything(
            self, _no_early_exit):
        rd = self._eroded(kept=5)

        out = _render(messaging.report_SCG_genome_filtering_update, rd)

        assert "1 of the 4 remaining SCG-sets" in out
        assert "50%" in out
        assert "55 retained" in out
        # flagged, not filtered -- the set is still in the run
        assert len(rd.get_all_SCG_targets_remaining()) == 4
        assert _no_early_exit == []

    def test_losing_every_set_exits_here_rather_than_at_the_align_stage(
            self, _no_early_exit):
        """
        The align stage would eventually notice and exit, but by then nothing has said
        why, and SCG-info.tsv hasn't been rewritten with the genome-filter reasons.
        """
        rd = self._eroded(kept=55)
        for scg in rd.SCG_targets:
            scg.mark_removed("no hits remaining after genome filtering (`-G`)",
                             SCGRemovalStage.GENOME_FILTER)

        _render(messaging.report_SCG_genome_filtering_update, rd)

        assert _no_early_exit == [rd]

    def test_earlier_stage_drops_are_not_attributed_to_genome_filtering(
            self, _no_early_exit):
        rd = self._eroded(kept=55)
        rd.SCG_targets[0].mark_removed("no hits in any genome", SCGRemovalStage.NO_HITS)
        rd.SCG_targets[1].mark_removed("too few genomes with hits (1 < 3 required)",
                                       SCGRemovalStage.GENE_FILTER)

        out = _render(messaging.report_SCG_genome_filtering_update, rd)

        assert "lost every genome" not in out
        assert _no_early_exit == []

    def test_no_cutoff_recorded_means_no_erosion_claim(self, _no_early_exit):
        """
        run_data.gene_representation_cutoff is None for a run started before it was
        stored; the report has to stay quiet rather than divide by a missing cutoff.
        """
        rd = self._eroded(kept=5)
        rd.gene_representation_cutoff = None
        assert _render(messaging.report_SCG_genome_filtering_update, rd) == ""
