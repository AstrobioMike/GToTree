"""Unit tests for gtotree/main_stages/filtering_genes.py."""


import pytest  # type: ignore

from gtotree.main_stages import filtering_genes
from gtotree.main_stages.filtering_genes import filter_genes
from gtotree.utils.misc.general import GenomeData, RunData, SCGset
from gtotree.utils.misc.general import required_count
from gtotree.utils.misc.summary_info import (write_SCG_info_table,
                                             SCG_INFO_COLUMNS,
                                             SCG_sets_below_representation_cutoff)
from gtotree.utils.misc.stages import (GenomeRemovalStage, PipelineStage,
                                  SCGRemovalStage)


class _Args:
    gene_representation_cutoff = 0.5
    num_jobs = 1


def _run_data(tmp_path, n_genomes=10, scg_ids=("A", "B")):
    rd = RunData()
    rd.run_files_dir = str(tmp_path)
    rd.SCG_targets = [SCGset.from_id(s) for s in scg_ids]
    for i in range(n_genomes):
        gd = GenomeData.from_acc(f"G{i}")
        gd.processing_done = True
        gd.hmm_search_done = True
        rd.ncbi_accs.append(gd)
    rd.update_all_input_genomes()
    return rd


class TestMinGenomesRequired:
    """
    The threshold is `int(n * cutoff + 0.5)` -- half-up rounding, chosen deliberately
    over round(), which rounds half to even and would make the threshold depend on
    whether the genome count happens to be odd.
    """

    @pytest.mark.parametrize("n,cutoff,expected", [
        (10, 0.5, 5),
        (100, 0.1, 10),
        (3, 0.5, 2),      # 1.5 rounds up, not to even
        (5, 0.5, 3),      # 2.5 rounds up, not to even
        (1, 1.0, 1),
        (10, 0.0, 0),
    ])
    def test_rounds_half_up(self, n, cutoff, expected):
        assert required_count(n, cutoff) == expected

    def test_differs_from_bankers_rounding_where_it_matters(self):
        # round(2.5) == 2 in Python; the threshold must be 3
        assert required_count(5, 0.5) == 3
        assert round(5 * 0.5) == 2


def _rows(tmp_path):
    text = (tmp_path / "SCG-info.tsv").read_text().splitlines()
    header = text[0].split("\t")
    return header, [dict(zip(header, line.split("\t"))) for line in text[1:]]


class TestWriteSCGInfoTable:
    """
    SCG-info.tsv replaced target-SCGs-dropped-from-analysis.tsv, which carried only the
    three removal columns and only for sets that had already left. The reason for the
    swap is that the interesting case isn't a removed set at all -- it's a retained one
    whose representation collapsed under `-G`, which the old file had no way to show.
    """

    def test_every_target_gets_a_row_not_just_the_removed_ones(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=("A", "B", "C"))
        rd.SCG_targets[1].mark_removed("no hits in any genome",
                                       SCGRemovalStage.NO_HITS)

        write_SCG_info_table(rd)

        header, rows = _rows(tmp_path)
        assert header == list(SCG_INFO_COLUMNS)
        assert [r["target_SCG"] for r in rows] == ["A", "B", "C"]
        assert [r["retained"] for r in rows] == ["yes", "no", "yes"]

    def test_removal_stage_and_reason_are_carried_through(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=("A", "B"))
        rd.SCG_targets[1].mark_removed("trimal failed", SCGRemovalStage.ALIGNMENT)

        write_SCG_info_table(rd)

        _header, rows = _rows(tmp_path)
        assert rows[1]["stage_removed"] == SCGRemovalStage.ALIGNMENT
        assert rows[1]["reason_removed"] == "trimal failed"
        # a retained set has no removal to report, and NA is not an empty cell
        assert rows[0]["stage_removed"] == "NA"
        assert rows[0]["reason_removed"] == "NA"

    def test_stages_that_have_not_run_yet_are_NA_not_zero(self, tmp_path):
        """
        A zero would read as "represented in no genomes", which is a real and very
        different outcome from "genome filtering hasn't happened yet".
        """
        rd = _run_data(tmp_path, scg_ids=("A",))
        rd.SCG_targets[0].num_genomes_after_copy_filtering = 8

        write_SCG_info_table(rd)

        _header, rows = _rows(tmp_path)
        assert rows[0]["num_genomes_after_copy_filtering"] == "8"
        assert rows[0]["num_genomes_after_genome_filtering"] == "NA"
        assert rows[0]["perc_genomes_after_genome_filtering"] == "NA"

    def test_the_two_percent_columns_use_their_own_denominators(self, tmp_path):
        """
        The pre columns are against the genomes alive out of the HMM search -- the same
        pool `-r` was evaluated against -- and the post column against the genomes left
        after `-G`. Sharing one denominator would make the two sides incomparable, which
        is the entire point of the table.
        """
        rd = _run_data(tmp_path, n_genomes=10, scg_ids=("A",))
        for gd in rd.all_input_genomes[:5]:
            gd.mark_removed("too few unique SCG hits",
                            GenomeRemovalStage.SCG_HIT_FILTER)

        scg = rd.SCG_targets[0]
        scg.num_genomes_after_copy_filtering = 10
        scg.num_genomes_after_length_filtering = 10
        scg.num_genomes_after_genome_filtering = 5

        write_SCG_info_table(rd)

        _header, rows = _rows(tmp_path)
        # 10 of the 10 that were searched
        assert rows[0]["perc_genomes_after_copy_filtering"] == "100.0"
        assert rows[0]["perc_genomes_after_length_filtering"] == "100.0"
        # 5 of the 5 retained -- not 5 of 10
        assert rows[0]["perc_genomes_after_genome_filtering"] == "100.0"

    def test_erosion_under_genome_filtering_is_visible(self, tmp_path):
        """
        The case that motivated the table: a set that passed `-r` against the pre-`-G`
        pool and then lost most of its genomes to `-G`.
        """
        rd = _run_data(tmp_path, n_genomes=100, scg_ids=("A",))
        for gd in rd.all_input_genomes[:45]:
            gd.mark_removed("too few unique SCG hits",
                            GenomeRemovalStage.SCG_HIT_FILTER)

        scg = rd.SCG_targets[0]
        scg.num_genomes_after_length_filtering = 50
        scg.num_genomes_after_genome_filtering = 5

        write_SCG_info_table(rd)

        _header, rows = _rows(tmp_path)
        assert rows[0]["perc_genomes_after_length_filtering"] == "50.0"   # met `-r 0.5`
        assert rows[0]["perc_genomes_after_genome_filtering"] == "9.1"    # 5 of 55
        assert rows[0]["retained"] == "yes"

    def test_rewriting_after_a_later_stage_does_not_duplicate_rows(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=("A", "B"))
        rd.SCG_targets[0].mark_removed("no hits in any genome",
                                       SCGRemovalStage.NO_HITS)
        write_SCG_info_table(rd)

        rd.SCG_targets[1].mark_removed("trimal failed", SCGRemovalStage.ALIGNMENT)
        write_SCG_info_table(rd)

        _header, rows = _rows(tmp_path)
        assert [r["target_SCG"] for r in rows] == ["A", "B"]

    def test_writes_nothing_when_there_are_no_targets_at_all(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=())
        write_SCG_info_table(rd)
        assert not (tmp_path / "SCG-info.tsv").exists()


class TestSCGSetsBelowRepresentationCutoff:
    """
    `-r` is applied once, in filter_genes, against the genomes alive at that point, and
    never re-checked after `-G` removes more. We deliberately don't re-filter -- that
    would change the gene set out from under the user -- so this is what makes the
    erosion visible instead.
    """

    def _eroded(self, tmp_path, kept_after_genome_filter):
        rd = _run_data(tmp_path, n_genomes=100, scg_ids=("A",))
        rd.gene_representation_cutoff = 0.5
        for gd in rd.all_input_genomes[:45]:
            gd.mark_removed("too few unique SCG hits",
                            GenomeRemovalStage.SCG_HIT_FILTER)
        rd.SCG_targets[0].num_genomes_after_length_filtering = 50
        rd.SCG_targets[0].num_genomes_after_genome_filtering = kept_after_genome_filter
        return rd

    def test_flags_a_set_that_fell_below_the_cutoff(self, tmp_path):
        rd = self._eroded(tmp_path, kept_after_genome_filter=5)   # 5/55 = 9%
        assert [s.id for s in SCG_sets_below_representation_cutoff(rd)] == ["A"]

    def test_leaves_a_set_still_meeting_the_cutoff_alone(self, tmp_path):
        rd = self._eroded(tmp_path, kept_after_genome_filter=40)  # 40/55 = 73%
        assert SCG_sets_below_representation_cutoff(rd) == []

    def test_removed_sets_are_not_flagged(self, tmp_path):
        """
        A set that already left has its own row and reason; reporting it again as
        under-represented would double-count it in the summary.
        """
        rd = self._eroded(tmp_path, kept_after_genome_filter=5)
        rd.SCG_targets[0].mark_removed("no hits remaining after genome filtering (`-G`)",
                                       SCGRemovalStage.GENOME_FILTER)
        assert SCG_sets_below_representation_cutoff(rd) == []

    def test_sets_that_have_not_reached_genome_filtering_are_not_flagged(self, tmp_path):
        rd = self._eroded(tmp_path, kept_after_genome_filter=5)
        rd.SCG_targets[0].num_genomes_after_genome_filtering = None
        assert SCG_sets_below_representation_cutoff(rd) == []


class TestCompletedStageIsNotRecomputed:
    """
    filter_genes marked FILTER_GENES complete but, unlike FILTER_GENOMES and
    ALIGN_SCG_SETS, had no guard against re-running -- so a resume re-derived the `-r`
    threshold from a genome pool `filter_genomes` had since shrunk. That can't currently
    give a wrong answer (the threshold only ever drops, and every surviving set already
    cleared the higher one), but it's a live dependency on a coincidence, and these pin
    the skip instead.
    """

    def _completed(self, tmp_path):
        rd = _run_data(tmp_path, n_genomes=10, scg_ids=("A", "B"))
        rd.run_data_path = str(tmp_path / "run-data.json")
        rd.found_SCG_seqs_dir = str(tmp_path)
        rd.general_ext = ".faa"
        rd.seq_length_cutoff = 0.2
        rd.gene_representation_cutoff = 0.5
        for scg in rd.SCG_targets:
            scg.gene_length_filtered = True
            scg.num_genomes_after_length_filtering = 8
            scg.num_genomes_after_copy_filtering = 8
        rd.mark_stage_complete(PipelineStage.FILTER_GENES)
        return rd

    def test_a_completed_stage_removes_nothing_further(self, tmp_path, monkeypatch):
        rd = self._completed(tmp_path)
        # half the genomes gone to `-G`, as they would be on a resume past filter_genomes
        for gd in rd.all_input_genomes[:5]:
            gd.mark_removed("too few unique SCG hits",
                            GenomeRemovalStage.SCG_HIT_FILTER)

        monkeypatch.setattr(filtering_genes, "report_processing_stage",
                            lambda *a, **kw: None)
        monkeypatch.setattr(filtering_genes, "report_message", lambda *a, **kw: None)
        monkeypatch.setattr(filtering_genes, "report_SCG_set_filtering_update",
                            lambda *a, **kw: None)

        filter_genes(_Args(), rd)

        assert [s.id for s in rd.SCG_targets if s.removed] == []

    def test_the_info_table_is_still_refreshed_on_the_skip_path(self, tmp_path,
                                                                monkeypatch):
        """
        Skipping the work isn't the same as skipping the output -- SCG-info.tsv has to
        stay current across a resume, since it's what every report points at.
        """
        rd = self._completed(tmp_path)
        monkeypatch.setattr(filtering_genes, "report_processing_stage",
                            lambda *a, **kw: None)
        monkeypatch.setattr(filtering_genes, "report_message", lambda *a, **kw: None)
        monkeypatch.setattr(filtering_genes, "report_SCG_set_filtering_update",
                            lambda *a, **kw: None)

        filter_genes(_Args(), rd)

        assert (tmp_path / "SCG-info.tsv").exists()
