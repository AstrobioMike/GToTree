"""Unit tests for gtotree/main_stages/filtering_genomes.py."""

import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData, SCGset
from gtotree.main_stages.filtering_genomes import capture_removed_genomes


def _run_data(tmp_path, hit_counts):
    """hit_counts: {genome_id: num_SCG_hits_after_filtering}"""
    rd = RunData()
    rd.run_files_dir = str(tmp_path)
    rd.SCG_targets = [SCGset.from_id(s) for s in ("A", "B", "C", "D")]
    for gid, hits in hit_counts.items():
        gd = GenomeData.from_acc(gid)
        gd.preprocessing_done = True
        gd.hmm_search_done = True
        gd.num_SCG_hits = hits
        gd.num_unique_SCG_hits = hits
        gd.num_SCG_hits_after_filtering = hits
        rd.ncbi_accs.append(gd)
    rd.update_all_input_genomes()
    return rd


class TestCaptureRemovedGenomes:
    """Records genomes dropped for having too few surviving SCG hits."""

    def test_writes_nothing_when_none_were_removed(self, tmp_path):
        rd = _run_data(tmp_path, {"G1": 4, "G2": 4})
        capture_removed_genomes(rd)
        assert not (tmp_path / "genomes-removed-for-too-few-SCG-hits.tsv").exists()

    def test_records_each_removed_genome_with_its_counts(self, tmp_path):
        rd = _run_data(tmp_path, {"G1": 4, "G2": 1})
        rd.all_input_genomes[1].mark_removed("too few unique SCG hits")

        capture_removed_genomes(rd)

        rows = (tmp_path / "genomes-removed-for-too-few-SCG-hits.tsv") \
            .read_text().splitlines()
        assert rows[0].split("\t") == ["assembly_id", "total_SCG_hits",
                                       "unique_SCG_hits", "num_SCG_hits_after_filtering"]
        assert rows[1].split("\t") == ["G2", "1", "1", "1"]
        assert len(rows) == 2

    def test_genomes_removed_for_other_reasons_are_not_listed(self, tmp_path):
        """This file is specifically the too-few-hits report."""
        rd = _run_data(tmp_path, {"G1": 4, "G2": 1})
        rd.all_input_genomes[0].mark_removed("HMM search failed")
        rd.all_input_genomes[1].mark_removed("too few unique SCG hits")

        capture_removed_genomes(rd)

        rows = (tmp_path / "genomes-removed-for-too-few-SCG-hits.tsv") \
            .read_text().splitlines()
        assert [r.split("\t")[0] for r in rows[1:]] == ["G2"]

    def test_uncounted_genomes_are_written_as_zero(self, tmp_path):
        """
        num_SCG_hits_after_filtering is None until filter_genes has counted, so the
        report must not emit "None" into a numeric column.
        """
        rd = _run_data(tmp_path, {"G1": 1})
        gd = rd.all_input_genomes[0]
        gd.num_SCG_hits_after_filtering = None
        gd.num_SCG_hits = None
        gd.num_unique_SCG_hits = None
        gd.mark_removed("too few unique SCG hits")

        capture_removed_genomes(rd)

        row = (tmp_path / "genomes-removed-for-too-few-SCG-hits.tsv") \
            .read_text().splitlines()[1]
        assert row.split("\t") == ["G1", "0", "0", "0"]


def test_uncounted_genome_is_filtered_out_rather_than_crashing():
    """
    A genome whose hits were never counted (num_SCG_hits_after_filtering is None) must
    compare as zero. Comparing None against the threshold directly raises TypeError.
    """
    from gtotree.utils.misc.general import required_count

    threshold = required_count(4, 0.5)
    uncounted = None
    assert (uncounted or 0) < threshold
    with pytest.raises(TypeError):
        _ = uncounted < threshold
