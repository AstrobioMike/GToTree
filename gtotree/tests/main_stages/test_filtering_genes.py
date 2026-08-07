"""Unit tests for gtotree/main_stages/filtering_genes.py."""


import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData, SCGset
from gtotree.utils.misc.general import required_count
from gtotree.main_stages.filtering_genes import write_out_removed_SCG_targets


def _run_data(tmp_path, n_genomes=10, scg_ids=("A", "B")):
    rd = RunData()
    rd.run_files_dir = str(tmp_path)
    rd.SCG_targets = [SCGset.from_id(s) for s in scg_ids]
    for i in range(n_genomes):
        gd = GenomeData.from_acc(f"G{i}")
        gd.preprocessing_done = True
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


class TestWriteOutRemovedSCGTargets:

    def test_writes_nothing_when_no_targets_were_removed(self, tmp_path):
        rd = _run_data(tmp_path)
        write_out_removed_SCG_targets(rd)
        assert not (tmp_path / "target-SCGs-dropped-from-analysis.tsv").exists()

    def test_records_each_removed_target_with_its_reason(self, tmp_path):
        rd = _run_data(tmp_path, scg_ids=("A", "B", "C"))
        rd.SCG_targets[1].mark_removed("too few genomes with hits (2 < 5 required)")

        write_out_removed_SCG_targets(rd)

        rows = (tmp_path / "target-SCGs-dropped-from-analysis.tsv").read_text().splitlines()
        assert rows[0] == "target_SCG\treason_removed"
        assert rows[1].split("\t")[0] == "B"
        assert "too few genomes" in rows[1]
        assert len(rows) == 2
