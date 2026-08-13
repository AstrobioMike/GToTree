"""
Unit tests for the Pfam/KO search columns in the main driver's genomes-summary-info.tsv
(gtotree/utils/misc/summary_info.py).

A full GToTree run reports nothing else about the additional searches -- they happen
inside regular genome processing, with no banner of their own -- so this table is the
only place a failed Pfam or KO search is visible. Before it, a genome whose search
failed still satisfied the counts writer's `search_done and not removed` filter and
landed in the hit-counts matrix as a row of zeros, indistinguishable from a genome that
really was searched and hit nothing.
"""

import types
import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData
from gtotree.utils.misc.stages import GenomeRemovalStage
from gtotree.utils.misc.summary_info import (generate_primary_summary_table,
                                             search_completed_value)


def _args():
    return types.SimpleNamespace(add_gtdb_tax=False, add_ncbi_tax=False)


def _run_data(tmp_path, **fields):
    rd = RunData()
    rd.output_dir = str(tmp_path)
    for key, value in fields.items():
        setattr(rd, key, value)
    return rd


def _add(rd, acc):
    gd = GenomeData.from_acc(acc)
    rd.ncbi_accs.append(gd)
    rd.update_all_input_genomes()
    return gd


def _table(tmp_path):
    lines = (tmp_path / "genomes-summary-info.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]
    return header, {row["genome_id"]: row for row in rows}


################################################################################
# when the columns appear at all
################################################################################

class TestColumnsAreConditional:
    """
    Same rule the taxonomy columns follow: without the flag every row would read NA,
    which is noise in a table people load into pandas or R.
    """

    def test_neither_column_without_the_flags(self, tmp_path):
        rd = _run_data(tmp_path)
        _add(rd, "G1")

        generate_primary_summary_table(_args(), rd)

        header, _rows = _table(tmp_path)
        assert "pfam_search_completed" not in header
        assert "ko_search_completed" not in header

    def test_only_the_pfam_column_under_dash_p(self, tmp_path):
        rd = _run_data(tmp_path, target_pfams_file="targets.txt")
        _add(rd, "G1")

        generate_primary_summary_table(_args(), rd)

        header, _rows = _table(tmp_path)
        assert "pfam_search_completed" in header
        assert "ko_search_completed" not in header

    def test_only_the_ko_column_under_dash_K(self, tmp_path):
        rd = _run_data(tmp_path, target_kos_file="targets.txt")
        _add(rd, "G1")

        generate_primary_summary_table(_args(), rd)

        header, _rows = _table(tmp_path)
        assert "ko_search_completed" in header
        assert "pfam_search_completed" not in header

    def test_both_columns_when_both_were_asked_for(self, tmp_path):
        rd = _run_data(tmp_path, target_pfams_file="p.txt", target_kos_file="k.txt")
        _add(rd, "G1")

        generate_primary_summary_table(_args(), rd)

        header, _rows = _table(tmp_path)
        assert header.index("pfam_search_completed") < header.index("ko_search_completed")
        assert header.index("ko_search_completed") < header.index("in_final_tree")


################################################################################
# what they say
################################################################################

class TestSearchCompletedValues:

    def test_a_completed_search_is_yes(self, tmp_path):
        rd = _run_data(tmp_path, target_pfams_file="p.txt")
        gd = _add(rd, "G1")
        gd.pfam_search_done = True

        generate_primary_summary_table(_args(), rd)

        _header, rows = _table(tmp_path)
        assert rows["G1"]["pfam_search_completed"] == "Yes"

    def test_a_failed_search_is_no_even_though_the_done_flag_is_set(self, tmp_path):
        """
        `mark_pfam_search_failed` sets the done flag too -- it records that the search
        was attempted. Reading the done flag alone as success is what let a failed
        search through into the counts matrix.
        """
        rd = _run_data(tmp_path, target_pfams_file="p.txt")
        gd = _add(rd, "G1")
        gd.pfam_search_done = True
        gd.pfam_search_failed = True

        generate_primary_summary_table(_args(), rd)

        _header, rows = _table(tmp_path)
        assert rows["G1"]["pfam_search_completed"] == "No"

    def test_a_genome_lost_in_preprocessing_is_na_not_no(self, tmp_path):
        """
        NA is not No. It never reached the search, so calling it a failed search blames
        this stage for a loss that happened before it.
        """
        rd = _run_data(tmp_path, target_pfams_file="p.txt")
        gd = _add(rd, "G1")
        gd.mark_removed("accession not found at NCBI",
                        GenomeRemovalStage.NCBI_LOOKUP)

        generate_primary_summary_table(_args(), rd)

        _header, rows = _table(tmp_path)
        assert rows["G1"]["pfam_search_completed"] == "NA"

    def test_the_two_columns_are_answered_independently(self, tmp_path):
        """
        Why these aren't one merged column: the searches are separate code paths
        (pyhmmer in-process vs. shelling out to exec_annotation), so a genome can fail
        one and pass the other, and one column couldn't say which.
        """
        rd = _run_data(tmp_path, target_pfams_file="p.txt", target_kos_file="k.txt")
        gd = _add(rd, "G1")
        gd.pfam_search_done = True
        gd.ko_search_done = True
        gd.ko_search_failed = True

        generate_primary_summary_table(_args(), rd)

        _header, rows = _table(tmp_path)
        assert rows["G1"]["pfam_search_completed"] == "Yes"
        assert rows["G1"]["ko_search_completed"] == "No"


class TestRemovedLaterThanTheSearch:
    """
    The fused worker runs the Pfam and KO searches *ahead* of the SCG search, so a
    genome removed after them still has a real answer to report.
    """

    @pytest.mark.parametrize("stage", [GenomeRemovalStage.HMM_SEARCH,
                                       GenomeRemovalStage.SCG_HIT_FILTER])
    def test_yes_survives_a_later_removal(self, tmp_path, stage):
        rd = _run_data(tmp_path, target_pfams_file="p.txt")
        gd = _add(rd, "G1")
        gd.pfam_search_done = True
        gd.mark_removed("dropped later", stage)

        generate_primary_summary_table(_args(), rd)

        _header, rows = _table(tmp_path)
        # the search really did finish; in_final_tree is what explains its absence
        # from the hit-counts table
        assert rows["G1"]["pfam_search_completed"] == "Yes"
        assert rows["G1"]["in_final_tree"] == "No"


################################################################################
# the helper the search subcommands share
################################################################################

def test_helper_is_the_one_the_subcommands_use():
    """
    `search_completed` in gtt search-pfams/search-kos comes from this same helper, so
    the column means the same thing in both tools.
    """
    from gtotree.utils.target_search import target_search_outputs

    assert target_search_outputs.search_completed_value is search_completed_value
