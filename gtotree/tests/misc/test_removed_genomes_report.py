"""
Unit tests for the consolidated removed-genomes report in
gtotree/utils/misc/summary_info.py.

This one file replaces the seven per-stage failure files that used to be scattered
across run-files/. It's the live, mid-run record of why inputs left the run, written
from the stage-scoped `removed_at`/`reason_removed` state rather than accumulated
separately, so it can't drift from the counts the OVERALL SUMMARY reports.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData
from gtotree.utils.misc.stages import GenomeRemovalStage
from gtotree.utils.misc.summary_info import (REMOVED_GENOMES_COLUMNS,
                                             REMOVED_GENOMES_FILENAME,
                                             write_removed_genomes_report)


def _run_data(tmp_path):
    rd = RunData()
    rd.run_files_dir = str(tmp_path)
    return rd


def _add_acc(rd, acc):
    gd = GenomeData.from_acc(acc)
    rd.ncbi_accs.append(gd)
    rd.update_all_input_genomes()
    return gd


def _add_file(rd, path, source_field, source):
    gd = GenomeData.from_path(path, source)
    getattr(rd, source_field).append(gd)
    rd.update_all_input_genomes()
    return gd


def _rows(tmp_path):
    return (tmp_path / REMOVED_GENOMES_FILENAME).read_text().splitlines()


class TestWriteRemovedGenomesReport:

    def test_writes_nothing_when_nothing_was_removed(self, tmp_path):
        rd = _run_data(tmp_path)
        _add_acc(rd, "GCF_000005845.2")

        write_removed_genomes_report(rd)

        assert not (tmp_path / REMOVED_GENOMES_FILENAME).exists()

    def test_header_matches_the_declared_columns(self, tmp_path):
        rd = _run_data(tmp_path)
        _add_acc(rd, "G1").mark_removed("nope", GenomeRemovalStage.NCBI_LOOKUP)

        write_removed_genomes_report(rd)

        assert _rows(tmp_path)[0].split("\t") == list(REMOVED_GENOMES_COLUMNS)

    def test_carries_the_stage_and_the_reason(self, tmp_path):
        rd = _run_data(tmp_path)
        _add_acc(rd, "G1").mark_removed("accession not found at NCBI",
                                        GenomeRemovalStage.NCBI_LOOKUP)

        write_removed_genomes_report(rd)

        row = _rows(tmp_path)[1].split("\t")
        assert row[0] == "G1"
        assert row[3] == GenomeRemovalStage.NCBI_LOOKUP
        assert row[4] == "accession not found at NCBI"

    def test_file_inputs_report_the_path_the_user_provided(self, tmp_path):
        """
        The whole reason this file isn't just a pointer at genomes-summary-info.tsv:
        that table keys on `assembly_id`, which for a file input is the basename minus
        its extension. Two `genome_1.fa` files in different directories collapse to the
        same id there, and neither tells you which file to go fix.
        """
        rd = _run_data(tmp_path)
        gd = _add_file(rd, "some/nested/dir/genome_1.fa.gz",
                       "fasta_files", "nt-fasta-file")
        gd.mark_removed("fasta-file processing failed",
                        GenomeRemovalStage.FASTA_PREP)

        write_removed_genomes_report(rd)

        row = _rows(tmp_path)[1].split("\t")
        assert row[0] == "genome_1"
        assert row[1] == "some/nested/dir/genome_1.fa.gz"
        # the human-facing label, not the internal GenomeData.source value; every
        # table that reports an input source now goes through genome_source_label
        assert row[2] == "nucleotide-fasta"

    def test_accessions_fall_back_to_the_id_for_the_input_column(self, tmp_path):
        """`provided_path` is None for accessions, so the accession itself stands in."""
        rd = _run_data(tmp_path)
        _add_acc(rd, "GCF_999999999.1").mark_removed(
            "accession not found at NCBI", GenomeRemovalStage.NCBI_LOOKUP)

        write_removed_genomes_report(rd)

        assert _rows(tmp_path)[1].split("\t")[1] == "GCF_999999999.1"

    def test_every_source_lands_in_the_one_file(self, tmp_path):
        rd = _run_data(tmp_path)
        _add_acc(rd, "A1").mark_removed("r", GenomeRemovalStage.NCBI_LOOKUP)
        _add_file(rd, "g.gb", "genbank_files", "genbank-file").mark_removed(
            "r", GenomeRemovalStage.GENBANK_PREP)
        _add_file(rd, "f.fa", "fasta_files", "nt-fasta-file").mark_removed(
            "r", GenomeRemovalStage.FASTA_PREP)
        _add_file(rd, "a.faa", "amino_acid_files", "aa-fasta-file").mark_removed(
            "r", GenomeRemovalStage.AMINO_ACID_PREP)

        write_removed_genomes_report(rd)

        assert [r.split("\t")[0] for r in _rows(tmp_path)[1:]] == ["A1", "g", "f", "a"]

    def test_rows_come_out_in_pipeline_order(self, tmp_path):
        """
        Reading top-to-bottom should walk the run forwards, regardless of the order
        removals happened to be recorded in.
        """
        rd = _run_data(tmp_path)
        _add_acc(rd, "late").mark_removed("r", GenomeRemovalStage.SCG_HIT_FILTER)
        _add_acc(rd, "early").mark_removed("r", GenomeRemovalStage.NCBI_LOOKUP)
        _add_acc(rd, "middle").mark_removed("r", GenomeRemovalStage.HMM_SEARCH)

        write_removed_genomes_report(rd)

        assert [r.split("\t")[3] for r in _rows(tmp_path)[1:]] == [
            GenomeRemovalStage.NCBI_LOOKUP,
            GenomeRemovalStage.HMM_SEARCH,
            GenomeRemovalStage.SCG_HIT_FILTER,
        ]

    def test_genomes_still_in_the_run_are_not_listed(self, tmp_path):
        rd = _run_data(tmp_path)
        _add_acc(rd, "kept")
        _add_acc(rd, "dropped").mark_removed("r", GenomeRemovalStage.HMM_SEARCH)

        write_removed_genomes_report(rd)

        assert [r.split("\t")[0] for r in _rows(tmp_path)[1:]] == ["dropped"]

    def test_rewriting_replaces_rather_than_appends(self, tmp_path):
        """
        Called once per removing stage, so it has to be idempotent. Removals are
        monotonic (nothing clears `removed_at`), which is what makes regenerating from
        state the right move -- appending would double up rows on the second call and
        again on every resume.
        """
        rd = _run_data(tmp_path)
        gd = _add_acc(rd, "G1")
        gd.mark_removed("r", GenomeRemovalStage.NCBI_LOOKUP)

        write_removed_genomes_report(rd)
        first = _rows(tmp_path)
        write_removed_genomes_report(rd)

        assert _rows(tmp_path) == first
        assert len(first) == 2

    def test_later_removals_are_added_to_the_existing_rows(self, tmp_path):
        rd = _run_data(tmp_path)
        _add_acc(rd, "G1").mark_removed("r", GenomeRemovalStage.NCBI_LOOKUP)
        write_removed_genomes_report(rd)

        _add_acc(rd, "G2").mark_removed("r", GenomeRemovalStage.HMM_SEARCH)
        write_removed_genomes_report(rd)

        assert [r.split("\t")[0] for r in _rows(tmp_path)[1:]] == ["G1", "G2"]

    def test_a_missing_reason_is_written_as_NA(self, tmp_path):
        """`reason_removed` is typed as optional, so the column can't emit "None"."""
        rd = _run_data(tmp_path)
        gd = _add_acc(rd, "G1")
        gd.mark_removed("r", GenomeRemovalStage.NCBI_LOOKUP)
        gd.reason_removed = None

        write_removed_genomes_report(rd)

        assert _rows(tmp_path)[1].split("\t")[4] == "NA"

    def test_no_partial_file_is_left_behind_when_writing_fails(self, tmp_path):
        """Atomic-write contract: the `.part` is cleaned up rather than orphaned."""
        rd = _run_data(tmp_path)
        gd = _add_acc(rd, "G1")
        gd.mark_removed("r", GenomeRemovalStage.NCBI_LOOKUP)
        gd.id = object()  # not joinable, so the write blows up mid-flight

        with pytest.raises(TypeError):
            write_removed_genomes_report(rd)

        assert not (tmp_path / f"{REMOVED_GENOMES_FILENAME}.part").exists()
        assert not (tmp_path / REMOVED_GENOMES_FILENAME).exists()
