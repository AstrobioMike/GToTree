"""Unit tests for gtotree/utils/seqs.py."""

import statistics

import pytest  # type: ignore
from Bio import SeqIO  # type: ignore

from gtotree.utils.misc.seqs import (_fasta_seq_lengths,
                                add_needed_gap_seqs,
                                filter_seqs_by_genome_ids,
                                filter_seqs_by_length)


def _write_fasta(path, records, wrap=None):
    with open(path, "w") as f:
        for name, seq in records:
            body = ("\n".join(seq[i:i + wrap] for i in range(0, len(seq), wrap))
                    if wrap else seq)
            f.write(f">{name}\n{body}\n")
    return str(path)


class TestFastaSeqLengths:
    """Sequence lengths, read as text without constructing SeqRecords."""

    def test_returns_one_length_per_record(self, tmp_path):
        path = _write_fasta(tmp_path / "a.faa",
                            [("g1", "A" * 10), ("g2", "C" * 250), ("g3", "D" * 3)])
        assert _fasta_seq_lengths(path) == [10, 250, 3]

    def test_sums_across_wrapped_lines(self, tmp_path):
        path = _write_fasta(tmp_path / "a.faa", [("g1", "M" * 187)], wrap=60)
        assert _fasta_seq_lengths(path) == [187]

    def test_empty_file_has_no_lengths(self, tmp_path):
        (tmp_path / "empty.faa").write_text("")
        assert _fasta_seq_lengths(str(tmp_path / "empty.faa")) == []


class TestFilterSeqsByLength:
    """Keeps records within `cutoff` of the median length."""

    def test_keeps_records_inside_the_cutoff_and_drops_those_outside(self, tmp_path):
        path = _write_fasta(tmp_path / "scg.faa", [
            ("keep1", "M" * 100), ("keep2", "M" * 105),
            ("short", "M" * 5), ("long", "M" * 900),
        ])
        kept = filter_seqs_by_length(path, 0.2)

        out = {r.id: str(r.seq) for r in
               SeqIO.parse(str(tmp_path / "scg-gene-filtered.faa"), "fasta")}
        assert kept == ["keep1", "keep2"]
        assert sorted(out) == ["keep1", "keep2"]
        assert out["keep1"] == "M" * 100

    def test_cutoff_is_relative_to_the_median(self, tmp_path):
        # median 100, cutoff 0.1 -> keep 90..110
        recs = [(f"g{i}", "M" * n) for i, n in
                enumerate([80, 95, 100, 105, 130])]
        path = _write_fasta(tmp_path / "scg.faa", recs)
        assert statistics.median([80, 95, 100, 105, 130]) == 100
        assert filter_seqs_by_length(path, 0.1) == ["g1", "g2", "g3"]

    def test_returned_ids_are_in_input_order(self, tmp_path):
        recs = [(f"g{i}", "M" * 100) for i in range(20)]
        path = _write_fasta(tmp_path / "scg.faa", recs)
        assert filter_seqs_by_length(path, 0.2) == [f"g{i}" for i in range(20)]

    def test_wrapped_input_is_measured_and_written_whole(self, tmp_path):
        path = _write_fasta(tmp_path / "scg.faa",
                            [("a", "M" * 300), ("b", "M" * 300)], wrap=60)
        assert sorted(filter_seqs_by_length(path, 0.2)) == ["a", "b"]
        out = {r.id: str(r.seq) for r in
               SeqIO.parse(str(tmp_path / "scg-gene-filtered.faa"), "fasta")}
        assert out["a"] == "M" * 300

    def test_empty_input_raises_valueerror(self, tmp_path):
        (tmp_path / "empty.faa").write_text("")
        with pytest.raises(ValueError):
            filter_seqs_by_length(str(tmp_path / "empty.faa"), 0.2)


class TestFilterSeqsByGenomeIds:

    def test_writes_everything_except_the_named_genomes(self, tmp_path):
        path = _write_fasta(tmp_path / "in.faa",
                            [("g1", "AAA"), ("g2", "CCC"), ("g3", "DDD")])
        out_path = str(tmp_path / "out.faa")

        filter_seqs_by_genome_ids(path, {"g2"}, out_path)

        assert [(r.id, str(r.seq)) for r in SeqIO.parse(out_path, "fasta")] == \
            [("g1", "AAA"), ("g3", "DDD")]

    def test_empty_removal_set_copies_everything(self, tmp_path):
        path = _write_fasta(tmp_path / "in.faa", [("g1", "AAA"), ("g2", "CCC")])
        out_path = str(tmp_path / "out.faa")
        filter_seqs_by_genome_ids(path, set(), out_path)
        assert len(list(SeqIO.parse(out_path, "fasta"))) == 2


class _RunData:
    def __init__(self, ids):
        self._ids = ids

    def get_all_remaining_input_genome_ids(self):
        return self._ids


class TestAddNeededGapSeqs:
    """Emits one record per remaining genome, gap-filling any without a hit."""

    def test_output_follows_genome_order_with_gaps_for_missing(self, tmp_path):
        # order matters: concatenation needs every per-SCG file in the same order
        path = _write_fasta(tmp_path / "trimmed.faa",
                            [("g3", "MMMM"), ("g1", "KKKK")])
        out_path = str(tmp_path / "final.faa")

        add_needed_gap_seqs(_RunData(["g1", "g2", "g3", "g4"]), path, out_path)

        assert [(r.id, str(r.seq)) for r in SeqIO.parse(out_path, "fasta")] == \
            [("g1", "KKKK"), ("g2", "----"), ("g3", "MMMM"), ("g4", "----")]

    def test_gap_row_matches_the_alignment_length(self, tmp_path):
        path = _write_fasta(tmp_path / "trimmed.faa", [("g1", "M" * 17)])
        out_path = str(tmp_path / "final.faa")

        add_needed_gap_seqs(_RunData(["g1", "g2"]), path, out_path)

        out = {r.id: str(r.seq) for r in SeqIO.parse(out_path, "fasta")}
        assert out["g2"] == "-" * 17

    def test_empty_alignment_raises_valueerror(self, tmp_path):
        (tmp_path / "empty.faa").write_text("")
        with pytest.raises(ValueError):
            add_needed_gap_seqs(_RunData(["g1"]), str(tmp_path / "empty.faa"),
                                str(tmp_path / "o.faa"))
