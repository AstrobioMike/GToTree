"""Unit tests for input validation in gtotree/utils/preflight_checks.py."""

import pytest  # type: ignore

from gtotree.utils.misc.preflight_checks import (report_duplicate_entries,
                                            check_for_whitespace,
                                            check_mapping_file_problem_chars_and_fields)
from gtotree.utils.misc.general import read_single_column_file


class TestMappingFileValidation:
    """
    Rejects characters that break downstream output. '/' is allowed in column 1, which
    holds input genome paths; everything else is rejected in every column.
    """

    def _write(self, tmp_path, line):
        p = tmp_path / "map.tsv"
        p.write_text(line + "\n")
        return str(p)

    @pytest.mark.parametrize("char", [";", "(", ")", "&"])
    def test_rejects_problematic_chars_in_the_label_column(self, tmp_path, char):
        assert check_mapping_file_problem_chars_and_fields(
            self._write(tmp_path, f"genome1\tbad{char}label"))

    def test_rejects_problematic_chars_in_the_first_column(self, tmp_path):
        assert check_mapping_file_problem_chars_and_fields(
            self._write(tmp_path, "gen;ome1\tgoodlabel"))

    def test_rejects_slash_outside_the_first_column(self, tmp_path):
        assert check_mapping_file_problem_chars_and_fields(
            self._write(tmp_path, "genome1\tbad/label"))

    def test_allows_slash_in_the_first_column(self, tmp_path):
        assert check_mapping_file_problem_chars_and_fields(
            self._write(tmp_path, "path/to/genome1.faa\tgoodlabel")) == []

    def test_accepts_a_clean_file(self, tmp_path):
        assert check_mapping_file_problem_chars_and_fields(
            self._write(tmp_path, "genome1\tNice Label")) == []


def test_duplicate_entries_are_dropped_in_input_order_without_touching_the_file(tmp_path,
                                                                                capsys):
    """
    De-duplication must be order-preserving; genome order flows into the alignment.

    It also has to happen in memory. This used to write a `<path>-unique` sibling next
    to the user's input and hand back that new path, which then flowed into the resume
    fingerprint -- a leftover from when these files were parsed by shell tooling.
    """

    src = tmp_path / "accs.txt"
    ordered = [f"GCF_{i:09d}.1" for i in range(40)]
    src.write_text("\n".join(ordered + [ordered[3], ordered[17]]) + "\n")
    before = src.read_text()

    report_duplicate_entries(str(src), "-a")

    assert read_single_column_file(str(src)) == ordered
    # the user's file is untouched, and no sibling was created
    assert src.read_text() == before
    assert [p.name for p in tmp_path.iterdir()] == ["accs.txt"]
    # ...but they were told about it
    out = capsys.readouterr().out
    assert "2 repeated entries" in out
    assert ordered[3] in out


def test_a_clean_file_produces_no_duplicate_notice(tmp_path, capsys):
    src = tmp_path / "accs.txt"
    src.write_text("GCF_000000001.1\nGCF_000000002.1\n")

    report_duplicate_entries(str(src), "-a")

    assert capsys.readouterr().out == ""


def test_crlf_input_needs_no_rewrite_to_be_read_correctly(tmp_path):
    """
    The `-unix` rewrite is gone; `.strip()` in the shared reader absorbs the `\\r`.
    """
    src = tmp_path / "accs.txt"
    src.write_bytes(b"GCF_000000001.1\r\nGCF_000000002.1\r\n")

    assert read_single_column_file(str(src)) == ["GCF_000000001.1", "GCF_000000002.1"]
    assert [p.name for p in tmp_path.iterdir()] == ["accs.txt"]


def test_blank_lines_are_dropped_by_the_shared_reader(tmp_path):
    """
    Previously only cleaned as a side effect of the duplicate rewrite, so a file with
    a blank line but no duplicates kept the blank and became an empty genome id.
    """
    src = tmp_path / "accs.txt"
    src.write_text("GCF_000000001.1\n\n   \nGCF_000000002.1\n")

    assert read_single_column_file(str(src)) == ["GCF_000000001.1", "GCF_000000002.1"]


def test_the_mapping_file_still_gets_its_crlf_rewrite(tmp_path):
    """
    `-m` is the one file that still needs it: pd.read_csv leaves the `\\r` inside the
    last field rather than at the end of a line, silently corrupting output labels.
    """
    from gtotree.utils.misc.preflight_checks import check_line_endings

    src = tmp_path / "map.tsv"
    src.write_bytes(b"genome1\tLabel One\r\ngenome2\tLabel Two\r\n")

    new_path = check_line_endings(str(src), "-m")

    assert new_path != str(src)
    assert b"\r" not in open(new_path, "rb").read()


def test_whitespace_in_an_input_file_exits_and_says_which_line(tmp_path, capsys):
    src = tmp_path / "accs.txt"
    src.write_text("GCF_000000001.1\nGCF 000000002.1\n")

    with pytest.raises(SystemExit):
        check_for_whitespace(str(src), "-a")

    out = capsys.readouterr().out
    assert "spaces" in out
    assert "GCF 000000002.1" in out


class TestSharedSingleColumnReader:
    """
    One definition of "what's in a single-column input file", used by the preflight
    checks, populate_run_data, the search subcommands, and the Pfam/KO target readers.
    They used to each do their own thing, and only agreed because preflight rewrote
    the file out from under them.
    """

    def _write(self, tmp_path, text):
        p = tmp_path / "in.txt"
        p.write_text(text)
        return str(p)

    def test_preserves_first_seen_order(self, tmp_path):
        path = self._write(tmp_path, "c\na\nb\na\nc\n")
        assert read_single_column_file(path) == ["c", "a", "b"]

    def test_strips_surrounding_whitespace(self, tmp_path):
        path = self._write(tmp_path, "  a  \n\tb\t\n")
        assert read_single_column_file(path) == ["a", "b"]

    def test_a_file_with_no_trailing_newline_keeps_its_last_entry(self, tmp_path):
        path = self._write(tmp_path, "a\nb")
        assert read_single_column_file(path) == ["a", "b"]

    def test_does_not_treat_hashes_as_comments(self, tmp_path):
        """
        Deliberate: nothing in the main pipeline ever has, and quietly starting to
        would put this reader and the target-ID shape check at odds about what an
        entry is.
        """
        path = self._write(tmp_path, "#a\nb\n")
        assert read_single_column_file(path) == ["#a", "b"]

    def test_populate_run_data_no_longer_makes_empty_genome_ids(self, tmp_path):
        """
        A blank line used to survive whenever the file had no duplicates to trigger a
        rewrite, giving GenomeData an empty id.
        """
        import argparse
        from gtotree.utils.misc.general import populate_run_data

        path = self._write(tmp_path, "GCF_000000001.1\n\nGCF_000000002.1\n")
        args = argparse.Namespace(ncbi_accessions=path, genbank_files=None,
                                  fasta_files=None, amino_acid_files=None,
                                  run_files_dir=str(tmp_path),
                                  run_files_dir_rel=str(tmp_path))

        ids = [gd.id for gd in populate_run_data(args).ncbi_accs]
        assert ids == ["GCF_000000001.1", "GCF_000000002.1"]

    def test_duplicate_accessions_are_collapsed_by_populate_run_data(self, tmp_path):
        import argparse
        from gtotree.utils.misc.general import populate_run_data

        path = self._write(tmp_path, "GCF_1\nGCF_2\nGCF_1\n")
        args = argparse.Namespace(ncbi_accessions=path, genbank_files=None,
                                  fasta_files=None, amino_acid_files=None,
                                  run_files_dir=str(tmp_path),
                                  run_files_dir_rel=str(tmp_path))

        assert [gd.id for gd in populate_run_data(args).ncbi_accs] == ["GCF_1", "GCF_2"]

    def test_check_expected_single_column_input_returns_the_path_unchanged(self,
                                                                            tmp_path):
        from gtotree.utils.misc.preflight_checks import (
            check_expected_single_column_input)

        path = self._write(tmp_path, "GCF_1\r\nGCF_2\r\nGCF_1\r\n")
        returned, count = check_expected_single_column_input(path, "-a", get_count=True)

        assert returned == path
        assert count == 2
        assert [p.name for p in tmp_path.iterdir()] == ["in.txt"]
