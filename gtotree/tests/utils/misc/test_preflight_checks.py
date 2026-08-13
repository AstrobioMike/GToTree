"""Unit tests for input validation in gtotree/utils/preflight_checks.py."""

import pytest  # type: ignore

from gtotree.utils.misc.preflight_checks import (check_for_duplicates,
                                            check_for_whitespace,
                                            check_mapping_file_problem_chars_and_fields)


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


def test_duplicate_entries_are_removed_in_input_order(tmp_path, monkeypatch):
    """De-duplication must be order-preserving; genome order flows into the alignment."""
    import gtotree.utils.misc.preflight_checks as pc
    monkeypatch.setattr(pc.time, "sleep", lambda *_a: None)

    src = tmp_path / "accs.txt"
    ordered = [f"GCF_{i:09d}.1" for i in range(40)]
    src.write_text("\n".join(ordered + [ordered[3], ordered[17]]) + "\n")

    new_path = check_for_duplicates(str(src), "-a")

    assert new_path != str(src)
    assert [ln for ln in open(new_path).read().splitlines() if ln.strip()] == ordered


def test_whitespace_in_an_input_file_exits_and_says_which_line(tmp_path, capsys):
    src = tmp_path / "accs.txt"
    src.write_text("GCF_000000001.1\nGCF 000000002.1\n")

    with pytest.raises(SystemExit):
        check_for_whitespace(str(src), "-a")

    out = capsys.readouterr().out
    assert "spaces" in out
    assert "GCF 000000002.1" in out
