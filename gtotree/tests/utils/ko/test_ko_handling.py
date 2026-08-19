"""
Tests for reading the `-K` targets file.

`get_wanted_KOs` had no coverage, and the gap it left was an invisible one: a blank
line in the targets file became an empty-string "wanted" KO, which is never in
`ko_list`, so it ended up in `failed_ko_targets`. The found count still matched
`total_ko_targets` (which skips blanks), so the run reported every target found in
green *and* wrote a `failed-ko-targets.txt` holding a single empty line. Nothing about
that is visible from the happy path, hence the regression tests.
"""

from types import SimpleNamespace

import pytest  # type: ignore

from gtotree.utils.ko.ko_handling import get_wanted_KOs, get_target_KOs_tab


@pytest.fixture
def ko_targets_file(tmp_path):
    """Returns a callable writing a `-K` file and handing back a run_data stub."""

    def _write(entries, newline="\n"):
        path = tmp_path / "target-kos.txt"
        path.write_bytes(newline.join(entries).encode() + newline.encode())
        return SimpleNamespace(target_kos_file=str(path))

    return _write


def test_blank_lines_are_not_read_as_wanted_kos(ko_targets_file):
    run_data = ko_targets_file(["K01601", "", "K00927", "   "])
    assert get_wanted_KOs(run_data) == ["K01601", "K00927"]


def test_a_trailing_blank_line_does_not_manufacture_a_failed_target(ko_targets_file,
                                                                    tmp_path):
    """
    The whole point: an otherwise-perfect targets file that happens to end with a blank
    line must not report a failure.
    """
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("knum\tdefinition\nK01601\tRuBisCO\nK00927\tPGK\n")

    run_data = ko_targets_file(["K01601", "K00927", ""])
    wanted = get_wanted_KOs(run_data)

    found, missing = get_target_KOs_tab(str(ko_list), wanted,
                                        str(tmp_path / "target-kos.tsv"))

    assert found == ["K01601", "K00927"]
    assert missing == []


def test_genuinely_missing_kos_are_still_reported(ko_targets_file, tmp_path):
    """The blank-line fix must not swallow real misses."""
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("knum\tdefinition\nK01601\tRuBisCO\n")

    run_data = ko_targets_file(["K01601", "K99999", ""])
    found, missing = get_target_KOs_tab(str(ko_list), get_wanted_KOs(run_data),
                                        str(tmp_path / "target-kos.tsv"))

    assert found == ["K01601"]
    assert missing == ["K99999"]


def test_crlf_targets_files_are_read_cleanly(ko_targets_file):
    run_data = ko_targets_file(["K01601", "K00927"], newline="\r\n")
    assert get_wanted_KOs(run_data) == ["K01601", "K00927"]


def test_the_wanted_list_matches_what_the_shared_line_counter_counts(ko_targets_file):
    """
    `total_ko_targets` comes from check_expected_single_column_input and drives the
    "N of M found" reporting, so the two readers have to agree on what a target is.
    """
    from gtotree.utils.misc.preflight_checks import check_expected_single_column_input

    run_data = ko_targets_file(["K01601", "", "K00927"])
    _path, total = check_expected_single_column_input(run_data.target_kos_file, "-K",
                                                      get_count=True)

    assert len(get_wanted_KOs(run_data)) == total


def test_failed_targets_come_out_in_input_order(ko_targets_file, tmp_path):
    """
    `missing_KOs` was a set difference, so failed-ko-targets.txt came out in a
    different (arbitrary) row order on every run. Ordering is asserted against a
    ko_list whose own order disagrees with the input's, so a passing result can't be
    an accident of the two happening to match.
    """
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("knum\tdefinition\nK00003\tc\nK00001\ta\n")

    wanted = ["K09999", "K00001", "K07777", "K00003", "K01111"]
    run_data = ko_targets_file(wanted)

    _found, missing = get_target_KOs_tab(str(ko_list), get_wanted_KOs(run_data),
                                         str(tmp_path / "target-kos.tsv"))

    assert missing == ["K09999", "K07777", "K01111"]


def test_duplicate_kos_are_collapsed(ko_targets_file):
    """
    De-duplication used to happen by rewriting the file during preflight; it now
    happens in the reader, so this has to hold without any preflight involvement.
    """
    run_data = ko_targets_file(["K01601", "K00927", "K01601"])
    assert get_wanted_KOs(run_data) == ["K01601", "K00927"]
