"""
Tests for the `-p` / `-K` target-ID shape check.

Two things matter here and neither is visible from the happy path:

  * the check has to fire *early* -- its whole reason for existing is that a swapped
    targets file is otherwise indistinguishable from "none of your targets were found",
    reported after the managed dataset has been downloaded and every genome searched.
    So the wiring tests pin down that it runs before `check_for_required_dbs()` and
    before `check_output_dir()` acts on `-F`.
  * a wholesale swap and a few stray bad IDs are different mistakes and get different
    messages; a regression that collapsed them would print a 4,000-line list of a
    user's KO file back at them.
"""

import argparse

import pytest  # type: ignore

from gtotree.utils.misc.target_id_checks import (TargetIDFormatError,
                                                 PFAM_FORMAT, KO_FORMAT,
                                                 WRONG_CASE, MALFORMED,
                                                 check_target_id_file,
                                                 classify_entry,
                                                 read_target_ids,
                                                 wholesale_swap)


@pytest.fixture
def targets(tmp_path):
    """Returns a callable writing a targets file and handing back its path."""

    def _write(name, entries, newline="\n"):
        path = tmp_path / name
        path.write_bytes(newline.join(entries).encode() + newline.encode())
        return str(path)

    return _write


################################################################################
# classification
################################################################################

@pytest.mark.parametrize("entry", ["PF00789", "PF00001", "PF23000", "PF00789.19",
                                   "PF00001.27", "PF90001"])
def test_accepts_well_formed_pfam_ids(entry):
    assert classify_entry(entry, PFAM_FORMAT) is None


@pytest.mark.parametrize("entry", ["K01601", "K00927", "K23000"])
def test_accepts_well_formed_ko_ids(entry):
    assert classify_entry(entry, KO_FORMAT) is None


def test_a_ko_id_in_a_pfam_file_is_named_as_a_ko():
    assert classify_entry("K01601", PFAM_FORMAT) == "ko"


def test_a_pfam_id_in_a_ko_file_is_named_as_a_pfam():
    assert classify_entry("PF00789", KO_FORMAT) == "pfam"


def test_a_versioned_pfam_id_in_a_ko_file_is_still_named_as_a_pfam():
    assert classify_entry("PF00789.19", KO_FORMAT) == "pfam"


@pytest.mark.parametrize("entry", ["pf00789", "Pf00789", "k01601"])
def test_lowercase_ids_are_flagged_as_a_case_problem(entry):
    fmt = PFAM_FORMAT if entry.upper().startswith("PF") else KO_FORMAT
    assert classify_entry(entry, fmt) == WRONG_CASE


@pytest.mark.parametrize("entry", ["7tm_1", "GCF_000008865.2", "PF", "PFxxxxx",
                                   "PF00789.", "ribosomal protein S4"])
def test_junk_is_flagged_as_malformed(entry):
    assert classify_entry(entry, PFAM_FORMAT) == MALFORMED


def test_a_pfam_accession_is_not_mistaken_for_a_ko_by_its_leading_letters():
    """'PF...' must never classify as a KO, or the swap advice would point backwards."""
    assert classify_entry("PF00789", PFAM_FORMAT) is None
    assert classify_entry("PF00789", KO_FORMAT) == "pfam"


################################################################################
# reading
################################################################################

def test_blank_lines_are_skipped_and_line_numbers_still_refer_to_the_real_file(targets):
    path = targets("t.txt", ["PF00001", "", "   ", "K01601"])
    assert read_target_ids(path) == [(1, "PF00001"), (4, "K01601")]


def test_crlf_files_are_tolerated(targets):
    """
    The check runs before line-ending normalization in the main run, so a stray '\\r'
    must not turn every ID into a malformed one.
    """
    path = targets("t.txt", ["PF00001", "PF00002"], newline="\r\n")
    assert check_target_id_file(path, "pfam") == [(1, "PF00001"), (2, "PF00002")]


################################################################################
# check_target_id_file
################################################################################

def test_a_clean_pfam_file_passes(targets):
    path = targets("pfams.txt", ["PF00789", "PF00001.27"])
    assert check_target_id_file(path, "pfam") == [(1, "PF00789"), (2, "PF00001.27")]


def test_a_clean_ko_file_passes(targets):
    path = targets("kos.txt", ["K01601", "K00927"])
    assert check_target_id_file(path, "ko") == [(1, "K01601"), (2, "K00927")]


def test_an_empty_file_is_rejected(targets):
    path = targets("pfams.txt", ["", "  "])
    with pytest.raises(TargetIDFormatError, match="is empty"):
        check_target_id_file(path, "pfam")


def test_a_ko_file_passed_to_p_points_at_the_K_flag(targets):
    path = targets("kos.txt", ["K01601", "K00927", "K02358"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    message = str(excinfo.value)
    assert "look like KO IDs" in message
    assert "`-K`/`--target-kos`" in message


def test_a_pfam_file_passed_to_K_points_at_the_p_flag(targets):
    path = targets("pfams.txt", ["PF00789", "PF00001"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "ko")

    message = str(excinfo.value)
    assert "look like Pfam IDs" in message
    assert "`-p`/`--target-pfams`" in message


def test_a_wholesale_swap_does_not_print_line_numbers(targets):
    """
    The swap message names the mistake instead of itemizing it; line-by-line output
    here would just be the user's whole file read back to them.
    """
    path = targets("kos.txt", ["K01601", "K00927"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    assert "line 1:" not in str(excinfo.value)


def test_a_long_swapped_file_is_truncated(targets):
    path = targets("kos.txt", [f"K{i:05d}" for i in range(1, 51)])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    message = str(excinfo.value)
    assert "All 50 of its entries" in message
    assert "...and 40 more." in message
    assert "K00050" not in message


def test_a_few_stray_bad_ids_are_itemized_with_line_numbers(targets):
    path = targets("pfams.txt", ["PF00789", "7tm_1", "PF00001", "K01601"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    message = str(excinfo.value)
    assert "2 entries" in message
    assert "line 2: 7tm_1" in message
    assert "line 4: K01601" in message
    # the good ones aren't listed (checked by line number: the IDs themselves also
    # turn up in the message's "e.g. 'PF00789'" prose)
    assert "line 1:" not in message
    assert "line 3:" not in message


def test_a_mixed_file_still_names_the_other_flag_for_the_stray_ids(targets):
    path = targets("pfams.txt", ["PF00789", "K01601"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    assert "Those go to `-K`/`--target-kos`" in str(excinfo.value)


def test_lowercase_ids_get_the_case_hint(targets):
    path = targets("pfams.txt", ["PF00789", "pf00001"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    assert "uppercased" in str(excinfo.value)


def test_a_single_bad_entry_reads_in_the_singular(targets):
    path = targets("pfams.txt", ["PF00789", "7tm_1"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    assert "1 entry that doesn't look like" in str(excinfo.value)


def test_an_all_lowercase_file_is_not_reported_as_a_swap(targets):
    """
    Every entry is bad, but they're all *this* type in the wrong case -- telling the
    user to pass the file to the other flag would be actively wrong.
    """
    path = targets("pfams.txt", ["pf00789", "pf00001"])
    with pytest.raises(TargetIDFormatError) as excinfo:
        check_target_id_file(path, "pfam")

    message = str(excinfo.value)
    assert "--target-kos" not in message
    assert "uppercased" in message


def test_wholesale_swap_needs_every_entry_to_agree():
    bad = [(1, "K01601", "ko"), (2, "7tm_1", MALFORMED)]
    assert wholesale_swap(bad, total=2) is None
    assert wholesale_swap([(1, "K01601", "ko")], total=1) is KO_FORMAT
    # a bad entry among good ones is never a wholesale swap
    assert wholesale_swap([(1, "K01601", "ko")], total=5) is None


################################################################################
# wiring: the main GToTree run
################################################################################

def _preflight_args(**overrides):
    base = {
        "target_pfams_file": None,
        "target_kos_file": None,
        "output_dir": "should-not-be-touched",
        "force_overwrite": True,
        "resume": False,
    }
    base.update(overrides)
    return argparse.Namespace(**base)


def test_preflight_exits_on_a_swapped_pfam_file(targets, capsys):
    from gtotree.utils.misc.preflight_checks import check_target_id_formats

    path = targets("kos.txt", ["K01601", "K00927"])
    with pytest.raises(SystemExit):
        check_target_id_formats(_preflight_args(target_pfams_file=path))

    assert "--target-kos" in capsys.readouterr().out


def test_preflight_exits_on_a_swapped_ko_file(targets, capsys):
    from gtotree.utils.misc.preflight_checks import check_target_id_formats

    path = targets("pfams.txt", ["PF00789", "PF00001"])
    with pytest.raises(SystemExit):
        check_target_id_formats(_preflight_args(target_kos_file=path))

    assert "--target-pfams" in capsys.readouterr().out


def test_preflight_is_a_no_op_without_target_files():
    from gtotree.utils.misc.preflight_checks import check_target_id_formats

    assert check_target_id_formats(_preflight_args()) is None


def test_preflight_reports_a_missing_targets_file_before_anything_else(targets):
    """
    A `-p` path that doesn't exist currently isn't noticed until after the Pfam
    dataset has been fetched; the early check catches it too.
    """
    from gtotree.utils.misc.preflight_checks import check_target_id_formats

    with pytest.raises(SystemExit):
        check_target_id_formats(_preflight_args(target_pfams_file="nope.txt"))


def test_the_check_runs_before_dbs_are_fetched_and_before_F_removes_the_output(targets,
                                                                               monkeypatch):
    """
    Placement is the whole point, so it's asserted rather than assumed: nothing
    expensive or destructive may run first.
    """
    import gtotree.utils.misc.preflight_checks as pc

    called = []
    monkeypatch.setattr(pc, "check_output_dir",
                        lambda args: called.append("output_dir") or args)
    monkeypatch.setattr(pc, "check_input_genome_files",
                        lambda args: called.append("input_genomes") or args)
    monkeypatch.setattr(pc, "check_for_required_dbs",
                        lambda *a, **k: called.append("dbs"))
    for name in ("check_for_minimum_args", "check_optional_deps", "check_set_values",
                 "check_lineage", "check_tree_program", "checks_for_nucleotide_mode",
                 "check_wanted_ref_tax_args"):
        monkeypatch.setattr(pc, name, lambda *a, **k: None)

    path = targets("kos.txt", ["K01601", "K00927"])
    args = _preflight_args(target_pfams_file=path)

    with pytest.raises(SystemExit):
        pc.primary_args_validation(args)

    assert called == []


def test_input_genome_files_are_validated_before_F_removes_the_output(monkeypatch):
    """
    `check_output_dir()` is what acts on `-F` and rmtree's the previous run, so a
    malformed `-g`/`-f`/`-A` listing has to be caught ahead of it -- otherwise a typo
    in an input list costs you the previous run's output before anything reads it.
    """
    import gtotree.utils.misc.preflight_checks as pc

    called = []
    monkeypatch.setattr(pc, "check_output_dir",
                        lambda args: called.append("output_dir") or args)
    monkeypatch.setattr(pc, "check_input_genome_files",
                        lambda args: called.append("input_genomes") or args)
    monkeypatch.setattr(pc, "check_target_id_formats",
                        lambda args: called.append("target_ids"))
    for name in ("check_for_minimum_args", "check_optional_deps", "check_set_values",
                 "check_lineage", "check_tree_program", "checks_for_nucleotide_mode",
                 "check_wanted_ref_tax_args"):
        monkeypatch.setattr(pc, name, lambda *a, **k: None)

    pc.primary_args_validation(_preflight_args())

    assert called == ["target_ids", "input_genomes", "output_dir"]


################################################################################
# wiring: target-type validation (gtt search-annotations -p / -K)
################################################################################

def test_pfam_targets_flag_rejects_a_ko_targets_file(targets):
    from gtotree.utils.target_search.target_search_setup import (TargetSearchError,
                                                                 validate_input_files,
                                                                 make_args)
    from gtotree.utils.target_search.target_search_spec import get_spec

    path = targets("kos.txt", ["K01601", "K00927"])
    args = make_args(target_pfams_file=path, amino_acid_files=None)

    with pytest.raises(TargetSearchError, match="--target-kos"):
        validate_input_files(args, get_spec("pfam"))


def test_ko_targets_flag_rejects_a_pfam_targets_file(targets):
    from gtotree.utils.target_search.target_search_setup import (TargetSearchError,
                                                                 validate_input_files,
                                                                 make_args)
    from gtotree.utils.target_search.target_search_spec import get_spec

    path = targets("pfams.txt", ["PF00789", "PF00001"])
    args = make_args(target_kos_file=path)

    with pytest.raises(TargetSearchError, match="--target-pfams"):
        validate_input_files(args, get_spec("ko"))


def test_pfam_targets_flag_accepts_a_real_pfam_targets_file(targets):
    from gtotree.utils.target_search.target_search_setup import (validate_input_files,
                                                                 make_args)
    from gtotree.utils.target_search.target_search_spec import get_spec

    path = targets("pfams.txt", ["PF00789", "PF00001.27"])
    args = validate_input_files(make_args(target_pfams_file=path), get_spec("pfam"))

    assert args.total_targets == 2


def test_both_specs_name_a_real_target_id_format():
    from gtotree.utils.misc.target_id_checks import FORMATS
    from gtotree.utils.target_search.target_search_spec import get_spec

    for name in ("pfam", "ko"):
        assert get_spec(name).target_id_format in FORMATS
