"""
The updating-headers stage, and what it does when it can't do its job.

Swapping in fuller labels sits on top of an alignment that's already finished, so a
failure here costs the nicer names and nothing else. The stage falls back to the
original genome IDs and lets the run go on to build its tree -- except when the
alignment itself has gone missing, which nothing downstream can survive.
"""

import argparse
import os
from datetime import datetime
import pytest  # type: ignore

from gtotree.main_stages import updating_headers as U
from gtotree.main_stages.updating_headers import update_headers
from gtotree.utils.misc.general import RunData
from gtotree.utils.misc.seqs import swap_labels_in_alignment as real_swap
from gtotree.utils.misc.stages import PipelineStage


ALIGNMENT = ">G0\nMMMMM\n>G1\nKKKKK\n"


@pytest.fixture(autouse=True)
def quiet_side_effects(monkeypatch):
    """Keep the stage banner and run-data writing out of these."""
    monkeypatch.setattr(U, "report_processing_stage", lambda *a, **kw: None)
    monkeypatch.setattr(U, "write_run_data", lambda *a, **kw: None)


def make_run_data(tmp_path, with_alignment="output"):
    out_dir = tmp_path / "output"
    run_files_dir = out_dir / "run-files"
    run_files_dir.mkdir(parents=True)

    # enough of a real run's furniture for the early-exit path to archive its log
    logs_dir = out_dir / "logs"
    (logs_dir / "gtotree-logs").mkdir(parents=True)
    (out_dir / "gtotree-runlog.txt").write_text("")

    rd = RunData()
    rd.output_dir = str(out_dir)
    rd.run_files_dir = str(run_files_dir)
    rd.logs_dir = str(logs_dir)
    rd.start_time = datetime.now()
    rd.general_ext = ".faa"
    rd.updating_headers = True
    rd.mapping_dict = {"G0": "Genus species G0", "G1": "Genus species G1"}
    rd.concatenated_alignment_path = str(out_dir / "aligned-SCGs.faa")

    if with_alignment == "output":
        (out_dir / "aligned-SCGs.faa").write_text(ALIGNMENT)
    elif with_alignment == "run-files":
        (run_files_dir / "aligned-SCGs.faa").write_text(ALIGNMENT)

    return rd


def make_args(**overrides):
    defaults = dict(add_gtdb_tax=False, add_ncbi_tax=False)
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def explode(message):
    def boom(*a, **kw):
        raise RuntimeError(message)
    return boom


################################################################################
# the happy path, and opting out
################################################################################

def test_headers_get_swapped_in(tmp_path):
    rd = make_run_data(tmp_path)

    rd = update_headers(make_args(), rd)

    assert rd.final_alignment_path.endswith("aligned-SCGs-mod-names.faa")
    assert ">Genus species G0" in open(rd.final_alignment_path).read()
    assert rd.stage_is_complete(PipelineStage.UPDATE_HEADERS)
    assert not rd.header_update_error


def test_a_run_not_updating_headers_just_takes_the_concatenated_alignment(tmp_path):
    rd = make_run_data(tmp_path)
    rd.updating_headers = False

    rd = update_headers(make_args(), rd)

    assert rd.final_alignment_path == rd.concatenated_alignment_path
    assert not rd.header_update_error


################################################################################
# failing without losing the run
################################################################################

def test_a_failed_taxonomy_lookup_carries_on_with_the_original_ids(tmp_path,
                                                                   monkeypatch):
    monkeypatch.setattr(U, "update_mapping_dict_with_gtdb_tax_info",
                        explode("GTDB is unreachable"))
    rd = make_run_data(tmp_path)

    rd = update_headers(make_args(add_gtdb_tax=True), rd)

    assert rd.final_alignment_path == os.path.join(rd.output_dir, "aligned-SCGs.faa")
    assert "GTDB is unreachable" in rd.header_update_error


def test_a_failed_swap_carries_on_with_the_original_ids(tmp_path, monkeypatch):
    monkeypatch.setattr(U, "swap_labels_in_alignment", explode("the swap blew up"))
    rd = make_run_data(tmp_path)

    rd = update_headers(make_args(), rd)

    assert rd.final_alignment_path == os.path.join(rd.output_dir, "aligned-SCGs.faa")
    assert ">G0" in open(rd.final_alignment_path).read()


def test_the_fallback_finds_the_alignment_wherever_it_sits(tmp_path, monkeypatch):
    """
    A run interrupted after a previous swap moved the original aside re-enters with
    the alignment in run-files/ rather than the output dir.
    """
    monkeypatch.setattr(U, "swap_labels_in_alignment", explode("nope"))
    rd = make_run_data(tmp_path, with_alignment="run-files")

    rd = update_headers(make_args(), rd)

    assert rd.final_alignment_path == os.path.join(rd.run_files_dir, "aligned-SCGs.faa")


def test_the_failure_is_reported_to_the_user(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(U, "swap_labels_in_alignment", explode("the swap blew up"))
    rd = make_run_data(tmp_path)

    update_headers(make_args(), rd)

    # the message is wrapped to the terminal width, so compare on collapsed whitespace
    out = " ".join(capsys.readouterr().out.split())
    assert "the swap blew up" in out
    assert "carrying on with the original genome IDs" in out


def test_a_failed_swap_leaves_the_stage_incomplete_so_a_resume_retries_it(tmp_path,
                                                                          monkeypatch):
    monkeypatch.setattr(U, "swap_labels_in_alignment", explode("transient"))
    rd = make_run_data(tmp_path)

    rd = update_headers(make_args(), rd)

    assert not rd.stage_is_complete(PipelineStage.UPDATE_HEADERS)

    # and a resume, with whatever went wrong now fixed, gets the fuller labels
    monkeypatch.setattr(U, "swap_labels_in_alignment", real_swap)
    rd = update_headers(make_args(), rd)

    assert rd.final_alignment_path.endswith("aligned-SCGs-mod-names.faa")
    assert rd.stage_is_complete(PipelineStage.UPDATE_HEADERS)
    assert not rd.header_update_error


################################################################################
# the one failure that can't be salvaged
################################################################################

def test_a_missing_alignment_still_exits(tmp_path):
    """
    With no alignment in either place there's nothing to tree from under any set of
    labels, so continuing would only defer the failure to the tree builder.
    """
    rd = make_run_data(tmp_path, with_alignment=None)

    with pytest.raises(SystemExit):
        update_headers(make_args(), rd)
