import os

import pytest # type: ignore

from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms_cli import (
    DEFAULT_PERCENT_SINGLE_COPY,
    DEFAULT_THREADS,
    build_parser,
    check_args,
    setup_output_dir,
)


def _parse(*argv):
    return build_parser().parse_args(list(argv))


################################################################################
# input sources
################################################################################

def test_requires_at_least_one_input():
    with pytest.raises(GenSCGHMMsError, match="need some target genomes"):
        check_args(_parse())


@pytest.mark.parametrize("argv", [
    ["-a", "accs.txt"],
    ["-g", "genbanks.txt"],
    ["-f", "fastas.txt"],
    ["-A", "aa.txt"],
    ["-W", "Nitrospirota"],
])
def test_each_input_source_accepted_alone(argv):
    """All five genome sources the main GToTree driver takes are valid on their own."""
    assert check_args(_parse(*argv)) is not None


def test_input_sources_can_be_combined():
    args = _parse("-a", "accs.txt", "-W", "Nitrospirota",
                  "-g", "gb.txt", "-f", "fa.txt", "-A", "aa.txt")
    check_args(args)
    assert args.target_accessions == "accs.txt"
    assert args.wanted_ref_tax == "Nitrospirota"
    assert args.genbank_files == "gb.txt"
    assert args.fasta_files == "fa.txt"
    assert args.amino_acid_files == "aa.txt"


################################################################################
# numeric validation
################################################################################

@pytest.mark.parametrize("percent", [0, -1, 101, 200])
def test_percent_single_copy_bounds(percent):
    with pytest.raises(GenSCGHMMsError, match="percent-single-copy"):
        check_args(_parse("-a", "x.txt", "-p", str(percent)))


@pytest.mark.parametrize("percent", [1, 50, 90, 100])
def test_percent_single_copy_valid_values(percent):
    args = check_args(_parse("-a", "x.txt", "-p", str(percent)))
    assert args.percent_single_copy == percent


def test_threads_must_be_positive():
    with pytest.raises(GenSCGHMMsError, match="num-threads"):
        check_args(_parse("-a", "x.txt", "-t", "0"))


def test_num_jobs_must_be_positive():
    with pytest.raises(GenSCGHMMsError, match="num-jobs"):
        check_args(_parse("-a", "x.txt", "-j", "0"))


################################################################################
# output directory / resume interactions
################################################################################

def test_creates_fresh_output_dir(tmp_path):
    out = tmp_path / "out"
    args = _parse("-a", "x.txt", "-o", str(out))
    out_dir, work_dir = setup_output_dir(args)

    assert os.path.isdir(work_dir)
    assert work_dir.startswith(out_dir)


def test_refuses_existing_dir_without_flags(tmp_path):
    out = tmp_path / "out"
    args = _parse("-a", "x.txt", "-o", str(out))
    setup_output_dir(args)

    with pytest.raises(GenSCGHMMsError, match="already exists"):
        setup_output_dir(_parse("-a", "x.txt", "-o", str(out)))


def test_force_overwrite_wipes_previous_run(tmp_path):
    out = tmp_path / "out"
    _, work_dir = setup_output_dir(_parse("-a", "x.txt", "-o", str(out)))
    marker = os.path.join(work_dir, "marker.txt")
    open(marker, "w").write("x")

    setup_output_dir(_parse("-a", "x.txt", "-o", str(out), "-F"))
    assert not os.path.exists(marker)


def test_resume_preserves_working_dir(tmp_path):
    out = tmp_path / "out"
    _, work_dir = setup_output_dir(_parse("-a", "x.txt", "-o", str(out)))
    marker = os.path.join(work_dir, "marker.txt")
    open(marker, "w").write("x")

    setup_output_dir(_parse("-a", "x.txt", "-o", str(out), "--resume"))
    assert os.path.exists(marker)


def test_resume_and_force_are_mutually_exclusive(tmp_path):
    """One reuses the previous run, the other deletes it."""
    out = tmp_path / "out"
    with pytest.raises(GenSCGHMMsError, match="can't be used together"):
        setup_output_dir(_parse("-a", "x.txt", "-o", str(out), "--resume", "-F"))


def test_resume_without_previous_run_starts_fresh(tmp_path, capsys):
    out = tmp_path / "out"
    out_dir, work_dir = setup_output_dir(
        _parse("-a", "x.txt", "-o", str(out), "--resume"))

    assert os.path.isdir(work_dir)
    assert "start fresh" in capsys.readouterr().out


def test_output_dir_trailing_slash_normalized(tmp_path):
    out = str(tmp_path / "out") + "/"
    out_dir, _ = setup_output_dir(_parse("-a", "x.txt", "-o", out))
    assert not out_dir.endswith("/")


################################################################################
# parser wiring
################################################################################

def test_parser_sets_func_default():
    assert _parse("-a", "x.txt").func == "gen_scg_hmms"


def test_source_accepts_lowercase():
    assert _parse("-W", "X", "-S", "ncbi").source == "ncbi"


def test_derep_rank_accepts_auto():
    assert _parse("-W", "X", "--derep-rank", "auto").derep_rank == "auto"
