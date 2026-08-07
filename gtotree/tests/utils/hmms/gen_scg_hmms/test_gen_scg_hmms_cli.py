import os
import pytest # type: ignore
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as cli
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import (
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
    ["-w", "Nitrospirota"],
])
def test_each_input_source_accepted_alone(argv):
    """All five genome sources the main GToTree driver takes are valid on their own."""
    assert check_args(_parse(*argv)) is not None


def test_input_sources_can_be_combined():
    args = _parse("-a", "accs.txt", "-w", "Nitrospirota",
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


def test_resume_without_working_dir_but_no_final_output_starts_fresh(tmp_path, capsys):
    """
    Output dir exists (e.g. left behind by a run that died before writing outputs)
    but has no final table and no working dir -> genuinely nothing to resume, so
    starting fresh is right and the message should say so.
    """
    out = tmp_path / "out"
    out_dir, work_dir = setup_output_dir(_parse("-a", "x.txt", "-o", str(out)))
    import shutil
    shutil.rmtree(work_dir)  # no final table was ever written

    _, work_dir2 = setup_output_dir(_parse("-a", "x.txt", "-o", str(out), "--resume"))
    assert os.path.isdir(work_dir2)
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
    assert _parse("-w", "X", "-S", "ncbi").source == "ncbi"


def test_derep_rank_accepts_auto():
    assert _parse("-w", "X", "--derep-rank", "auto").derep_rank == "auto"


################################################################################
# reference data
################################################################################

def _resolve_args(**overrides):
    import argparse

    base = dict(target_accessions=None, wanted_ref_tax=None, source="gtdb",
                target_rank=None, derep_rank="off", genbank_files=None,
                fasta_files=None, amino_acid_files=None)
    base.update(overrides)
    return argparse.Namespace(**base)


def test_phase_one_fetches_the_ncbi_table_for_target_accessions(reference_data_fetches,
                                                                monkeypatch, tmp_path):
    """
    Regression test. This module reaches the NCBI table through `ncbi_data_table_path`,
    which only resolves a path -- nothing here ever downloaded the asset. On a machine
    that had never run the main GToTree program that's a missing-file failure, and on
    one that had, it silently used a stale copy forever.
    """
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000008865.2\n")

    monkeypatch.setattr(cli, "read_accessions_file", lambda p: ["GCF_000008865.2"])
    monkeypatch.setattr(cli, "build_local_genomes", lambda args: ([], []))

    cli.phase_resolve_genomes(_resolve_args(target_accessions=str(accs)))

    assert reference_data_fetches == ["ncbi"]


def test_phase_one_fetches_both_tables_for_a_gtdb_wanted_ref_tax(reference_data_fetches,
                                                                 monkeypatch):
    """
    A GTDB selection needs both: GTDB resolves the taxon, but the assemblies it names
    are NCBI accessions, screened against and downloaded using the NCBI summary.
    """
    monkeypatch.setattr(cli, "build_local_genomes", lambda args: ([], []))
    monkeypatch.setattr(cli, "describe_source_version", lambda source: None)

    class _Selection:
        canonical = "Nitrospirota"
        resolved_rank = "phylum"
        effective_derep_rank = None
        warnings = []

    monkeypatch.setattr(cli, "resolve_wanted_ref_tax_accessions",
                        lambda *a, **kw: (["GCF_000008865.2"], _Selection()))

    cli.phase_resolve_genomes(_resolve_args(wanted_ref_tax="Nitrospirota",
                                            source="gtdb"))

    assert sorted(reference_data_fetches) == ["gtdb", "ncbi"]


def test_phase_one_fetches_nothing_for_a_local_files_only_run(reference_data_fetches,
                                                              monkeypatch):
    """Local genome files never touch either table, so neither should be fetched."""
    from gtotree.utils.misc.general import GenomeData

    gd = GenomeData.from_path("/tmp/g1.faa", "amino-acid")
    monkeypatch.setattr(cli, "build_local_genomes", lambda args: ([gd], []))

    cli.phase_resolve_genomes(_resolve_args(amino_acid_files="aa.txt"))

    assert reference_data_fetches == []
