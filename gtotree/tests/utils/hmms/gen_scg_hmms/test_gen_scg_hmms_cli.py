import os
import pytest # type: ignore
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as cli
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import GenSCGHMMsError
from gtotree.utils.misc.general import OutputDirExistsError
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import (
    build_parser,
    check_args,
    setup_output_dir,
)
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import build_run_data


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
    assert args.ncbi_accessions == "accs.txt"
    assert args.wanted_ref_tax == ["Nitrospirota"]
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

    with pytest.raises(OutputDirExistsError, match="already exists"):
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
    assert _parse("-w", "X", "--source", "ncbi").source == "ncbi"


def test_derep_rank_accepts_auto():
    assert _parse("-w", "X", "--derep-rank", "auto").derep_rank == "auto"


################################################################################
# reference data
################################################################################

def _resolve_args(**overrides):
    """
    An args namespace for phase_resolve_genomes tests.

    Defaults are taken from the real parser so that adding a CLI flag can't silently
    break these with an AttributeError -- phase 1 reads args attributes directly, and a
    hand-maintained dict here is exactly the kind of parallel list that drifts out of
    sync.
    """
    import argparse

    base = vars(cli.build_parser().parse_args([]))
    # these tests drive phase 1 directly and set their own inputs, so start from
    # "nothing requested" rather than the parser's placeholder values
    base.update(ncbi_accessions=None, wanted_ref_tax=None, source="gtdb",
                derep_rank="off")
    base.update(overrides)
    return argparse.Namespace(**base)


def _run_phase_one(args, tmp_path):
    """Build the RunData the driver would and run phase 1 over it."""
    out_dir = str(tmp_path / "out")
    work_dir = str(tmp_path / "work")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)
    run_data = build_run_data(args, out_dir, work_dir)
    return cli.phase_resolve_genomes(args, run_data)


def _patch_wanted_ref_tax(monkeypatch, by_taxon):
    """
    Stand in for the taxonomy core.

    Patched on `wanted_ref_tax` itself rather than on this CLI module, because phase 1
    is `general.resolve_input_genomes` now and reaches the taxonomy layer through its
    own import.
    """
    from gtotree.utils.taxonomy import wanted_ref_tax as wrt

    monkeypatch.setattr(wrt, "describe_source_version", lambda source: None)
    monkeypatch.setattr(wrt, "resolve_wanted_ref_tax_accessions",
                        lambda source, taxon, **kw: by_taxon[taxon])


class _Sel:
    def __init__(self, canonical, accessions=(), resolved_rank="phylum",
                 effective_derep_rank=None):
        self.canonical = canonical
        self.accessions = list(accessions)
        self.resolved_rank = resolved_rank
        self.effective_derep_rank = effective_derep_rank
        self.warnings = []


def test_phase_one_fetches_the_ncbi_table_for_accessions(reference_data_fetches,
                                                         tmp_path):
    """
    Regression test. This module reaches the NCBI table through `ncbi_data_table_path`,
    which only resolves a path -- nothing here ever downloaded the asset. On a machine
    that had never run the main GToTree program that's a missing-file failure, and on
    one that had, it silently used a stale copy forever.
    """
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000008865.2\n")

    _run_phase_one(_resolve_args(ncbi_accessions=str(accs)), tmp_path)

    assert reference_data_fetches == ["ncbi"]


def test_phase_one_fetches_both_tables_for_a_gtdb_wanted_ref_tax(reference_data_fetches,
                                                                 monkeypatch, tmp_path):
    """
    A GTDB selection needs both: GTDB resolves the taxon, but the assemblies it names
    are NCBI accessions, screened against and downloaded using the NCBI summary.
    """
    _patch_wanted_ref_tax(monkeypatch, {
        "Nitrospirota": (["GCF_000008865.2"],
                         _Sel("Nitrospirota", ["GCF_000008865.2"]))})

    _run_phase_one(_resolve_args(wanted_ref_tax=["Nitrospirota"], source="gtdb"),
                   tmp_path)

    assert sorted(reference_data_fetches) == ["gtdb", "ncbi"]


def test_repeated_wanted_ref_tax_pools_every_taxon(reference_data_fetches,
                                                   monkeypatch, tmp_path):
    """`-w Bacteria -w Archaea` resolves each taxon and merges their genomes."""
    _patch_wanted_ref_tax(monkeypatch, {
        "Bacteria": (["GCF_000000001.1", "GCF_000000002.1"],
                     _Sel("Bacteria", ["GCF_000000001.1", "GCF_000000002.1"],
                          "domain", "family")),
        "Archaea": (["GCF_000000003.1"],
                    _Sel("Archaea", ["GCF_000000003.1"], "domain", "family")),
    })

    run_data, selections = _run_phase_one(
        _resolve_args(wanted_ref_tax=["Bacteria", "Archaea"], source="gtdb"), tmp_path)

    assert sorted(run_data.get_input_ncbi_accs()) == [
        "GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1"]
    assert [s.canonical for s in selections] == ["Bacteria", "Archaea"]

    # provenance records which taxon each genome came from
    by_id = {gd.id: gd for gd in run_data.ncbi_accs}
    assert by_id["GCF_000000001.1"].wanted_ref_tax_taxon == "Bacteria"
    assert by_id["GCF_000000003.1"].wanted_ref_tax_taxon == "Archaea"


def test_repeated_wanted_ref_tax_dedups_overlap(reference_data_fetches,
                                                monkeypatch, tmp_path):
    """A genome selected by two overlapping taxa is kept once, attributed to the
    first taxon that introduced it."""
    shared = "GCF_000000009.1"
    _patch_wanted_ref_tax(monkeypatch, {
        "TaxonA": ([shared, "GCF_000000010.1"],
                   _Sel("TaxonA", [shared, "GCF_000000010.1"])),
        "TaxonB": ([shared, "GCF_000000011.1"],
                   _Sel("TaxonB", [shared, "GCF_000000011.1"])),
    })

    run_data, _ = _run_phase_one(
        _resolve_args(wanted_ref_tax=["TaxonA", "TaxonB"], source="gtdb"), tmp_path)

    accessions = run_data.get_input_ncbi_accs()
    assert accessions.count(shared) == 1
    assert len(accessions) == 3

    by_id = {gd.id: gd for gd in run_data.ncbi_accs}
    assert by_id[shared].wanted_ref_tax_taxon == "TaxonA"

    # and the overlap is accounted for rather than silently swallowed
    assert run_data.wanted_ref_tax_selections[1]["num_selected"] == 2
    assert run_data.wanted_ref_tax_selections[1]["num_added"] == 1


def test_phase_one_fetches_nothing_for_a_local_files_only_run(reference_data_fetches,
                                                              tmp_path):
    """Local genome files never touch either table, so neither should be fetched."""
    genome = tmp_path / "g1.faa"
    genome.write_text(">p1\nMKVLAAAL\n")
    listing = tmp_path / "aa.txt"
    listing.write_text(f"{genome}\n")

    run_data, _ = _run_phase_one(_resolve_args(amino_acid_files=str(listing)),
                                 tmp_path)

    assert reference_data_fetches == []
    assert [gd.id for gd in run_data.all_input_genomes] == ["g1"]


def test_phase_one_raises_when_nothing_resolves(tmp_path):
    args = _resolve_args()
    out_dir = str(tmp_path / "out"); os.makedirs(out_dir, exist_ok=True)
    work_dir = str(tmp_path / "work"); os.makedirs(work_dir, exist_ok=True)
    run_data = build_run_data(args, out_dir, work_dir)

    with pytest.raises(GenSCGHMMsError, match="No input genomes"):
        cli.phase_resolve_genomes(args, run_data)
