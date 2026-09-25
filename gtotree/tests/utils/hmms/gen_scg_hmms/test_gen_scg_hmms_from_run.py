"""
`--from-run` builds a new SCG set from a finished run's output tables with different
selection parameters, without redoing the genome retrieval or the search.

The end-to-end tests here drive the real `gen_scg_hmms` entry point, with the managed
Pfam data pointed at the mock profiles, so a full scan and the reselections from it
both run exactly as they would from the command line.
"""

import os
import shutil

import pytest  # type: ignore

import gtotree.utils.pfam.get_pfam_data as pfam_data_mod
from gtotree.tests.paths import DATA_DIR
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as cli
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_outputs as outputs
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import (
    GenSCGHMMsError,
    PfamProfileInfo,
    cap_targets,
    read_hmm_accessions,
    select_scg_targets,
    shared_pair_key,
)
from gtotree.utils.misc.messaging import REMOVED_GENOMES_FILENAME


MOTIFS = {
    "A": "MKVLAAAL",    # PF90001.3, coverage 75.5
    "B": "MARTKQTA",    # PF90002.7, coverage 60.2
    "C": "MSDKIIHL",    # PF90003.1, coverage 50.0
    "D": "MAHHWWGS",    # PF90004.2, coverage 12.4
}
A, B, C, D = "PF90001.3", "PF90002.7", "PF90003.1", "PF90004.2"

MOCK_PFAM_VERSION = "38.2-mock"


def _parse(*argv):
    return cli.build_parser().parse_args(list(argv))


def _run(*argv):
    cli.gen_scg_hmms(_parse(*argv))


def _targets(run_dir):
    with open(os.path.join(run_dir, outputs.HMM_INFO_FILENAME)) as f:
        return [line.split("\t")[0] for line in f.read().splitlines()[1:]]


def _params(run_dir):
    with open(os.path.join(run_dir, outputs.SELECTION_PARAMS_FILENAME)) as f:
        return dict(line.split("\t") for line in f.read().splitlines()[1:])


def _snapshot(run_dir):
    return {name: open(os.path.join(run_dir, name), "rb").read()
            for name in sorted(os.listdir(run_dir))}


@pytest.fixture
def mock_pfam(tmp_path, monkeypatch):
    """Point the managed Pfam data at the mock profiles."""
    pfam_dir = tmp_path / "pfam"
    pfam_dir.mkdir()
    shutil.copy(DATA_DIR / "mock-pfams.hmm", pfam_dir / pfam_data_mod.HMM_FILENAME)
    shutil.copy(DATA_DIR / "mock-pfamA.txt", pfam_dir / pfam_data_mod.INFO_FILENAME)

    state = {"version": MOCK_PFAM_VERSION}
    monkeypatch.setattr(pfam_data_mod, "get_pfam_data", lambda *a, **kw: str(pfam_dir))
    monkeypatch.setattr(pfam_data_mod, "get_stored_pfam_version",
                        lambda location: state["version"])
    return state


@pytest.fixture
def scan_dir(tmp_path, mock_pfam, monkeypatch):
    """
    A finished full run over six genomes:

      g0-g4:  A, B, and a fused C+D protein
      g5:     two copies of A, B, C, and no D

    so A and D are single-copy in 5/6 (83%) and B and C in 6/6, and C and D share a
    protein in 5/6. At the default -p 90, only B and C make it, and since D isn't
    retained, the C/D pair never becomes a conflict in this run.
    """
    monkeypatch.chdir(tmp_path)
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()

    genomes = {f"g{i}": [["A"], ["B"], ["C", "D"]] for i in range(5)}
    genomes["g5"] = [["A"], ["A"], ["B"], ["C"]]

    paths = []
    for genome_id, proteins in genomes.items():
        path = genome_dir / f"{genome_id}.faa"
        path.write_text("".join(f">p{i}\n{''.join(MOTIFS[x] for x in prot)}\n"
                                for i, prot in enumerate(proteins, 1)))
        paths.append(str(path))
    (tmp_path / "aa.txt").write_text("\n".join(paths) + "\n")

    # the mock C sits at exactly 50% coverage, which the default filter drops
    _run("-A", "aa.txt", "-o", "scan", "-t", "1", "-j", "1",
         "--min-pfam-coverage", "10")
    return "scan"


################################################################################
# end to end
################################################################################

def test_full_run_writes_the_scan_record(scan_dir):
    assert _targets(scan_dir) == [B, C]
    for name in outputs.SCAN_RECORD_FILENAMES:
        assert os.path.isfile(os.path.join(scan_dir, name))
    assert _params(scan_dir) == {"percent_single_copy": "90", "max_hmms": "250",
                                 "max_shared_protein_percent": "10", "from_run": "NA"}


def test_from_run_selects_anew_without_touching_the_original(scan_dir):
    before = _snapshot(scan_dir)

    _run("--from-run", scan_dir, "-o", "reselected", "-p", "80")

    assert _snapshot(scan_dir) == before

    # A and D now pass, but D shares a protein with C in 5/6 genomes. That conflict
    # never came up in the original run, so this only works because the full pair
    # record was carried, not just that run's exclusions.
    assert _targets("reselected") == [A, B, C]
    assert sorted(read_hmm_accessions("reselected/reselected.hmm")) == [A, B, C]
    with open("reselected/" + outputs.SHARED_PROTEIN_EXCLUSIONS_FILENAME) as f:
        rows = [line.split("\t") for line in f.read().splitlines()[1:]]
    assert [(r[0], r[2]) for r in rows] == [(D, C)]

    params = _params("reselected")
    assert params["percent_single_copy"] == "80"
    assert params["from_run"] == os.path.abspath(scan_dir)

    # the scan record came along unchanged, so the new dir is a complete run
    for name in outputs.SCAN_RECORD_FILENAMES:
        assert before[name] == open(os.path.join("reselected", name), "rb").read()


def test_from_run_matches_what_a_full_run_would_select(scan_dir, tmp_path):
    _run("--from-run", scan_dir, "-o", "reselected", "-p", "80", "--max-hmms", "2")
    _run("-A", "aa.txt", "-o", "rescanned", "-t", "1", "-j", "1",
         "--min-pfam-coverage", "10", "-p", "80", "--max-hmms", "2")

    assert _targets("reselected") == _targets("rescanned") == [B, C]
    for name in (outputs.HMM_INFO_FILENAME, outputs.SHARED_PROTEIN_EXCLUSIONS_FILENAME):
        assert (open(os.path.join("reselected", name)).read()
                == open(os.path.join("rescanned", name)).read())


def test_from_run_of_a_from_run_output(scan_dir):
    _run("--from-run", scan_dir, "-o", "first", "-p", "80")
    _run("--from-run", "first", "-o", "second", "-p", "90")
    assert _targets("second") == _targets(scan_dir)


def test_cap_ties_go_to_higher_coverage_end_to_end(scan_dir):
    """
    With the shared-protein check out of the way, all four pass at -p 80; a cap of 3
    keeps B and C (6/6), and of A and D (both 5/6) the higher-coverage A.
    """
    _run("--from-run", scan_dir, "-o", "capped", "-p", "80", "--max-hmms", "3",
         "--max-shared-protein-percent", "100")
    assert _targets("capped") == [A, B, C]


def test_from_run_carries_the_removed_genomes_report(scan_dir):
    report = os.path.join(scan_dir, REMOVED_GENOMES_FILENAME)
    with open(report, "w") as f:
        f.write("genome_id\treason\ngX\tsomething\n")
    _run("--from-run", scan_dir, "-o", "reselected", "-p", "80")
    assert open(os.path.join("reselected", REMOVED_GENOMES_FILENAME)).read() \
        == open(report).read()


def test_from_run_refuses_a_different_pfam_version(scan_dir, mock_pfam):
    mock_pfam["version"] = "39.0"
    with pytest.raises(GenSCGHMMsError, match="built with Pfam 38.2-mock"):
        _run("--from-run", scan_dir, "-o", "reselected")
    assert not os.path.exists("reselected")


def test_from_run_refuses_an_existing_output_dir_without_force(scan_dir):
    os.makedirs("taken")
    with pytest.raises(GenSCGHMMsError, match="already exists"):
        _run("--from-run", scan_dir, "-o", "taken")

    _run("--from-run", scan_dir, "-o", "taken", "-F")
    assert _targets("taken") == [B, C]


def test_from_run_reports_nothing_passing(scan_dir):
    """-p 100 still keeps B and C; with every Pfam multi-copy somewhere, nothing would."""
    with open(os.path.join(scan_dir, outputs.HIT_COUNTS_FILENAME)) as f:
        lines = f.read().splitlines()
    lines[1] = lines[1].split("\t")[0] + "\t2\t2\t2\t2"
    with open(os.path.join(scan_dir, outputs.HIT_COUNTS_FILENAME), "w") as f:
        f.write("\n".join(lines) + "\n")

    with pytest.raises(GenSCGHMMsError, match="No Pfams were found"):
        _run("--from-run", scan_dir, "-o", "reselected", "-p", "100")


################################################################################
# argument checks
################################################################################

def test_from_run_alone_is_enough_input(tmp_path):
    (tmp_path / "prev").mkdir()
    args = cli.check_args(_parse("--from-run", str(tmp_path / "prev"),
                                 "-o", str(tmp_path / "new")))
    assert args.from_run == str(tmp_path / "prev")


@pytest.mark.parametrize("extra,flag", [
    (["-w", "Nitrospirota"], "-w"),
    (["-a", "accs.txt"], "-a"),
    (["-A", "aa.txt"], "-A"),
    (["--source", "ncbi"], "--source"),
    (["--derep-rank", "genus"], "--derep-rank"),
    (["--min-pfam-coverage", "60"], "--min-pfam-coverage"),
    (["-R"], "-R/--resume"),
])
def test_from_run_refuses_scan_parameters(tmp_path, extra, flag):
    with pytest.raises(GenSCGHMMsError, match="can't be combined with") as err:
        cli.check_args(_parse("--from-run", str(tmp_path), "-o", str(tmp_path / "new"),
                              *extra))
    assert flag in str(err.value)


def test_from_run_names_every_refused_flag_at_once(tmp_path):
    with pytest.raises(GenSCGHMMsError) as err:
        cli.check_args(_parse("--from-run", str(tmp_path), "-o", str(tmp_path / "new"),
                              "-w", "Nitrospirota", "--derep-rank", "genus"))
    assert "-w" in str(err.value) and "--derep-rank" in str(err.value)


def test_from_run_ignores_execution_only_parameters(tmp_path):
    cli.check_args(_parse("--from-run", str(tmp_path), "-o", str(tmp_path / "new"),
                          "-t", "4", "-j", "2", "--keep-working-dir",
                          "-p", "80", "--max-hmms", "100"))


def test_from_run_refuses_writing_over_itself(tmp_path):
    with pytest.raises(GenSCGHMMsError, match="different directory"):
        cli.check_args(_parse("--from-run", str(tmp_path), "-o", str(tmp_path) + "/"))


def test_max_hmms_default_and_bounds():
    assert _parse("-w", "X").max_hmms == 250
    cli.check_args(_parse("-w", "X", "--max-hmms", "0"))
    with pytest.raises(GenSCGHMMsError, match="--max-hmms"):
        cli.check_args(_parse("-w", "X", "--max-hmms", "-1"))


################################################################################
# reading a previous run
################################################################################

def test_check_previous_run_needs_the_directory(tmp_path):
    with pytest.raises(GenSCGHMMsError, match="can't be found"):
        outputs.check_previous_run(str(tmp_path / "nope"))


def test_check_previous_run_needs_a_finished_run(tmp_path):
    with pytest.raises(GenSCGHMMsError, match="doesn't hold a finished"):
        outputs.check_previous_run(str(tmp_path))


def test_check_previous_run_names_missing_tables(tmp_path):
    (tmp_path / outputs.HMM_INFO_FILENAME).write_text("pfam_id\n")
    for name in outputs.SCAN_RECORD_FILENAMES:
        if name != outputs.SHARED_PROTEIN_PAIRS_FILENAME:
            (tmp_path / name).write_text("x\n")
    with pytest.raises(GenSCGHMMsError, match=outputs.SHARED_PROTEIN_PAIRS_FILENAME):
        outputs.check_previous_run(str(tmp_path))


def test_hit_counts_round_trip(tmp_path):
    counts = {"g1": {A: 1, B: 2}, "g2": {B: 1}, "g3": {}}
    outputs.write_hit_counts(str(tmp_path), ["g1", "g2", "g3"], [A, B], counts)

    genome_ids, accs, hits = outputs.read_hit_counts(
        str(tmp_path / outputs.HIT_COUNTS_FILENAME))
    assert genome_ids == ["g1", "g2", "g3"]
    assert accs == [A, B]
    assert hits == {"g1": {A: 1, B: 2}, "g2": {B: 1}, "g3": {}}


@pytest.mark.parametrize("content,match", [
    ("not\ta\theader\n", "doesn't look like"),
    ("genome\tPF1.1\tPF2.1\ng1\t1\n", "columns"),
    ("genome\tPF1.1\ng1\tx\n", "non-integer"),
    ("genome\tPF1.1\n", "no genomes"),
])
def test_hit_counts_malformed(tmp_path, content, match):
    path = tmp_path / outputs.HIT_COUNTS_FILENAME
    path.write_text(content)
    with pytest.raises(GenSCGHMMsError, match=match):
        outputs.read_hit_counts(str(path))


def test_shared_protein_pairs_round_trip(tmp_path):
    sharing = {shared_pair_key(C, D): 5, shared_pair_key(A, B): 2,
               shared_pair_key(A, C): 0}
    info = {A: PfamProfileInfo(A, "MockA", "d", 75.5)}
    path = outputs.write_shared_protein_pairs(str(tmp_path), sharing, 6, info)

    lines = open(path).read().splitlines()
    assert lines[1].split("\t") == [C, "NA", D, "NA", "5", "83.33"]
    assert lines[2].split("\t") == [A, "MockA", B, "NA", "2", "33.33"]
    assert len(lines) == 3      # the zero-genome pair isn't written

    assert outputs.read_shared_protein_pairs(path) == {
        shared_pair_key(C, D): 5, shared_pair_key(A, B): 2}


def test_shared_protein_pairs_written_even_when_empty(tmp_path):
    path = outputs.write_shared_protein_pairs(str(tmp_path), {}, 3, {})
    assert outputs.read_shared_protein_pairs(path) == {}


################################################################################
# the cap
################################################################################

def _info(**coverages):
    return {acc: PfamProfileInfo(acc, acc, "d", cov) for acc, cov in coverages.items()}


def test_cap_keeps_the_most_widely_single_copy():
    kept, out = cap_targets(["p1", "p2", "p3"], {"p1": 7, "p2": 9, "p3": 8},
                            _info(p1=99, p2=10, p3=10), 2)
    assert kept == ["p2", "p3"]
    assert out == ["p1"]


def test_cap_breaks_ties_on_coverage_then_accession():
    single = {"p1": 9, "p2": 9, "p3": 9, "p4": 9}
    info = _info(p1=60, p2=90, p3=90, p4=70)
    kept, out = cap_targets(["p4", "p3", "p2", "p1"], single, info, 2)
    assert kept == ["p3", "p2"]         # original order kept
    assert out == ["p4", "p1"]          # in rank order


def test_cap_ranks_missing_info_as_zero_coverage():
    kept, _ = cap_targets(["p1", "p2"], {"p1": 5, "p2": 5}, _info(p2=1), 1)
    assert kept == ["p2"]


@pytest.mark.parametrize("max_hmms", [0, None, 3, 10])
def test_cap_not_applied_when_off_or_not_reached(max_hmms):
    kept, out = cap_targets(["p1", "p2", "p3"], {}, {}, max_hmms)
    assert kept == ["p1", "p2", "p3"]
    assert out == []


def test_cap_applies_after_shared_protein_conflicts():
    """
    Capping first would pick p1 and p2 and then lose p2 to its conflict with p1,
    leaving one target under a cap of two.
    """
    genomes = [f"g{i}" for i in range(10)]
    hits = {g: {"p1": 1, "p2": 1, "p3": 1} for g in genomes}
    hits["g0"]["p3"] = 2        # p3 single-copy in 9, the others in 10
    sharing = {shared_pair_key("p1", "p2"): 10}

    selection = select_scg_targets(hits, genomes, ["p1", "p2", "p3"], sharing, {},
                                   90, max_hmms=2)
    assert selection.targets == ["p1", "p3"]
    assert [ex.dropped_acc for ex in selection.shared_exclusions] == ["p2"]
    assert selection.capped_out == []
    assert selection.num_single_copy == 3
    assert selection.multi_copy_counts == {"p3": 1}
