import argparse
import os
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore
import gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli as cli
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import (
    TargetGenomeError,
    build_run_data,
    resolve_download_info,
)
from gtotree.utils.misc.messaging import REMOVED_GENOMES_FILENAME
from gtotree.utils.misc.stages import GenomeRemovalStage
from gtotree.tests.paths import DATA_DIR


MOTIFS = {
    "PF90001.3": "MKVLAAAL",
    "PF90002.7": "MARTKQTA",
}

_COLS = ["assembly_accession", "asm_name", "ftp_path", "organism_name"]


def _write_ncbi_table(path, rows):
    cols = {c: [r.get(c, "") for r in rows] for c in _COLS}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
        str(path))


def _phase_args(**kw):
    base = dict(ncbi_accessions=None, wanted_ref_tax=None, genbank_files=None,
                fasta_files=None, amino_acid_files=None, source="GTDB",
                target_rank=None, derep_rank="off", num_jobs=4, num_threads=1,
                min_completeness=None, max_contamination=None)
    base.update(kw)
    return argparse.Namespace(**base)


def _listing(tmp_path, name, paths):
    p = tmp_path / name
    p.write_text("\n".join(str(x) for x in paths) + "\n")
    return str(p)


################################################################################
# resolve_download_info
################################################################################

def test_resolve_download_info_finds_accessions(tmp_path):
    table = tmp_path / "ncbi.parquet"
    _write_ncbi_table(table, [
        {"assembly_accession": "GCF_000000001.1", "asm_name": "ASM1",
         "ftp_path": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/001/GCF_000000001.1_ASM1",
         "organism_name": "Testus one"},
    ])

    info, not_found = resolve_download_info(["GCF_000000001.1"], str(table))

    assert not_found == []
    entry = info["GCF_000000001.1"]
    assert entry["organism_name"] == "Testus one"
    # ftp:// is normalized to https://
    assert entry["base_link"].startswith("https://")


def test_resolve_download_info_reports_not_found(tmp_path):
    table = tmp_path / "ncbi.parquet"
    _write_ncbi_table(table, [
        {"assembly_accession": "GCF_000000001.1", "asm_name": "ASM1",
         "ftp_path": "https://example.org/GCF_000000001.1_ASM1",
         "organism_name": "Testus one"},
    ])

    info, not_found = resolve_download_info(
        ["GCF_000000001.1", "GCF_999999999.9"], str(table))

    assert list(info) == ["GCF_000000001.1"]
    assert not_found == ["GCF_999999999.9"]


def test_resolve_download_info_raises_without_table(tmp_path):
    with pytest.raises(TargetGenomeError, match="gtt data get ncbi"):
        resolve_download_info(["GCF_000000001.1"], str(tmp_path / "missing.parquet"))


################################################################################
# phase_resolve_genomes
################################################################################

def _run_data(args, tmp_path):
    out_dir = str(tmp_path / "out"); os.makedirs(out_dir, exist_ok=True)
    work_dir = str(tmp_path / "work"); os.makedirs(work_dir, exist_ok=True)
    return build_run_data(args, out_dir, work_dir)


def _resolved(args, tmp_path):
    """Phase 1's output: the RunData carrying the whole input genome set."""
    run_data, _selections = cli.phase_resolve_genomes(args, _run_data(args, tmp_path))
    return run_data


def test_phase_resolve_genomes_from_accessions_file(tmp_path):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000000001.1\nGCF_000000002.1\n")

    run_data = _resolved(_phase_args(ncbi_accessions=str(accs)), tmp_path)

    assert run_data.get_input_ncbi_accs() == ["GCF_000000001.1", "GCF_000000002.1"]
    assert run_data.get_wanted_ref_tax_accs() == []


def test_phase_resolve_genomes_from_local_files(tmp_path):
    aa = tmp_path / "g1.faa"
    aa.write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    run_data = _resolved(args, tmp_path)

    assert run_data.get_input_ncbi_accs() == []
    assert [g.id for g in run_data.amino_acid_files] == ["g1"]


def test_phase_resolve_genomes_combines_sources(tmp_path):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000000001.1\n")
    aa = tmp_path / "g1.faa"
    aa.write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(ncbi_accessions=str(accs),
                       amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    run_data = _resolved(args, tmp_path)

    assert [g.id for g in run_data.all_input_genomes] == ["GCF_000000001.1", "g1"]


def test_phase_resolve_genomes_raises_when_nothing_resolves(tmp_path):
    args = _phase_args()
    with pytest.raises(GenSCGHMMsError, match="No input genomes"):
        cli.phase_resolve_genomes(args, _run_data(args, tmp_path))


def test_phase_resolve_genomes_keeps_missing_local_files_as_removed(tmp_path):
    """
    A listed-but-absent file is a genome that left the run, not a nothing. Keeping it
    is what puts it in removed-genomes.tsv with a stage and a reason.
    """
    real = tmp_path / "g1.faa"
    real.write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [real, tmp_path / "ghost.faa"]))
    run_data = _resolved(args, tmp_path)

    ghost = {g.id: g for g in run_data.all_input_genomes}["ghost"]
    assert ghost.removed
    assert ghost.removed_at == GenomeRemovalStage.AMINO_ACID_PREP


################################################################################
# phase_get_amino_acids
################################################################################

def _amino_acids(args, tmp_path, work=None):
    """Phase 1 then phase 2, the way the driver runs them."""
    run_data = _resolved(args, tmp_path)
    work = work or str(tmp_path / "work")
    os.makedirs(work, exist_ok=True)
    # the accession lookup is part of phase 1 now, and phase 2 takes what it produced
    to_fetch, download_info = cli._resolve_accession_downloads(run_data)
    combined, kept = cli.phase_get_amino_acids(run_data, work, args,
                                               to_fetch, download_info)
    return run_data, combined, kept


def _removed(run_data):
    return {g.id: g for g in run_data.all_input_genomes if g.removed}


def test_phase_get_amino_acids_local_only(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for name in ("g1", "g2"):
        (genome_dir / f"{name}.faa").write_text(
            ">p1\n" + MOTIFS["PF90001.3"] + "\n>p2\n" + MOTIFS["PF90002.7"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", sorted(genome_dir.glob("*.faa"))))

    run_data, combined, kept = _amino_acids(args, tmp_path)

    assert kept == ["g1", "g2"]
    assert _removed(run_data) == {}

    # The combined fasta is written in COMPLETION order, not input order -- the pool
    # applies results as they land. That's by design and harmless, because genome
    # identity travels in each header and the output tables are built from the
    # input-ordered `kept` list (asserted above). So check content, not sequence.
    headers = [l.strip() for l in open(combined) if l.startswith(">")]
    assert sorted(headers) == [">g1_1", ">g1_2", ">g2_1", ">g2_2"]

    # each genome's own proteins must still be contiguous and correctly numbered
    for genome_id in ("g1", "g2"):
        own = [h for h in headers if h.startswith(f">{genome_id}_")]
        assert own == [f">{genome_id}_1", f">{genome_id}_2"]


def test_phase_get_amino_acids_preserves_input_order(tmp_path):
    """
    The pool completes out of order, so kept_ids must come back in input order --
    otherwise the output tables would vary run to run.
    """
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    names = [f"g{i:02d}" for i in range(12)]
    for name in names:
        (genome_dir / f"{name}.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / f"{n}.faa" for n in names]))

    _run, _combined, kept = _amino_acids(args, tmp_path)

    assert kept == names


def test_phase_get_amino_acids_records_failures(tmp_path):
    """A genome that can't be processed is reported, and the run continues."""
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "good.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")
    (genome_dir / "empty.faa").write_text("")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "good.faa", genome_dir / "empty.faa"]))

    run_data, _combined, kept = _amino_acids(args, tmp_path)

    assert kept == ["good"]
    assert list(_removed(run_data)) == ["empty"]
    assert _removed(run_data)["empty"].removed_at == \
        GenomeRemovalStage.AMINO_ACID_PREP


def test_phase_get_amino_acids_writes_the_removed_genomes_report(tmp_path):
    """
    The losses have to reach the shared report file, not just the in-memory RunData --
    that file is the only place the user ever sees them.
    """
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "good.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")
    (genome_dir / "empty.faa").write_text("")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "good.faa", genome_dir / "empty.faa",
                             tmp_path / "ghost.faa"]))

    run_data, _combined, _kept = _amino_acids(args, tmp_path)

    report = os.path.join(run_data.run_files_dir, REMOVED_GENOMES_FILENAME)
    rows = [line.rstrip("\n").split("\t") for line in open(report)]

    assert rows[0] == ["genome_id", "input", "source", "stage_removed",
                       "reason_removed"]
    assert {r[0] for r in rows[1:]} == {"empty", "ghost"}
    assert {r[3] for r in rows[1:]} == {GenomeRemovalStage.AMINO_ACID_PREP}


def test_phase_get_amino_acids_raises_when_all_fail(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "empty.faa").write_text("")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "empty.faa"]))

    with pytest.raises(GenSCGHMMsError, match="No amino-acid sequences"):
        _amino_acids(args, tmp_path)


def test_phase_get_amino_acids_raises_when_nothing_resolvable(tmp_path):
    """Every input genome absent up front leaves nothing to even attempt."""
    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [tmp_path / "ghost.faa"]))

    with pytest.raises(GenSCGHMMsError, match="None of the target genomes"):
        _amino_acids(args, tmp_path)


################################################################################
# phase_search
################################################################################

def test_phase_search_returns_hit_counts(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "g1.faa").write_text(
        ">p1\n" + MOTIFS["PF90001.3"] + "\n>p2\n" + MOTIFS["PF90002.7"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "g1.faa"]))

    work = tmp_path / "work"
    _run, combined, kept = _amino_acids(args, tmp_path, work=str(work))

    # phase_search takes the genome count to size a per-genome progress bar; the caller
    # already knows it, so no pre-pass over the combined fasta is needed
    hits = cli.phase_search(str(DATA_DIR / "mock-pfams.hmm"), combined, len(kept), args)

    assert hits["g1"]["PF90001.3"] == 1
    assert hits["g1"]["PF90002.7"] == 1


def test_phase_search_checkpoints_into_the_work_dir(tmp_path):
    """
    The phase has to actually place a checkpoint where a resumed run will look for it;
    the library-level machinery is useless if the CLI never wires a path into it.
    """
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "g1.faa").write_text(
        ">p1\n" + MOTIFS["PF90001.3"] + "\n>p2\n" + MOTIFS["PF90002.7"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "g1.faa"]))

    work = tmp_path / "work"
    _run, combined, kept = _amino_acids(args, tmp_path, work=str(work))

    cli.phase_search(str(DATA_DIR / "mock-pfams.hmm"), combined, len(kept), args,
                     work_dir=str(work))

    assert (work / cli.SEARCH_CHECKPOINT_FILENAME).is_file()


def test_phase_search_without_a_work_dir_writes_no_checkpoint(tmp_path):
    """Checkpointing is opt-in on a work dir being supplied, and must degrade quietly."""
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "g1.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "g1.faa"]))

    work = tmp_path / "work"
    _run, combined, kept = _amino_acids(args, tmp_path, work=str(work))

    hits = cli.phase_search(str(DATA_DIR / "mock-pfams.hmm"), combined, len(kept), args)
    assert hits["g1"]["PF90001.3"] == 1
    assert not (work / cli.SEARCH_CHECKPOINT_FILENAME).exists()


################################################################################
# NCBI lookup removals
################################################################################

def test_accessions_missing_from_ncbi_are_removed_at_the_lookup_stage(tmp_path,
                                                                     monkeypatch,
                                                                     capsys):
    """
    An accession NCBI no longer lists leaves at `ncbi-lookup`, not at `ncbi-download`:
    nothing was ever attempted for it. The distinction is the whole point of the
    stage_removed column, and it's the same one the search subcommands draw.
    """
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000000001.1\nGCF_999999999.9\n")

    monkeypatch.setattr(
        cli, "resolve_download_info",
        lambda ids: ({"GCF_000000001.1": {"base_link": "https://example.org/x",
                                          "organism_name": "Testus one"}},
                     ["GCF_999999999.9"]))

    args = _phase_args(ncbi_accessions=str(accs))
    run_data = _resolved(args, tmp_path)

    to_fetch, info = cli._resolve_accession_downloads(run_data)

    assert [gd.id for gd in to_fetch] == ["GCF_000000001.1"]
    assert list(info) == ["GCF_000000001.1"]

    missing = {gd.id: gd for gd in run_data.ncbi_accs}["GCF_999999999.9"]
    assert missing.removed_at == GenomeRemovalStage.NCBI_LOOKUP
    assert missing.acc_was_found is False

    # organism names ride on the genome now rather than in a parallel dict
    found = {gd.id: gd for gd in run_data.ncbi_accs}["GCF_000000001.1"]
    assert found.organism_name == "Testus one"

    # the message matches what `gtt search-annotations` print for the
    # same lookup, since it's the same phase-1 step in all three
    out = capsys.readouterr().out
    assert "1 accession(s) not found at NCBI" in out
    assert REMOVED_GENOMES_FILENAME in out
