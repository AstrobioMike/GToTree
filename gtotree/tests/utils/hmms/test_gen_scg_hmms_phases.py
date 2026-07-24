import argparse
import os
from pathlib import Path

import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest # type: ignore

import gtotree.utils.hmms.gen_scg_hmms_cli as cli
import gtotree.utils.hmms.gen_scg_hmms_genomes as genomes_mod
from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms_genomes import (
    TargetGenomeError,
    resolve_download_info,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"

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
    base = dict(target_accessions=None, wanted_ref_tax=None, genbank_files=None,
                fasta_files=None, amino_acid_files=None, source="GTDB",
                target_rank=None, derep_rank="off", num_jobs=4, num_threads=1)
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

def test_phase_resolve_genomes_from_accessions_file(tmp_path, capsys):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000000001.1\nGCF_000000002.1\n")

    args = _phase_args(target_accessions=str(accs))
    accessions, sources, local, missing = cli.phase_resolve_genomes(args)

    assert accessions == ["GCF_000000001.1", "GCF_000000002.1"]
    assert set(sources.values()) == {"input-accessions"}
    assert local == [] and missing == []


def test_phase_resolve_genomes_from_local_files(tmp_path):
    aa = tmp_path / "g1.faa"
    aa.write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    accessions, sources, local, missing = cli.phase_resolve_genomes(args)

    assert accessions == []
    assert [g.id for g in local] == ["g1"]
    assert missing == []


def test_phase_resolve_genomes_combines_sources(tmp_path):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCF_000000001.1\n")
    aa = tmp_path / "g1.faa"
    aa.write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(target_accessions=str(accs),
                       amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    accessions, sources, local, missing = cli.phase_resolve_genomes(args)

    assert accessions == ["GCF_000000001.1"]
    assert [g.id for g in local] == ["g1"]


def test_phase_resolve_genomes_raises_when_nothing_resolves(tmp_path):
    aa_listing = tmp_path / "aa.txt"
    aa_listing.write_text(str(tmp_path / "ghost.faa") + "\n")

    args = _phase_args(amino_acid_files=str(aa_listing))
    with pytest.raises(GenSCGHMMsError, match="No target genomes"):
        cli.phase_resolve_genomes(args)


def test_phase_resolve_genomes_reports_missing_local_files(tmp_path):
    real = tmp_path / "g1.faa"
    real.write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [real, tmp_path / "ghost.faa"]))
    _, _, local, missing = cli.phase_resolve_genomes(args)

    assert [g.id for g in local] == ["g1"]
    assert [m[0] for m in missing] == ["ghost"]


################################################################################
# phase_get_amino_acids
################################################################################

def test_phase_get_amino_acids_local_only(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for name in ("g1", "g2"):
        (genome_dir / f"{name}.faa").write_text(
            ">p1\n" + MOTIFS["PF90001.3"] + "\n>p2\n" + MOTIFS["PF90002.7"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", sorted(genome_dir.glob("*.faa"))))
    _, _, local, missing = cli.phase_resolve_genomes(args)

    work = tmp_path / "work"
    work.mkdir()
    combined, kept, missed, organisms, sources_extra = cli.phase_get_amino_acids(
        [], local, missing, str(work), args)

    assert kept == ["g1", "g2"]
    assert missed == []
    assert sources_extra == {"g1": "amino-acid", "g2": "amino-acid"}

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
    The pool completes out of order, so kept_ids must be re-sorted to input order --
    otherwise the output tables would vary run to run.
    """
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    names = [f"g{i:02d}" for i in range(12)]
    for name in names:
        (genome_dir / f"{name}.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    listing = _listing(tmp_path, "aa.txt",
                       [genome_dir / f"{n}.faa" for n in names])
    args = _phase_args(amino_acid_files=listing)
    _, _, local, missing = cli.phase_resolve_genomes(args)

    work = tmp_path / "work"
    work.mkdir()
    _, kept, _, _, _ = cli.phase_get_amino_acids([], local, missing, str(work), args)

    assert kept == names


def test_phase_get_amino_acids_records_failures(tmp_path):
    """A genome that can't be processed is reported, and the run continues."""
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "good.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")
    (genome_dir / "empty.faa").write_text("")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "good.faa", genome_dir / "empty.faa"]))
    _, _, local, missing = cli.phase_resolve_genomes(args)

    work = tmp_path / "work"
    work.mkdir()
    _, kept, missed, _, _ = cli.phase_get_amino_acids(
        [], local, missing, str(work), args)

    assert kept == ["good"]
    assert [m[0] for m in missed] == ["empty"]


def test_phase_get_amino_acids_raises_when_all_fail(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "empty.faa").write_text("")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "empty.faa"]))
    _, _, local, missing = cli.phase_resolve_genomes(args)

    work = tmp_path / "work"
    work.mkdir()
    with pytest.raises(GenSCGHMMsError, match="No amino-acid sequences"):
        cli.phase_get_amino_acids([], local, missing, str(work), args)


def test_phase_get_amino_acids_carries_missing_into_report(tmp_path):
    """Files listed but absent must reach missed-accessions.tsv, not vanish."""
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    (genome_dir / "good.faa").write_text(">p1\n" + MOTIFS["PF90001.3"] + "\n")

    args = _phase_args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [genome_dir / "good.faa", tmp_path / "ghost.faa"]))
    _, _, local, missing = cli.phase_resolve_genomes(args)

    work = tmp_path / "work"
    work.mkdir()
    _, kept, missed, _, _ = cli.phase_get_amino_acids(
        [], local, missing, str(work), args)

    assert kept == ["good"]
    assert [m[0] for m in missed] == ["ghost"]


def test_phase_get_amino_acids_raises_when_nothing_resolvable(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    with pytest.raises(GenSCGHMMsError, match="None of the target genomes"):
        cli.phase_get_amino_acids([], [], [], str(work), _phase_args())


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
    _, _, local, missing = cli.phase_resolve_genomes(args)

    work = tmp_path / "work"
    work.mkdir()
    combined, kept, _, _, _ = cli.phase_get_amino_acids(
        [], local, missing, str(work), args)

    hits = cli.phase_search(str(DATA_DIR / "mock-pfams.hmm"), combined,
                            ["PF90001.3", "PF90002.7"], args)

    assert hits["g1"]["PF90001.3"] == 1
    assert hits["g1"]["PF90002.7"] == 1
