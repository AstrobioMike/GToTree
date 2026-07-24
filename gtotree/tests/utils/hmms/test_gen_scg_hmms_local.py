import argparse
import gzip
import hashlib
import os

import pytest # type: ignore

import gtotree.utils.hmms.gen_scg_hmms_local as mod
from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms_genomes import TargetGenomeError
from gtotree.utils.hmms.gen_scg_hmms_local import (
    SOURCE_AMINO_ACID,
    SOURCE_FASTA,
    SOURCE_GENBANK,
    build_local_genomes,
    process_local_genome,
    read_paths_file,
)


# a GenBank record with three CDS features, one of them flagged frameshifted -- the
# reused main-run extractor is expected to skip that one
GENBANK = """LOCUS       TESTCTG                  300 bp    DNA     linear   BCT 01-JAN-2020
DEFINITION  test contig.
ACCESSION   TESTCTG
FEATURES             Location/Qualifiers
     source          1..300
     CDS             1..90
                     /translation="MKVLAAALRAELMKVLAAALLA"
     CDS             100..180
                     /translation="MARTKQTALATKAARKS"
     CDS             200..260
                     /note="frameshifted"
                     /translation="MBADSHOULDBESKIPPED"
ORIGIN
        1 atgaaagttc tggccgccgc actgctggcg ctgggcgcct gcagcgcgcc ggaaccggtg
       61 cgcgcggaac tgatgaaagt tctggccgcc gcactgctgg cgctgggcgc ctgcagcgcg
      121 ccggaaccgg tgcgcgcgga actgatgaaa gttctggccg ccgcactgct ggcgctgggc
      181 gcctgcagcg cgccggaacc ggtgcgcgcg gaactgatga aagttctggc cgccgcactg
      241 ctggcgctgg gcgcctgcag cgcgccggaa ccggtgcgcg cggaactgat gaaagttctg
//
"""

AMINO_ACIDS = ">p1 some description\nMKVLAAALRAEL\n>p2\nMARTKQTALATKAARKS\n"


def _listing(tmp_path, name, paths):
    path = tmp_path / name
    path.write_text("\n".join(str(p) for p in paths) + "\n")
    return str(path)


def _args(**kw):
    base = dict(genbank_files=None, fasta_files=None, amino_acid_files=None)
    base.update(kw)
    return argparse.Namespace(**base)


################################################################################
# paths listing
################################################################################

def test_read_paths_file(tmp_path):
    listing = _listing(tmp_path, "l.txt", ["/a/x.gb", "/a/y.gb"])
    assert read_paths_file(listing, "genbank-files") == ["/a/x.gb", "/a/y.gb"]


def test_read_paths_file_skips_blanks_and_comments(tmp_path):
    path = tmp_path / "l.txt"
    path.write_text("# note\n\n/a/x.gb\n  \n/a/y.gb\n")
    assert read_paths_file(str(path), "genbank-files") == ["/a/x.gb", "/a/y.gb"]


def test_read_paths_file_dedupes(tmp_path):
    path = tmp_path / "l.txt"
    path.write_text("/a/x.gb\n/a/y.gb\n/a/x.gb\n")
    assert read_paths_file(str(path), "genbank-files") == ["/a/x.gb", "/a/y.gb"]


def test_read_paths_file_missing(tmp_path):
    with pytest.raises(GenSCGHMMsError, match="can't be found"):
        read_paths_file(str(tmp_path / "nope.txt"), "fasta-files")


def test_read_paths_file_empty(tmp_path):
    path = tmp_path / "l.txt"
    path.write_text("# nothing\n")
    with pytest.raises(GenSCGHMMsError, match="no paths"):
        read_paths_file(str(path), "fasta-files")


################################################################################
# building GenomeData
################################################################################

def test_build_local_genomes_tags_sources(tmp_path):
    gb = tmp_path / "g1.gb"; gb.write_text(GENBANK)
    fa = tmp_path / "g2.fna"; fa.write_text(">c\nATGAAA\n")
    aa = tmp_path / "g3.faa"; aa.write_text(AMINO_ACIDS)

    args = _args(genbank_files=_listing(tmp_path, "gb.txt", [gb]),
                 fasta_files=_listing(tmp_path, "fa.txt", [fa]),
                 amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))

    genomes, missing = build_local_genomes(args)

    assert missing == []
    assert [(g.id, g.source) for g in genomes] == [
        ("g1", SOURCE_GENBANK), ("g2", SOURCE_FASTA), ("g3", SOURCE_AMINO_ACID)]


def test_build_local_genomes_reports_missing_up_front(tmp_path):
    """
    A typo in a long listing should be reported before any expensive work starts,
    not midway through downloading and gene-calling.
    """
    real = tmp_path / "g1.faa"; real.write_text(AMINO_ACIDS)
    args = _args(amino_acid_files=_listing(
        tmp_path, "aa.txt", [real, tmp_path / "ghost.faa"]))

    genomes, missing = build_local_genomes(args)

    assert [g.id for g in genomes] == ["g1"]
    assert len(missing) == 1
    assert missing[0][0] == "ghost"
    assert "not found" in missing[0][1]


def test_build_local_genomes_with_no_inputs():
    genomes, missing = build_local_genomes(_args())
    assert genomes == [] and missing == []


def test_build_local_genomes_strips_gz_for_id(tmp_path):
    path = tmp_path / "g1.gb.gz"
    with gzip.open(path, "wt") as f:
        f.write(GENBANK)
    args = _args(genbank_files=_listing(tmp_path, "gb.txt", [path]))

    genomes, _ = build_local_genomes(args)
    assert [g.id for g in genomes] == ["g1"]


################################################################################
# processing each source type
################################################################################

def test_process_amino_acid_file(tmp_path):
    aa = tmp_path / "g3.faa"; aa.write_text(AMINO_ACIDS)
    args = _args(amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    gd = build_local_genomes(args)[0][0]

    work = tmp_path / "work"; work.mkdir()
    out_path, used_prodigal = process_local_genome(gd, str(work))

    assert used_prodigal is False
    headers = [l.strip() for l in open(out_path) if l.startswith(">")]
    assert headers == [">g3_1", ">g3_2"]


def test_process_genbank_uses_existing_translations(tmp_path):
    """
    GenBank CDS translations are the submitter's own gene calls, so they're preferred
    over re-calling with prodigal. The frameshifted CDS is dropped by the reused
    main-run extractor, leaving 2 of 3.
    """
    gb = tmp_path / "g1.gb"; gb.write_text(GENBANK)
    args = _args(genbank_files=_listing(tmp_path, "gb.txt", [gb]))
    gd = build_local_genomes(args)[0][0]

    work = tmp_path / "work"; work.mkdir()
    out_path, used_prodigal = process_local_genome(gd, str(work))

    assert used_prodigal is False
    headers = [l.strip() for l in open(out_path) if l.startswith(">")]
    assert headers == [">g1_1", ">g1_2"]


def test_process_amino_acid_gz_is_handled(tmp_path):
    path = tmp_path / "g3.faa.gz"
    with gzip.open(path, "wt") as f:
        f.write(AMINO_ACIDS)
    args = _args(amino_acid_files=_listing(tmp_path, "aa.txt", [path]))
    gd = build_local_genomes(args)[0][0]

    work = tmp_path / "work"; work.mkdir()
    out_path, _ = process_local_genome(gd, str(work))

    assert len([l for l in open(out_path) if l.startswith(">")]) == 2


def test_gzipped_input_does_not_clobber_user_file(tmp_path):
    """
    Regression test. `general.gunzip_if_needed` decompresses ALONGSIDE the original
    (stripping `.gz`), so if a user has both `x.faa.gz` and `x.faa`, decompressing and
    then cleaning up the temp would destroy their uncompressed file. This helper stages
    into the working dir instead and must never touch the input directory.
    """
    gz_path = tmp_path / "g3.faa.gz"
    with gzip.open(gz_path, "wt") as f:
        f.write(AMINO_ACIDS)

    # a DIFFERENT file that happens to share the decompressed name
    sibling = tmp_path / "g3.faa"
    sibling.write_text(">mine\nMMMMM\n")
    before = hashlib.md5(sibling.read_bytes()).hexdigest()

    args = _args(amino_acid_files=_listing(tmp_path, "aa.txt", [gz_path]))
    gd = build_local_genomes(args)[0][0]

    work = tmp_path / "work"; work.mkdir()
    process_local_genome(gd, str(work))

    assert sibling.exists(), "user's uncompressed file was deleted"
    assert hashlib.md5(sibling.read_bytes()).hexdigest() == before, \
        "user's uncompressed file was overwritten"


def test_gzip_staging_is_cleaned_up(tmp_path):
    path = tmp_path / "g3.faa.gz"
    with gzip.open(path, "wt") as f:
        f.write(AMINO_ACIDS)
    args = _args(amino_acid_files=_listing(tmp_path, "aa.txt", [path]))
    gd = build_local_genomes(args)[0][0]

    work = tmp_path / "work"; work.mkdir()
    process_local_genome(gd, str(work))

    leftovers = [f for f in os.listdir(work) if f.startswith("_unzipped_")]
    assert leftovers == []


def test_process_empty_amino_acid_file_raises(tmp_path):
    aa = tmp_path / "g3.faa"; aa.write_text("")
    args = _args(amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    gd = build_local_genomes(args)[0][0]

    work = tmp_path / "work"; work.mkdir()
    with pytest.raises(TargetGenomeError):
        process_local_genome(gd, str(work))


def test_process_fasta_requires_prodigal(tmp_path, monkeypatch):
    """Nucleotide input always needs gene calling."""
    fa = tmp_path / "g2.fna"; fa.write_text(">c1\n" + "ATGAAAGTTCTGGCC" * 40 + "\n")
    args = _args(fasta_files=_listing(tmp_path, "fa.txt", [fa]))
    gd = build_local_genomes(args)[0][0]

    called = {}

    def fake_prodigal(nt_path, out_path):
        called["nt"] = nt_path
        with open(out_path, "w") as f:
            f.write(">pred1\nMKVLAAALLA\n>pred2\nMARTKQTARK\n")

    monkeypatch.setattr(mod, "run_prodigal", fake_prodigal)

    work = tmp_path / "work"; work.mkdir()
    out_path, used_prodigal = process_local_genome(gd, str(work))

    assert used_prodigal is True
    assert called["nt"].endswith("g2.fna")
    headers = [l.strip() for l in open(out_path) if l.startswith(">")]
    assert headers == [">g2_1", ">g2_2"]


def test_unrecognized_source_raises(tmp_path):
    aa = tmp_path / "g3.faa"; aa.write_text(AMINO_ACIDS)
    args = _args(amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))
    gd = build_local_genomes(args)[0][0]
    gd.source = "not-a-real-source"

    work = tmp_path / "work"; work.mkdir()
    with pytest.raises(TargetGenomeError, match="unrecognized genome source"):
        process_local_genome(gd, str(work))
