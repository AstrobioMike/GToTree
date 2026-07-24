import argparse
import io
import os
from pathlib import Path

import pytest # type: ignore

from gtotree.utils.hmms.gen_scg_hmms import (
    count_single_copy_hits,
    load_coverage_filtered_pfams,
    read_hmm_accessions,
    write_filtered_pfam_hmms,
)
from gtotree.utils.hmms.gen_scg_hmms_genomes import relabel_and_append
from gtotree.utils.hmms.gen_scg_hmms_local import build_local_genomes, process_local_genome
from gtotree.utils.hmms.gen_scg_hmms_search import search_profiles
from gtotree.utils.hmms import gen_scg_hmms_outputs as outputs

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
MOCK_PFAM_HMM = DATA_DIR / "mock-pfams.hmm"
MOCK_PFAM_INFO = DATA_DIR / "mock-pfamA.txt"

MOTIFS = {
    "PF90001.3": "MKVLAAAL",
    "PF90002.7": "MARTKQTA",
    "PF90003.1": "MSDKIIHL",
    "PF90004.2": "MAHHWWGS",
}


def _write_aa_file(path, accs):
    with open(path, "w") as f:
        for i, acc in enumerate(accs, 1):
            f.write(f">p{i}\n{MOTIFS[acc]}\n")
    return path


def _listing(tmp_path, name, paths):
    p = tmp_path / name
    p.write_text("\n".join(str(x) for x in paths) + "\n")
    return str(p)


def test_full_pipeline_from_local_amino_acid_files(tmp_path):
    """
    Local files -> combined fasta -> search -> single-copy -> outputs.

    g1/g2/g3 each carry one copy of PF90001/2/3. g3 additionally has a SECOND copy of
    PF90001 and no PF90004, so at 90%:
      PF90001 -> single-copy in 2/3 (67%)  -> dropped
      PF90002 -> single-copy in 3/3        -> kept
      PF90003 -> single-copy in 3/3        -> kept
      PF90004 -> single-copy in 2/3 (67%)  -> dropped
    """
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    _write_aa_file(genome_dir / "g1.faa", list(MOTIFS))
    _write_aa_file(genome_dir / "g2.faa", list(MOTIFS))
    _write_aa_file(genome_dir / "g3.faa",
                   ["PF90001.3", "PF90002.7", "PF90003.1", "PF90001.3"])

    args = argparse.Namespace(
        genbank_files=None, fasta_files=None,
        amino_acid_files=_listing(tmp_path, "aa.txt",
                                  sorted(genome_dir.glob("*.faa"))))

    work = tmp_path / "work"
    work.mkdir()

    # --- genome stage -------------------------------------------------------
    genomes, missing = build_local_genomes(args)
    assert missing == []
    assert [g.id for g in genomes] == ["g1", "g2", "g3"]

    combined_path = work / "all-target-proteins.faa"
    kept_ids = []
    with open(combined_path, "w") as combined:
        for gd in genomes:
            aa_path, _ = process_local_genome(gd, str(work))
            relabel_and_append(gd.id, aa_path, combined)
            kept_ids.append(gd.id)

    # --- pfam stage ---------------------------------------------------------
    pfam_info = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO), min_coverage=10.0)
    filtered_hmm = work / "filtered.hmm"
    filtered_accs = write_filtered_pfam_hmms(
        str(MOCK_PFAM_HMM), set(pfam_info), str(filtered_hmm))

    # --- search stage -------------------------------------------------------
    hits = search_profiles(str(filtered_hmm), str(combined_path), threads=1)
    assert set(hits) == {"g1", "g2", "g3"}
    assert hits["g3"]["PF90001.3"] == 2

    # --- single-copy --------------------------------------------------------
    wanted, per_genome = count_single_copy_hits(hits, kept_ids, filtered_accs, 90)
    assert set(wanted) == {"PF90002.7", "PF90003.1"}

    # every genome survives this stage, even the one that lost its Pfams
    assert sorted(per_genome) == ["g1", "g2", "g3"]

    # --- outputs ------------------------------------------------------------
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    final_hmm = out_dir / "Test-Set.hmm"
    write_filtered_pfam_hmms(str(filtered_hmm), wanted, str(final_hmm))
    assert sorted(read_hmm_accessions(str(final_hmm))) == sorted(wanted)

    outputs.write_scg_targets_info(str(out_dir), wanted, pfam_info)
    outputs.write_hit_counts(str(out_dir), kept_ids, filtered_accs, per_genome)
    outputs.write_target_genomes(str(out_dir), kept_ids,
                                 {g.id: g.source for g in genomes})
    outputs.write_pfam_version(str(out_dir), "38.2-mock")
    assert outputs.write_missed_accessions(str(out_dir), []) is None

    produced = sorted(p.name for p in out_dir.iterdir())
    assert produced == [
        "Pfam-hit-counts.tsv",
        "SCG-targets-info.tsv",
        "Test-Set.hmm",
        "pfam-version-used.txt",
        "target-genomes.tsv",
    ]


def test_pipeline_lowering_threshold_keeps_more(tmp_path):
    """The same data at a looser threshold retains the Pfams g3 spoiled."""
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    _write_aa_file(genome_dir / "g1.faa", list(MOTIFS))
    _write_aa_file(genome_dir / "g2.faa", list(MOTIFS))
    _write_aa_file(genome_dir / "g3.faa",
                   ["PF90001.3", "PF90002.7", "PF90003.1", "PF90001.3"])

    args = argparse.Namespace(
        genbank_files=None, fasta_files=None,
        amino_acid_files=_listing(tmp_path, "aa.txt",
                                  sorted(genome_dir.glob("*.faa"))))
    work = tmp_path / "work"; work.mkdir()

    genomes, _ = build_local_genomes(args)
    combined_path = work / "combined.faa"
    with open(combined_path, "w") as combined:
        for gd in genomes:
            aa_path, _ = process_local_genome(gd, str(work))
            relabel_and_append(gd.id, aa_path, combined)

    pfam_info = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO), min_coverage=10.0)
    filtered_hmm = work / "filtered.hmm"
    filtered_accs = write_filtered_pfam_hmms(
        str(MOCK_PFAM_HMM), set(pfam_info), str(filtered_hmm))
    hits = search_profiles(str(filtered_hmm), str(combined_path), threads=1)

    ids = [g.id for g in genomes]
    loose, _ = count_single_copy_hits(hits, ids, filtered_accs, 66)
    assert set(loose) == set(MOTIFS)


def test_mixed_sources_share_one_header_convention(tmp_path):
    """
    A genbank, a fasta, and an amino-acid input must all end up with `{id}_{n}` headers
    so genome ids recover uniformly downstream, regardless of where they came from.
    """
    genbank = tmp_path / "gb1.gb"
    genbank.write_text(
        "LOCUS       TESTCTG                   60 bp    DNA     linear   BCT 01-JAN-2020\n"
        "DEFINITION  t.\n"
        "ACCESSION   TESTCTG\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..60\n"
        "     CDS             1..60\n"
        '                     /translation="MKVLAAALRAEL"\n'
        "ORIGIN\n"
        "        1 atgaaagttc tggccgccgc actgctggcg ctgggcgcct gcagcgcgcc ggaaccggtg\n"
        "//\n")
    aa = _write_aa_file(tmp_path / "aa1.faa", ["PF90002.7", "PF90003.1"])

    args = argparse.Namespace(
        genbank_files=_listing(tmp_path, "gb.txt", [genbank]),
        fasta_files=None,
        amino_acid_files=_listing(tmp_path, "aa.txt", [aa]))

    work = tmp_path / "work"; work.mkdir()
    genomes, _ = build_local_genomes(args)

    combined = io.StringIO()
    for gd in genomes:
        aa_path, _ = process_local_genome(gd, str(work))
        relabel_and_append(gd.id, aa_path, combined)

    headers = [l[1:] for l in combined.getvalue().splitlines() if l.startswith(">")]
    assert headers == ["gb1_1", "aa1_1", "aa1_2"]

    # and no prodigal stop characters survive into the combined file
    assert "*" not in combined.getvalue()
