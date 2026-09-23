"""
One gene can be the chosen hit for two SCG targets (e.g., domains of one fused protein).
Pulling it for both would put the same sequence into the tree twice, so it's kept only
for the target it scores best against, and the conflict is reported.
"""

import pytest # type: ignore

from gtotree.utils.hmms.hmm_searching import (
    COMBINED_SHARED_HITS_FILENAME,
    GENOME_SHARED_HITS_FILENAME,
    get_seqs,
    no_hits_reason,
    parse_hmmer_results,
    parse_hmmer_results_with_shared,
    read_genome_shared_gene_hits,
    rebuild_combined_SCG_outputs,
    resolve_shared_gene_hits,
    write_genome_hit_counts,
    write_genome_shared_gene_hits,
)
from gtotree.utils.misc.general import RunData, SCGset, GenomeData


def _run_data(targets, best_hit_mode=False):
    rd = RunData()
    rd.SCG_targets = [SCGset.from_id(t) for t in targets]
    rd.best_hit_mode = best_hit_mode
    return rd


def _tblout(tmp_path, rows, name="SCG-hits-hmm.txt"):
    path = tmp_path / name
    with open(path, "w") as f:
        f.write("# target name  accession  query name  accession  E-value  score\n")
        for gene, scg, score in rows:
            f.write(f"{gene}\t-\t{scg}\tACC\t1e-20\t{score}\n")
    return str(path)


################################################################################
# resolving
################################################################################

def test_gene_goes_to_the_better_scoring_target(tmp_path):
    rd = _run_data(["A", "B", "C"])
    path = _tblout(tmp_path, [("g1", "A", 50.0), ("g1", "B", 300.0), ("g2", "C", 80.0)])

    counts, gene_ids, n_hits, n_unique, shared = parse_hmmer_results_with_shared(path, rd)

    assert gene_ids == {"A": None, "B": "g1", "C": "g2"}
    assert len(shared) == 1
    hit = shared[0]
    assert (hit.gene_id, hit.kept_SCG, hit.dropped_SCG) == ("g1", "B", "A")
    assert (hit.kept_score, hit.dropped_score) == (300.0, 50.0)
    # the counts still describe what the search found
    assert counts == {"A": 1, "B": 1, "C": 1}
    assert (n_hits, n_unique) == (3, 3)


def test_tie_goes_to_the_target_first_in_the_hmm_file(tmp_path):
    rd = _run_data(["B", "A"])
    path = _tblout(tmp_path, [("g1", "A", 100.0), ("g1", "B", 100.0)])

    _, gene_ids, _, _, shared = parse_hmmer_results_with_shared(path, rd)
    assert gene_ids == {"B": "g1", "A": None}
    assert shared[0].kept_SCG == "B"


def test_three_targets_on_one_gene_leave_one(tmp_path):
    rd = _run_data(["A", "B", "C"])
    path = _tblout(tmp_path, [("g1", "A", 10.0), ("g1", "B", 30.0), ("g1", "C", 20.0)])

    _, gene_ids, _, _, shared = parse_hmmer_results_with_shared(path, rd)
    assert gene_ids == {"A": None, "B": "g1", "C": None}
    assert sorted(h.dropped_SCG for h in shared) == ["A", "C"]
    assert {h.kept_SCG for h in shared} == {"B"}


def test_no_gene_is_ever_assigned_twice(tmp_path):
    rd = _run_data(["A", "B", "C", "D"])
    path = _tblout(tmp_path, [("g1", "A", 5.0), ("g1", "B", 6.0),
                              ("g2", "C", 7.0), ("g2", "D", 7.5)])
    _, gene_ids, _, _, _ = parse_hmmer_results_with_shared(path, rd)
    picked = [g for g in gene_ids.values() if g is not None]
    assert len(picked) == len(set(picked))


def test_best_hit_mode_conflicts_are_resolved_too(tmp_path):
    """In -B mode a multi-copy target takes its best hit, which can collide too."""
    rd = _run_data(["A", "B"], best_hit_mode=True)
    path = _tblout(tmp_path, [("g1", "A", 200.0), ("g2", "A", 150.0),
                              ("g1", "B", 90.0)])
    _, gene_ids, _, _, shared = parse_hmmer_results_with_shared(path, rd)
    assert gene_ids == {"A": "g1", "B": None}
    assert shared[0].dropped_SCG == "B"


def test_no_conflicts_is_a_no_op(tmp_path):
    rd = _run_data(["A", "B"])
    path = _tblout(tmp_path, [("g1", "A", 5.0), ("g2", "B", 6.0)])
    _, gene_ids, _, _, shared = parse_hmmer_results_with_shared(path, rd)
    assert gene_ids == {"A": "g1", "B": "g2"}
    assert shared == []


def test_plain_parse_is_also_conflict_free(tmp_path):
    """Callers of the 4-tuple version get the same protection."""
    rd = _run_data(["A", "B"])
    path = _tblout(tmp_path, [("g1", "A", 5.0), ("g1", "B", 6.0)])
    _, gene_ids, _, _ = parse_hmmer_results(path, rd)
    assert gene_ids == {"A": None, "B": "g1"}


def test_resolve_handles_a_missing_score_by_losing():
    gene_ids = {"A": "g1", "B": "g1"}
    shared = resolve_shared_gene_hits(gene_ids, {("g1", "B"): 1.0}, ["A", "B"])
    assert gene_ids == {"A": None, "B": "g1"}
    assert shared[0].dropped_score is None


################################################################################
# get_seqs never silently collapses a double assignment
################################################################################

def test_get_seqs_refuses_a_gene_assigned_twice(tmp_path):
    faa = tmp_path / "g.faa"
    faa.write_text(">g_1\nMKTAYIAKQR\n")
    seqs, failed = get_seqs({"A": "g_1", "B": "g_1"}, str(faa))
    assert failed is True
    assert seqs is None


def test_get_seqs_normal_case_still_works(tmp_path):
    faa = tmp_path / "g.faa"
    faa.write_text(">g_1\nMKTAYIAKQR\n>g_2\nMSSHEGGKKK\n")
    seqs, failed = get_seqs({"A": "g_1", "B": "g_2", "C": None}, str(faa))
    assert not failed
    assert seqs == {"A": "MKTAYIAKQR", "B": "MSSHEGGKKK", "C": None}


################################################################################
# reporting
################################################################################

def test_per_genome_report_round_trips(tmp_path):
    rd = _run_data(["A", "B"])
    path = _tblout(tmp_path, [("g1", "A", 50.0), ("g1", "B", 300.0)])
    *_, shared = parse_hmmer_results_with_shared(path, rd)

    out = tmp_path / GENOME_SHARED_HITS_FILENAME
    write_genome_shared_gene_hits(str(out), shared)
    assert read_genome_shared_gene_hits(str(out)) == [["g1", "B", "300", "A", "50"]]


def test_per_genome_report_is_removed_when_there_are_no_conflicts(tmp_path):
    out = tmp_path / GENOME_SHARED_HITS_FILENAME
    out.write_text("stale\n")
    write_genome_shared_gene_hits(str(out), [])
    assert not out.exists()


def _genome(gid):
    gd = GenomeData(id=gid, source="amino_acid_files", full_path=gid,
                    provided_path=gid, basename=gid)
    gd.hmm_search_done = True
    return gd


def test_rebuild_combines_reports_and_counts_losses(tmp_path):
    rd = _run_data(["A", "B"])
    rd.hmm_results_dir = str(tmp_path / "hmm")
    rd.found_SCG_seqs_dir = str(tmp_path / "found")
    rd.output_dir = str(tmp_path / "out")
    rd.run_files_dir = str(tmp_path / "run-files")
    rd.general_ext = ".faa"
    for d in (rd.found_SCG_seqs_dir, rd.output_dir, rd.run_files_dir):
        (tmp_path / d).mkdir(parents=True, exist_ok=True)

    genomes = [_genome("G1"), _genome("G2")]
    rd.all_input_genomes = genomes

    for gd in genomes:
        gdir = tmp_path / "hmm" / gd.id
        gdir.mkdir(parents=True)
        write_genome_hit_counts(str(gdir / "SCG-hit-counts.txt"), {"A": 1, "B": 1})
        path = _tblout(gdir, [(f"{gd.id}_1", "A", 50.0), (f"{gd.id}_1", "B", 300.0)])
        *_, shared = parse_hmmer_results_with_shared(path, rd)
        write_genome_shared_gene_hits(str(gdir / GENOME_SHARED_HITS_FILENAME), shared)
        (gdir / "SCG-hits.faa").write_text(">B\nMKTAYIAKQR\n")

    rebuild_combined_SCG_outputs(rd)

    assert rd.num_shared_gene_hits == 2
    by_id = {scg.id: scg for scg in rd.SCG_targets}
    assert by_id["A"].num_genomes_lost_to_shared_gene == 2
    assert by_id["B"].num_genomes_lost_to_shared_gene == 0
    assert by_id["A"].num_genomes_after_copy_filtering == 0

    combined = (tmp_path / "run-files" / COMBINED_SHARED_HITS_FILENAME).read_text()
    assert combined.splitlines()[0].startswith("genome_id\tgene_id\tkept_SCG")
    assert len(combined.splitlines()) == 3

    # and the removal reason says what actually happened, not "never single-copy"
    reason = no_hits_reason(by_id["A"], best_hit_mode=False)
    assert "scored better against another SCG target" in reason
