"""
Pfams that sit on the same protein (e.g., domains of a fused, multi-functional protein)
can each look single-copy while really being one gene. `gen-scg-hmms` records which
profiles share proteins during the search and keeps only one of any pair that commonly
does, since main GToTree extracts whole proteins and would otherwise put the same
sequence into two SCG sets.
"""

import json

from gtotree.tests.paths import DATA_DIR
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import (
    count_genomes_sharing,
    count_single_copy_hits,
    resolve_shared_protein_conflicts,
    shared_pair_key,
    split_shared_pair_key,
)
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_search import (
    search_profiles,
    shared_pairs_from_protein_hits,
)
from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_outputs as outputs
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import PfamProfileInfo


MOCK_PFAM_HMM = DATA_DIR / "mock-pfams.hmm"

MOTIFS = {
    "PF90001.3": "MKVLAAAL",
    "PF90002.7": "MARTKQTA",
    "PF90003.1": "MSDKIIHL",
    "PF90004.2": "MAHHWWGS",
}
A, B, C, D = MOTIFS


def _write(path, genomes):
    """
    genomes: {genome_id: [protein, ...]} where each protein is a list of accs whose
    motifs are concatenated into one sequence -- so ["A", "B"] is a fused protein
    """
    with open(path, "w") as f:
        for genome_id, proteins in genomes.items():
            for i, accs in enumerate(proteins, 1):
                f.write(f">{genome_id}_{i}\n{''.join(MOTIFS[a] for a in accs)}\n")
    return str(path)


################################################################################
# pair keys
################################################################################

def test_pair_key_is_order_independent_and_round_trips():
    assert shared_pair_key(B, A) == shared_pair_key(A, B)
    assert split_shared_pair_key(shared_pair_key(B, A)) == (A, B)


def test_pairs_from_protein_hits_counts_every_pair_of_a_multi_hit_protein():
    shared = shared_pairs_from_protein_hits({
        "g1_1": {A, B, C},
        "g1_2": {D},
        "g2_1": {A, B},
        "g2_7": {A, B},
    })
    assert shared == {
        "g1": {shared_pair_key(A, B): 1, shared_pair_key(A, C): 1,
               shared_pair_key(B, C): 1},
        "g2": {shared_pair_key(A, B): 2},
    }


################################################################################
# search
################################################################################

def test_search_records_fused_proteins(tmp_path):
    faa = _write(tmp_path / "t.faa", {
        "g1": [[A, B], [C], [D]],
        "g2": [[A], [B], [C], [D]],
    })
    shared = {}
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, shared_proteins=shared)

    # the per-profile hit counts are unchanged by fusion: each is still hit once
    assert hits["g1"] == {A: 1, B: 1, C: 1, D: 1}
    assert shared == {"g1": {shared_pair_key(A, B): 1}}


def test_search_without_shared_arg_behaves_as_before(tmp_path):
    faa = _write(tmp_path / "t.faa", {"g1": [[A, B], [C]]})
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)
    assert hits == {"g1": {A: 1, B: 1, C: 1}}


def test_shared_pairs_are_chunk_independent(tmp_path):
    faa = _write(tmp_path / "t.faa", {
        f"g{i}": [[A, B], [C, D], [A]] for i in range(6)
    })
    results = []
    for budget in (8, 16, 40, 10 ** 9):
        shared = {}
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=budget,
                        shared_proteins=shared)
        results.append(shared)
    assert all(r == results[0] for r in results)
    assert results[0]["g0"] == {shared_pair_key(A, B): 1, shared_pair_key(C, D): 1}


def test_resumed_search_restores_shared_pairs(tmp_path):
    faa = _write(tmp_path / "t.faa", {f"g{i}": [[A, B], [C]] for i in range(5)})

    expected_shared = {}
    expected = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                               shared_proteins=expected_shared)

    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))
    lines = [l for l in ckpt.read_text().splitlines() if l.strip()]
    assert len(lines) > 2, "need more than one chunk for this to be a real test"
    ckpt.write_text("\n".join(lines[:2]) + "\n")

    resumed_shared = {}
    resumed = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                              checkpoint_path=str(ckpt), resume=True,
                              shared_proteins=resumed_shared)
    assert resumed == expected
    assert resumed_shared == expected_shared


def test_old_format_checkpoint_is_not_resumed_from(tmp_path):
    """
    A checkpoint from before shared pairs were recorded has none, so resuming from it
    would silently lose them. It must be treated as unusable and the search redone.
    """
    faa = _write(tmp_path / "t.faa", {f"g{i}": [[A, B]] for i in range(5)})
    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))

    lines = [l for l in ckpt.read_text().splitlines() if l.strip()]
    header = json.loads(lines[0])
    header.pop("format")
    old_records = []
    for line in lines[1:]:
        rec = json.loads(line)
        rec.pop("shared")
        old_records.append(json.dumps(rec))
    ckpt.write_text("\n".join([json.dumps(header)] + old_records) + "\n")

    shared = {}
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt), resume=True, shared_proteins=shared)
    assert all(shared[g] == {shared_pair_key(A, B): 1} for g in [f"g{i}" for i in range(5)])


################################################################################
# resolving conflicts
################################################################################

def _resolve(wanted, shared_by_genome, genome_ids, hits, **kw):
    """The per-genome search record goes through the same aggregation a run does."""
    return resolve_shared_protein_conflicts(
        wanted, count_genomes_sharing(shared_by_genome, genome_ids), genome_ids, hits,
        **kw)


def _hits(genome_ids, accs):
    return {g: {a: 1 for a in accs} for g in genome_ids}


def test_pair_always_fused_keeps_only_one():
    genomes = [f"g{i}" for i in range(10)]
    hits = _hits(genomes, [A, B, C])
    shared = {g: {shared_pair_key(A, B): 1} for g in genomes}

    kept, exclusions = _resolve([A, B, C], shared, genomes, hits)

    assert kept == [A, C]           # tie on single-copy count -> lower acc kept
    assert len(exclusions) == 1
    ex = exclusions[0]
    assert (ex.dropped_acc, ex.kept_acc) == (B, A)
    assert ex.num_genomes_shared == 10
    assert ex.percent_genomes_shared == 100


def test_the_more_widely_single_copy_pfam_is_kept():
    genomes = [f"g{i}" for i in range(10)]
    hits = _hits(genomes, [A, B])
    hits["g0"][A] = 2               # A is single-copy in 9, B in 10
    shared = {g: {shared_pair_key(A, B): 1} for g in genomes}

    kept, exclusions = _resolve([A, B], shared, genomes, hits)
    assert kept == [B]
    assert exclusions[0].dropped_acc == A


def test_rare_sharing_below_the_cutoff_is_left_alone():
    genomes = [f"g{i}" for i in range(20)]
    hits = _hits(genomes, [A, B])
    shared = {"g0": {shared_pair_key(A, B): 1}}     # 5% of genomes

    kept, exclusions = _resolve(
        [A, B], shared, genomes, hits, max_shared_percent=10)
    assert kept == [A, B]
    assert exclusions == []

    kept, exclusions = _resolve(
        [A, B], shared, genomes, hits, max_shared_percent=5)
    assert len(kept) == 1 and len(exclusions) == 1


def test_three_domains_on_one_protein_leave_one():
    genomes = [f"g{i}" for i in range(5)]
    hits = _hits(genomes, [A, B, C, D])
    shared = {g: {shared_pair_key(A, B): 1, shared_pair_key(A, C): 1,
                  shared_pair_key(B, C): 1} for g in genomes}

    kept, exclusions = _resolve([A, B, C, D], shared, genomes, hits)
    assert kept == [A, D]
    assert {ex.dropped_acc for ex in exclusions} == {B, C}


def test_pairs_involving_a_non_retained_pfam_are_ignored():
    genomes = [f"g{i}" for i in range(5)]
    hits = _hits(genomes, [A, B])
    shared = {g: {shared_pair_key(A, D): 1} for g in genomes}   # D wasn't retained

    kept, exclusions = _resolve([A, B], shared, genomes, hits)
    assert kept == [A, B]
    assert exclusions == []


def test_genomes_outside_the_kept_set_do_not_count():
    genomes = [f"g{i}" for i in range(10)]
    hits = _hits(genomes + ["dropped"], [A, B])
    shared = {"dropped": {shared_pair_key(A, B): 1}}

    kept, exclusions = _resolve(
        [A, B], shared, genomes, hits, max_shared_percent=0)
    assert kept == [A, B]
    assert exclusions == []


def test_count_genomes_sharing_counts_genomes_not_proteins():
    genomes = ["g1", "g2", "g3"]
    shared = {
        "g1": {shared_pair_key(A, B): 3},       # three fused proteins, one genome
        "g2": {shared_pair_key(A, B): 1, shared_pair_key(C, D): 1},
        "g3": {shared_pair_key(C, D): 0},       # recorded but empty
    }
    assert count_genomes_sharing(shared, genomes) == {
        shared_pair_key(A, B): 2, shared_pair_key(C, D): 1}


def test_count_genomes_sharing_ignores_genomes_outside_the_set():
    shared = {"g1": {shared_pair_key(A, B): 1}, "dropped": {shared_pair_key(A, B): 1}}
    assert count_genomes_sharing(shared, ["g1", "g2"]) == {shared_pair_key(A, B): 1}


def test_resolution_is_deterministic_and_preserves_order():
    genomes = [f"g{i}" for i in range(4)]
    hits = _hits(genomes, [D, C, B, A])
    shared = {g: {shared_pair_key(C, D): 1} for g in genomes}
    kept, _ = _resolve([D, C, B, A], shared, genomes, hits)
    assert kept == [C, B, A]


def test_end_to_end_fused_pair_yields_one_target(tmp_path):
    """
    The situation that motivated all this: two domains fused into one protein in every
    genome, each hit exactly once, so both pass the single-copy cutoff.
    """
    genomes = {f"g{i}": [[A, B], [C], [D]] for i in range(6)}
    faa = _write(tmp_path / "t.faa", genomes)
    shared = {}
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, shared_proteins=shared)

    ids = list(genomes)
    wanted, _ = count_single_copy_hits(hits, ids, [A, B, C, D], 90)
    assert wanted == [A, B, C, D]

    kept, exclusions = _resolve(wanted, shared, ids, hits)
    assert kept == [A, C, D]
    assert [(e.dropped_acc, e.kept_acc) for e in exclusions] == [(B, A)]


################################################################################
# output table
################################################################################

def test_exclusions_table_is_written_with_names(tmp_path):
    genomes = [f"g{i}" for i in range(4)]
    hits = _hits(genomes, [A, B])
    shared = {g: {shared_pair_key(A, B): 1} for g in genomes[:3]}
    _, exclusions = _resolve([A, B], shared, genomes, hits)

    info = {A: PfamProfileInfo(A, "MockA", "desc", 90.0),
            B: PfamProfileInfo(B, "MockB", "desc", 90.0)}
    path = outputs.write_shared_protein_exclusions(str(tmp_path), exclusions, info)

    rows = [l.split("\t") for l in open(path).read().splitlines()]
    assert rows[0][0] == "dropped_pfam_id"
    assert rows[1] == [B, "MockB", A, "MockA", "3", "75"]


def test_exclusions_table_is_written_even_when_empty(tmp_path):
    path = outputs.write_shared_protein_exclusions(str(tmp_path), [], {})
    assert open(path).read().splitlines() == [
        "dropped_pfam_id\tdropped_name\tkept_pfam_id\tkept_name\t"
        "num_genomes_with_both_in_one_protein\tperc_genomes_with_both_in_one_protein"]
