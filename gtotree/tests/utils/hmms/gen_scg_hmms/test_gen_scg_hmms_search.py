import json
import pytest # type: ignore
import gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_search as mod
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import count_single_copy_hits
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_search import (
    HmmSearchError,
    search_profiles,
    CHUNK_ENV_VAR,
    DEFAULT_MAX_RESIDUES_PER_CHUNK,
    ChunkSizeError,
    CheckpointError,
    resolve_residue_budget,
    estimated_chunk_bytes,
    _estimate_total_residues,
    _read_env_override,
)
from gtotree.tests.paths import DATA_DIR


MOCK_PFAM_HMM = DATA_DIR / "mock-pfams.hmm"

# the four motifs the mock HMMs were built from; each matches exactly one profile
MOTIFS = {
    "PF90001.3": "MKVLAAAL",
    "PF90002.7": "MARTKQTA",
    "PF90003.1": "MSDKIIHL",
    "PF90004.2": "MAHHWWGS",
}


def _write_proteome(path, genomes):
    """
    genomes: {genome_id: [pfam_acc, ...]} -- one protein written per listed acc,
    so repeating an acc gives that genome two copies of that domain.
    """
    with open(path, "w") as f:
        for genome_id, accs in genomes.items():
            for i, acc in enumerate(accs, 1):
                f.write(f">{genome_id}_{i}\n{MOTIFS[acc]}\n")
    return str(path)


################################################################################
# opening / streaming
################################################################################

def test_open_target_proteins_yields_a_streaming_file(tmp_path):
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS)})
    alphabet, seq_file = mod.open_target_proteins(faa)
    with seq_file:
        block = seq_file.read_block(residues=10 ** 9)
        assert len(block) == 4


def test_search_raises_on_empty_file(tmp_path):
    """
    An empty file fails at pyhmmer's format detection; either way it must surface as
    HmmSearchError so the CLI can translate it rather than leaking a raw pyhmmer error.
    """
    empty = tmp_path / "empty.faa"
    empty.write_text("")
    with pytest.raises(HmmSearchError):
        search_profiles(str(MOCK_PFAM_HMM), str(empty), threads=1)


def test_search_raises_on_whitespace_only(tmp_path):
    blank = tmp_path / "blank.faa"
    blank.write_text("\n\n")
    with pytest.raises(HmmSearchError):
        search_profiles(str(MOCK_PFAM_HMM), str(blank), threads=1)


def test_search_raises_on_missing_file(tmp_path):
    with pytest.raises(HmmSearchError, match="failed to read"):
        search_profiles(str(MOCK_PFAM_HMM), str(tmp_path / "nope.faa"), threads=1)


def test_chunks_are_read_lazily_not_all_at_once(tmp_path, monkeypatch):
    """
    The memory bound depends on reading a chunk per iteration rather than reading the
    whole file up front. Assert that directly: with a budget forcing several chunks,
    read_block must be called more than once (a single read_block() of everything --
    the thing that would silently undo the bound -- would call it once).
    """
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS), "g2": list(MOTIFS), "g3": list(MOTIFS),
    })

    calls = {"n": 0}
    real_read_chunk = mod._read_chunk

    def counting_read_chunk(seq_file, budget):
        calls["n"] += 1
        return real_read_chunk(seq_file, budget)

    monkeypatch.setattr(mod, "_read_chunk", counting_read_chunk)
    # motifs are 8 residues; a 16-residue budget gives 2 seqs/chunk over 12 seqs
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=16)

    # 6 chunks + 1 final empty read that terminates the loop
    assert calls["n"] > 2


################################################################################
# searching
################################################################################

def test_search_counts_hits_per_genome(tmp_path):
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS),
        "g2": list(MOTIFS),
    })
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)

    assert set(hits) == {"g1", "g2"}
    assert hits["g1"] == dict.fromkeys(MOTIFS, 1)


def test_search_counts_duplicates(tmp_path):
    """A genome with two copies of a domain must register 2, not 1."""
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS) + ["PF90001.3"],
    })
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)
    assert hits["g1"]["PF90001.3"] == 2


def test_search_keys_are_versioned_accession_strings(tmp_path):
    """
    Keys must be the versioned accession as `str`. On pyhmmer 0.11.0 the underlying
    values come back as bytes, so a missing _decode would produce b'...' keys that
    silently fail to join against the Pfam info table.
    """
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS)})
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)

    for acc in hits["g1"]:
        assert isinstance(acc, str)
        assert acc in MOTIFS


def test_search_recovers_genome_ids_with_underscores(tmp_path):
    """
    Genome ids contain underscores (GCF_000091665.1) and headers append `_<n>`, so
    only the final segment may be stripped.
    """
    faa = _write_proteome(tmp_path / "t.faa", {"GCF_000091665.1": list(MOTIFS)})
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)
    assert set(hits) == {"GCF_000091665.1"}


def test_search_progress_reports_genomes_completed(tmp_path):
    """
    The callback reports genomes *finished*, so the reported total must equal the
    genome count -- not the sequence count, and not the chunk count.
    """
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS),   # 4 seqs
        "g2": list(MOTIFS),   # 4 seqs
        "g3": list(MOTIFS),   # 4 seqs
    })  # 12 sequences, 3 genomes, each seq 8 residues

    counts = []
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=16,
                    progress_callback=lambda n: counts.append(n))
    assert sum(counts) == 3
    assert all(c >= 1 for c in counts)


def test_search_progress_counts_each_genome_exactly_once(tmp_path):
    """
    Chunks don't respect genome boundaries, so a genome can straddle a boundary. It
    must be reported once and only once -- never double counted by the chunk that ends
    mid-genome and the chunk that finishes it, and never dropped.

    Genomes have differing protein counts so boundaries land mid-genome at some budgets.
    """
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": ["PF90001.3"] * 5,
        "g2": ["PF90002.7"] * 3,
        "g3": ["PF90003.1"] * 7,
        "g4": ["PF90004.2"] * 1,
        "g5": ["PF90001.3"] * 4,
    })

    for budget in (8, 16, 24, 40, 10 ** 9):
        counts = []
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=budget,
                        progress_callback=lambda n: counts.append(n))
        assert sum(counts) == 5, f"budget {budget} reported {sum(counts)} genomes"


def test_search_progress_total_is_chunk_independent(tmp_path):
    """The reported genome total must not change with the residue budget."""
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS), "g2": list(MOTIFS), "g3": list(MOTIFS),
    })

    def total_reported(budget):
        counts = []
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=budget,
                        progress_callback=lambda n: counts.append(n))
        return sum(counts)

    assert total_reported(1) == 3           # every sequence its own chunk
    assert total_reported(10 ** 9) == 3     # one chunk


def test_search_progress_reaches_total_even_in_one_chunk(tmp_path):
    """
    The single-chunk case is the one that looked broken in practice (bar stuck at 0,
    then jumping to 100%). Even with one chunk every genome must still be reported, so
    the bar finishes filled rather than short.
    """
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS), "g2": list(MOTIFS),
    })
    counts = []
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=10 ** 9,
                    progress_callback=lambda n: counts.append(n))
    assert sum(counts) == 2


def test_search_progress_is_reported_during_a_chunk_not_just_after(tmp_path):
    """
    The bar has to move *while* a chunk is being searched. Reporting only at chunk
    boundaries is what made a multi-hour phase look hung, since a large run is only a
    handful of chunks.

    Asserted by relative ordering: with one chunk covering many genomes, the callback
    must fire before the search of that chunk has finished, i.e. interleaved with the
    per-profile passes rather than all afterwards.
    """
    faa = _write_proteome(tmp_path / "t.faa",
                          {f"g{i}": list(MOTIFS) for i in range(12)})

    events = []
    real_search_one_chunk = mod._search_one_chunk

    def tracking_search(pressed_base, seq_block, threads, hits_by_genome,
                        on_profile=None):
        def wrapped():
            events.append("profile")
            if on_profile is not None:
                on_profile()
        real_search_one_chunk(pressed_base, seq_block, threads, hits_by_genome,
                              on_profile=wrapped)
        events.append("chunk_searched")

    monkeypatch_target = mod
    original = monkeypatch_target._search_one_chunk
    monkeypatch_target._search_one_chunk = tracking_search
    try:
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=10 ** 9,
                        progress_callback=lambda n: events.append("progress"))
    finally:
        monkeypatch_target._search_one_chunk = original

    assert "progress" in events
    assert events.index("progress") < events.index("chunk_searched")


def test_search_progress_updates_more_than_once_per_chunk(tmp_path):
    """
    A single chunk covering many genomes should produce many updates, not one. This is
    the actual user-visible complaint -- a bar that ticks once per chunk.
    """
    faa = _write_proteome(tmp_path / "t.faa",
                          {f"g{i}": list(MOTIFS) for i in range(12)})
    counts = []
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=10 ** 9,
                    progress_callback=lambda n: counts.append(n))

    assert sum(counts) == 12
    assert len(counts) > 1


def test_press_profiles_reports_the_profile_count(tmp_path):
    """The count is the denominator for spreading progress across a chunk."""
    base, n = mod._press_profiles(str(MOCK_PFAM_HMM), str(tmp_path))
    assert n == len(MOTIFS)
    assert base.startswith(str(tmp_path))


################################################################################
# _ChunkProgress
################################################################################

def test_chunk_progress_awards_whole_genomes_only():
    """
    The callback contract is whole genomes; fractions are accumulated internally so a
    caller driving an integer-unit bar never sees a partial genome.
    """
    seen = []
    p = mod._ChunkProgress(seen.append, genomes_pending=3, num_profiles=10)
    for _ in range(10):
        p.on_profile()
    p.finish()
    assert all(isinstance(n, int) and n >= 1 for n in seen)
    assert sum(seen) == 3


def test_chunk_progress_spreads_credit_across_profiles():
    """Credit arrives gradually rather than all at the end."""
    seen = []
    p = mod._ChunkProgress(seen.append, genomes_pending=100, num_profiles=100)
    for _ in range(50):
        p.on_profile()
    assert sum(seen) == 50, "half the profiles done should mean about half the credit"


def test_chunk_progress_finish_settles_rounding():
    """
    Rounding down at each step would leave the bar permanently short of its total, so
    finish() has to hand over the remainder.
    """
    seen = []
    p = mod._ChunkProgress(seen.append, genomes_pending=7, num_profiles=3)
    for _ in range(3):
        p.on_profile()
    p.finish()
    assert sum(seen) == 7
    assert p.awarded == 7


def test_chunk_progress_never_overshoots():
    """Extra profile callbacks must not award more than the chunk actually earned."""
    seen = []
    p = mod._ChunkProgress(seen.append, genomes_pending=2, num_profiles=4)
    for _ in range(40):
        p.on_profile()
    p.finish()
    assert sum(seen) == 2


def test_chunk_progress_is_inert_with_no_callback():
    p = mod._ChunkProgress(None, genomes_pending=5, num_profiles=5)
    p.on_profile()
    p.finish()
    assert p.awarded == 0


def test_chunk_progress_handles_a_chunk_that_completes_nothing():
    """
    A single genome spanning several chunks earns no credit until it ends; those chunks
    must report nothing rather than dividing by zero or awarding phantom genomes.
    """
    seen = []
    p = mod._ChunkProgress(seen.append, genomes_pending=0, num_profiles=10)
    for _ in range(10):
        p.on_profile()
    p.finish()
    assert seen == []


def test_search_wraps_lazy_iteration_errors(tmp_path, monkeypatch):
    """
    pyhmmer.hmmsearch returns a LAZY generator -- the pipeline is only built (and its
    options validated) when results are consumed. A try/except around just the call
    would let real failures escape as raw pyhmmer errors past the library/CLI seam.
    """
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS)})

    def fake_hmmsearch(*args, **kwargs):
        def gen():
            yield from ()
            raise RuntimeError("simulated mid-iteration failure")
        return gen()

    monkeypatch.setattr(mod.pyhmmer, "hmmsearch", fake_hmmsearch)

    with pytest.raises(HmmSearchError, match="hmmsearch step failed"):
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)


################################################################################
# chunking: results must be identical regardless of residue budget
################################################################################

# genomes differ from each other, so a mis-merged accumulator would actually corrupt
# a per-genome Counter rather than being masked by identical genomes. g2 and g4 carry
# duplicate domains, which is what the single-copy math is sensitive to.
_MIXED = {
    "g1": list(MOTIFS),
    "g2": ["PF90001.3", "PF90002.7", "PF90001.3"],   # two copies of PF90001
    "g3": ["PF90003.1", "PF90004.2"],
    "g4": list(MOTIFS) + ["PF90002.7"],              # two copies of PF90002
    "g5": ["PF90001.3"],
}


def test_chunked_equals_unchunked(tmp_path):
    """
    Core invariant: for any residue budget, accumulated hits are identical to a
    single-pass search. Budget 1 (every sequence its own chunk) and a huge budget
    (one chunk) must agree exactly.
    """
    faa = _write_proteome(tmp_path / "t.faa", _MIXED)

    one_pass = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1,
                               residue_budget=10 ** 9)      # one chunk
    per_seq = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1,
                              residue_budget=1)              # 1 seq/chunk
    assert one_pass == per_seq


def test_all_budgets_agree(tmp_path):
    """Every intermediate budget must give the same result, not just the extremes."""
    faa = _write_proteome(tmp_path / "t.faa", _MIXED)
    # motifs are 8 residues; budgets from below one seq to well above the whole set
    results = [
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=b)
        for b in (1, 8, 9, 16, 24, 100, 10 ** 9)
    ]
    first = results[0]
    for r in results[1:]:
        assert r == first


def test_genome_split_across_chunks_still_counts_correctly(tmp_path):
    """
    The reason genome boundaries can be dropped: a genome whose proteins straddle a
    chunk boundary must still tally correctly, because hits merge by genome id in a
    shared accumulator. g2 has three proteins (two of them PF90001); a budget that
    forces a cut *between* them must still yield PF90001 == 2 for g2.
    """
    faa = _write_proteome(tmp_path / "t.faa", _MIXED)
    # 8 residues per motif; budget 16 fits two proteins per chunk, so any genome with
    # 3+ proteins (g2, g1, g4) is guaranteed to split across a boundary
    split = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=16)
    whole = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=10 ** 9)

    assert split == whole
    assert split["g2"]["PF90001.3"] == 2      # the duplicate survived the split
    assert split["g4"]["PF90002.7"] == 2


def test_env_override_drives_chunking(tmp_path):
    """The GTT_SCG_CHUNK env override (residues) must produce the same result too."""
    import os as _os
    faa = _write_proteome(tmp_path / "t.faa", _MIXED)
    baseline = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=10 ** 9)

    _os.environ["GTT_SCG_CHUNK"] = "8"     # one motif per chunk
    try:
        via_env = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)
    finally:
        del _os.environ["GTT_SCG_CHUNK"]
    assert via_env == baseline


def test_genome_with_no_hits_still_appears(tmp_path):
    """
    A genome whose proteins hit nothing must still be a key in the result (empty dict),
    because the single-copy denominator counts it. The up-front seeding of every genome
    id guarantees this even when hmmsearch reports nothing for it.
    """
    faa = tmp_path / "t.faa"
    with open(faa, "w") as f:
        f.write(">g1_1\n" + MOTIFS["PF90001.3"] + "\n")
        f.write(">g2_1\nWWWWWWWW\n")   # not a motif; hits nothing
    hits = search_profiles(str(MOCK_PFAM_HMM), str(faa), threads=1)
    assert "g2" in hits
    assert hits["g2"] == {}


def test_oversized_sequence_runs_as_its_own_chunk(tmp_path):
    """
    A single protein longer than the budget can't be deferred; it must still be
    searched (as a lone chunk) and its hit counted, not dropped or hung on.
    """
    faa = _write_proteome(tmp_path / "t.faa", {"g1": ["PF90001.3"]})
    # budget far below the 8-residue motif: the one sequence exceeds it alone
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=2)
    assert hits["g1"]["PF90001.3"] == 1


################################################################################
# search -> single-copy, end to end
################################################################################

def test_search_feeds_single_copy_determination(tmp_path):
    """
    g1/g2 carry one copy of everything; g3 has PF90001 twice and no PF90004.
    At 90% both of g3's oddities drop their Pfams; at 66% they survive.
    """
    faa = _write_proteome(tmp_path / "t.faa", {
        "g1": list(MOTIFS),
        "g2": list(MOTIFS),
        "g3": ["PF90001.3", "PF90002.7", "PF90003.1", "PF90001.3"],
    })
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1)
    genomes = ["g1", "g2", "g3"]
    accs = sorted(MOTIFS)

    strict, _ = count_single_copy_hits(hits, genomes, accs, 90)
    assert strict == ["PF90002.7", "PF90003.1"]

    loose, _ = count_single_copy_hits(hits, genomes, accs, 66)
    assert loose == sorted(MOTIFS)


################################################################################
# env override parsing
################################################################################

def test_env_override_unset_is_none():
    assert _read_env_override(env={}) is None
    assert _read_env_override(env={CHUNK_ENV_VAR: ""}) is None


def test_env_override_parses_int():
    assert _read_env_override(env={CHUNK_ENV_VAR: "1000000"}) == 1_000_000


def test_env_override_rejects_garbage():
    """
    A set-but-unparseable value must fail loudly: silently falling back to the default
    would hide the user's mistake behind a run that looks like it worked.
    """
    with pytest.raises(ChunkSizeError):
        _read_env_override(env={CHUNK_ENV_VAR: "fiftymillion"})


def test_env_override_rejects_nonpositive():
    with pytest.raises(ChunkSizeError):
        _read_env_override(env={CHUNK_ENV_VAR: "0"})
    with pytest.raises(ChunkSizeError):
        _read_env_override(env={CHUNK_ENV_VAR: "-3"})


################################################################################
# budget resolution
################################################################################

def test_resolve_prefers_override():
    assert resolve_residue_budget(env={CHUNK_ENV_VAR: "500"}) == 500


def test_budget_does_not_depend_on_run_size():
    """
    The budget is a flat cap, deliberately independent of how big the run is. Chunk
    boundaries have to be reproducible between a run and its `--resume`, so the budget
    must not vary with anything that could differ between the two.
    """
    assert resolve_residue_budget(env={}) == DEFAULT_MAX_RESIDUES_PER_CHUNK


def test_budget_leaves_room_for_the_pressed_profile_set():
    """
    The chunk is sized so the ~1.4 GB pressed profile set -- re-read every chunk -- can
    stay in page cache alongside it. A block big enough to evict it turns every chunk's
    profile read into real disk I/O.
    """
    assert estimated_chunk_bytes(DEFAULT_MAX_RESIDUES_PER_CHUNK) < 300_000_000


def test_estimate_total_residues_is_roughly_right(tmp_path):
    """
    The estimate only picks chunk granularity, so it needs to be in the right ballpark
    rather than exact -- but a wildly wrong figure would defeat the purpose.
    """
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS), "g2": list(MOTIFS)})
    actual = 8 * 8   # 8 sequences of 8 residues
    est = _estimate_total_residues(faa)
    assert est is not None
    assert 0.2 * actual < est < 5 * actual


def test_estimate_returns_none_for_missing_file(tmp_path):
    assert _estimate_total_residues(str(tmp_path / "nope.faa")) is None


def test_resolve_falls_back_to_default():
    assert resolve_residue_budget(env={}) == DEFAULT_MAX_RESIDUES_PER_CHUNK


def test_resolve_reads_real_environment(monkeypatch):
    """With no env argument it must consult the actual environment."""
    monkeypatch.setenv(CHUNK_ENV_VAR, "1234")
    assert resolve_residue_budget() == 1234
    monkeypatch.delenv(CHUNK_ENV_VAR)
    assert resolve_residue_budget() == DEFAULT_MAX_RESIDUES_PER_CHUNK


def test_default_is_positive():
    assert DEFAULT_MAX_RESIDUES_PER_CHUNK >= 1


################################################################################
# the measured-footprint helper
################################################################################

def test_estimated_bytes_scales_with_budget():
    small = estimated_chunk_bytes(1_000_000)
    big = estimated_chunk_bytes(10_000_000)
    assert big == pytest.approx(small * 10, rel=0.01)


def test_estimated_bytes_is_in_a_sane_range():
    """
    Guards the measured constant against a typo (e.g. losing a factor of 1000). At the
    measured ~2.9 bytes/residue, 10M residues should be tens of MB, not KB or GB.
    """
    est = estimated_chunk_bytes(10_000_000)
    assert 10 * 1024 ** 2 < est < 200 * 1024 ** 2


################################################################################
# chunk-level checkpointing
################################################################################

def _five_genome_proteome(tmp_path, name="t.faa"):
    """A proteome big enough to span several chunks at a small budget."""
    return _write_proteome(tmp_path / name, {
        "g1": list(MOTIFS),
        "g2": list(MOTIFS),
        "g3": list(MOTIFS),
        "g4": list(MOTIFS),
        "g5": list(MOTIFS),
    })


def test_checkpoint_is_written_per_chunk(tmp_path):
    """
    One record per finished chunk, not one at the end -- an all-at-the-end write would
    be worthless for the crash it's meant to survive.
    """
    faa = _five_genome_proteome(tmp_path)
    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))

    lines = [l for l in ckpt.read_text().splitlines() if l.strip()]
    header = json.loads(lines[0])
    assert header["budget"] == 32
    assert header["num_profiles"] == len(MOTIFS)
    records = [json.loads(l) for l in lines[1:]]
    assert len(records) > 1
    assert [r["chunk"] for r in records] == list(range(len(records)))


def test_resume_reproduces_the_uninterrupted_result(tmp_path):
    """
    The invariant that matters: resuming from a partial checkpoint must give exactly
    what an uninterrupted run gives. Simulated by truncating a complete checkpoint to
    its first chunk, which is indistinguishable from having died after that chunk.
    """
    faa = _five_genome_proteome(tmp_path)
    expected = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32)

    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))
    lines = [l for l in ckpt.read_text().splitlines() if l.strip()]
    assert len(lines) > 2, "need more than one chunk for this to be a real test"
    ckpt.write_text("\n".join(lines[:2]) + "\n")

    resumed = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                              checkpoint_path=str(ckpt), resume=True)
    assert resumed == expected


def test_resume_from_a_complete_checkpoint_searches_nothing(tmp_path):
    """
    Every chunk already recorded means no chunk needs searching again -- that's the
    whole point. Result must still be complete.
    """
    faa = _five_genome_proteome(tmp_path)
    ckpt = tmp_path / "ckpt.jsonl"
    expected = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                               checkpoint_path=str(ckpt))

    searched = {"n": 0}
    real = mod._search_one_chunk

    def counting(*a, **kw):
        searched["n"] += 1
        return real(*a, **kw)

    mod._search_one_chunk = counting
    try:
        resumed = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                                  checkpoint_path=str(ckpt), resume=True)
    finally:
        mod._search_one_chunk = real

    assert searched["n"] == 0
    assert resumed == expected


def test_resume_reports_every_genome_on_the_progress_bar(tmp_path):
    """
    The bar has to refill over the replayed chunks, otherwise a resumed run finishes
    showing a fraction of the genomes it actually processed.
    """
    faa = _five_genome_proteome(tmp_path)
    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))

    counts = []
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt), resume=True,
                    progress_callback=lambda n: counts.append(n))
    assert sum(counts) == 5


def test_torn_trailing_record_costs_only_that_chunk(tmp_path):
    """
    A process killed mid-write leaves a half-written last line. Replay must stop there
    and redo that chunk, not refuse the whole checkpoint and not ingest garbage.
    """
    faa = _five_genome_proteome(tmp_path)
    expected = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32)

    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))
    lines = [l for l in ckpt.read_text().splitlines() if l.strip()]
    truncated = lines[1][: len(lines[1]) // 2]
    ckpt.write_text("\n".join(lines[:2]) + "\n" + truncated + "\n")

    resumed = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                              checkpoint_path=str(ckpt), resume=True)
    assert resumed == expected


def test_checkpoint_from_a_different_budget_is_discarded(tmp_path):
    """
    Chunk boundaries are a function of the budget, so replaying a checkpoint written
    under a different one would double-count or drop hits. It must be ignored and the
    work redone rather than trusted.
    """
    faa = _five_genome_proteome(tmp_path)
    expected = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32)

    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))

    resumed = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=24,
                              checkpoint_path=str(ckpt), resume=True)
    assert resumed == expected


def test_checkpoint_from_a_different_fasta_is_discarded(tmp_path):
    """Same reasoning: different targets means different boundaries."""
    faa = _five_genome_proteome(tmp_path)
    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))

    other = _write_proteome(tmp_path / "other.faa", {
        "h1": list(MOTIFS), "h2": list(MOTIFS),
    })
    resumed = search_profiles(str(MOCK_PFAM_HMM), other, threads=1, residue_budget=32,
                              checkpoint_path=str(ckpt), resume=True)
    assert set(resumed) == {"h1", "h2"}, "stale genomes must not leak in from the old run"


def test_missing_checkpoint_just_starts_over(tmp_path):
    """Resuming with nothing to resume from is a full run, not an error."""
    faa = _five_genome_proteome(tmp_path)
    expected = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32)
    resumed = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                              checkpoint_path=str(tmp_path / "nope.jsonl"), resume=True)
    assert resumed == expected


def test_resume_does_not_truncate_the_checkpoint_it_read(tmp_path):
    """
    A resumed run appends to the existing checkpoint. Starting a fresh one would throw
    away the replayed chunks, so a second crash would restart from zero.
    """
    faa = _five_genome_proteome(tmp_path)
    ckpt = tmp_path / "ckpt.jsonl"
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt))
    lines = [l for l in ckpt.read_text().splitlines() if l.strip()]
    ckpt.write_text("\n".join(lines[:2]) + "\n")

    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                    checkpoint_path=str(ckpt), resume=True)

    after = [json.loads(l) for l in ckpt.read_text().splitlines() if l.strip()][1:]
    assert [r["chunk"] for r in after] == list(range(len(lines) - 1))


def test_no_checkpoint_path_means_no_checkpointing(tmp_path):
    """The stage must work unchanged when no checkpoint location is supplied."""
    faa = _five_genome_proteome(tmp_path)
    hits = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                           checkpoint_path=None)
    assert set(hits) == {"g1", "g2", "g3", "g4", "g5"}
    assert not list(tmp_path.glob("*.jsonl"))


def test_unwritable_checkpoint_surfaces_as_a_named_error(tmp_path):
    """
    Library code raises a named exception for the CLI to translate, rather than letting
    a raw OSError escape the seam.
    """
    faa = _five_genome_proteome(tmp_path)
    unwritable = tmp_path / "missing-dir" / "ckpt.jsonl"
    with pytest.raises(CheckpointError):
        search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, residue_budget=32,
                        checkpoint_path=str(unwritable))


def test_merge_hits_sums_across_chunks():
    """
    A genome straddling a chunk boundary contributes to two deltas; replay has to add
    them, not let the later one win.
    """
    target = {"g1": {"PF90001.3": 1}}
    mod._merge_hits(target, {"g1": {"PF90001.3": 2, "PF90002.7": 1}, "g2": {}})
    assert target == {"g1": {"PF90001.3": 3, "PF90002.7": 1}, "g2": {}}
