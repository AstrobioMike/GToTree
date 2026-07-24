import pytest # type: ignore
from pathlib import Path

import pyhmmer # type: ignore

import gtotree.utils.hmms.gen_scg_hmms_search as mod
from gtotree.utils.hmms.gen_scg_hmms import count_single_copy_hits
from gtotree.utils.hmms.gen_scg_hmms_search import (
    HmmSearchError,
    load_target_proteins,
    search_profiles,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
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
# loading
################################################################################

def test_load_target_proteins(tmp_path):
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS)})
    alphabet, seqs = load_target_proteins(faa)
    assert len(seqs) == 4


def test_load_target_proteins_raises_on_empty(tmp_path):
    """
    An empty file fails at pyhmmer's format detection rather than reaching the
    empty-block check, but either way it must surface as HmmSearchError so the CLI can
    translate it instead of leaking a raw pyhmmer error.
    """
    empty = tmp_path / "empty.faa"
    empty.write_text("")
    with pytest.raises(HmmSearchError, match="failed to read"):
        load_target_proteins(str(empty))


def test_load_target_proteins_raises_on_whitespace_only(tmp_path):
    blank = tmp_path / "blank.faa"
    blank.write_text("\n\n")
    with pytest.raises(HmmSearchError):
        load_target_proteins(str(blank))


def test_load_target_proteins_raises_on_unreadable(tmp_path):
    with pytest.raises(HmmSearchError, match="failed to read"):
        load_target_proteins(str(tmp_path / "does-not-exist.faa"))


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
    assert hits["g1"] == {acc: 1 for acc in MOTIFS}


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


def test_search_progress_callback_fires_per_profile(tmp_path):
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS)})
    calls = []
    search_profiles(str(MOCK_PFAM_HMM), faa, threads=1,
                    progress_callback=lambda: calls.append(1))
    assert len(calls) == 4


def test_search_respects_block_size(tmp_path):
    """Blocking is an internal detail; results must not depend on it."""
    faa = _write_proteome(tmp_path / "t.faa", {"g1": list(MOTIFS)})
    a = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, block_size=1)
    b = search_profiles(str(MOCK_PFAM_HMM), faa, threads=1, block_size=100)
    assert a == b


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
