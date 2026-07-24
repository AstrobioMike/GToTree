import os
import pytest # type: ignore
from collections import Counter
from pathlib import Path

import pyhmmer # type: ignore

from gtotree.utils.hmms.gen_scg_hmms import (
    DEFAULT_MIN_PFAM_COVERAGE,
    GenSCGHMMsError,
    PfamDataError,
    _decode,
    count_single_copy_hits,
    load_coverage_filtered_pfams,
    pfam_data_paths,
    read_hmm_accessions,
    write_filtered_pfam_hmms,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
MOCK_PFAM_INFO = DATA_DIR / "mock-pfamA.txt"
MOCK_PFAM_HMM = DATA_DIR / "mock-pfams.hmm"

# what the fixture holds, for reference in the tests below:
#   PF90001.3  coverage 75.5   -> kept
#   PF90002.7  coverage 60.2   -> kept
#   PF90003.1  coverage 50.0   -> DROPPED (filter is strictly >, not >=)
#   PF90004.2  coverage 12.4   -> dropped
#   PF90005.9  coverage 104.7  -> kept (real Pfam coverages can exceed 100; see below)


################################################################################
# _decode
################################################################################

def test_decode_handles_bytes_and_str():
    """
    pyhmmer returns name/accession as bytes on 0.11.x and str on later versions, so
    every call site goes through _decode. Both must normalize to str.
    """
    assert _decode(b"PF00001.27") == "PF00001.27"
    assert _decode("PF00001.27") == "PF00001.27"
    assert _decode(None) is None


################################################################################
# coverage filtering
################################################################################

def test_load_coverage_filtered_pfams_applies_cutoff():
    kept = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO))
    assert set(kept) == {"PF90001.3", "PF90002.7", "PF90005.9"}


def test_coverage_cutoff_is_strictly_greater_than():
    """
    PF90003 sits exactly at the default cutoff of 50.0 and must be excluded -- v1 used
    `> 50`, not `>= 50`, and that boundary decides whether a marker is considered.
    """
    kept = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO))
    assert "PF90003.1" not in kept

    # lowering the cutoff below 50 lets it through, confirming it's the cutoff and not
    # some other property of that row
    kept_lower = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO), min_coverage=49.9)
    assert "PF90003.1" in kept_lower


def test_coverage_above_100_is_kept():
    """
    Coverage over 100 is real in Pfam (envelope coordinates can overhang; ~15 of 30k
    rows at Pfam 38.2, mostly viral/secreted families that are essentially all domain).
    Those rows must not be treated as invalid and silently dropped.
    """
    kept = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO))
    assert "PF90005.9" in kept
    assert kept["PF90005.9"].coverage == pytest.approx(104.7)


def test_accession_key_is_versioned():
    """
    The key must be `{acc}.{version}` because that's what the master Pfam-A.hmm carries
    in its ACC lines -- it's the join key the whole pipeline depends on.
    """
    kept = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO))
    assert "PF90001.3" in kept
    assert "PF90001" not in kept
    assert kept["PF90001.3"].name == "MockA"
    assert kept["PF90001.3"].description == "Mock domain A"


def test_load_coverage_filtered_pfams_custom_cutoff():
    kept = load_coverage_filtered_pfams(str(MOCK_PFAM_INFO), min_coverage=70.0)
    assert set(kept) == {"PF90001.3", "PF90005.9"}


def test_load_coverage_filtered_pfams_skips_short_rows(tmp_path):
    """Malformed/truncated rows are skipped rather than raising."""
    path = tmp_path / "partial.txt"
    good = MOCK_PFAM_INFO.read_text().splitlines()[0]
    path.write_text("too\tfew\tcolumns\n" + good + "\n")

    kept = load_coverage_filtered_pfams(str(path))
    assert set(kept) == {"PF90001.3"}


def test_load_coverage_filtered_pfams_raises_when_unparsable(tmp_path):
    path = tmp_path / "junk.txt"
    path.write_text("nothing\tusable\there\n")

    with pytest.raises(PfamDataError, match="no usable rows"):
        load_coverage_filtered_pfams(str(path))


def test_load_coverage_filtered_pfams_raises_when_none_pass(tmp_path):
    path = tmp_path / "info.txt"
    path.write_text(MOCK_PFAM_INFO.read_text())

    with pytest.raises(PfamDataError, match="no Pfam profiles passed"):
        load_coverage_filtered_pfams(str(path), min_coverage=200.0)


def test_non_numeric_coverage_row_is_skipped(tmp_path):
    row = MOCK_PFAM_INFO.read_text().splitlines()[0].split("\t")
    row[33] = "NA"
    path = tmp_path / "info.txt"
    path.write_text("\t".join(row) + "\n" + MOCK_PFAM_INFO.read_text().splitlines()[1] + "\n")

    kept = load_coverage_filtered_pfams(str(path))
    assert set(kept) == {"PF90002.7"}


################################################################################
# pfam_data_paths
################################################################################

def test_pfam_data_paths_raises_on_missing(tmp_path):
    with pytest.raises(PfamDataError, match="was not found"):
        pfam_data_paths(str(tmp_path))


def test_pfam_data_paths_raises_on_empty_file(tmp_path):
    from gtotree.utils.pfam.get_pfam_data import HMM_FILENAME, INFO_FILENAME
    (tmp_path / HMM_FILENAME).write_text("")
    (tmp_path / INFO_FILENAME).write_text("x")

    with pytest.raises(PfamDataError):
        pfam_data_paths(str(tmp_path))


################################################################################
# profile extraction
################################################################################

def test_write_filtered_pfam_hmms_selects_wanted(tmp_path):
    out = tmp_path / "subset.hmm"
    found = write_filtered_pfam_hmms(str(MOCK_PFAM_HMM),
                                     {"PF90001.3", "PF90003.1"}, str(out))

    assert sorted(found) == ["PF90001.3", "PF90003.1"]
    assert sorted(read_hmm_accessions(str(out))) == ["PF90001.3", "PF90003.1"]


def test_write_filtered_pfam_hmms_round_trips_exactly(tmp_path):
    """
    Order must be preserved: `filtered_accs` becomes the column order of the
    hit-count matrix, so a reordering here would silently scramble that table.
    """
    out = tmp_path / "subset.hmm"
    wanted = set(read_hmm_accessions(str(MOCK_PFAM_HMM)))
    found = write_filtered_pfam_hmms(str(MOCK_PFAM_HMM), wanted, str(out))

    assert read_hmm_accessions(str(out)) == found


def test_write_filtered_pfam_hmms_progress_callback(tmp_path):
    calls = []
    write_filtered_pfam_hmms(str(MOCK_PFAM_HMM), {"PF90001.3"},
                             str(tmp_path / "s.hmm"),
                             progress_callback=lambda: calls.append(1))
    # called once per profile scanned, not per profile matched
    assert len(calls) == 4


def test_write_filtered_pfam_hmms_raises_when_nothing_matches(tmp_path):
    with pytest.raises(PfamDataError, match="none of the coverage-filtered"):
        write_filtered_pfam_hmms(str(MOCK_PFAM_HMM), {"PF00000.1"},
                                 str(tmp_path / "s.hmm"))


def test_write_filtered_pfam_hmms_is_atomic_on_failure(tmp_path):
    """
    A failed extraction must leave neither a `.part` nor a truncated file at the
    destination -- a later run (or --resume) would otherwise trust the partial file.
    """
    out = tmp_path / "subset.hmm"
    boom = RuntimeError("interrupted mid-write")

    def explode():
        raise boom

    with pytest.raises(RuntimeError):
        write_filtered_pfam_hmms(str(MOCK_PFAM_HMM), {"PF90001.3"}, str(out),
                                 progress_callback=explode)

    assert not out.exists()
    assert not Path(str(out) + ".part").exists()


def test_read_hmm_accessions_raises_on_empty(tmp_path):
    empty = tmp_path / "empty.hmm"
    empty.write_bytes(b"")
    with pytest.raises((PfamDataError, Exception)):
        read_hmm_accessions(str(empty))


################################################################################
# single-copy determination
################################################################################

def _hits(pattern, n_genomes=10):
    """Build hits_by_genome from {acc: [count per genome]}."""
    genomes = [f"G{i}" for i in range(n_genomes)]
    out = {}
    for i, g in enumerate(genomes):
        counts = {acc: counts_list[i] for acc, counts_list in pattern.items()
                  if counts_list[i] > 0}
        out[g] = Counter(counts)
    return genomes, out


def test_single_copy_requires_exactly_one_copy():
    """
    The criterion is 'hit exactly once', not 'present'. A Pfam in EVERY genome but at
    two copies is a bad phylogenetic marker (you can't tell which copy is the ortholog)
    and must be dropped -- this is the whole point of a single-copy set.
    """
    genomes, hits = _hits({
        "single_in_all": [1] * 10,
        "duplicated_in_all": [2] * 10,
    })
    wanted, _ = count_single_copy_hits(hits, genomes,
                                       ["single_in_all", "duplicated_in_all"], 90)
    assert wanted == ["single_in_all"]


def test_duplication_and_absence_penalized_identically():
    """Both simply fail to be 'exactly one', so both cost the same."""
    genomes, hits = _hits({
        "nine_then_dup":     [1] * 9 + [2],
        "nine_then_absent":  [1] * 9 + [0],
    })
    wanted, _ = count_single_copy_hits(
        hits, genomes, ["nine_then_dup", "nine_then_absent"], 90)
    assert wanted == ["nine_then_dup", "nine_then_absent"]


def test_threshold_boundary_is_inclusive():
    """90% must pass at exactly 9/10; 95% must not."""
    genomes, hits = _hits({"nine_of_ten": [1] * 9 + [2]})

    wanted, _ = count_single_copy_hits(hits, genomes, ["nine_of_ten"], 90)
    assert wanted == ["nine_of_ten"]

    wanted, _ = count_single_copy_hits(hits, genomes, ["nine_of_ten"], 95)
    assert wanted == []


def test_matches_v1_criterion_across_thresholds():
    """
    Cross-check against v1's criterion transcribed directly:
        sum(count == 1) / n_genomes * 100 >= percent_single_copy
    """
    pattern = {
        "all_single":    [1] * 10,
        "nine_single":   [1] * 9 + [2],
        "five_single":   [1] * 5 + [0] * 5,
        "never_single":  [2] * 10,
    }
    genomes, hits = _hits(pattern)
    accs = list(pattern)

    def v1(percent):
        keep = []
        for acc in accs:
            n_one = sum(1 for g in genomes if hits[g].get(acc, 0) == 1)
            if n_one / len(genomes) * 100 >= percent:
                keep.append(acc)
        return keep

    for percent in (50, 66, 90, 95, 100):
        mine, _ = count_single_copy_hits(hits, genomes, accs, percent)
        assert mine == v1(percent), f"mismatch at -p {percent}"


def test_genomes_are_never_dropped():
    """
    This stage excludes Pfams, never genomes. Every genome must still appear in the
    returned count matrix even if it contributed no single-copy hits at all.
    """
    genomes, hits = _hits({"x": [1] * 9 + [0]})
    _, per_genome = count_single_copy_hits(hits, genomes, ["x"], 90)

    assert sorted(per_genome) == sorted(genomes)
    assert len(per_genome) == 10


def test_count_matrix_includes_all_searched_profiles():
    genomes, hits = _hits({"a": [1] * 10, "b": [0] * 10})
    _, per_genome = count_single_copy_hits(hits, genomes, ["a", "b"], 90)
    # 'b' was searched but never hit -- it should simply be absent per genome,
    # and the output writer fills 0
    assert per_genome["G0"] == {"a": 1}


def test_count_single_copy_hits_raises_without_genomes():
    with pytest.raises(GenSCGHMMsError, match="no genomes"):
        count_single_copy_hits({}, [], ["a"], 90)


def test_wanted_order_follows_filtered_accs():
    """Output ordering must follow the searched-profile order, not hit counts."""
    genomes, hits = _hits({"z": [1] * 10, "a": [1] * 10, "m": [1] * 10})
    wanted, _ = count_single_copy_hits(hits, genomes, ["z", "a", "m"], 90)
    assert wanted == ["z", "a", "m"]
