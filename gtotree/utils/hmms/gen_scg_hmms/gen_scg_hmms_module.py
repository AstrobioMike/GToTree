"""
Core library for generating single-copy-gene (SCG) HMM sets

  1. Resolve a set of target genomes (accessions file and/or `--wanted-ref-tax`)
  2. Get amino acids for each (download the NCBI protein file, or call genes with
     prodigal off the nucleotide file when no protein file exists)
  3. Take the Pfam profiles whose underlying proteins are well-covered by the model
     (average coverage > 50%), since partial-domain models make poor SCG markers
  4. Hmmsearch all target proteins against that filtered Pfam set
  5. Keep the Pfams hit exactly once in >= `percent_single_copy` of the genomes
  6. Of any retained Pfams that commonly sit on the *same* protein, keep only one
  7. If more than `max_hmms` remain, keep the best-ranked ones (see `cap_targets`)
  8. Write those profiles out as a new SCG-HMM set

Steps 5-7 only need the search results, which is what lets `--from-run` build a new
set from a finished run's output tables without redoing 1-4 (see `select_scg_targets`).
"""

import os
from collections import Counter

import pyhmmer  # type: ignore

from gtotree.utils.misc.general import decode_pyhmmer_text
from gtotree.utils.pfam.get_pfam_data import HMM_FILENAME, INFO_FILENAME


# Pfam profiles whose average coverage of the underlying proteins is at or below this
# are dropped before searching: this is to try to avoid multi-domain proteins from being used
DEFAULT_MIN_PFAM_COVERAGE = 50.0

# two retained Pfams hitting the same protein in at least this percent of the genomes
# are treated as one gene, and only one of them is kept (see
# `resolve_shared_protein_conflicts`)
DEFAULT_MAX_SHARED_PROTEIN_PERCENT = 10.0

DEFAULT_MAX_HMMS = 250

# column positions (0-based) in pfamA.txt that we depend on
#   col 0  -> pfamA_acc            e.g. "PF00001"
#   col 1  -> pfamA_id (name)      e.g. "7tm_1"
#   col 3  -> description
#   col 27 -> version integer      e.g. 27, giving the full acc "PF00001.27"
#   col 33 -> average coverage     e.g. 66.58
_PFAM_ACC_COL = 0
_PFAM_NAME_COL = 1
_PFAM_DESC_COL = 3
_PFAM_VERSION_COL = 27
_PFAM_COVERAGE_COL = 33
_PFAM_MIN_COLS = 34


class GenSCGHMMsError(Exception):
    """Something went wrong that the CLI should report and exit on."""


class PfamDataError(GenSCGHMMsError):
    """The managed Pfam data is missing or not usable."""


class PfamProfileInfo:
    """Minimal metadata for one Pfam profile, for the output info table."""

    __slots__ = ("acc", "name", "description", "coverage")

    def __init__(self, acc, name, description, coverage):
        self.acc = acc
        self.name = name
        self.description = description
        self.coverage = coverage


def pfam_data_paths(pfam_data_dir):
    """
    Return (master_hmm_path, info_path) inside the managed Pfam data dir, raising if
    either is missing so callers fail before doing any expensive work.
    """
    hmm_path = os.path.join(pfam_data_dir, HMM_FILENAME)
    info_path = os.path.join(pfam_data_dir, INFO_FILENAME)

    for path in (hmm_path, info_path):
        if not (os.path.isfile(path) and os.path.getsize(path) > 0):
            raise PfamDataError(
                f"the required Pfam file '{os.path.basename(path)}' was not found in "
                f"'{pfam_data_dir}'.")

    return hmm_path, info_path


def load_coverage_filtered_pfams(info_path, min_coverage=DEFAULT_MIN_PFAM_COVERAGE):
    """
    Parse pfamA.txt and return (kept, total_profiles) where `kept` is
    {versioned_acc: PfamProfileInfo} for profiles whose average coverage of their
    underlying proteins exceeds `min_coverage`, and `total_profiles` is the count
    of usable rows (i.e. the number of profiles in the master set, before the
    coverage filter).

    `min_coverage=None` keeps every profile; `--from-run` uses that, since the
    previous run's hit-count table already records which profiles were searched.

    `total_profiles` is returned so callers can size a progress bar for the
    subsequent streaming pass over the master HMM without a second pass to count.

    Keyed by versioned accession (e.g. "PF00001.27") because that is what the master
    Pfam-A.hmm carries in its ACC lines, making it the reliable join key.
    """
    kept = {}
    usable_rows = 0

    with open(info_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < _PFAM_MIN_COLS:
                continue

            acc = parts[_PFAM_ACC_COL].strip()
            if not acc:
                continue

            try:
                coverage = float(parts[_PFAM_COVERAGE_COL])
            except (TypeError, ValueError):
                continue

            usable_rows += 1

            if min_coverage is not None and coverage <= min_coverage:
                continue

            version = parts[_PFAM_VERSION_COL].strip()
            versioned_acc = f"{acc}.{version}" if version else acc

            kept[versioned_acc] = PfamProfileInfo(
                acc=versioned_acc,
                name=parts[_PFAM_NAME_COL].strip(),
                description=parts[_PFAM_DESC_COL].strip(),
                coverage=coverage,
            )

    if usable_rows == 0:
        raise PfamDataError(
            f"no usable rows were parsed out of the Pfam info table ('{info_path}'). "
            "The Pfam data may be corrupt; re-downloading it with "
            "`gtt data get pfam --force-update` may resolve it.")

    if not kept:
        raise PfamDataError(
            f"no Pfam profiles passed the average-coverage filter (> {min_coverage}%).")

    return kept, usable_rows


def write_filtered_pfam_hmms(master_hmm_path, wanted_accs, out_path,
                             progress_callback=None):
    """
    Stream the master Pfam-A.hmm and write out only the profiles in `wanted_accs`.

    `progress_callback`, if given, is called once per profile *scanned* (not per
    profile matched), so a caller can drive a bar whose total is the number of
    profiles in the master file -- giving a steady fill as the file streams.

    Returns the list of versioned accessions actually found and written.
    """
    wanted = set(wanted_accs)
    found = []
    tmp_path = out_path + ".part"

    try:
        with pyhmmer.plan7.HMMFile(master_hmm_path) as hmm_file, \
                open(tmp_path, "wb") as out:
            for hmm in hmm_file:
                acc = decode_pyhmmer_text(hmm.accession)
                if acc in wanted:
                    hmm.write(out)
                    found.append(acc)
                if progress_callback is not None:
                    progress_callback()
        os.replace(tmp_path, out_path)
    except BaseException:
        _remove_quietly(tmp_path)
        raise

    if not found:
        raise PfamDataError(
            "none of the coverage-filtered Pfam profiles were found in the master "
            "Pfam HMM file; the Pfam data may be out of sync.")

    return found


def count_single_copy_hits(hits_by_genome, genome_ids, filtered_accs,
                           percent_single_copy):
    """
    Decide which Pfams qualify as single-copy markers.

    `hits_by_genome` maps genome_id -> Counter of versioned Pfam acc -> hit count.
    A Pfam qualifies when it is hit EXACTLY once in at least `percent_single_copy`
    percent of `genome_ids`.

    Returns (wanted_accs, per_genome_counts) where per_genome_counts is the full
    genome x pfam count matrix as {genome_id: {acc: count}}, for the output table.
    """
    if not genome_ids:
        raise GenSCGHMMsError("no genomes were available to determine single-copy genes from.")

    total = len(genome_ids)
    threshold = percent_single_copy / 100.0 * total

    single_copy_counts = count_single_copy_genomes(hits_by_genome, genome_ids)

    wanted = [acc for acc in filtered_accs if single_copy_counts.get(acc, 0) >= threshold]

    per_genome_counts = {
        genome_id: dict(hits_by_genome.get(genome_id, {})) for genome_id in genome_ids
    }

    return wanted, per_genome_counts


SHARED_PAIR_SEP = "|"


def shared_pair_key(acc_a, acc_b):
    """
    Order-independent key for a pair of profile accessions
    """
    a, b = sorted((acc_a, acc_b))
    return f"{a}{SHARED_PAIR_SEP}{b}"


def split_shared_pair_key(key):
    a, _sep, b = key.partition(SHARED_PAIR_SEP)
    return a, b


def count_single_copy_genomes(hits_by_genome, genome_ids):
    """
    Counter of acc -> number of `genome_ids` in which it was hit exactly once
    """
    single_copy_counts = Counter()
    for genome_id in genome_ids:
        counts = hits_by_genome.get(genome_id, {})
        for acc, n in counts.items():
            if n == 1:
                single_copy_counts[acc] += 1
    return single_copy_counts


def count_multi_copy_genomes(hits_by_genome, genome_ids):
    """
    Counter of acc -> number of `genome_ids` in which it was hit more than once
    """
    multi_copy_counts = Counter()
    for genome_id in genome_ids:
        counts = hits_by_genome.get(genome_id, {})
        for acc, n in counts.items():
            if n > 1:
                multi_copy_counts[acc] += 1
    return multi_copy_counts


def count_genomes_sharing(shared_by_genome, genome_ids):
    """
    Collapse the per-genome shared-protein record into {pair_key: number of
    `genome_ids` in which at least one protein was hit by both Pfams of the pair}.

    `shared_by_genome` is {genome_id: {pair_key: n proteins}} from the search. This
    aggregate is all the conflict resolution needs, and it's what gets written to
    `shared-protein-pairs.tsv`, so a finished run carries it for `--from-run`.
    """
    genomes_sharing = Counter()
    for genome_id in genome_ids:
        for key, n in (shared_by_genome.get(genome_id) or {}).items():
            if n > 0:
                genomes_sharing[key] += 1
    return genomes_sharing


class SharedProteinExclusion:
    """
    One Pfam dropped for sharing proteins with a Pfam that was kept
    """

    __slots__ = ("dropped_acc", "kept_acc", "num_genomes_shared",
                 "percent_genomes_shared")

    def __init__(self, dropped_acc, kept_acc, num_genomes_shared, percent_genomes_shared):
        self.dropped_acc = dropped_acc
        self.kept_acc = kept_acc
        self.num_genomes_shared = num_genomes_shared
        self.percent_genomes_shared = percent_genomes_shared


def resolve_shared_protein_conflicts(wanted_accs, genomes_sharing, genome_ids,
                                     hits_by_genome,
                                     max_shared_percent=DEFAULT_MAX_SHARED_PROTEIN_PERCENT):
    """
    Of the retained Pfams, keep only one of any pair that commonly hits the same protein.

    Two Pfams can each be hit exactly once per genome and still be the same gene, e.g.,
    two domains of a fused, multi-functional protein. Main GToTree extracts whole
    proteins, so keeping both would put the same sequence into two SCG sets
    (double-weighting that gene in the tree), and where the fusion is lineage-specific
    the fused proteins run about twice the usual length and get dropped by the
    length filter.

    `genomes_sharing` is {pair_key: number of genomes with a protein hit by both}, as
    built by `count_genomes_sharing` over `genome_ids`.

    A pair conflicts when both Pfams are in `wanted_accs` and at least one protein is
    hit by both in >= `max_shared_percent`% of `genome_ids`. Conflicts are settled
    most-shared first; each time, the Pfam that is single-copy in more genomes is kept
    (ties go to the lower accession, just so it's deterministic). A Pfam already dropped
    doesn't knock anything else out, so three domains on one protein leave one behind.

    Returns (kept_accs in their original order, [SharedProteinExclusion, ...]).
    """
    total = len(genome_ids)
    if not total or not wanted_accs:
        return list(wanted_accs), []

    wanted_set = set(wanted_accs)
    threshold = max_shared_percent / 100.0 * total

    conflicts = []
    for key, n in genomes_sharing.items():
        if n <= 0 or n < threshold:
            continue
        a, b = split_shared_pair_key(key)
        if a in wanted_set and b in wanted_set:
            conflicts.append((key, n))

    if not conflicts:
        return list(wanted_accs), []

    single_copy_counts = count_single_copy_genomes(hits_by_genome, genome_ids)

    # most-shared first, then by key so the order is fully deterministic
    conflicts.sort(key=lambda item: (-item[1], item[0]))

    dropped = {}
    exclusions = []
    for key, n in conflicts:
        a, b = split_shared_pair_key(key)
        if a in dropped or b in dropped:
            continue
        keep, drop = sorted((a, b), key=lambda acc: (-single_copy_counts.get(acc, 0), acc))
        dropped[drop] = keep
        exclusions.append(SharedProteinExclusion(
            dropped_acc=drop, kept_acc=keep, num_genomes_shared=n,
            percent_genomes_shared=round(n / total * 100, 2)))

    kept = [acc for acc in wanted_accs if acc not in dropped]
    return kept, exclusions


def cap_targets(wanted_accs, single_copy_counts, pfam_info, max_hmms):
    """
    Keep at most `max_hmms` of `wanted_accs`, the best-ranked ones.

    Ranked by the number of genomes the Pfam is single-copy in (most first), then by
    its average coverage of the underlying proteins (highest first, since a model
    spanning more of the protein is closer to a whole-gene marker), then by
    accession so the result is deterministic. A Pfam with no info entry ranks as
    zero coverage.

    `max_hmms` of 0 or None means no cap.

    Returns (kept in their original order, capped-out in rank order).
    """
    if not max_hmms or len(wanted_accs) <= max_hmms:
        return list(wanted_accs), []

    def _coverage(acc):
        info = pfam_info.get(acc)
        return info.coverage if info is not None else 0.0

    ranked = sorted(wanted_accs,
                    key=lambda acc: (-single_copy_counts.get(acc, 0), -_coverage(acc), acc))
    keep = set(ranked[:max_hmms])

    return [acc for acc in wanted_accs if acc in keep], ranked[max_hmms:]


class SCGSelection:
    """
    The outcome of `select_scg_targets`, with what was dropped at each step
    """

    __slots__ = ("targets", "num_single_copy", "shared_exclusions", "capped_out",
                 "single_copy_counts", "multi_copy_counts", "num_genomes")

    def __init__(self, targets, num_single_copy, shared_exclusions, capped_out,
                 single_copy_counts, multi_copy_counts, num_genomes):
        self.targets = targets
        self.num_single_copy = num_single_copy
        self.shared_exclusions = shared_exclusions
        self.capped_out = capped_out
        self.single_copy_counts = single_copy_counts
        self.multi_copy_counts = multi_copy_counts
        self.num_genomes = num_genomes


def select_scg_targets(hits_by_genome, genome_ids, filtered_accs, genomes_sharing,
                       pfam_info, percent_single_copy,
                       max_shared_percent=DEFAULT_MAX_SHARED_PROTEIN_PERCENT,
                       max_hmms=DEFAULT_MAX_HMMS):
    """
    Everything from the search results to the final target list: the single-copy
    cutoff, then shared-protein conflicts, then the `max_hmms` cap.

    The cap comes last so that Pfams dropped for sharing a protein don't leave the
    set short of the cap. Used by both a full run and `--from-run`, so the two can't
    drift apart in how they select.

    Returns an SCGSelection; `targets` may be empty, which the caller reports.
    """
    wanted, _per_genome = count_single_copy_hits(
        hits_by_genome, genome_ids, filtered_accs, percent_single_copy)
    num_single_copy = len(wanted)

    wanted, shared_exclusions = resolve_shared_protein_conflicts(
        wanted, genomes_sharing, genome_ids, hits_by_genome,
        max_shared_percent=max_shared_percent)

    single_copy_counts = count_single_copy_genomes(hits_by_genome, genome_ids)
    targets, capped_out = cap_targets(wanted, single_copy_counts, pfam_info, max_hmms)

    return SCGSelection(
        targets=targets,
        num_single_copy=num_single_copy,
        shared_exclusions=shared_exclusions,
        capped_out=capped_out,
        single_copy_counts=single_copy_counts,
        multi_copy_counts=count_multi_copy_genomes(hits_by_genome, genome_ids),
        num_genomes=len(genome_ids),
    )


def read_hmm_accessions(hmm_path):
    """
    Read just the versioned accessions out of an existing HMM file.

    Used when resuming, to recover the searched-profile list from an already-written
    filtered HMM without redoing the expensive extraction pass.
    """
    import pyhmmer  # type: ignore

    accs = []
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        for hmm in hmm_file:
            acc = decode_pyhmmer_text(hmm.accession)
            if acc is None:
                acc = decode_pyhmmer_text(hmm.name)
            accs.append(acc)

    if not accs:
        raise PfamDataError(
            f"no profiles could be read out of '{hmm_path}'.")

    return accs


def _remove_quietly(path):
    try:
        os.remove(path)
    except OSError:
        pass
