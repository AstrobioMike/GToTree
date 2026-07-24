"""
Core library for generating single-copy-gene (SCG) HMM sets

  1. resolve a set of target genomes (accessions file and/or `--wanted-ref-tax`)
  2. get amino acids for each (download the NCBI protein file, or call genes with
     prodigal off the nucleotide file when no protein file exists)
  3. take the Pfam profiles whose underlying proteins are well-covered by the model
     (average coverage > 50%), since partial-domain models make poor SCG markers
  4. hmmsearch all target proteins against that filtered Pfam set
  5. keep the Pfams hit exactly once in >= `percent_single_copy`% of the genomes
  6. write those profiles out as a new SCG-HMM set
"""

import os
from collections import Counter

import pyhmmer  # type: ignore

from gtotree.utils.pfam.get_pfam_data import HMM_FILENAME, INFO_FILENAME


# Pfam profiles whose average coverage of the underlying proteins is at or below this
# are dropped before searching: this is to try to avoid multi-domain proteins from being used
DEFAULT_MIN_PFAM_COVERAGE = 50.0

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


def _decode(value):
    """
    pyhmmer has returned `name`/`accession` as bytes in some versions and str in
    others, so normalize rather than assuming either.
    """
    if value is None:
        return None
    if isinstance(value, (bytes, bytearray)):
        return value.decode()
    return str(value)


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
    Parse pfamA.txt and return {versioned_acc: PfamProfileInfo} for profiles whose
    average coverage of their underlying proteins exceeds `min_coverage`.

    Keyed by versioned accession (e.g. "PF00001.27") because that is what the master
    Pfam-A.hmm carries in its ACC lines, making it the reliable join key.
    """
    kept = {}
    usable_rows = 0

    with open(info_path, "r", encoding="utf-8", errors="ignore") as f:
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

            if coverage <= min_coverage:
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

    return kept


def write_filtered_pfam_hmms(master_hmm_path, wanted_accs, out_path,
                             progress_callback=None):
    """
    Stream the master Pfam-A.hmm and write out only the profiles in `wanted_accs`.

    Returns the list of versioned accessions actually found and written.
    """
    wanted = set(wanted_accs)
    found = []
    tmp_path = out_path + ".part"

    try:
        with pyhmmer.plan7.HMMFile(master_hmm_path) as hmm_file, \
                open(tmp_path, "wb") as out:
            for hmm in hmm_file:
                acc = _decode(hmm.accession)
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

    single_copy_counts = Counter()
    for genome_id in genome_ids:
        counts = hits_by_genome.get(genome_id, {})
        for acc, n in counts.items():
            if n == 1:
                single_copy_counts[acc] += 1

    wanted = [acc for acc in filtered_accs if single_copy_counts.get(acc, 0) >= threshold]

    per_genome_counts = {
        genome_id: dict(hits_by_genome.get(genome_id, {})) for genome_id in genome_ids
    }

    return wanted, per_genome_counts


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
            acc = _decode(hmm.accession)
            if acc is None:
                acc = _decode(hmm.name)
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
