"""
In-process HMM searching with pyhmmer for the main GToTree path

GATHERING THRESHOLDS
--------------------
Searches run with `bit_cutoffs="gathering"`, the equivalent of what i used to
do with `--cut_ga` at the command-line invocation
"""

import os

import pyhmmer  # type: ignore
import pyhmmer.easel as easel  # type: ignore
import pyhmmer.plan7 as plan7  # type: ignore


class MissingGatheringCutoffs(Exception):
    """
    A profile set has no gathering thresholds, so --cut_ga-equivalent can't be used
    """


def profiles_missing_gathering_cutoffs(hmm_path):
    """
    Names of profiles in `hmm_path` that carry no GA line.

    Empty list means the set is usable with gathering thresholds.
    """
    missing = []
    with plan7.HMMFile(hmm_path) as hmm_file:
        for hmm in hmm_file:
            if not hmm.cutoffs.gathering_available():
                name = hmm.name
                missing.append(name.decode() if isinstance(name, bytes) else str(name))
    return missing


def press_profiles(hmm_path, press_dir, basename):
    """
    `hmmpress` a profile set into `press_dir` once, returning the pressed base path.

    Returns None if pressing fails, so callers fall back to reading the plain HMM file
    rather than failing the run over an optimization.

    The HMMFile is handed to hmmpress directly rather than being materialized with
    list(): for large sets (the Pfam path especially) holding every profile as a live
    object at once is a needless memory spike.
    """
    base = os.path.join(press_dir, basename)
    try:
        with plan7.HMMFile(hmm_path) as hmm_file:
            pyhmmer.hmmer.hmmpress(hmm_file, base)
    except Exception:
        return None
    return base


def _open_profiles(pressed_base, hmm_path):
    """
    Open the pressed profile set if one was prepared, else the plain HMM file.

    The fallback keeps every entry point working when no pressing step ran -- the
    standalone search helpers, and any caller that didn't go through the fused stage.
    """
    if pressed_base:
        try:
            return plan7.HMMPressedFile(pressed_base)
        except Exception:
            pass
    return plan7.HMMFile(hmm_path)


def search_one_genome(hmm_path, sequence_path, tblout_path,
                      pressed_base=None, cpus=1):
    """
    Search one genome's proteins against a profile set, writing a --tblout table
    """
    alphabet = easel.Alphabet.amino()

    with easel.SequenceFile(sequence_path, digital=True, alphabet=alphabet) as seq_file:
        sequences = seq_file.read_block()

    results = []

    with _open_profiles(pressed_base, hmm_path) as profiles, \
            open(tblout_path, "wb") as out:
        first = True
        for hits in pyhmmer.hmmsearch(profiles, sequences, cpus=cpus,
                                      bit_cutoffs="gathering"):
            hits.write(out, format="targets", header=first)
            first = False
            results.append(hits)

    return results
