"""
The hmmsearch stage of `gtt gen-scg-hmms`

Searches all target proteins against the coverage-filtered Pfam profiles and tallies,
per genome, how many times each profile was hit.

Gathering thresholds (`--cut_ga`) are used just like in main GToTree

MEMORY SHAPE
------------
The target proteins are the dominant resident cost of this phase and grow linearly
with the run size, so they are never all held at once. Instead the combined fasta is
streamed in residue-budgeted chunks (`SequenceFile.read_block(residues=...)`): each
chunk is searched against the whole profile set, folded into the hit accumulator, then
released before the next chunk is read. Resident sequence memory is therefore bounded
by the budget rather than by genome count.

The profile set is `hmmpress`ed once up front so the repeated per-chunk passes read
optimized profiles instead of re-parsing the HMM text each time; the press streams
the profiles.
"""

import os
import tempfile
import pyhmmer  # type: ignore
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import GenSCGHMMsError, _decode
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import genome_id_from_protein_name


class HmmSearchError(GenSCGHMMsError):
    """The hmmsearch stage failed."""


def open_target_proteins(fasta_path):
    """
    Open the combined protein fasta for streaming reads, returning (alphabet, file).

    Failures here are translated so a bad/missing/empty fasta surfaces as
    HmmSearchError at the library/CLI seam rather than as a raw pyhmmer error.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    try:
        seq_file = pyhmmer.easel.SequenceFile(
            fasta_path, digital=True, alphabet=alphabet)
    except Exception as e:
        raise HmmSearchError(f"failed to read the combined protein fasta: {e}")

    return alphabet, seq_file


def _read_chunk(seq_file, residue_budget):
    """
    Read the next residue-budgeted chunk, or None at end of file.

    pyhmmer fills the block with whole sequences until the residue limit is reached, so
    a chunk never splits a protein, and a protein longer than the budget comes back as
    a one-sequence chunk (rather than looping forever on an empty read).
    """
    try:
        block = seq_file.read_block(residues=residue_budget)
    except Exception as e:
        raise HmmSearchError(f"failed to read the combined protein fasta: {e}")

    if block is None or len(block) == 0:
        return None
    return block


def _press_profiles(filtered_hmm_path, press_dir):
    """
    hmmpress the filtered profile set into `press_dir`, returning the pressed base path.

    Pressing once means every chunk pass reads optimized profiles rather than
    re-parsing the HMM text. The HMMFile is handed to hmmpress directly rather than
    being materialized with list(): the filtered set runs to many thousands of
    profiles, and holding them all as live objects would be a large spike in exactly
    the code meant to reduce peak memory.
    """
    base = os.path.join(press_dir, "filtered-pfams-pressed")
    try:
        with pyhmmer.plan7.HMMFile(filtered_hmm_path) as hmm_file:
            pressed = pyhmmer.hmmer.hmmpress(hmm_file, base)
    except Exception as e:
        raise HmmSearchError(f"failed to prepare (hmmpress) the filtered profiles: {e}")

    if not pressed:
        raise HmmSearchError(
            "the filtered Pfam HMM held no profiles to search with.")

    return base


def _search_one_chunk(pressed_base, seq_block, threads, hits_by_genome):
    """
    Search one sequence chunk against the whole pressed profile set, folding hits into
    the shared `hits_by_genome` accumulator.

    Opening the pressed file fresh per chunk keeps each pass reading optimized profiles
    off disk without holding the whole profile set resident between chunks.
    """
    try:
        with pyhmmer.plan7.HMMPressedFile(pressed_base) as profiles:

            results = pyhmmer.hmmsearch(
                profiles, seq_block, cpus=threads, bit_cutoffs="gathering",
            )

            for top_hits in results:

                query = getattr(top_hits, "query", None)
                acc = _decode(getattr(query, "accession", None))
                if acc is None:
                    acc = _decode(getattr(query, "name", None))

                for hit in top_hits:
                    if not hit.included:
                        continue
                    genome_id = genome_id_from_protein_name(_decode(hit.name))
                    counts = hits_by_genome.setdefault(genome_id, {})
                    counts[acc] = counts.get(acc, 0) + 1

    except KeyboardInterrupt:
        raise
    except HmmSearchError:
        raise
    except Exception as e:
        raise HmmSearchError(f"the hmmsearch step failed: {e}")


def search_profiles(filtered_hmm_path, fasta_path, threads=1,
                    progress_callback=None, residue_budget=None):
    """
    Search the filtered Pfam profiles against the combined target proteins.

    Returns `hits_by_genome`: {genome_id: {versioned_pfam_acc: hit_count}}, including
    an empty dict for any genome that was searched but hit nothing.

    Proteins are streamed in residue-budgeted chunks

    `progress_callback`, if given, is called with the number of genomes *finished* since
    the last call, so a caller can drive a bar totalled by genome count. A genome counts
    as finished once a protein from a later genome has been seen (or the file ends):
    chunks don't respect genome boundaries, so the genome a chunk ends on is usually
    still mid-read and isn't reported until the chunk that completes it.

    `residue_budget` overrides the chunk policy for deterministic testing; left None it
    uses `resolve_residue_budget`, which also takes the size of this fasta into account
    so that a small run still gets several chunks (and so a steadily moving bar) rather
    than a single one that jumps from 0 to 100% at the end.
    """
    if residue_budget is not None:
        budget = residue_budget
    else:
        budget = resolve_residue_budget(total_residues=_estimate_total_residues(fasta_path))

    hits_by_genome = {}
    total_seqs = 0

    # genome-completion tracking; proteins of a genome are contiguous in the combined
    # fasta (see relabel_and_append), so a change of genome id means the previous one
    # is finished
    genomes_started = 0
    genomes_reported = 0
    current_genome = None

    alphabet, seq_file = open_target_proteins(fasta_path)

    with seq_file, tempfile.TemporaryDirectory(prefix="gtt-scg-press-") as press_dir:
        pressed_base = _press_profiles(filtered_hmm_path, press_dir)

        while True:
            chunk = _read_chunk(seq_file, budget)
            if chunk is None:
                break

            # register every genome in this chunk before searching, so genomes that hit
            # nothing still end up as keys; the same pass tracks genome transitions
            for seq in chunk:
                genome_id = genome_id_from_protein_name(_decode(seq.name))
                hits_by_genome.setdefault(genome_id, {})
                if genome_id != current_genome:
                    genomes_started += 1
                    current_genome = genome_id

            _search_one_chunk(pressed_base, chunk, threads, hits_by_genome)

            total_seqs += len(chunk)

            # all but the genome this chunk ended on are definitely complete; that last
            # one may continue into the next chunk, so it waits
            if progress_callback is not None:
                finished = max(0, genomes_started - 1)
                if finished > genomes_reported:
                    progress_callback(finished - genomes_reported)
                    genomes_reported = finished

            # drop the chunk before reading the next one -- this is what bounds the
            # resident sequence memory
            del chunk

    # end of file: whatever genome we ended on is complete too
    if progress_callback is not None and genomes_started > genomes_reported:
        progress_callback(genomes_started - genomes_reported)

    if total_seqs == 0:
        raise HmmSearchError(
            "no protein sequences were available to search; none of the target genomes "
            "yielded usable amino-acid sequences.")

    return hits_by_genome


"""
Chunk-size policy details

The search stage streams the target proteins past the filtered Pfam set. Read all at
once, those proteins are the dominant resident cost of the phase and grow linearly
with the run size, so instead they're read in *residue-budgeted chunks*: each chunk is
searched and discarded before the next is read, which bounds resident sequence memory
at roughly the budget no matter how many genomes the run covers.

This module answers one question: *what is the residue budget?* The actual chunking is
done by pyhmmer's `SequenceFile.read_block(residues=...)`, which fills a block up to a
residue limit while streaming from disk (see `gen_scg_hmms_search`).

CALIBRATION
-----------
Resident cost of a chunk was measured directly (240k proteins / 77M residues, digital
amino blocks, pyhmmer 0.12.1): peak RSS tracks the residue budget linearly at

    ~2.8 - 3.0 bytes resident per residue

(3.5 B/res at a 1M-residue budget, converging to 2.80 at 64M as fixed per-block
overhead amortizes). That figure covers per-residue storage plus per-sequence object
overhead at realistic protein lengths (~150-500 aa); very short sequences shift it up.
So, to size the budget for a target sequence-block footprint:

    budget_residues  ~=  target_bytes / 2.9
"""

# Env var to force a fixed residues-per-chunk, bypassing the default. Two audiences:
# (1) calibration: pin the budget to read a peak off each of several runs;
# (2) an escape hatch for a user hitting an edge the default gets wrong.
# I made this an env var rather than a CLI flag to avoid cluttering `--help`
CHUNK_ENV_VAR = "GTT_SCG_CHUNK"

# Measured resident cost per residue for a digital amino-acid sequence block; see the
# CALIBRATION note. Used only to document/derive the default, not at runtime.
MEASURED_BYTES_PER_RESIDUE = 2.9

# Default residues per chunk. At the measured ~2.9 B/residue this is a ~1 GB
# sequence block. Raising it means fewer, larger chunk passes; the pressed profile set
# makes each pass cheap, so the cost of a smaller budget is modest.

# This only ever binds on very large runs, the granularity rule below tightens the
# budget so a smaller run is still split into ~TARGET_CHUNKS pieces (so progress bar is
# informative, as there is minimal impact on speed or memory use).
DEFAULT_MAX_RESIDUES_PER_CHUNK = 350_000_000


class ChunkSizeError(ValueError):
    """A chunk-size override was given but wasn't a usable positive integer."""


def _read_env_override(env=None):
    """
    Return the fixed residues-per-chunk from the environment, or None if unset.

    Raises ChunkSizeError on a set-but-garbage value rather than silently ignoring it:
    a user who exported `GTT_SCG_CHUNK=fiftymillion` meant to pin the budget, and
    quietly falling back to the default would hide the mistake behind a run that looks
    like it worked.
    """
    env = os.environ if env is None else env
    raw = env.get(CHUNK_ENV_VAR)
    if raw is None or raw == "":
        return None

    try:
        value = int(raw)
    except (TypeError, ValueError):
        raise ChunkSizeError(
            f"{CHUNK_ENV_VAR} must be a positive integer (residues) if set, got {raw!r}.")

    if value < 1:
        raise ChunkSizeError(
            f"{CHUNK_ENV_VAR} must be >= 1 if set, got {value}.")

    return value


# Aim for at least roughly this many chunks so the progress bar actually moves during
# the phase instead of jumping from 0 to 100% at the end
TARGET_CHUNKS = 10

# Bytes of fasta per residue, used only to estimate a file's residue count from its size
# for the granularity calculation above. Protein fasta carries headers and newlines on
# top of the residues themselves, so this is deliberately approximate. It never affects
# the memory cap, which is applied separately and exactly.
_APPROX_FASTA_BYTES_PER_RESIDUE = 1.1


def _estimate_total_residues(fasta_path):
    """
    Cheap O(1) estimate of the residues in a fasta, from its size on disk.

    Used only to pick a chunk granularity; being off by 10-20% just shifts the chunk
    count slightly. Returns None if the size can't be read, in which case the caller
    falls back to the plain memory-capped budget.
    """
    try:
        size = os.path.getsize(fasta_path)
    except OSError:
        return None
    return int(size / _APPROX_FASTA_BYTES_PER_RESIDUE)


def resolve_residue_budget(env=None, total_residues=None, target_chunks=TARGET_CHUNKS):
    """
    The residues-per-chunk budget for this run.

    Precedence:
      1. `GTT_SCG_CHUNK`, if set, wins outright.
      2. otherwise the memory cap (`DEFAULT_MAX_RESIDUES_PER_CHUNK`), tightened if
         needed so a run smaller than the cap is still split into `target_chunks`
         chunks rather than one.

    The memory cap is only ever *lowered* by the granularity rule, never raised, so a
    large run stays bounded by memory and simply gets more than `target_chunks` chunks.
    """
    override = _read_env_override(env)
    if override is not None:
        return override

    budget = DEFAULT_MAX_RESIDUES_PER_CHUNK

    if total_residues and target_chunks and target_chunks > 1:
        # ceil division, so target_chunks is a floor on the chunk count rather than a
        # value we can undershoot by rounding
        per_chunk = -(-total_residues // target_chunks)
        budget = min(budget, per_chunk)

    return max(1, budget)


def estimated_chunk_bytes(residue_budget):
    """
    Rough resident bytes for a chunk of `residue_budget` residues, from the measured
    per-residue figure. For reporting and for sanity-checking a chosen budget; the
    hmmsearch pipeline's own working memory is *not* included (see CALIBRATION).
    """
    return int(residue_budget * MEASURED_BYTES_PER_RESIDUE)
