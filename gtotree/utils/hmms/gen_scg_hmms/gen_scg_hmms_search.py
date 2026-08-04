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

    `progress_callback`, if given, is called once per chunk with that chunk's sequence
    count, so a caller can drive a bar totalled by sequence count.

    `residue_budget` overrides the chunk policy for deterministic testing; left None it
    uses `resolve_residue_budget`
    """
    budget = residue_budget if residue_budget is not None else resolve_residue_budget()

    hits_by_genome = {}
    total_seqs = 0

    alphabet, seq_file = open_target_proteins(fasta_path)

    with seq_file, tempfile.TemporaryDirectory(prefix="gtt-scg-press-") as press_dir:
        pressed_base = _press_profiles(filtered_hmm_path, press_dir)

        while True:
            chunk = _read_chunk(seq_file, budget)
            if chunk is None:
                break

            # register every genome in this chunk before searching, so genomes that hit
            # nothing still end up as keys
            for seq in chunk:
                hits_by_genome.setdefault(
                    genome_id_from_protein_name(_decode(seq.name)), {})

            _search_one_chunk(pressed_base, chunk, threads, hits_by_genome)

            total_seqs += len(chunk)
            if progress_callback is not None:
                progress_callback(len(chunk))

            # drop the chunk before reading the next one -- this is what bounds the
            # resident sequence memory
            del chunk

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

What is NOT yet measured is the *hmmsearch pipeline's* own working memory per chunk
(DP matrices, per-thread scratch, accumulated TopHits), which sits on top of the block
and scales with profile count and thread count rather than with genome count. The
default below therefore leaves generous headroom for that term. Reading the real
combined figure off a `GTT_DEBUG_TIMING=1` run against the true ~10-20k-profile
filtered Pfam set is what would let the default be tightened.
"""

# Env var to force a fixed residues-per-chunk, bypassing the default. Two audiences:
# (1) calibration: pin the budget to read a peak off each of several runs;
# (2) an escape hatch for a user hitting an edge the default gets wrong. An env var
# rather than a CLI flag so it avoids cluttering `--help`
CHUNK_ENV_VAR = "GTT_SCG_CHUNK"

# Measured resident cost per residue for a digital amino-acid sequence block; see the
# CALIBRATION note. Used only to document/derive the default, not at runtime.
MEASURED_BYTES_PER_RESIDUE = 2.9

# Default residues per chunk. At the measured ~2.9 B/residue this is a ~140 MB
# sequence block, chosen to leave room for the (still unmeasured) hmmsearch pipeline
# term on top. Raising it means fewer, larger chunk passes; the pressed profile set
# makes each pass cheap, so the cost of a smaller budget is modest.
DEFAULT_MAX_RESIDUES_PER_CHUNK = 50_000_000


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


def resolve_residue_budget(env=None):
    """
    The residues-per-chunk budget for this run: the env override if set, else the
    default. Always >= 1.
    """
    override = _read_env_override(env)
    return override if override is not None else DEFAULT_MAX_RESIDUES_PER_CHUNK


def estimated_chunk_bytes(residue_budget):
    """
    Rough resident bytes for a chunk of `residue_budget` residues, from the measured
    per-residue figure. For reporting and for sanity-checking a chosen budget; the
    hmmsearch pipeline's own working memory is *not* included (see CALIBRATION).
    """
    return int(residue_budget * MEASURED_BYTES_PER_RESIDUE)
