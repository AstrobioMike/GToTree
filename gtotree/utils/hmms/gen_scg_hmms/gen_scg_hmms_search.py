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

CRASH RECOVERY
--------------
Chunks are also the checkpoint unit. This phase runs for hours at realistic sizes, and
each finished chunk appends what it found to a JSONL checkpoint, so an interrupted run
resumes at the first chunk that never completed instead of starting over. See
`_SearchCheckpoint`.
"""

import os
import json
import tempfile
import pyhmmer  # type: ignore
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import GenSCGHMMsError
from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import genome_id_from_protein_name
from gtotree.utils.misc.general import decode_pyhmmer_text


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
        raise HmmSearchError(f"failed to read the combined protein fasta: {e}") from e

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
        raise HmmSearchError(f"failed to read the combined protein fasta: {e}") from e

    if block is None or len(block) == 0:
        return None
    return block


def _press_profiles(filtered_hmm_path, press_dir):
    """
    hmmpress the filtered profile set into `press_dir`, returning
    (pressed base path, number of profiles pressed).

    Pressing once means every chunk pass reads optimized profiles rather than
    re-parsing the HMM text. The HMMFile is handed to hmmpress directly rather than
    being materialized with list(): the filtered set runs to many thousands of
    profiles, and holding them all as live objects would be a large spike in exactly
    the code meant to reduce peak memory.

    The profile count comes back because it's the denominator for progress reporting:
    each chunk pass yields exactly one TopHits per profile, so the count says how many
    progress steps a chunk is worth (see `_ChunkProgress`).
    """
    base = os.path.join(press_dir, "filtered-pfams-pressed")
    try:
        with pyhmmer.plan7.HMMFile(filtered_hmm_path) as hmm_file:
            pressed = pyhmmer.hmmer.hmmpress(hmm_file, base)
    except Exception as e:
        raise HmmSearchError(f"failed to prepare (hmmpress) the filtered profiles: {e}") from e

    if not pressed:
        raise HmmSearchError(
            "the filtered Pfam HMM held no profiles to search with.")

    return base, pressed


class _ChunkProgress:
    """
    Spreads the genome-completion credit one chunk earns across the profile passes
    inside that chunk (so progress bar is more informative)

    Without this the callback only fires at chunk boundaries, so a run with ~10 chunks
    gets ~10 bar updates across a multi-hour phase and looks hung. `pyhmmer.hmmsearch`
    yields one TopHits per profile and (measured) yields them near-linearly in time, so
    profile completions are a good clock to interpolate against.
    """

    __slots__ = ("_callback", "_pending", "_step", "_credit", "awarded")

    def __init__(self, callback, genomes_pending, num_profiles):
        self._callback = callback
        self._pending = genomes_pending
        self._credit = 0.0
        self.awarded = 0
        if callback is None or genomes_pending <= 0 or not num_profiles:
            self._step = 0.0
        else:
            self._step = genomes_pending / num_profiles

    def on_profile(self):
        """
        One profile finished against this chunk.

        The award is clamped to what the chunk actually earned. In normal operation
        this can't bind, as hmmsearch yields exactly one TopHits per profile, but the
        clamp means a miscounted denominator can only make the bar lag, never make it
        overshoot its total or double-count genomes.
        """
        if not self._step:
            return
        self._credit += self._step
        whole = min(int(self._credit), self._pending) - self.awarded
        if whole > 0:
            self.awarded += whole
            self._callback(whole)

    def finish(self):
        """
        Hand over any genomes this chunk earned but rounding hasn't released yet, so
        the cumulative reported count stays exact rather than drifting low.
        """
        if self._callback is None:
            return
        remaining = self._pending - self.awarded
        if remaining > 0:
            self.awarded += remaining
            self._callback(remaining)


def _search_one_chunk(pressed_base, seq_block, threads, hits_by_genome,
                      on_profile=None):
    """
    Search one sequence chunk against the whole pressed profile set, folding hits into
    the shared `hits_by_genome` accumulator.

    Opening the pressed file fresh per chunk keeps each pass reading optimized profiles
    off disk without holding the whole profile set resident between chunks. (Measured:
    re-opening and streaming a 13k-profile pressed set costs well under a second, so
    the per-chunk cost is negligible against the search itself.)

    `on_profile`, if given, is called once per profile finished against this chunk
    """
    try:
        with pyhmmer.plan7.HMMPressedFile(pressed_base) as profiles:

            results = pyhmmer.hmmsearch(
                profiles, seq_block, cpus=threads, bit_cutoffs="gathering",
            )

            for top_hits in results:

                if on_profile is not None:
                    on_profile()

                query = getattr(top_hits, "query", None)
                acc = decode_pyhmmer_text(getattr(query, "accession", None))
                if acc is None:
                    acc = decode_pyhmmer_text(getattr(query, "name", None))

                for hit in top_hits:
                    if not hit.included:
                        continue
                    genome_id = genome_id_from_protein_name(decode_pyhmmer_text(hit.name))
                    counts = hits_by_genome.setdefault(genome_id, {})
                    counts[acc] = counts.get(acc, 0) + 1

    except KeyboardInterrupt:
        raise
    except HmmSearchError:
        raise
    except Exception as e:
        raise HmmSearchError(f"the hmmsearch step failed: {e}") from e


class CheckpointError(GenSCGHMMsError):
    """The search checkpoint exists but can't be used."""


class _SearchCheckpoint:
    """
    Chunk-level crash recovery for the search.
    """

    def __init__(self, path, budget, num_profiles, fasta_path):
        self.path = path
        self._header = self._build_header(budget, num_profiles, fasta_path)
        self._started = False

    @staticmethod
    def _build_header(budget, num_profiles, fasta_path):
        try:
            st = os.stat(fasta_path)
            size, mtime = st.st_size, int(st.st_mtime)
        except OSError:
            size, mtime = None, None
        return {"budget": budget, "num_profiles": num_profiles,
                "fasta_size": size, "fasta_mtime": mtime}

    def load(self):
        """
        Replay the checkpoint, returning (hits_by_genome, chunks_done).

        Returns ({}, 0) when there's nothing usable to resume from, no file, an
        unreadable one, or a header describing a different run. That's a fallback to
        redoing the work, never a wrong answer, so it isn't an error.
        """
        if not self.path or not os.path.isfile(self.path):
            return {}, 0

        hits = {}
        chunks_done = 0
        try:
            with open(self.path, encoding="utf-8") as f:
                first = f.readline()
                if not first.strip():
                    return {}, 0
                try:
                    header = json.loads(first)
                except ValueError:
                    return {}, 0
                if header != self._header:
                    return {}, 0

                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        record = json.loads(line)
                    except ValueError:
                        # torn trailing write: everything before it is still good
                        break
                    if record.get("chunk") != chunks_done:
                        # out of order means the file isn't what we think it is
                        break
                    _merge_hits(hits, record.get("hits") or {})
                    chunks_done += 1
        except OSError as e:
            raise CheckpointError(
                f"the search checkpoint at '{self.path}' couldn't be read: {e}") from e

        return hits, chunks_done

    def record(self, chunk_index, delta_hits):
        """
        Append one finished chunk. Flushed and fsynced so a checkpoint survives the
        process being killed, which is the entire point of having one.
        """
        if not self.path:
            return
        try:
            if not self._started:
                self._begin()
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(json.dumps({"chunk": chunk_index, "hits": delta_hits},
                                   separators=(",", ":")) + "\n")
                f.flush()
                os.fsync(f.fileno())
        except OSError as e:
            raise CheckpointError(
                f"the search checkpoint at '{self.path}' couldn't be written: {e}") from e

    def _begin(self):
        """Write the header, truncating anything already there that we didn't accept."""
        with open(self.path, "w", encoding="utf-8") as f:
            f.write(json.dumps(self._header, separators=(",", ":")) + "\n")
            f.flush()
            os.fsync(f.fileno())
        self._started = True

    def adopt(self, chunks_done):
        """
        Continue appending to an existing checkpoint rather than starting a new one.

        Called after a successful `load()`; without it the first `record()` would
        truncate the file being resumed from.
        """
        if chunks_done:
            self._started = True


def _merge_hits(target, delta):
    """Fold one chunk's hit counts into an accumulator, summing where they overlap."""
    for genome_id, counts in delta.items():
        into = target.setdefault(genome_id, {})
        for acc, n in counts.items():
            into[acc] = into.get(acc, 0) + n


def search_profiles(filtered_hmm_path, fasta_path, threads=1,
                    progress_callback=None, residue_budget=None,
                    checkpoint_path=None, resume=False):
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

    Those genomes aren't all reported at once at the end of the chunk, though. A chunk
    can take a long time. On a large run the whole phase is only a handful of chunks,
    and reporting per chunk left the bar apparently frozen for tens of minutes. So the
    genomes a chunk is going to complete are counted up front and then paid out across
    the chunk's profile passes as they finish (see `_ChunkProgress`), which is a fine
    enough clock to keep the bar moving without changing what the numbers mean.

    `checkpoint_path`, if given, is where each finished chunk is recorded, and with
    `resume=True` it's replayed so an interrupted search picks up at the first chunk
    that never completed (see `_SearchCheckpoint`). Skipped chunks are still *read*,
    but they aren't searched

    `residue_budget` overrides the chunk policy for deterministic testing; left None it
    uses `resolve_residue_budget`.
    """
    if residue_budget is not None:
        budget = residue_budget
    else:
        budget = resolve_residue_budget()

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
        pressed_base, num_profiles = _press_profiles(filtered_hmm_path, press_dir)

        checkpoint = _SearchCheckpoint(
            checkpoint_path, budget, num_profiles, fasta_path)

        if resume:
            hits_by_genome, chunks_done = checkpoint.load()
            checkpoint.adopt(chunks_done)
        else:
            chunks_done = 0

        chunk_index = 0

        while True:
            chunk = _read_chunk(seq_file, budget)
            if chunk is None:
                break

            # register every genome in this chunk before searching, so genomes that hit
            # nothing still end up as keys; the same pass tracks genome transitions.
            for seq in chunk:
                genome_id = genome_id_from_protein_name(decode_pyhmmer_text(seq.name))
                hits_by_genome.setdefault(genome_id, {})
                if genome_id != current_genome:
                    genomes_started += 1
                    current_genome = genome_id

            # all but the genome this chunk ends on are complete once the chunk is
            # searched; that last one may continue into the next chunk, so it waits.
            # The credit is known before the search runs, which is what lets it be paid
            # out gradually *during* the search instead of in one jump afterwards.
            finished = max(0, genomes_started - 1)
            progress = _ChunkProgress(
                progress_callback, finished - genomes_reported, num_profiles)

            if chunk_index < chunks_done:
                # already searched in the interrupted run; its hits came back with the
                # checkpoint, so just settle the bar and move on
                progress.finish()
            else:
                # a fresh dict per chunk, so what this chunk found can be checkpointed
                # on its own; it's folded into the running totals straight after
                chunk_hits = {}
                _search_one_chunk(pressed_base, chunk, threads, chunk_hits,
                                  on_profile=progress.on_profile)
                progress.finish()
                _merge_hits(hits_by_genome, chunk_hits)
                checkpoint.record(chunk_index, chunk_hits)

            genomes_reported += progress.awarded

            total_seqs += len(chunk)
            chunk_index += 1

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

MEASURED_BYTES_PER_RESIDUE = 2.9
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
    except (TypeError, ValueError) as e:
        raise ChunkSizeError(
            f"{CHUNK_ENV_VAR} must be a positive integer (residues) if set, got {raw!r}.") from e

    if value < 1:
        raise ChunkSizeError(
            f"{CHUNK_ENV_VAR} must be >= 1 if set, got {value}.")

    return value


_APPROX_FASTA_BYTES_PER_RESIDUE = 1.1


def _estimate_total_residues(fasta_path):
    """
    Cheap O(1) estimate of the residues in a fasta, from its size on disk.

    Only used for reporting an expected chunk count; nothing depends on it being exact.
    Returns None if the size can't be read.
    """
    try:
        size = os.path.getsize(fasta_path)
    except OSError:
        return None
    return int(size / _APPROX_FASTA_BYTES_PER_RESIDUE)


def resolve_residue_budget(env=None):
    """
    The residues-per-chunk budget for this run: `GTT_SCG_CHUNK` if set, otherwise the
    memory cap.

    Deliberately independent of the run's size. The budget is what makes chunk
    boundaries reproducible between a run and its `--resume`, so it must not drift with
    anything that could differ between the two.
    """
    override = _read_env_override(env)
    if override is not None:
        return override

    return max(1, DEFAULT_MAX_RESIDUES_PER_CHUNK)


def estimated_chunk_bytes(residue_budget):
    """
    Rough resident bytes for a chunk of `residue_budget` residues, from the measured
    per-residue figure. For reporting and for sanity-checking a chosen budget; the
    hmmsearch pipeline's own working memory is *not* included (see CALIBRATION).
    """
    return int(residue_budget * MEASURED_BYTES_PER_RESIDUE)
