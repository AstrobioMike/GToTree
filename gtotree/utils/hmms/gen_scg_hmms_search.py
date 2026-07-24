"""
The hmmsearch stage of `gtt gen-scg-hmms`

Searches all target proteins against the coverage-filtered Pfam profiles and tallies,
per genome, how many times each profile was hit. That tally is what the single-copy
determination in `gen_scg_hmms.count_single_copy_hits` uses

Gathering thresholds (`--cut_ga`) are used just like in main GToTree
"""

import pyhmmer  # type: ignore

from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError, _decode
from gtotree.utils.hmms.gen_scg_hmms_genomes import genome_id_from_protein_name


# how many profiles to hand pyhmmer per batch
DEFAULT_HMM_BLOCK_SIZE = 200


class HmmSearchError(GenSCGHMMsError):
    """The hmmsearch stage failed."""


def load_target_proteins(fasta_path):
    """
    Read the combined protein fasta into a pyhmmer digital sequence block.

    The proteins are searched repeatedly (once per profile block), so they're loaded
    once up front and reused rather than re-read per block.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    try:
        with pyhmmer.easel.SequenceFile(fasta_path, digital=True,
                                        alphabet=alphabet) as seq_file:
            sequences = seq_file.read_block()
    except Exception as e:
        raise HmmSearchError(f"failed to read the combined protein fasta: {e}")

    if sequences is None or len(sequences) == 0:
        raise HmmSearchError(
            "no protein sequences were available to search; none of the target genomes "
            "yielded usable amino-acid sequences.")

    return alphabet, sequences


def _iter_hmm_blocks(hmm_path, block_size):
    """ Yield lists of profiles from `hmm_path`, `block_size` at a time. """
    block = []
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        for hmm in hmm_file:
            block.append(hmm)
            if len(block) >= block_size:
                yield block
                block = []
    if block:
        yield block


def search_profiles(filtered_hmm_path, fasta_path, threads=1,
                    block_size=DEFAULT_HMM_BLOCK_SIZE, progress_callback=None):
    """
    Search the filtered Pfam profiles against the combined target proteins.

    Returns `hits_by_genome`: {genome_id: {versioned_pfam_acc: hit_count}}.

    Only hits above each profile's gathering threshold are counted (bit_cutoffs="gathering",
    equivalent to hmmsearch's `--cut_ga`). Every reported hit increments the count for its
    source genome, so a genome with two copies of a domain counts 2
    """
    alphabet, sequences = load_target_proteins(fasta_path)

    hits_by_genome = {}

    for block in _iter_hmm_blocks(filtered_hmm_path, block_size):
        # NOTE: pyhmmer.hmmsearch returns a lazy generator, the pipeline is only
        # constructed (and its options validated) when results are consumed, so the
        # try/except has to wrap the ITERATION, not just the call. Wrapping only the
        # call lets real search errors escape as raw pyhmmer exceptions.
        try:
            results = pyhmmer.hmmsearch(
                block, sequences, cpus=threads, bit_cutoffs="gathering",
            )

            for top_hits in results:
                # pyhmmer exposes the profile that produced these hits as `.query`; take
                # its versioned accession (PF00001.27) since that's the key everything
                # else in this pipeline joins on, falling back to the profile name if
                # absent. In pyhmmer 0.11.0 these come back as bytes and in later
                # versions as str, so both go through _decode for an attempt at future-proofing
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

                if progress_callback is not None:
                    progress_callback()

        except KeyboardInterrupt:
            raise
        except HmmSearchError:
            raise
        except Exception as e:
            raise HmmSearchError(f"the hmmsearch step failed: {e}")

    return hits_by_genome
