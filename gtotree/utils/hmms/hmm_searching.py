import os
import contextlib
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
import pyhmmer.easel as easel #type: ignore
from gtotree.utils.misc.general import (search_threads_per_genome,
                                   atomic_write_text)
from gtotree.utils.hmms.hmm_searching_engine import search_one_genome
from gtotree.utils.misc.stages import GenomeRemovalStage


MAX_OPEN_SCG_HANDLES = 64

def _hmm_search_worker(genome, run_data, aa_path=None, nt_path=None, pressed_base=None):
    """
    Per-genome SCG search.

    `aa_path`/`nt_path` override the paths recorded on the GenomeData. The fused
    processing stage searches a genome inside the same worker that just produced its
    FASTA, before apply_result has copied those paths onto the GenomeData, so it passes
    them explicitly; the standalone path leaves them None and reads the object.

    Wrapped so it cannot raise: run_pooled_stage aborts the whole stage on a worker
    exception, leaving some items applied and some not.
    """
    try:
        return _hmm_search_worker_inner(genome, run_data, aa_path, nt_path, pressed_base)
    except BaseException as e:
        return {
            "hmm_search_failed": True,
            "extract_seqs_failed": False,
            "num_SCG_hits": 0,
            "num_unique_SCG_hits": 0,
            "error": f"{type(e).__name__}: {e}",
        }


def _hmm_search_worker_inner(genome, run_data, aa_path=None, nt_path=None, pressed_base=None):

    ID = genome.id
    AA_path = aa_path if aa_path is not None else genome.final_AA_path
    out_dir = f"{run_data.hmm_results_dir}/{ID}"
    os.makedirs(out_dir, exist_ok=True)
    hmm_out_path = f"{out_dir}/SCG-hits-hmm.txt"

    hmm_search_failed = run_hmm_search(ID, run_data, AA_path, hmm_out_path,
                                       pressed_base=pressed_base)

    num_SCG_hits = 0
    num_unique_SCG_hits = 0
    extract_seqs_failed = False

    if not hmm_search_failed:
        (dict_of_hit_counts, dict_of_hit_gene_ids,
         num_SCG_hits, num_unique_SCG_hits,
         shared_gene_hits) = parse_hmmer_results_with_shared(hmm_out_path, run_data)

        # These per-genome files are the source of truth the combined outputs are
        # rebuilt from at the end of the stage, so they're written atomically
        write_genome_hit_counts(f"{out_dir}/SCG-hit-counts.txt", dict_of_hit_counts)
        write_genome_shared_gene_hits(f"{out_dir}/{GENOME_SHARED_HITS_FILENAME}",
                                      shared_gene_hits)

        AA_hit_seqs_dict, extract_seqs_failed = get_seqs(dict_of_hit_gene_ids, AA_path)

        if not extract_seqs_failed:
            def _write_aa(f):
                for gene_id, seq in AA_hit_seqs_dict.items():
                    if seq is not None:
                        f.write(f">{gene_id}\n{seq}\n")

            atomic_write_text(f"{out_dir}/SCG-hits.faa", _write_aa)

            if run_data.nucleotide_mode:
                genome_nt_path = nt_path if nt_path is not None else genome.final_nt_path
                nt_hit_seqs_dict, _ = get_seqs(dict_of_hit_gene_ids, genome_nt_path)

                def _write_nt(f):
                    for gene_id, seq in nt_hit_seqs_dict.items():
                        if seq is not None:
                            f.write(f">{gene_id}\n{seq}\n")

                atomic_write_text(f"{out_dir}/SCG-hits.fasta", _write_nt)

    return {
        "hmm_search_failed": bool(hmm_search_failed),
        "extract_seqs_failed": bool(extract_seqs_failed),
        "num_SCG_hits": int(num_SCG_hits),
        "num_unique_SCG_hits": int(num_unique_SCG_hits),
    }


def _apply_hmm_search_result(genome, status, run_data):

    hmm_search_failed = bool(status.get("hmm_search_failed", False))
    extract_seqs_failed = bool(status.get("extract_seqs_failed", False))

    if hmm_search_failed:
        genome.mark_hmm_search_failed()
        genome.mark_removed("HMM search failed", GenomeRemovalStage.HMM_SEARCH)
        genome.num_SCG_hits = 0
    else:
        if extract_seqs_failed:
            genome.mark_extract_seqs_failed()
            genome.num_SCG_hits = 0
            genome.mark_removed("extracting sequences after HMM search failed",
                                GenomeRemovalStage.HMM_SEARCH)
        else:
            genome.mark_hmm_search_done()
            genome.num_SCG_hits = int(status.get("num_SCG_hits", 0))
            genome.num_unique_SCG_hits = int(status.get("num_unique_SCG_hits", 0))


def run_hmm_search(id, run_data, inpath, outpath, args=None, pressed_base=None):
    """
    Search one genome against the SCG profiles, in-process via pyhmmer.

    Was a subprocess `hmmsearch --cut_ga` call, which spawned a process and re-parsed
    the whole profile set from text for every genome. The tblout written here is
    field-for-field identical to what the CLI produced, so everything downstream is
    unchanged.

    `pressed_base` points at the once-per-run hmmpress output; without it the plain HMM
    file is read instead (slower, but keeps standalone callers working).
    """
    try:
        search_one_genome(run_data.hmm_path, inpath, outpath,
                          pressed_base=pressed_base,
                          cpus=search_threads_per_genome(args))
        hmm_search_failed = False
    except Exception as e:
        print(f"[SCG search failed for {id}] {type(e).__name__}: {e}")
        hmm_search_failed = True

    return hmm_search_failed


def read_hmmer_results(inpath):
    # col 5 is the full-sequence bit score, used to settle a gene hit by >1 SCG target
    colnames = ["gene_id", "target_SCG", "accession", "evalue", "score"]
    df = pd.read_csv(inpath, sep=r'\s+',
                     comment="#", header=None,
                     usecols=[0,2,3,4,5],
                     names=colnames)

    return df


def parse_hmmer_results(inpath, run_data):
    """
    Tally per-SCG hits for one genome and pick the gene id to pull for each

    A gene is never assigned to more than one SCG target (see
    `parse_hmmer_results_with_shared`, which also reports any such conflicts)
    """
    return parse_hmmer_results_with_shared(inpath, run_data)[:4]


class SharedGeneHit:
    """One gene picked for more than one SCG target, and how it was settled."""

    __slots__ = ("gene_id", "kept_SCG", "kept_score", "dropped_SCG", "dropped_score")

    def __init__(self, gene_id, kept_SCG, kept_score, dropped_SCG, dropped_score):
        self.gene_id = gene_id
        self.kept_SCG = kept_SCG
        self.kept_score = kept_score
        self.dropped_SCG = dropped_SCG
        self.dropped_score = dropped_score


def resolve_shared_gene_hits(dict_of_hit_gene_ids, scores, target_order):
    """
    Make sure no gene is pulled for more than one SCG target.

    A gene can be the chosen hit for two targets, e.g., when they're domains of one
    fused protein. Pulling it for both would put the same sequence into the tree twice,
    so the gene goes to the target it scores best against, and the other target gets no
    sequence from this genome. Ties go to whichever target comes first in the HMM file.

    `scores` is {(gene_id, target): full-sequence bit score}. Mutates
    `dict_of_hit_gene_ids` in place and returns [SharedGeneHit, ...].
    """
    rank = {target: i for i, target in enumerate(target_order)}

    targets_by_gene = {}
    for target, gene_id in dict_of_hit_gene_ids.items():
        if gene_id is not None:
            targets_by_gene.setdefault(gene_id, []).append(target)

    shared = []
    for gene_id, targets in targets_by_gene.items():
        if len(targets) < 2:
            continue
        ordered = sorted(targets, key=lambda t: (-scores.get((gene_id, t), float("-inf")),
                                                 rank.get(t, len(rank))))
        kept = ordered[0]
        for dropped in ordered[1:]:
            dict_of_hit_gene_ids[dropped] = None
            shared.append(SharedGeneHit(gene_id, kept, scores.get((gene_id, kept)),
                                        dropped, scores.get((gene_id, dropped))))
    return shared


def parse_hmmer_results_with_shared(inpath, run_data):
    """
    `parse_hmmer_results`, plus the list of SharedGeneHit conflicts it settled
    """
    df = read_hmmer_results(inpath)

    remaining_SCG_targets = [SCG_target.id for SCG_target in run_data.get_all_SCG_targets_remaining()]

    counts = df["target_SCG"].value_counts()
    first_gene_by_scg = df.groupby("target_SCG", sort=False)["gene_id"].first()

    dict_of_hit_counts = {scg: int(counts.get(scg, 0)) for scg in remaining_SCG_targets}
    dict_of_hit_gene_ids = dict.fromkeys(remaining_SCG_targets, None)

    num_SCG_hits = 0
    num_unique_SCG_hits = 0
    for scg, count in dict_of_hit_counts.items():

        if count == 1:
            num_SCG_hits += 1
            num_unique_SCG_hits += 1
            dict_of_hit_gene_ids[scg] = first_gene_by_scg[scg]

        elif count > 1:
            num_SCG_hits += 1
            if run_data.best_hit_mode:
                dict_of_hit_gene_ids[scg] = first_gene_by_scg[scg]

    scores = {}
    for gene_id, target, score in zip(df["gene_id"], df["target_SCG"], df["score"]):
        key = (gene_id, target)
        # hmmsearch emits one row per gene-target pair; keep the first if ever repeated
        if key not in scores:
            scores[key] = float(score)

    shared = resolve_shared_gene_hits(dict_of_hit_gene_ids, scores,
                                      remaining_SCG_targets)

    # the hit counts are left as found; they describe the search, and the table built
    # from them should still show that the losing target was hit

    return (dict_of_hit_counts, dict_of_hit_gene_ids, num_SCG_hits,
            num_unique_SCG_hits, shared)


def get_seqs(dict_of_hit_gene_ids, path):

    try:
        hit_seqs_dict = dict.fromkeys(dict_of_hit_gene_ids.keys(), None)
        reverse_lookup = {}
        for target, seq_id in dict_of_hit_gene_ids.items():
            if seq_id is None:
                continue
            if seq_id in reverse_lookup:
                # parse_hmmer_results settles these; getting here means a caller
                # skipped that, and silently keeping one would lose a target unseen
                raise ValueError(f"gene '{seq_id}' was assigned to both "
                                 f"'{reverse_lookup[seq_id]}' and '{target}'")
            reverse_lookup[seq_id] = target
        easel_alphabet = easel.Alphabet.amino() if path.endswith(".faa") else easel.Alphabet.dna()

        with easel.SequenceFile(path, digital=True, alphabet=easel_alphabet) as seq_file:
            for seq in seq_file:
                # pyhmmer changed seq.name from bytes to str across some versions i was testing on, but i didn't pin it down
                # so just handling either
                seq_id = seq.name.decode("utf8") if isinstance(seq.name, bytes) else seq.name
                if seq_id in reverse_lookup:
                    target_scg = reverse_lookup[seq_id]
                    hit_seqs_dict[target_scg] = easel_alphabet.decode(seq)
        extract_seqs_failed = False
    except Exception:
        extract_seqs_failed = True
        hit_seqs_dict = None

    return hit_seqs_dict, extract_seqs_failed


GENOME_HIT_COUNTS_HEADER = "target_SCG\tnum_hits"

GENOME_SHARED_HITS_FILENAME = "SCG-shared-gene-hits.tsv"
SHARED_HITS_COLUMNS = ["gene_id", "kept_SCG", "kept_score", "dropped_SCG",
                       "dropped_score"]
COMBINED_SHARED_HITS_FILENAME = "SCG-shared-gene-hits.tsv"


def _score_str(score):
    return "NA" if score is None else f"{score:g}"


def write_genome_shared_gene_hits(path, shared):
    """
    Write one genome's SharedGeneHit conflicts, or remove a stale file if there were
    none (so a re-search on resume can't leave an old report behind)
    """
    if not shared:
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        return

    def _write(f):
        f.write("\t".join(SHARED_HITS_COLUMNS) + "\n")
        for hit in shared:
            f.write("\t".join([hit.gene_id, hit.kept_SCG, _score_str(hit.kept_score),
                               hit.dropped_SCG, _score_str(hit.dropped_score)]) + "\n")

    atomic_write_text(path, _write)


def read_genome_shared_gene_hits(path):
    """
    One genome's shared-gene rows as lists of fields, [] if there are none
    """
    try:
        with open(path) as f:
            lines = [line.rstrip("\n") for line in f if line.strip()]
    except OSError:
        return []
    if not lines or lines[0] != "\t".join(SHARED_HITS_COLUMNS):
        return []
    return [line.split("\t") for line in lines[1:]
            if len(line.split("\t")) == len(SHARED_HITS_COLUMNS)]


def write_genome_hit_counts(path, dict_of_hit_counts):
    """
    Write one genome's per-SCG hit counts, keyed by target name
    """
    def _write(f):
        f.write(GENOME_HIT_COUNTS_HEADER + "\n")
        for scg, count in dict_of_hit_counts.items():
            f.write(f"{scg}\t{count}\n")

    atomic_write_text(path, _write)


def read_genome_hit_counts(path, target_SCG_ids):
    """
    One genome's {target SCG id -> hit count}, or None if it can't be trusted
    """
    try:
        with open(path) as f:
            lines = [line.rstrip("\n") for line in f if line.strip()]
    except OSError:
        return None

    if not lines:
        return None

    if lines[0] == GENOME_HIT_COUNTS_HEADER:
        counts = {}
        for line in lines[1:]:
            name, _tab, value = line.partition("\t")
            try:
                counts[name] = int(value)
            except ValueError:
                continue
        return counts

    fields = lines[0].split("\t")[1:]  # [0] was the assembly id
    if len(fields) != len(target_SCG_ids):
        return None

    counts = {}
    for scg_id, field in zip(target_SCG_ids, fields):
        try:
            counts[scg_id] = int(field)
        except ValueError:
            continue
    return counts


def rebuild_combined_SCG_outputs(run_data):
    """
    Rebuild the combined SCG outputs from the per-genome artifacts
    """
    genomes = [gd for gd in run_data.all_input_genomes
               if gd.hmm_search_done and not gd.removed]

    scgs = run_data.get_all_SCG_targets_remaining()
    target_SCG_ids = [SCG.id for SCG in scgs]

    genomes_with_any_hit = dict.fromkeys(target_SCG_ids, 0)
    genomes_with_usable_seq = dict.fromkeys(target_SCG_ids, 0)

    # --- the combined per-genome hit-count table ---
    def _write_table(out):
        out.write("genome_id\t" + "\t".join(target_SCG_ids) + "\n")
        for gd in genomes:
            counts = read_genome_hit_counts(
                f"{run_data.hmm_results_dir}/{gd.id}/SCG-hit-counts.txt",
                target_SCG_ids)
            if counts is None:
                continue

            row = []
            for scg_id in target_SCG_ids:
                count = counts.get(scg_id, 0)
                if count > 0:
                    genomes_with_any_hit[scg_id] += 1
                row.append(str(count))

            out.write(gd.id + "\t" + "\t".join(row) + "\n")

    atomic_write_text(f"{run_data.output_dir}/SCG-hit-counts.tsv", _write_table)

    ext = run_data.general_ext
    out_paths = {t: f"{run_data.found_SCG_seqs_dir}/{t}{ext}" for t in target_SCG_ids}

    hits_paths = [(gd.id, f"{run_data.hmm_results_dir}/{gd.id}/SCG-hits{ext}")
                  for gd in genomes]
    hits_paths = [(gid, p) for gid, p in hits_paths if os.path.isfile(p)]

    try:
        for i in range(0, len(target_SCG_ids), MAX_OPEN_SCG_HANDLES):
            batch = target_SCG_ids[i:i + MAX_OPEN_SCG_HANDLES]
            with contextlib.ExitStack() as stack:
                handles = {
                    target: stack.enter_context(open(f"{out_paths[target]}.part", "w"))
                    for target in batch
                }
                for gid, hits_path in hits_paths:
                    with open(hits_path) as infile:
                        for record in SeqIO.parse(infile, "fasta"):
                            handle = handles.get(record.id)
                            if handle is not None:
                                handle.write(f">{gid}\n{record.seq}\n")
                                genomes_with_usable_seq[record.id] += 1
    except BaseException:
        for path in out_paths.values():
            try:
                os.remove(f"{path}.part")
            except FileNotFoundError:
                pass
        raise

    for path in out_paths.values():
        os.replace(f"{path}.part", path)

    for scg in scgs:
        scg.num_genomes_with_hit = genomes_with_any_hit[scg.id]
        scg.num_genomes_after_copy_filtering = genomes_with_usable_seq[scg.id]

    rebuild_combined_shared_gene_hits(run_data, genomes, scgs)

    return run_data


def rebuild_combined_shared_gene_hits(run_data, genomes, scgs):
    """
    Combine the per-genome shared-gene reports into run-files/, and record on each
    SCG-set how many genomes' hits it lost to another target that way. Returns the
    number of conflicts (0 means no combined file was written).
    """
    lost = {scg.id: 0 for scg in scgs}
    rows = []
    for gd in genomes:
        path = f"{run_data.hmm_results_dir}/{gd.id}/{GENOME_SHARED_HITS_FILENAME}"
        for fields in read_genome_shared_gene_hits(path):
            rows.append([gd.id] + fields)
            dropped = fields[3]
            if dropped in lost:
                lost[dropped] += 1

    for scg in scgs:
        scg.num_genomes_lost_to_shared_gene = lost[scg.id]

    run_data.num_shared_gene_hits = len(rows)

    if not run_data.run_files_dir:
        return len(rows)

    out_path = f"{run_data.run_files_dir}/{COMBINED_SHARED_HITS_FILENAME}"
    if not rows:
        try:
            os.remove(out_path)
        except FileNotFoundError:
            pass
        return 0

    def _write(f):
        f.write("genome_id\t" + "\t".join(SHARED_HITS_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    atomic_write_text(out_path, _write)
    return len(rows)


def no_hits_reason(scg, best_hit_mode):
    """
    Why an SCG-set came out of the search with no usable sequences
    """
    num_with_hits = scg.num_genomes_with_hit or 0

    if not num_with_hits:
        return "no hits in any genome"

    plural = "" if num_with_hits == 1 else "s"

    num_lost = getattr(scg, "num_genomes_lost_to_shared_gene", None) or 0
    if num_lost:
        if num_lost >= num_with_hits:
            which = "that one" if num_with_hits == 1 else "all of them"
        else:
            which = f"{num_lost} of them"
        return (f"hits in {num_with_hits} genome{plural}, but in {which} the gene hit "
                "scored better against another SCG target and was kept for that one "
                f"instead (see {COMBINED_SHARED_HITS_FILENAME})")

    if not best_hit_mode:
        return (f"hits in {num_with_hits} genome{plural}, but never as a single copy "
                "(`-B`/`--best-hit-mode` would retain it)")

    return (f"hits in {num_with_hits} genome{plural}, but no sequences could be "
            "extracted for them")



