import os
import contextlib
import pandas as pd # type: ignore
from Bio import SeqIO # type: ignore
import pyhmmer.easel as easel #type: ignore
from gtotree.utils.misc.general import (search_threads_per_genome,
                                   atomic_write_text)
from gtotree.utils.hmms.hmm_searching_engine import search_one_genome
from gtotree.utils.misc.stages import GenomeRemovalStage


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
         num_SCG_hits, num_unique_SCG_hits) = parse_hmmer_results(hmm_out_path, run_data)

        # These per-genome files are the source of truth the combined outputs are
        # rebuilt from at the end of the stage, so they're written atomically
        def _write_counts(f):
            f.write(f"{ID}")
            for count in dict_of_hit_counts.values():
                f.write(f"\t{count}")
            f.write("\n")

        atomic_write_text(f"{out_dir}/SCG-hit-counts.txt", _write_counts)

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
    colnames = ["gene_id", "target_SCG", "accession", "evalue"]
    df = pd.read_csv(inpath, sep=r'\s+',
                     comment="#", header=None,
                     usecols=[0,2,3,4],
                     names=colnames)

    return df


def parse_hmmer_results(inpath, run_data):
    """
    Tally per-SCG hits for one genome and pick the gene id to pull for each
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

    return dict_of_hit_counts, dict_of_hit_gene_ids, num_SCG_hits, num_unique_SCG_hits


def get_seqs(dict_of_hit_gene_ids, path):

    try:
        hit_seqs_dict = dict.fromkeys(dict_of_hit_gene_ids.keys(), None)
        reverse_lookup = {seq_id: target for target, seq_id in dict_of_hit_gene_ids.items()}
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
    except:
        extract_seqs_failed = True
        hit_seqs_dict = None

    return hit_seqs_dict, extract_seqs_failed


def rebuild_combined_SCG_outputs(run_data):
    """
    Rebuild the combined SCG outputs from the per-genome artifacts
    """
    genomes = [gd for gd in run_data.all_input_genomes
               if gd.hmm_search_done and not gd.removed]

    target_SCG_ids = [SCG.id for SCG in run_data.get_all_SCG_targets_remaining()]

    # --- the combined per-genome hit-count table ---
    def _write_table(out):
        out.write("assembly_id\t" + "\t".join(target_SCG_ids) + "\n")
        for gd in genomes:
            row_path = f"{run_data.hmm_results_dir}/{gd.id}/SCG-hit-counts.txt"
            try:
                with open(row_path) as row_file:
                    contents = row_file.read()
            except OSError:
                continue
            if contents and not contents.endswith("\n"):
                contents += "\n"
            out.write(contents)

    atomic_write_text(f"{run_data.output_dir}/SCG-hit-counts.tsv", _write_table)

    ext = run_data.general_ext
    out_paths = {t: f"{run_data.found_SCG_seqs_dir}/{t}{ext}" for t in target_SCG_ids}

    with contextlib.ExitStack() as stack:
        handles = {
            target: stack.enter_context(open(f"{path}.part", "w"))
            for target, path in out_paths.items()
        }
        for gd in genomes:
            hits_path = f"{run_data.hmm_results_dir}/{gd.id}/SCG-hits{ext}"
            if not os.path.isfile(hits_path):
                continue
            with open(hits_path) as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    handle = handles.get(record.id)
                    if handle is not None:
                        handle.write(f">{gd.id}\n{record.seq}\n")

    for path in out_paths.values():
        os.replace(f"{path}.part", path)

    return run_data


def capture_hmm_search_failures(run_data):
    if len(run_data.get_failed_hmm_search_paths()) > 0:
        with open(run_data.run_files_dir + "/inputs-that-failed-at-the-hmm-search.txt", "w") as fail_file:
            for genome_id in run_data.get_failed_hmm_search_paths():
                fail_file.write(genome_id + "\n")
