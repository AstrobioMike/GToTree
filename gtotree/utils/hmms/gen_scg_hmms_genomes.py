"""
Target-genome handling for `gtt gen-scg-hmms`

Resolves the requested target genomes (an accessions file and/or `--wanted-ref-tax`)
and gets an amino-acid fasta for each, reusing GToTree's existing NCBI machinery:

  * `select_ref_genomes` (via `resolve_wanted_ref_tax_accessions`) for `-W`, the same
    shared selection core the driver and get-accs-from-* helpers use
  * the hosted NCBI Parquet table for accession -> download-directory resolution
  * `download_and_unzip_accession` for the actual fetch (atomic writes, retry with
    backoff, fast-fail on 404)
"""

import os
import subprocess
from types import SimpleNamespace

import pyarrow.parquet as pq  # type: ignore

from gtotree.utils.general import run_pooled_stage
from gtotree.utils.ncbi.get_ncbi_assembly_data import ncbi_data_table_path
from gtotree.utils.ncbi.parse_ncbi_assembly_summary import resolve_base_link
from gtotree.utils.hmms.gen_scg_hmms import GenSCGHMMsError, _remove_quietly


MAX_DOWNLOAD_THREADS = 20

# reasons an accession can drop out, recorded in the missed-accessions report
MISSED_NOT_FOUND = "not found in NCBI assembly data"
MISSED_NO_LINK = "no resolvable download location"
MISSED_DOWNLOAD_FAILED = "download failed"
MISSED_PRODIGAL_FAILED = "prodigal failed"
MISSED_NO_PROTEINS = "no protein sequences found"


class TargetGenomeError(GenSCGHMMsError):
    """The requested target genomes could not be resolved."""


def _download_and_unzip():
    """
    Lazily grab the shared atomic downloader from the preprocessing stage.

    Imported on first use rather than at module scope because
    `main_stages.preprocessing_genomes` pulls in a bunch of other crap,
    none of which this standalone helper otherwise needs. The
    download behavior itself (atomic write, retry/backoff, fast-fail on 404) is
    deliberately shared rather than reimplemented.
    """
    from gtotree.main_stages.preprocessing_genomes import download_and_unzip_accession
    return download_and_unzip_accession


def read_accessions_file(path):
    """
    Read an accessions file (one per line, blank lines and #-comments ignored),
    preserving first-seen order while removing duplicates.
    """
    if not os.path.isfile(path):
        raise TargetGenomeError(
            f"the specified target-accessions file '{path}' can't be found.")

    accessions = []
    seen = set()
    with open(path) as f:
        for line in f:
            acc = line.strip()
            if not acc or acc.startswith("#"):
                continue
            if acc not in seen:
                seen.add(acc)
                accessions.append(acc)

    if not accessions:
        raise TargetGenomeError(
            f"no accessions were read out of '{path}'.")

    return accessions


def resolve_download_info(accessions, table_path=None):
    """
    Look the wanted accessions up in the hosted NCBI Parquet table.

    Returns (info, not_found) where `info` maps accession -> dict with the
    download base link and assembly name, and `not_found` lists accessions absent
    from the table.

    Matching is done on the accession as given. GCF/GCA are both present in the
    table, so no cross-prefix translation is done here -- an accession is either
    in NCBI's data or it isn't.
    """
    if table_path is None:
        table_path = ncbi_data_table_path()

    if not os.path.isfile(table_path):
        raise TargetGenomeError(
            "the NCBI assembly data table wasn't found; it can be retrieved with "
            "`gtt data get ncbi`.")

    wanted = set(accessions)

    table = pq.read_table(
        table_path,
        columns=["assembly_accession", "asm_name", "ftp_path", "organism_name"],
        filters=[("assembly_accession", "in", wanted)],
    )

    info = {}
    for row in table.to_pylist():
        acc = row.get("assembly_accession")
        if not acc:
            continue
        base_link = resolve_base_link(row.get("ftp_path"), acc, row.get("asm_name"))
        info[acc] = {
            "base_link": base_link,
            "assembly_name": row.get("asm_name"),
            "organism_name": row.get("organism_name"),
        }

    not_found = [acc for acc in accessions if acc not in info]

    return info, not_found


def fetch_amino_acids_pooled(to_fetch, info, work_dir, args=None, num_jobs=1,
                             nucleotide_fallback=True, on_result=None,
                             progress_callback=None):
    """
    Fetch amino acids for many accessions concurrently.

    Thin wrapper over the shared `run_pooled_stage` helper so this program and the main
    GToTree run use one pooling implementation rather than two. That
    helper's contract is the one that matters here: the worker runs in a thread and
    touches only per-accession paths, while `apply_result` runs single-threaded on the
    main thread, letting each finished genome into the one combined fasta without records interleaving.

    Each worker leaves a per-genome faa on disk and `on_result(acc, aa_path,
    used_prodigal, error)` is invoked as it lands, so the caller can absorb and delete
    it immediately. That streaming is deliberate: accumulating every per-genome faa and
    combining at the end would mean tens of GB in the working dir at GToTree's 10k-30k
    genome scale. Only `num_jobs` files exist at once this way.

    Note combining order follows completion order, not input order. Genome identity is
    carried in each header and the output tables are written from the caller's ordered
    accession list, so results are unaffected.

    `args` supplies `num_jobs` for the shared helper; `num_jobs=` is accepted directly
    for callers that don't have an args namespace.
    """
    if not to_fetch:
        return

    if args is None:
        args = SimpleNamespace(num_jobs=num_jobs)

    def worker(acc, rd):
        """Thread-side: never raises, returns a plain status dict."""
        try:
            aa_path, used_prodigal = fetch_amino_acids(
                acc, rd["info"].get(acc), rd["work_dir"],
                nucleotide_fallback=rd["nucleotide_fallback"])
            return {"aa_path": aa_path, "prodigal": used_prodigal, "error": None}
        except TargetGenomeError as e:
            return {"aa_path": None, "prodigal": False, "error": str(e)}
        except BaseException as e:  # noqa: BLE001 - reported per-accession
            return {"aa_path": None, "prodigal": False, "error": f"unexpected failure: {e}"}

    def apply_result(acc, status, rd):
        """Main-thread side: safe to touch the shared combined fasta here."""
        if on_result is not None:
            on_result(acc, status["aa_path"], status["prodigal"], status["error"])
        if progress_callback is not None:
            progress_callback()

    run_pooled_stage(to_fetch, worker, apply_result, args,
                     {"info": info, "work_dir": work_dir,
                      "nucleotide_fallback": nucleotide_fallback},
                     max_workers_cap=MAX_DOWNLOAD_THREADS)


def fetch_amino_acids(accession, entry, work_dir, nucleotide_fallback=True):
    """
    Get an amino-acid fasta for one accession into `work_dir`.

    Tries the NCBI protein file first; on failure (commonly a 404) falls back to
    downloading the nucleotide fasta and calling genes with prodigal.

    Returns (out_path, used_prodigal) on success, or raises TargetGenomeError with a
    reason string suitable for the missed-accessions report.
    """
    base_link = (entry or {}).get("base_link")
    if not base_link or str(base_link).lower() == "na":
        raise TargetGenomeError(MISSED_NO_LINK)

    base_link = str(base_link).replace(" ", "_").rstrip("/")
    assembly_str = base_link.split("/")[-1]

    aa_path = os.path.join(work_dir, f"{accession}_protein.faa")

    download_and_unzip_accession = _download_and_unzip()

    try:
        download_and_unzip_accession(
            f"{base_link}/{assembly_str}_protein.faa.gz", aa_path)
        return aa_path, False
    except KeyboardInterrupt:
        raise
    except Exception:
        _remove_quietly(aa_path)
        if not nucleotide_fallback:
            raise TargetGenomeError(MISSED_DOWNLOAD_FAILED)

    # no protein file -> get nucleotides and call genes
    nt_path = os.path.join(work_dir, f"{accession}_genomic.fna")
    try:
        download_and_unzip_accession(
            f"{base_link}/{assembly_str}_genomic.fna.gz", nt_path)
    except KeyboardInterrupt:
        raise
    except Exception:
        _remove_quietly(nt_path)
        raise TargetGenomeError(MISSED_DOWNLOAD_FAILED)

    try:
        run_prodigal(nt_path, aa_path)
    finally:
        _remove_quietly(nt_path)

    return aa_path, True


def run_prodigal(nucleotide_path, out_aa_path):
    """
    Call genes on a nucleotide fasta with prodigal. Writes atomically so a failed or
    interrupted call never leaves a partial faa behind.
    """
    tmp_path = out_aa_path + ".part"
    try:
        result = subprocess.run(
            ["prodigal", "-c", "-q", "-i", nucleotide_path, "-a", tmp_path],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
        if result.returncode != 0:
            raise TargetGenomeError(MISSED_PRODIGAL_FAILED)
        os.replace(tmp_path, out_aa_path)
    except FileNotFoundError:
        _remove_quietly(tmp_path)
        raise TargetGenomeError(
            "prodigal doesn't seem to be available, but it's needed for genomes that "
            "have no protein file at NCBI.")
    except BaseException:
        _remove_quietly(tmp_path)
        raise


def relabel_and_append(accession, aa_path, combined_handle):
    """
    Append one genome's proteins to the combined fasta, renaming each header to
    `<accession>_<n>` so hits can be traced back to their genome, and stripping the
    `*` stop characters prodigal emits.

    Returns the number of sequences written.
    """
    n = 0
    seq_chunks = []

    def flush():
        if seq_chunks:
            combined_handle.write("".join(seq_chunks).replace("*", ""))
            combined_handle.write("\n")
            seq_chunks.clear()

    with open(aa_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                flush()
                n += 1
                combined_handle.write(f">{accession}_{n}\n")
            else:
                seq_chunks.append(line.strip())
        flush()

    if n == 0:
        raise TargetGenomeError(MISSED_NO_PROTEINS)

    return n


def genome_id_from_protein_name(protein_name):
    """
    Get the source accession from a relabeled protein header
    """
    return protein_name.rsplit("_", 1)[0]
