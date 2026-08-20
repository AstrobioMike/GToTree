"""
Per-source genome processing: turn each kind of input into a filtered, renamed
amino-acid FASTA (plus nucleotides in nucleotide mode).
"""

import pandas as pd # type: ignore
import urllib.request
import urllib.error
import gzip
import shutil
import os
import time
import random
import subprocess
import threading
from gtotree.utils.misc.seqs import (filter_and_rename_fasta,
                                extract_filter_and_rename_cds_amino_acids_from_gb,
                                extract_fasta_from_gb)
from gtotree.utils.misc.general import (gunzip_if_needed,
                                   remove_file_if_exists)
from gtotree.utils.misc.stages import GenomeRemovalStage


# Concurrent NCBI downloads are capped for politeness toward NCBI, NOT because of
# local resources. Capping the whole pool at this number would also throttle the
# CPU-bound following work to 20-wide no matter what `-j` said. So the cap is applied as a
# semaphore around the download call alone, the pool itself runs at --num-jobs
MAX_NCBI_DOWNLOAD_THREADS = 30
_NCBI_DOWNLOAD_SEMAPHORE = threading.Semaphore(MAX_NCBI_DOWNLOAD_THREADS)
NCBI_DOWNLOAD_MAX_RETRIES = 10

# Retry/backoff policy, mirroring bit's dl-ncbi-assemblies
NCBI_THROTTLE_STATUS = {429}

# ceiling on any single retry sleep (seconds), applied to the throttled path
NCBI_MAX_BACKOFF = 30

# separate, larger ceiling for a server-specified Retry-After
NCBI_MAX_RETRY_AFTER = 300

# Socket timeout for a single accession download, in seconds
NCBI_DOWNLOAD_TIMEOUT = 60

# length of the non-throttled sawtooth cycle, e.g., with 5 set here: 1, 2, 4, 8, 16, then start over
NCBI_SAWTOOTH_CYCLE = 5


class AssemblyFormatNotAvailable(Exception):
    """
    NCBI has no file of the requested format for this accession (HTTP 404).

    This is the ONE failure that legitimately means "try the other format". Everything
    else -- a dead network, a corrupt gzip, an unwritable directory, a bug in our own
    code -- means something is wrong, not that this assembly ships nucleotides instead
    of proteins.
    """


class AssemblyDownloadFailed(Exception):
    """
    Downloading this accession failed after exhausting retries: network, HTTP, or a
    corrupt archive. An accession-level failure, so it's reported per genome rather
    than aborting the run.
    """


def _process_one_ncbi_accession(acc_gd, run_data, base_link_map=None):
    """
    Worker: run one accession's download + (optional) prodigal + filter/rename. Runs in a worker
    thread, so it must not mutate shared run_data state, only local files and
    the returned dict. All GenomeData mutation happens back on the main thread in
    _apply_ncbi_accession_status.
    """
    try:
        # politeness cap toward NCBI, held only for the network work so the search
        # steps that follow aren't throttled by it (see _NCBI_DOWNLOAD_SEMAPHORE)
        with _NCBI_DOWNLOAD_SEMAPHORE:
            done, nt = prepare_accession(acc_gd.id, run_data, base_link_map=base_link_map)
        downloaded = bool(done)

        error = None
        if done and nt:
            done, error = run_prodigal(acc_gd.id, run_data, group="ncbi")
            prodigal_used = True
        else:
            prodigal_used = False

        if done:
            done, final_AA_path, num_genes, final_nt_path = \
                filter_and_rename_fasta(acc_gd.id, run_data, run_data.ncbi_processing_dir)
        else:
            final_AA_path = None
            final_nt_path = None
            num_genes = 0

        return {
            "done": bool(done),
            "downloaded": downloaded,
            "prodigal_used": prodigal_used,
            "final_AA_path": final_AA_path,
            "final_nt_path": final_nt_path,
            "num_genes": int(num_genes or 0),
            "error": error,
        }
    except BaseException as e:
        # a worker must never take down the whole pool; treat any unexpected
        # failure as a non-downloaded processing failure and let the bookkeeping
        # step mark it removed. The message is carried out with it
        return {
            "done": False,
            "downloaded": False,
            "prodigal_used": False,
            "final_AA_path": None,
            "final_nt_path": None,
            "num_genes": 0,
            "error": f"{type(e).__name__}: {e}",
        }


def _apply_ncbi_accession_status(acc_gd, status, run_data):
    """
    Apply one worker's result to its GenomeData. Called on the main thread only.
    """
    done = bool(status.get("done"))
    downloaded = bool(status.get("downloaded", False))
    prodigal_used = bool(status.get("prodigal_used", False))

    acc_gd.num_genes = int(status.get("num_genes", 0) or 0)

    if done:
        acc_gd.mark_processing_done()
        acc_gd.final_AA_path = status.get("final_AA_path")
        acc_gd.final_nt_path = status.get("final_nt_path")
        acc_gd.acc_was_downloaded = downloaded
    else:
        error = status.get("error")
        if downloaded:
            acc_gd.acc_was_downloaded = True
            reason = "acc processing failed after download"
        else:
            acc_gd.acc_was_downloaded = False
            reason = "acc download failed"
        # the reason lands in the summary table, so an unexpected failure says what it
        # was instead of being indistinguishable from a genuinely missing assembly
        acc_gd.mark_removed(f"{reason} ({error})" if error else reason,
                            GenomeRemovalStage.NCBI_DOWNLOAD)

    acc_gd.prodigal_used = prodigal_used
    if prodigal_used:
        run_data.tools_used.prodigal_used = True


def prepare_accession(acc, run_data, base_link_map=None):
    base_link, acc_assembly_str = get_base_link(acc, base_link_map)

    # an unresolvable download directory (no ftp_path and nothing to rebuild
    # from) comes through as "na" -> there's nothing to download
    if not base_link or base_link.lower() == "na":
        return False, False

    # Nucleotide mode goes straight to the genomic fasta. This used to be spelled
    # `raise Exception` inside the amino-acid try block, using an exception as a goto
    # into the fallback
    if not run_data.nucleotide_mode:
        amino_acid_link = base_link + "/" + acc_assembly_str + "_protein.faa.gz"
        amino_acid_filepath = run_data.ncbi_processing_dir + "/" + acc + "_protein.faa"
        try:
            download_and_unzip_accession(amino_acid_link, amino_acid_filepath)
            return True, False
        except (AssemblyFormatNotAvailable, AssemblyDownloadFailed):
            # no protein fasta for this assembly (or it wouldn't come down) -> fall
            # through and call genes off the nucleotides instead
            pass

    nucleotide_link = base_link + "/" + acc_assembly_str + "_genomic.fna.gz"
    nucleotide_file = run_data.ncbi_processing_dir + "/" + acc + "_genomic.fna"
    try:
        download_and_unzip_accession(nucleotide_link, nucleotide_file)
        return True, True
    except (AssemblyFormatNotAvailable, AssemblyDownloadFailed):
        return False, False


def download_and_unzip_accession(link, filepath, max_retries=NCBI_DOWNLOAD_MAX_RETRIES):
    # Transient failures (network blips, truncated/corrupt gzip, transient HTTP
    # 429/5xx) are retried, with the backoff policy chosen by why the attempt failed:
    # an explicit throttle (429, or any response carrying Retry-After) gets capped
    # exponential backoff, while a plain hiccup gets a sawtooth that never grows past
    # ~16s. A 404 is permanent (the requested format doesn't exist for this accession)
    # so we fail fast rather than burning retries; prepare_accession will then try the
    # other format.
    # A local filesystem problem is not a download problem, and must not be retried as
    # one: an unwritable destination would otherwise burn every retry and its backoff
    # before failing, per genome. Checked once, up front, and raised as a plain OSError,
    # deliberately NOT one of the typed exceptions above, so callers that handle
    # "this accession didn't work" can't accidentally swallow "our directory is wrong".
    dest_dir = os.path.dirname(os.path.abspath(filepath))
    if not os.path.isdir(dest_dir):
        raise OSError(f"the download directory '{dest_dir}' does not exist")
    if not os.access(dest_dir, os.W_OK):
        raise OSError(f"the download directory '{dest_dir}' is not writable")

    tmp_gzip = filepath + ".gz"
    tmp_out = filepath + ".part"

    def cleanup_partials():
        for tmp in (tmp_gzip, tmp_out):
            try:
                os.remove(tmp)
            except OSError:
                pass

    last_err = None
    throttled = False
    for attempt in range(1, max_retries + 1):
        try:
            _fetch_to_file(link, tmp_gzip)
            with gzip.open(tmp_gzip, 'rb') as f_in, open(tmp_out, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.replace(tmp_out, filepath)
            try:
                os.remove(tmp_gzip)
            except OSError:
                pass
            return
        except KeyboardInterrupt:
            cleanup_partials()
            raise
        except urllib.error.HTTPError as e:
            # NOTE: HTTPError subclasses URLError, so this must stay ahead of the
            # broader except below or status codes would never be seen.
            cleanup_partials()
            # 404 -> requested format not available; permanent, don't retry
            if e.code == 404:
                raise AssemblyFormatNotAvailable(
                    f"no file of this format at {link}") from e
            last_err = e
            throttled = _is_throttle(e)
        except (urllib.error.URLError, OSError, EOFError, gzip.BadGzipFile) as e:
            cleanup_partials()
            last_err = e
            # a connection reset or corrupt gzip carries no server instruction,
            # so it takes the sawtooth path
            throttled = False

        if attempt < max_retries:
            _sleep_backoff(attempt, err=last_err, throttled=throttled)

    # exhausted retries -> raise so prepare_accession falls back / marks failure
    if last_err is not None:
        raise AssemblyDownloadFailed(
            f"{type(last_err).__name__}: {last_err}") from last_err
    raise AssemblyDownloadFailed(f"failed to download {link}")


def _fetch_to_file(link, dest, timeout=None):
    """
    Stream `link` to `dest` with a socket timeout

    `urlopen(..., timeout=)` bounds the connect and every individual read, so a stalled
    server raises (and gets retried by the caller's backoff loop) rather than hanging.
    Note this is a per-operation timeout, not a total-transfer deadline: a server
    trickling bytes slowly will still be waited on, which is what we want, since large
    genomes over a slow link are legitimate.

    Errors propagate untouched so the caller keeps distinguishing HTTPError (status
    codes, Retry-After) from generic transport failures.
    """
    # resolved at call time, not bound as a default, so the module constant stays the
    # single source of truth (and remains patchable)
    if timeout is None:
        timeout = NCBI_DOWNLOAD_TIMEOUT

    req = urllib.request.Request(link, headers={"User-Agent": "curl/8.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out, length=1024 * 256)


def _sleep_backoff(attempt, err=None, throttled=False):
    """
    Sleep before a retry. Two policies, chosen by WHY the attempt failed (same split
    as bit's dl-ncbi-assemblies):

    throttled=True  -- the server is explicitly rate-limiting us (HTTP 429, or any
                       response carrying Retry-After). Honor Retry-After if given
                       (capped at NCBI_MAX_RETRY_AFTER), otherwise true exponential
                       backoff capped at NCBI_MAX_BACKOFF. Backing off progressively
                       is the point here: retrying faster accelerates into the thing
                       that's throttling us and tends to extend the penalty window.

    throttled=False -- a generic transient failure (5xx, connection reset, timeout,
                       truncated/corrupt gzip). These are blips or genuinely dead
                       URLs, and there's no politeness argument for stretching the
                       interval forever, so the wait SAWTOOTHS: 1, 2, 4, 8, 16, 1, 2,
                       4, 8, 16, ... That keeps a straggler from parking a pool thread
                       for many minutes while still spacing requests out, and we find
                       out a file is dead far sooner.

    `err` is the exception that caused the retry; its headers are consulted for
    Retry-After when throttled.
    """
    if throttled:
        retry_after = _retry_after_seconds(err)
        if retry_after is not None:
            time.sleep(min(retry_after, NCBI_MAX_RETRY_AFTER))
            return
        time.sleep(min((2 ** (attempt - 1)) + random.uniform(0, 1), NCBI_MAX_BACKOFF))
        return

    # sawtooth: exponential within a short cycle, then start over
    step = (attempt - 1) % NCBI_SAWTOOTH_CYCLE
    time.sleep((2 ** step) + random.uniform(0, 1))


def _retry_after_seconds(err):
    """
    Pull a numeric Retry-After (seconds) off an HTTPError, or None if absent/unusable.

    urllib's HTTPError carries the response headers, so unlike a bare urlretrieve call
    we can see the server's instruction here. A non-numeric value (HTTP allows an
    HTTP-date form) is treated as absent so the caller falls back to exponential
    backoff rather than guessing.
    """
    headers = getattr(err, "headers", None)
    if not headers:
        return None
    try:
        value = headers.get("Retry-After")
    except AttributeError:
        return None
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_throttle(err):
    """
    True when a failure is the server explicitly rate-limiting us, rather than a
    generic hiccup: an explicit throttle status, or any response that bothered to
    tell us when to come back.
    """
    code = getattr(err, "code", None)
    if code in NCBI_THROTTLE_STATUS:
        return True
    return _retry_after_seconds(err) is not None


def build_base_link_map(run_data):
    """
    {input accession -> its download base link}, from the NCBI sub-table
    """
    path = run_data.ncbi_sub_table_path
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(
            "the NCBI accession info table couldn't be found at "
            f'"{path or "<unset>"}", so there are no download links to work from')

    df = pd.read_csv(path, sep="\t",
                     usecols=["input_accession", "http_base_link"])
    return dict(zip(df["input_accession"], df["http_base_link"], strict=True))


def _normalize_base_link(raw_base_link):
    base_link = str(raw_base_link).replace(" ", "_")
    base_link = base_link.rstrip("/")
    acc_assembly_str = base_link.split("/")[-1]
    return base_link, acc_assembly_str


def get_base_link(acc, base_link_map):
    return _normalize_base_link(base_link_map[acc])


def _process_one_genbank_file(gb, run_data):
    """
    Worker: preprocess one genbank file. Tries CDS AA extraction first; falls back
    to extracting nucleotides + calling genes with prodigal. Runs in a worker
    thread, touches only per-file paths and returns a status dict; all GenomeData
    mutation happens on the main thread in _apply_genbank_status.
    """
    try:
        path, was_gzipped = gunzip_if_needed(gb.full_path,
                                             out_dir=run_data.genbank_processing_dir)

        prodigal_used = False
        final_nt_path = None
        final_AA_path = None
        num_genes = 0
        error = None

        done, final_AA_path, num_genes = extract_filter_and_rename_cds_amino_acids_from_gb(gb.id, path, run_data)

        if not done:
            extract_fasta_from_gb(gb.id, path, run_data)
            done, error = run_prodigal(gb.id, run_data, path, "genbank")
            prodigal_used = True
            if done:
                done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(gb.id, run_data, run_data.genbank_processing_dir)
            else:
                prodigal_used = False

        if was_gzipped:
            os.remove(path)

        return {
            "done": bool(done),
            "prodigal_used": prodigal_used,
            "was_gzipped": bool(was_gzipped),
            "final_AA_path": final_AA_path,
            "final_nt_path": final_nt_path,
            "num_genes": int(num_genes or 0),
            "error": error,
        }
    except BaseException as e:
        return {
            "done": False,
            "prodigal_used": False,
            "was_gzipped": False,
            "final_AA_path": None,
            "final_nt_path": None,
            "num_genes": 0,
            "error": f"{type(e).__name__}: {e}",
        }


def _apply_genbank_status(gb, status, run_data):
    """
    Apply one worker's result to its GenomeData. Main thread only.
    """
    gb.num_genes = int(status.get("num_genes", 0) or 0)

    if status.get("done"):
        if status.get("was_gzipped"):
            gb.mark_was_gzipped()
        gb.mark_processing_done()
        gb.final_AA_path = status.get("final_AA_path")
        gb.final_nt_path = status.get("final_nt_path")
    else:
        error = status.get("error")
        gb.mark_removed(f"genbank-file processing failed ({error})" if error
                        else "genbank-file processing failed",
                        GenomeRemovalStage.GENBANK_PREP)

    if status.get("prodigal_used"):
        gb.mark_prodigal_used()
        run_data.tools_used.prodigal_used = True


def _process_one_fasta_file(fasta, run_data):
    """
    Worker: preprocess one nucleotide FASTA -- call genes with prodigal, then
    filter/rename. Runs in a worker thread; per-file paths only, returns a status
    dict. GenomeData mutation happens on the main thread in _apply_fasta_status.
    """
    try:
        path, was_gzipped = gunzip_if_needed(fasta.full_path,
                                             out_dir=run_data.fasta_processing_dir)

        done, error = run_prodigal(fasta.id, run_data, path, "fasta")

        if was_gzipped:
            os.remove(path)

        if done:
            done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(fasta.id, run_data, run_data.fasta_processing_dir)
        else:
            final_AA_path = None
            final_nt_path = None
            num_genes = 0

        return {
            "done": bool(done),
            "was_gzipped": bool(was_gzipped),
            "prodigal_used": bool(done),
            "final_AA_path": final_AA_path,
            "final_nt_path": final_nt_path,
            "num_genes": int(num_genes or 0),
            "error": error,
        }
    except BaseException as e:
        return {
            "done": False,
            "was_gzipped": False,
            "prodigal_used": False,
            "final_AA_path": None,
            "final_nt_path": None,
            "num_genes": 0,
            "error": f"{type(e).__name__}: {e}",
        }


def _apply_fasta_status(fasta, status, run_data):
    """
    Apply one worker's result to its GenomeData. Main thread only.
    """
    fasta.num_genes = int(status.get("num_genes", 0) or 0)

    if status.get("done"):
        if status.get("was_gzipped"):
            fasta.mark_was_gzipped()
        fasta.mark_processing_done()
        fasta.final_AA_path = status.get("final_AA_path")
        fasta.final_nt_path = status.get("final_nt_path")
    else:
        error = status.get("error")
        fasta.mark_removed(f"fasta-file processing failed ({error})" if error
                           else "fasta-file processing failed",
                           GenomeRemovalStage.FASTA_PREP)

    if status.get("prodigal_used"):
        run_data.tools_used.prodigal_used = True


def _process_one_amino_acid_file(AA, run_data):
    """
    Worker: preprocess one amino-acid FASTA -- just filter/rename (no gene calling).
    Runs in a worker thread; per-file paths only, returns a status dict. GenomeData
    mutation happens on the main thread in _apply_amino_acid_status.
    """
    try:
        path, was_gzipped = gunzip_if_needed(AA.full_path,
                                             out_dir=run_data.AA_processing_dir)

        done, final_AA_path, num_genes, final_nt_path = filter_and_rename_fasta(AA.id, run_data, path, full_path=True)

        if was_gzipped:
            os.remove(path)

        return {
            "done": bool(done),
            "was_gzipped": bool(was_gzipped),
            "final_AA_path": final_AA_path,
            "final_nt_path": final_nt_path,
            "num_genes": int(num_genes or 0),
        }
    except BaseException as e:
        return {
            "done": False,
            "was_gzipped": False,
            "final_AA_path": None,
            "final_nt_path": None,
            "num_genes": 0,
            "error": f"{type(e).__name__}: {e}",
        }


def _apply_amino_acid_status(AA, status, run_data):
    """
    Apply one worker's result to its GenomeData. Main thread only.
    """
    AA.num_genes = int(status.get("num_genes", 0) or 0)

    if status.get("done"):
        if status.get("was_gzipped"):
            AA.mark_was_gzipped()
        AA.mark_processing_done()
        AA.final_AA_path = status.get("final_AA_path")
        AA.final_nt_path = status.get("final_nt_path")
    else:
        error = status.get("error")
        AA.mark_removed(f"amino-acid-file processing failed ({error})" if error
                        else "amino-acid-file processing failed",
                        GenomeRemovalStage.AMINO_ACID_PREP)


PRODIGAL_TOO_SHORT_EXIT = 10
PRODIGAL_TOO_SHORT_REASON = ("prodigal failed: genome too short for gene-calling")


def run_prodigal(id, run_data, full_inpath = None, group = None):
    """
    Call genes with prodigal and strip the trailing '*' stop characters

    Returns `(done, error)`. `error` is None on success, otherwise a short human-readable
    reason the caller can carry back to the main thread and hand to `mark_removed`, so it
    lands in `reason_removed` in genomes-summary-info.tsv
    """
    allowed_groups = ["ncbi", "fasta", "genbank"]
    if group not in allowed_groups:
        raise ValueError(f"Invalid group: {group}. Must be one of {', '.join(allowed_groups)}")

    if group == "ncbi":
        in_path = f"{run_data.ncbi_processing_dir}/{id}_genomic.fna"
        out_AA_path = f"{run_data.ncbi_processing_dir}/{id}_protein.faa"
        out_nt_path = f"{run_data.ncbi_processing_dir}/{id}_cds.fasta"
    elif group == "genbank":
        in_path = f"{run_data.genbank_processing_dir}/{id}.fasta"
        out_AA_path = f"{run_data.genbank_processing_dir}/{id}_protein.faa"
        out_nt_path = f"{run_data.genbank_processing_dir}/{id}_cds.fasta"
    else:  # fasta
        in_path = full_inpath
        out_AA_path = f"{run_data.fasta_processing_dir}/{id}_protein.faa"
        out_nt_path = f"{run_data.fasta_processing_dir}/{id}_cds.fasta"

    tmp_AA = f"{out_AA_path}.tmp"
    tmp_nt = f"{out_nt_path}.tmp"

    prodigal_cmd = [
        "prodigal", "-c", "-q",
        "-i", in_path,
        "-a", tmp_AA,
        "-d", tmp_nt,
    ]

    done = False
    error = None
    try:
        subprocess.run(prodigal_cmd,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE,
                       check=True)

        _strip_stop_chars(tmp_AA, out_AA_path)
        if run_data.nucleotide_mode:
            _strip_stop_chars(tmp_nt, out_nt_path)
        done = True

    except subprocess.CalledProcessError as e:
        if e.returncode == PRODIGAL_TOO_SHORT_EXIT:
            error = PRODIGAL_TOO_SHORT_REASON
        else:
            stderr = e.stderr.decode("utf-8", errors="replace") if e.stderr else ""
            error = f"prodigal failed (exit {e.returncode})"
            if stderr.strip():
                error += f": {_condense(stderr)}"
    except (OSError, ValueError) as e:
        # binary missing, unreadable input, or a failed strip/replace
        error = f"prodigal could not be run: {type(e).__name__}: {e}"
    finally:
        for tmp in (tmp_AA, tmp_nt):
            remove_file_if_exists(tmp)

    # if the output file doesn't exist or is empty, mark as not done
    if not os.path.exists(out_AA_path):
        done = False
    elif os.path.getsize(out_AA_path) == 0:
        os.remove(out_AA_path)
        done = False

    if not done and error is None:
        error = "prodigal produced no genes"

    return done, error


def _condense(text, limit=200):
    """
    Flatten a captured stderr blob to a single line for a summary-table cell
    """
    flat = " ".join(text.split())
    return flat if len(flat) <= limit else flat[:limit - 3] + "..."


def _strip_stop_chars(src, dest):
    """
    Copy `src` to `dest` with '*' characters removed, atomically
    """
    tmp = f"{dest}.part"
    try:
        with open(src, "r") as in_f, open(tmp, "w") as out_f:
            while True:
                chunk = in_f.read(1024 * 256)
                if not chunk:
                    break
                out_f.write(chunk.replace("*", ""))
            out_f.flush()
            os.fsync(out_f.fileno())
        os.replace(tmp, dest)
    except BaseException:
        remove_file_if_exists(tmp)
        raise
