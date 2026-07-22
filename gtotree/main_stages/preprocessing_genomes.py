import pandas as pd # type: ignore
import urllib.request
import urllib.error
import gzip
import shutil
import os
import time
import random
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm # type: ignore
from gtotree.utils.messaging import (report_early_exit,
                                     report_processing_stage,
                                     report_ncbi_update,
                                     report_genbank_update,
                                     report_fasta_update,
                                     report_AA_update,
                                     report_genome_preprocessing_update)
from gtotree.utils.ncbi.parse_ncbi_assembly_summary import parse_assembly_summary
from gtotree.utils.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_summary_tab
from gtotree.utils.seqs import filter_and_rename_fasta
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)

# capping concurrent NCBI download threads regardless of --num-jobs
# client to NCBI (mirrors bit's max_threads ceiling).
MAX_NCBI_DOWNLOAD_THREADS = 20
NCBI_DOWNLOAD_MAX_RETRIES = 5


def preprocess_genomes(args, run_data):

    run_data = preprocess_ncbi_genomes(args, run_data)
    run_data = preprocess_genbank_genomes(args, run_data)
    run_data = preprocess_fasta_genomes(args, run_data)
    run_data = preprocess_amino_acid_files(args, run_data)

    genome_preprocessing_update(run_data)

    return run_data


def preprocess_ncbi_genomes(args, run_data):

    if run_data.ncbi_accs:

        report_processing_stage("ncbi", run_data)

        run_data = parse_assembly_summary(get_ncbi_assembly_summary_tab(), run_data)

        accs_to_process = run_data.get_ncbi_accs_for_snakemake_preprocessing()

        if len(accs_to_process) > 0:
            run_data = download_and_preprocess_ncbi_accessions(accs_to_process, args, run_data)
            run_data = capture_ncbi_failed_downloads(run_data)

        report_ncbi_update(run_data)

    return run_data


def download_and_preprocess_ncbi_accessions(accs_to_process, args, run_data):

    num_workers = max(1, min(args.num_jobs, MAX_NCBI_DOWNLOAD_THREADS))

    base_link_map = build_base_link_map(run_data)

    bar_format = (
        "      {percentage:3.0f}%|{bar}| "
        "{n_fmt}/{total_fmt} "
        "[time elapsed: {elapsed} | est. remaining: {remaining}]"
    )

    print("")
    with ThreadPoolExecutor(max_workers=num_workers) as pool, \
         tqdm(total=len(accs_to_process), bar_format=bar_format, ncols=76) as pbar:

        futures = {
            pool.submit(_process_one_ncbi_accession, acc_gd, run_data, base_link_map): acc_gd
            for acc_gd in accs_to_process
        }

        for future in as_completed(futures):
            acc_gd = futures[future]
            status = future.result()
            _apply_ncbi_accession_status(acc_gd, status, run_data)
            pbar.update(1)

    return run_data


def _process_one_ncbi_accession(acc_gd, run_data, base_link_map=None):
    """
    Worker: run one accession's download + (optional) prodigal + filter/rename.
    Returns a status dict of the same shape the old Snakefile carried through its
    per-accession JSON, so the bookkeeping step is unchanged. Runs in a worker
    thread, so it must not mutate shared run_data state -- only local files and
    the returned dict. All GenomeData mutation happens back on the main thread in
    _apply_ncbi_accession_status.
    """
    try:
        done, nt = prepare_accession(acc_gd.id, run_data, base_link_map=base_link_map)
        downloaded = bool(done)

        if done and nt:
            done = run_prodigal(acc_gd.id, run_data, group="ncbi")
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
        }
    except BaseException:
        # a worker must never take down the whole pool; treat any unexpected
        # failure as a non-downloaded processing failure and let the bookkeeping
        # step mark it removed
        return {
            "done": False,
            "downloaded": False,
            "prodigal_used": False,
            "final_AA_path": None,
            "final_nt_path": None,
            "num_genes": 0,
        }


def _apply_ncbi_accession_status(acc_gd, status, run_data):
    """
    Apply one worker's result to its GenomeData. Mirrors the reduce logic in the
    old preprocess-ncbi-accessions.smk `rule all`. Called on the main thread only.
    """
    done = bool(status.get("done"))
    downloaded = bool(status.get("downloaded", False))
    prodigal_used = bool(status.get("prodigal_used", False))

    acc_gd.num_genes = int(status.get("num_genes", 0) or 0)

    if done:
        acc_gd.mark_preprocessing_done()
        acc_gd.final_AA_path = status.get("final_AA_path")
        acc_gd.final_nt_path = status.get("final_nt_path")
        acc_gd.acc_was_downloaded = downloaded
    else:
        if downloaded:
            acc_gd.acc_was_downloaded = True
            acc_gd.mark_removed("acc processing failed after download")
        else:
            acc_gd.acc_was_downloaded = False
            acc_gd.mark_removed("acc download failed")

    acc_gd.prodigal_used = prodigal_used
    if prodigal_used:
        run_data.tools_used.prodigal_used = True


def prepare_accession(acc, run_data, base_link_map=None):
    base_link, acc_assembly_str = get_base_link(acc, run_data, base_link_map=base_link_map)

    # an unresolvable download directory (no ftp_path and nothing to rebuild
    # from) comes through as "na" -> there's nothing to download, fail cleanly
    # rather than building a bogus URL.
    if not base_link or base_link.lower() == "na":
        return False, False

    # first trying amino acids
    try:
        # going directly to nucleotides if running in nucleotide mode
        # (so we can call our own genes and CDS and protein match up as expected)
        if run_data.nucleotide_mode:
            raise Exception
        amino_acid_link = base_link + "/" + acc_assembly_str + "_protein.faa.gz"
        amino_acid_filepath = run_data.ncbi_processing_dir + "/" + acc + "_protein.faa"
        download_and_unzip_accession(amino_acid_link, amino_acid_filepath)
        done = True
        nt = False
    except:
        # then trying nucleotides
        try:
            nucleotide_link = base_link + "/" + acc_assembly_str + "_genomic.fna.gz"
            nucleotide_file = run_data.ncbi_processing_dir + "/" + acc + "_genomic.fna"
            download_and_unzip_accession(nucleotide_link, nucleotide_file)
            done = True
            nt = True
        except:
            done = False
            nt = False

    return done, nt


def download_and_unzip_accession(link, filepath, max_retries=NCBI_DOWNLOAD_MAX_RETRIES):
    # Atomic write: fetch to a temp gzip and unzip into a `.part` temp alongside
    # the destination, then os.replace() into place only on a fully successful
    # unzip. An interrupted download or a truncated/corrupt gzip (process kill,
    # `-R` interrupted mid-run, disk full) therefore never leaves a truncated
    # file at `filepath` that a later os.path.isfile() / resume check would
    # mistake for a complete genome. Both temps are cleaned up on any failure.
    # The function still raises on final failure, which prepare_accession relies
    # on to fall back from amino-acid to nucleotide downloads.
    #
    # Transient failures (network blips, truncated/corrupt gzip, transient HTTP
    # 429/5xx) are retried with exponential backoff + jitter. A 404 is permanent
    # (the requested format doesn't exist for this accession) so we fail fast
    # rather than burning retries -- prepare_accession will then try the other format.
    tmp_gzip = filepath + ".gz"
    tmp_out = filepath + ".part"

    def cleanup_partials():
        for tmp in (tmp_gzip, tmp_out):
            try:
                os.remove(tmp)
            except OSError:
                pass

    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            urllib.request.urlretrieve(link, tmp_gzip)
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
            cleanup_partials()
            # 404 -> requested format not available; permanent, don't retry
            if e.code == 404:
                raise
            last_err = e
        except (urllib.error.URLError, OSError, EOFError, gzip.BadGzipFile) as e:
            cleanup_partials()
            last_err = e

        if attempt < max_retries:
            _sleep_backoff(attempt)

    # exhausted retries -> raise so prepare_accession falls back / marks failure
    raise last_err if last_err is not None else RuntimeError(f"failed to download {link}")


def _sleep_backoff(attempt):
    # exponential backoff with jitter (like i do in bit)
    time.sleep((2 ** (attempt - 1)) + random.uniform(0, 1))


def build_base_link_map(run_data):
    df = pd.read_csv(run_data.tmp_dir + "/ncbi-accessions-info.tsv", sep="\t",
                     usecols=["input_accession", "http_base_link"])
    return dict(zip(df["input_accession"], df["http_base_link"]))


def _normalize_base_link(raw_base_link):
    base_link = str(raw_base_link).replace(" ", "_")
    base_link = base_link.rstrip("/")
    acc_assembly_str = base_link.split("/")[-1]
    return base_link, acc_assembly_str


def get_base_link(acc, run_data, base_link_map=None):
    raw_base_link = base_link_map[acc]
    return _normalize_base_link(raw_base_link)


def capture_ncbi_failed_downloads(run_data):
    if len(run_data.get_ncbi_accs_not_downloaded()) > 0:
        with open(run_data.run_files_dir + "/ncbi-accessions-not-downloaded.txt", "w") as not_downloaded_file:
            for acc in run_data.get_ncbi_accs_not_downloaded():
                not_downloaded_file.write(acc + "\n")
    return run_data


def preprocess_genbank_genomes(args, run_data):
    if args.genbank_files:

        report_processing_stage("genbank", run_data)

        num_genbank_files_remaining = len([gd for gd in run_data.genbank_files if not gd.preprocessing_done and not gd.removed])

        if num_genbank_files_remaining > 0:
            # writing run_data to file so it can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-genbank-files.smk")
            description = "Preprocessing genbank files"

            run_snakemake(snakefile,
                          num_genbank_files_remaining,
                          args, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_failed_genbank_files(run_data)

        report_genbank_update(run_data)

    return run_data


def capture_failed_genbank_files(run_data):
    failed_genbank_files_list = run_data.get_failed_genbank_paths()
    if len(failed_genbank_files_list) > 0:
        with open(run_data.run_files_dir + "/failed-genbank-files.txt", "w") as not_parsed_file:
            for entry in failed_genbank_files_list:
                not_parsed_file.write(entry + "\n")


def preprocess_fasta_genomes(args, run_data):
    if args.fasta_files:

        report_processing_stage("fasta", run_data)

        num_fasta_files_remaining = len([fd for fd in run_data.fasta_files if not fd.preprocessing_done and not fd.removed])

        if num_fasta_files_remaining > 0:
            # writing run_data to file so it can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-fasta-files.smk")
            description = "Preprocessing fasta files"

            run_snakemake(snakefile,
                          num_fasta_files_remaining,
                          args, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_failed_fasta_files(run_data)

        report_fasta_update(run_data)

    return run_data


def capture_failed_fasta_files(run_data):
    failed_fasta_files_list = run_data.get_failed_fasta_paths()
    if len(failed_fasta_files_list) > 0:
        with open(run_data.run_files_dir + "/failed-fasta-files.txt", "w") as failed_fastas_file:
            for entry in failed_fasta_files_list:
                failed_fastas_file.write(entry + "\n")


def preprocess_amino_acid_files(args, run_data):
    if args.amino_acid_files:

        report_processing_stage("amino-acid", run_data)

        num_AA_files_remaining = len([fd for fd in run_data.amino_acid_files if not fd.preprocessing_done and not fd.removed])

        if num_AA_files_remaining > 0:
            # writing run_data to file so it can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-amino-acid-files.smk")
            description = "Preprocessing amino-acid files"

            run_snakemake(snakefile,
                          num_AA_files_remaining,
                          args, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_failed_amino_acid_files(run_data)

        report_AA_update(run_data)

    return run_data


def capture_failed_amino_acid_files(run_data):
    failed_amino_acid_files_list = run_data.get_failed_amino_acid_paths()
    if len(failed_amino_acid_files_list) > 0:
        with open(run_data.run_files_dir + "/failed-amino-acid-files.txt", "w") as failed_amino_acids_file:
            for entry in failed_amino_acid_files_list:
                failed_amino_acids_file.write(entry + "\n")


def genome_preprocessing_update(run_data):
    report_processing_stage("preprocessing-update", run_data)
    run_data.update_all_input_genomes()
    report_genome_preprocessing_update(run_data)


def run_prodigal(id, run_data, full_inpath = None, group = None):
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
    elif group == "fasta":
        in_path = full_inpath
        out_AA_path = f"{run_data.fasta_processing_dir}/{id}_protein.faa"
        out_nt_path = f"{run_data.fasta_processing_dir}/{id}_cds.fasta"
    else:
        report_early_exit(run_data, message = f"    Prodigal not yet implemented for \"{group}\".")

    prodigal_cmd = [
        "prodigal",
        "-c",
        "-q",
        "-i", f"{in_path}",
        "-a", f"{out_AA_path}.tmp",
        "-d", f"{out_nt_path}.tmp"
    ]

    remove_AA_ast_cmd = f"tr -d '*' < {out_AA_path}.tmp > {out_AA_path}"
    remove_nt_ast_cmd = f"tr -d '*' < {out_nt_path}.tmp > {out_nt_path}"

    try:
        subprocess.run(prodigal_cmd, stdout=subprocess.DEVNULL)
        subprocess.run(remove_AA_ast_cmd, shell=True)
        if run_data.nucleotide_mode:
            subprocess.run(remove_nt_ast_cmd, shell=True)
        os.remove(f"{out_AA_path}.tmp")
        os.remove(f"{out_nt_path}.tmp")
        done = True
    except:
        done = False

    # be defensive: if the output file doesn't exist or is empty, mark as not done
    if not os.path.exists(out_AA_path):
        done = False
    else:
        if os.path.getsize(out_AA_path) == 0:
            os.remove(out_AA_path)
            done = False

    return done
