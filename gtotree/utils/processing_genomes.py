import pandas as pd
import urllib.request
import gzip
import shutil
import os
import subprocess
from gtotree.utils.messaging import (report_ncbi_update,
                                     report_processing_stage)
from gtotree.utils.ncbi.parse_assembly_summary_file import parse_assembly_summary
from gtotree.utils.ncbi.get_ncbi_assembly_tables import NCBI_assembly_summary_tab
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   write_args,
                                   run_snakemake)

def process_genomes(args, run_data):
    process_ncbi_genomes(args, run_data)


def process_ncbi_genomes(args, run_data):
    if args.ncbi_accessions:

        report_processing_stage("ncbi")

        run_data = parse_assembly_summary(NCBI_assembly_summary_tab, run_data)

        if set(run_data.ncbi_accessions) != set(run_data.ncbi_accessions_done):
            # writing run_data and args objects to files so they can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("process-ncbi-accessions.smk")
            description = "Processing NCBI accessions"

            cmd = [
                "snakemake",
                "--snakefile", snakefile,
                "--cores", f"{args.num_jobs}",
                "--default-resources", f"tmpdir='{args.tmp_dir}'",
                "--config",
                f"run_data_path={run_data.run_data_path}"
            ]

            run_snakemake(cmd, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_ncbi_failed_downloads(run_data)

        report_ncbi_update(run_data)


def prepare_accession(acc, run_data):
    base_link, acc_assembly_str = get_base_link(acc, run_data)

    # first trying amino acids
    try:
        amino_acid_link = base_link + acc_assembly_str + "_protein.faa.gz"
        amino_acid_filepath = run_data.ncbi_downloads_dir + "/" + acc + "_protein.faa"
        download_and_unzip_accession(amino_acid_link, amino_acid_filepath)
        done = True
        nt = False
    except:
        # then trying nucleotides
        try:
            nucleotide_link = base_link + acc_assembly_str + "_genomic.fna.gz"
            nucleotide_file = run_data.ncbi_downloads_dir + "/" + acc + "_genomic.fna"
            download_and_unzip_accession(nucleotide_link, nucleotide_file)
            done = True
            nt = True
        except:
            done = False
            nt = False

    return done, nt


def download_and_unzip_accession(link, filepath):
    tmp_gzip = filepath + ".gz"
    urllib.request.urlretrieve(link, tmp_gzip)
    with gzip.open(tmp_gzip, 'rb') as f_in, open(filepath, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_gzip)


def get_base_link(acc, run_data):
    df = pd.read_csv(run_data.tmp_dir + "/ncbi-accessions-info.tsv", sep="\t",
                     usecols=["input_accession", "http_base_link"])
    base_link = df.loc[df['input_accession'] == acc, 'http_base_link'].values[0]
    acc_assembly_str = base_link.split("/")[-2]
    return base_link, acc_assembly_str


def capture_ncbi_failed_downloads(run_data):
    if len(run_data.ncbi_accs_not_downloaded) > 0:
        with open(run_data.run_files_dir + "/ncbi-accessions-not-downloaded.txt", "w") as not_downloaded_file:
            for acc in run_data.ncbi_accs_not_downloaded:
                not_downloaded_file.write(acc + "\n")
