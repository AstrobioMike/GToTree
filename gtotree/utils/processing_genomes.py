import pandas as pd
import urllib.request
import gzip
import shutil
import os
import subprocess
from gtotree.utils.messaging import (report_message,
                                     report_notice,
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
            run_data_path = write_run_data(run_data)
            snakefile = get_snakefile_path("process-ncbi-accessions.smk")
            description = "Processing NCBI accessions"

            cmd = [
                "snakemake",
                "--snakefile", snakefile,
                "--cores", f"{args.num_jobs}",
                "--default-resources", f"tmpdir='{args.tmp_dir}'",
                "--config",
                f"run_data_path={run_data_path}"
            ]

            run_snakemake(cmd, run_data, description)

            # reading in updated run_data object
            run_data = read_run_data(run_data_path)

        print(run_data)

##### MIGHT WANT TO HANDLE ALL CONVERSIONS IN HERE TOO (E.G., GTT-RENAME-FASTA, GTT-FILTER-SEQS-BY-LENGTH,
##### PRODIGAL WHEN NEEDED, ETC.) AND JUST GET TO FINAL AMINO-ACID FILES
        # if not run_data.ncbi_hmm_searches_done:
            # scan_genome() # this should be re-usable for all genomes, not just NCBI ones
            ## maybe i should get all genomes to amino-acid files first, then do the hmm search
            ## this way, i can do the hmm searches on all at once with one snakemake call

def prepare_accession(acc, run_data):
    base_link, acc_assembly_str = get_base_link(acc, run_data)

    # first trying amino acids
    try:
        amino_acid_link = base_link + acc_assembly_str + "_protein.faa.gz"
        amino_acid_filepath = run_data.ncbi_downloads_dir + "/" + acc + "_protein.faa"
        if acc == "GCF_000153765.1":
            amino_acid_link = "https://wrong"
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
