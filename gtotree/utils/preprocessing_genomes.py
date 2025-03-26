import pandas as pd
import urllib.request
import gzip
import shutil
import os
from gtotree.utils.messaging import (report_processing_stage,
                                     report_ncbi_update,
                                     report_genbank_update,
                                     report_fasta_update,
                                     report_AA_update,
                                     report_genome_preprocessing_update)
from gtotree.utils.ncbi.parse_assembly_summary_file import parse_assembly_summary
from gtotree.utils.ncbi.get_ncbi_assembly_tables import NCBI_assembly_summary_tab
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)


def preprocess_genomes(args, run_data):

    print(run_data)

    run_data = preprocess_ncbi_genomes(args, run_data)
    run_data = preprocess_genbank_genomes(args, run_data)
    run_data = preprocess_fasta_genomes(args, run_data)
    run_data = preprocess_amino_acid_files(args, run_data)

    genome_preprocessing_update(run_data)

    print(run_data)

    return args, run_data


def preprocess_ncbi_genomes(args, run_data):
    if args.ncbi_accessions:

        report_processing_stage("ncbi")

        run_data = parse_assembly_summary(NCBI_assembly_summary_tab, run_data)

        ncbi_accs_remaining = len(run_data.remaining_ncbi_accs())

        if ncbi_accs_remaining > 0:
            # writing run_data and args objects to files so they can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-ncbi-accessions.smk")
            description = "Preprocessing NCBI accessions"

            cmd = [
                "snakemake",
                "--snakefile", snakefile,
                "--cores", f"{args.num_jobs}",
                "--default-resources", f"tmpdir='{args.tmp_dir}'",
                "--config",
                f"run_data_path={run_data.run_data_path}"
            ]

            run_snakemake(cmd, ncbi_accs_remaining, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            run_data = capture_ncbi_failed_downloads(run_data)

        report_ncbi_update(run_data)

    return run_data


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


            if acc == "GCF_000153765.1":
                nucleotide_link = "https://wrong"


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
    if len(run_data.get_ncbi_accs_not_downloaded()) > 0:
        with open(run_data.run_files_dir + "/ncbi-accessions-not-downloaded.txt", "w") as not_downloaded_file:
            for acc in run_data.get_ncbi_accs_not_downloaded():
                not_downloaded_file.write(acc + "\n")
    return run_data


def preprocess_genbank_genomes(args, run_data):
    if args.genbank_files:

        report_processing_stage("genbank")

        if run_data.any_incomplete_genbank_files():
            # writing run_data and args objects to files so they can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-genbank-files.smk")
            description = "Preprocessing genbank files"

            cmd = [
                "snakemake",
                "--snakefile", snakefile,
                "--cores", f"{args.num_jobs}",
                "--default-resources", f"tmpdir='{args.tmp_dir}'",
                "--config",
                f"run_data_path={run_data.run_data_path}"
            ]

            run_snakemake(cmd, run_data.num_incomplete_genbank_files, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_failed_genbank_files(run_data)

        report_genbank_update(run_data)

    return run_data


def capture_failed_genbank_files(run_data):
    failed_genbank_files_list = run_data.failed_genbank_files()
    if len(failed_genbank_files_list) > 0:
        with open(run_data.run_files_dir + "/genbank-files-not-parsed.txt", "w") as not_parsed_file:
            for entry in failed_genbank_files_list:
                not_parsed_file.write(entry + "\n")
                run_data.remove_genbank_file(entry)

def preprocess_fasta_genomes(args, run_data):
    if args.fasta_files:

        report_processing_stage("fasta")

        if run_data.any_incomplete_fasta_files():
            # writing run_data and args objects to files so they can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-fasta-files.smk")
            description = "Preprocessing fasta files"

            cmd = [
                "snakemake",
                "--snakefile", snakefile,
                "--cores", f"{args.num_jobs}",
                "--default-resources", f"tmpdir='{args.tmp_dir}'",
                "--config",
                f"run_data_path={run_data.run_data_path}"
            ]

            run_snakemake(cmd, run_data.num_incomplete_fasta_files, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_failed_fasta_files(run_data)

        report_fasta_update(run_data)

    return run_data


def capture_failed_fasta_files(run_data):
    failed_fasta_files_list = run_data.failed_fasta_files()
    if len(failed_fasta_files_list) > 0:
        with open(run_data.run_files_dir + "/failed-fasta-files.txt", "w") as failed_fastas_file:
            for entry in failed_fasta_files_list:
                failed_fastas_file.write(entry + "\n")
                run_data.remove_fasta_file(entry)


def preprocess_amino_acid_files(args, run_data):
    if args.amino_acid_files:

        report_processing_stage("amino-acid")

        if run_data.any_incomplete_amino_acid_files():
            # writing run_data and args objects to files so they can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-amino-acid-files.smk")
            description = "Preprocessing amino-acid files"

            cmd = [
                "snakemake",
                "--snakefile", snakefile,
                "--cores", f"{args.num_jobs}",
                "--default-resources", f"tmpdir='{args.tmp_dir}'",
                "--config",
                f"run_data_path={run_data.run_data_path}"
            ]

            run_snakemake(cmd, run_data.num_incomplete_amino_acid_files, run_data, description)

            run_data = read_run_data(run_data.run_data_path)
            capture_failed_amino_acid_files(run_data)

        report_AA_update(run_data)

    return run_data


def capture_failed_amino_acid_files(run_data):
    failed_amino_acid_files_list = run_data.failed_amino_acid_files()
    if len(failed_amino_acid_files_list) > 0:
        with open(run_data.run_files_dir + "/failed-amino-acid-files.txt", "w") as failed_amino_acids_file:
            for entry in failed_amino_acid_files_list:
                failed_amino_acids_file.write(entry + "\n")
                run_data.remove_amino_acid_file(entry)


def genome_preprocessing_update(run_data):
    report_processing_stage("preprocessing-update")
    run_data.update_all_input_genomes()
    report_genome_preprocessing_update(run_data)

