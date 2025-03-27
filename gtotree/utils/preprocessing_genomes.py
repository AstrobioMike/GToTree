import pandas as pd
import urllib.request
import gzip
import shutil
import os
import subprocess
from gtotree.utils.messaging import (report_early_exit,
                                     report_processing_stage,
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

    run_data = preprocess_ncbi_genomes(args, run_data)
    run_data = preprocess_genbank_genomes(args, run_data)
    run_data = preprocess_fasta_genomes(args, run_data)
    run_data = preprocess_amino_acid_files(args, run_data)

    genome_preprocessing_update(run_data)

    return run_data


def preprocess_ncbi_genomes(args, run_data):
    if args.ncbi_accessions:

        report_processing_stage("ncbi")

        run_data = parse_assembly_summary(NCBI_assembly_summary_tab, run_data)

        num_ncbi_accs_remaining = len(run_data.remaining_ncbi_accs())

        if num_ncbi_accs_remaining > 0:
            # writing run_data to file so it can be accessed by snakemake
            write_run_data(run_data)
            snakefile = get_snakefile_path("preprocess-ncbi-accessions.smk")
            description = "Preprocessing NCBI accessions"

            run_snakemake(snakefile, num_ncbi_accs_remaining, args, run_data, description)

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
    if len(run_data.get_ncbi_accs_not_downloaded()) > 0:
        with open(run_data.run_files_dir + "/ncbi-accessions-not-downloaded.txt", "w") as not_downloaded_file:
            for acc in run_data.get_ncbi_accs_not_downloaded():
                not_downloaded_file.write(acc + "\n")
    return run_data


def preprocess_genbank_genomes(args, run_data):
    if args.genbank_files:

        report_processing_stage("genbank")

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
        with open(run_data.run_files_dir + "/genbank-files-not-parsed.txt", "w") as not_parsed_file:
            for entry in failed_genbank_files_list:
                not_parsed_file.write(entry + "\n")

def preprocess_fasta_genomes(args, run_data):
    if args.fasta_files:

        report_processing_stage("fasta")

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

        report_processing_stage("amino-acid")

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
    report_processing_stage("preprocessing-update")
    run_data.update_all_input_genomes()
    report_genome_preprocessing_update(run_data)


def run_prodigal(id, run_data, full_inpath = None, group = None):
    allowed_groups = ["ncbi", "fasta", "genbank"]
    if group not in allowed_groups:
        raise ValueError(f"Invalid group: {group}. Must be one of {', '.join(allowed_groups)}")

    if group == "ncbi":
        in_path = f"{run_data.ncbi_downloads_dir}/{id}_genomic.fna"
        out_path = f"{run_data.ncbi_downloads_dir}/{id}_protein.faa"
    elif group == "genbank":
        in_path = f"{run_data.genbank_processing_dir}/{id}.fasta"
        out_path = f"{run_data.genbank_processing_dir}/{id}_protein.faa"
        print(f"\n\n    {out_path}\n\n")
    elif group == "fasta":
        in_path = full_inpath
        out_path = f"{run_data.fasta_processing_dir}/{id}_protein.faa"
    else:
        report_early_exit(f"    Prodigal not yet implemented for \"{group}\".")

    prodigal_cmd = [
        "prodigal",
        "-c",
        "-q",
        "-i", f"{in_path}",
        "-a", f"{out_path}.tmp",
    ]
    print(f"\n\n    {prodigal_cmd}\n\n")

    remove_ast_cmd = f"tr -d '*' < {out_path}.tmp > {out_path}"

    try:
        subprocess.run(prodigal_cmd, stdout=subprocess.DEVNULL)
        subprocess.run(remove_ast_cmd, shell=True)
        os.remove(f"{out_path}.tmp")
        done = True
    except:
        done = False

    if os.path.getsize(out_path) == 0:
        os.remove(out_path)
        done = False

    return done
