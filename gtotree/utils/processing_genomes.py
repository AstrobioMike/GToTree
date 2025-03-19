import pandas as pd
import urllib.request
import gzip
import shutil
import os
from gtotree.utils.messaging import report_message, report_notice
from gtotree.utils.ncbi.parse_assembly_summary_file import parse_assembly_summary
from gtotree.utils.ncbi.get_ncbi_assembly_tables import NCBI_assembly_summary_tab

def process_genomes(args, genome_data):
    process_ncbi_genomes(args, genome_data)


def process_ncbi_genomes(args, genome_data):
    if args.ncbi_accessions:

        report_message("\n ##############################################################################"
                        " ####          Working on the genomes provided as NCBI accessions          ####"
                        " ##############################################################################"
                    , color = None)
        ncbi_accs_not_found = parse_assembly_summary(NCBI_assembly_summary_tab, genome_data.ncbi_accessions, args)
        if ncbi_accs_not_found:
            report_notice(f"    {len(ncbi_accs_not_found)} accession(s) not successfully found at NCBI.\n"
                          f"    Reported in {args.run_files_dir}/ncbi-accessions-not-found.txt")

            for acc in ncbi_accs_not_found:
                genome_data.remove_ncbi_accession(acc)

        for acc in genome_data.ncbi_accessions:

            downloaded = prepare_accession(acc, args, genome_data)
            if not downloaded:
                genome_data.remove_ncbi_accession(acc)
                continue

            # scan_genome()

##### MIGHT WANT TO HANDLE ALL CONVERSIONS IN HERE TOO (E.G., GTT-RENAME-FASTA, GTT-FILTER-SEQS-BY-LENGTH,
##### PRODIGAL WHEN NEEDED, ETC.) AND JUST GET TO FINAL AMINO-ACID FILES
def prepare_accession(acc, args, genome_data):
    base_link, acc_assembly_str = get_base_link(acc, args)

    # trying amino acids first
    try:
        amino_acid_link = base_link + acc_assembly_str + "_protein.faa.gz"
        amino_acid_filepath = args.ncbi_downloads_dir + "/" + acc + "_protein.faa"

        download_and_unzip_accession(amino_acid_link, amino_acid_filepath)
        downloaded = True
    except:
        # then trying nucleotides
        try:
            nucleotide_link = base_link + acc_assembly_str + "_genomic.fna.gz"
            nucleotide_file = args.ncbi_downloads_dir + "/" + acc + "_genomic.fna"

            download_and_unzip_accession(nucleotide_link, nucleotide_file)
            downloaded = True
        except:
            downloaded = False

    return downloaded


def download_and_unzip_accession(link, filepath):
    tmp_gzip = filepath + ".gz"
    urllib.request.urlretrieve(link, tmp_gzip)
    with gzip.open(tmp_gzip, 'rb') as f_in, open(filepath, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_gzip)





def get_base_link(acc, args):
    df = pd.read_csv(args.tmp_dir + "/ncbi-accessions-info.tsv", sep="\t",
                     usecols=["input_accession", "http_base_link"])
    base_link = df.loc[df['input_accession'] == acc, 'http_base_link'].values[0]
    acc_assembly_str = base_link.split("/")[-2]
    return base_link, acc_assembly_str
