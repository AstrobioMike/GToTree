from Bio import SeqIO

def filter_and_rename_fasta(prefix, run_data, max_length = 99999):

    num = 0
    infile = f"{run_data.ncbi_downloads_dir}/{prefix}_protein.faa"
    outfile = f"{run_data.all_input_genome_AA_files_dir}/{prefix}.faa"
    with open(outfile, "w") as outfile, open(infile, "r") as infile:

        for record in SeqIO.parse(infile, "fasta"):

            if len(record.seq) <= max_length:
                num += 1
                header = f">{prefix}_{num}"
                outfile.write(f"{header}\n{record.seq}\n")
