from gtotree.utils.messaging import report_message

def display_run_info(args, input_genome_data):

    print(" ---------------------------------  RUN INFO  ---------------------------------\n")

    report_message("  Input-genome sources include:")

    if args.ncbi_accessions:
        print(f"      - NCBI accessions listed in {args.ncbi_accessions} ({input_genome_data.num_ncbi_accessions} genomes)")
    if args.genbank_files:
        print(f"      - Genbank files listed in {args.genbank_files} ({input_genome_data.num_genbank_files} genomes)")
    if args.fasta_files:
        print(f"      - Fasta files listed in {args.fasta_files} ({input_genome_data.num_fasta_files} genomes)")
    if args.amino_acid_files:
        print(f"      - Amino-acid files listed in {args.amino_acid_files} ({input_genome_data.num_amino_acid_files} genomes)")

    total_input_genomes = input_genome_data.num_ncbi_accessions + input_genome_data.num_genbank_files + input_genome_data.num_fasta_files + input_genome_data.num_amino_acid_files
    report_message(f"                           Total input genomes: {total_input_genomes}", "green")

    report_message("  HMM source to be used:")
    if args.hmm == "Universal":
        print(f"      - Universal-Hug-et-al (<num_hmm_targets> targets)")

    ## need mechanism for getting number of targets in HMM file here
    ## see what we were doing before...