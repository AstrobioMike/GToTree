from gtotree.utils.messaging import report_message
from gtotree.utils.hmm_handling import get_number_of_targets
from gtotree.utils.preflight_checks import (check_input_genomes_amount,
                                            check_and_report_any_changed_default_behavior)

def display_initial_run_info(args, input_genome_data):

    print("\n ---------------------------------  RUN INFO  ---------------------------------\n")

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

    check_input_genomes_amount(total_input_genomes, args)

    report_message("  HMM source to be used:")
    print(f"      - {args.hmm} ({get_number_of_targets(args.hmm_path)} targets)")

    check_and_report_any_changed_default_behavior(args)
    ## need mechanism for getting number of targets in HMM file here
    ## see what we were doing before...