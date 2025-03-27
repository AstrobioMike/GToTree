import sys
import contextlib
from gtotree.utils.messaging import (report_message,
                                     check_and_report_any_changed_default_behavior,
                                     capture_stdout_to_log)
from gtotree.utils.context import log_file_var
from gtotree.utils.preflight_checks import check_input_genomes_amount


@capture_stdout_to_log(lambda: log_file_var.get())
def display_initial_run_info(args, run_data):

    # time.sleep(1)

    print("\n ----------------------------------- RUN INFO ----------------------------------- \n", flush=True)

    # time.sleep(1)

    report_message("  Input-genome sources include:")

    if args.ncbi_accessions:
        print(f"      - NCBI accessions listed in {args.ncbi_accessions} ({len(run_data.get_input_ncbi_accs())} genomes)", flush=True)
    if args.genbank_files:
        print(f"      - Genbank files listed in {args.genbank_files} ({len(run_data.get_input_genbank_ids())} genomes)", flush=True)
    if args.fasta_files:
        print(f"      - Fasta files listed in {args.fasta_files} ({len(run_data.get_input_fasta_ids())} genomes)", flush=True)
    if args.amino_acid_files:
        print(f"      - Amino-acid files listed in {args.amino_acid_files} ({len(run_data.get_input_amino_acid_ids())} genomes)", flush=True)

    report_message(f"                           Total input genomes: {len(run_data.all_input_genomes)}", "green")
    # time.sleep(1)

    with contextlib.redirect_stdout(sys.__stdout__):
        check_input_genomes_amount(len(run_data.all_input_genomes), args)

    report_message("  Single-copy gene HMM source to be used:")
    print(f"      - {args.hmm} ({len(run_data.initial_SCG_targets)} targets)", flush=True)
    # time.sleep(1)

    check_and_report_any_changed_default_behavior(args)
    # time.sleep(3)
