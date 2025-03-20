import re
import sys
import contextlib
import time
from gtotree.utils.messaging import report_message
from gtotree.utils.general import log_file_var
from gtotree.utils.hmm_handling import get_number_of_targets
from gtotree.utils.preflight_checks import (check_input_genomes_amount,
                                            check_and_report_any_changed_default_behavior)


class Tee:
    def __init__(self, *files):
        self.files = files

    ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')

    def write(self, s):
        for f in self.files:
            if hasattr(f, "isatty") and f.isatty():
                f.write(s)
            else:
                f.write(self.ansi_escape.sub("", s))

    def flush(self):
        for f in self.files:
            f.flush()

    def isatty(self):
        if self.files and hasattr(self.files[0], "isatty"):
            return self.files[0].isatty()
        return False


def capture_stdout_to_log(log_file):
    def decorator(func):
        def wrapper(*args, **kwargs):
            file_path = log_file() if callable(log_file) else log_file
            with open(file_path, 'a') as f:
                tee = Tee(sys.stdout, f)
                with contextlib.redirect_stdout(tee):
                    return func(*args, **kwargs)
        return wrapper
    return decorator


@capture_stdout_to_log(lambda: log_file_var.get())
def display_initial_run_info(args, genome_data):

    # time.sleep(1)

    print("\n ---------------------------------  RUN INFO  --------------------------------- \n")

    # time.sleep(1)

    report_message("  Input-genome sources include:")

    if args.ncbi_accessions:
        print(f"      - NCBI accessions listed in {args.ncbi_accessions} ({genome_data.num_ncbi_accessions} genomes)")
    if args.genbank_files:
        print(f"      - Genbank files listed in {args.genbank_files} ({genome_data.num_genbank_files} genomes)")
    if args.fasta_files:
        print(f"      - Fasta files listed in {args.fasta_files} ({genome_data.num_fasta_files} genomes)")
    if args.amino_acid_files:
        print(f"      - Amino-acid files listed in {args.amino_acid_files} ({genome_data.num_amino_acid_files} genomes)")

    report_message(f"                           Total input genomes: {genome_data.num_input_genomes}", "green")
    # time.sleep(1)

    check_input_genomes_amount(genome_data.num_input_genomes, args)

    report_message("  HMM source to be used:")
    print(f"      - {args.hmm} ({get_number_of_targets(args.hmm_path)} targets)")
    # time.sleep(1)

    check_and_report_any_changed_default_behavior(args)
    # time.sleep(1)
