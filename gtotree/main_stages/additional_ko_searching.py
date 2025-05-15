import subprocess
from gtotree.utils.ko.ko_handling import parse_kofamscan_targets
from gtotree.utils.messaging import (report_processing_stage,
                                     report_ko_searching_update)
from gtotree.utils.general import (write_run_data,
                                   read_run_data,
                                   get_snakefile_path,
                                   run_snakemake)


def search_kos(args, run_data):

    report_processing_stage("additional-ko-searching", run_data)

    if run_data.additional_ko_searching_done:
        report_ko_searching_update(run_data)
        return run_data

    # generating the subset of target KOs
    run_data = parse_kofamscan_targets(run_data)

    num_genomes_to_search = len(run_data.get_all_input_genomes_for_hmm_search())

    if num_genomes_to_search > 0:
        # writing run_data to file so it can be accessed by snakemake
        write_run_data(run_data)
        snakefile = get_snakefile_path("search-kos.smk")
        description = "Searching KOs"

        run_snakemake(snakefile, num_genomes_to_search, args, run_data, description)

        run_data = read_run_data(run_data.run_data_path)


    write_out_failed_ko_targets(run_data)

    report_ko_searching_update(run_data)

    return run_data


def run_ko_search(assembly_id, profiles_dir, ko_file, base_outpath, AA_file):

    outpath = f"{base_outpath}/{assembly_id}.kofamscan.tsv"
    tmp_path = f"{base_outpath}/{assembly_id}-tmp"
    cmd = [
        "exec_annotation",
        "-p", profiles_dir,
        "-k", ko_file,
        "--cpu", "1",
        "-f", "mapper",
        "--no-report-unannotated",
        "--tmp-dir", tmp_path,
        "-o", outpath,
        AA_file
    ]

    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
        kofamscan_failed = False
    except:
        kofamscan_failed = True

    return kofamscan_failed


def write_out_failed_ko_targets(run_data):
    if len(run_data.failed_ko_targets) > 0:
        with open(run_data.run_files_dir + "/failed-ko-targets.txt", "w") as fail_file:
            for KO in run_data.failed_ko_targets:
                fail_file.write(KO + "\n")
