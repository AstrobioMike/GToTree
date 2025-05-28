import os
import pandas as pd
import shutil

def parse_kofamscan_targets(run_data):

    target_KOs_tsv = run_data.ko_results_dir + "/target-kos.tsv"
    run_data.target_kos_tsv = target_KOs_tsv
    target_KO_profiles_dir = run_data.ko_results_dir + "/target-ko-profiles/"
    run_data.target_ko_profiles_dir = target_KO_profiles_dir

    KO_data_dir = os.environ["KO_data_dir"]
    full_KO_list_tsv = KO_data_dir + "/ko_list"
    full_KO_HMMs_dir = KO_data_dir + "/profiles"

    wanted_KOs = get_wanted_KOs(run_data)
    found_KOs, missing_KOs = get_target_KOs_tab(full_KO_list_tsv, wanted_KOs, target_KOs_tsv)
    run_data.wanted_ko_targets = wanted_KOs
    run_data.found_ko_targets = found_KOs
    run_data.failed_ko_targets = missing_KOs

    copy_over_target_ko_HMMs(found_KOs, full_KO_HMMs_dir, target_KO_profiles_dir)

    return(run_data)


def get_wanted_KOs(run_data):
    wanted_KOs = []
    with open(run_data.target_kos_file, "r") as target_KOs_file:
        for line in target_KOs_file:
            KO = line.strip()
            wanted_KOs.append(KO)
    return wanted_KOs


def get_target_KOs_tab(full_KO_list_tsv, wanted_KOs, target_KOs_tsv):
    full_KO_tab = pd.read_csv(full_KO_list_tsv, sep="\t", header=0)
    target_KOs_tab = full_KO_tab[full_KO_tab["knum"].isin(wanted_KOs)]
    found_KOs = target_KOs_tab["knum"].tolist()
    missing_KOs = list(set(wanted_KOs) - set(full_KO_tab["knum"].tolist()))

    target_KOs_tab.to_csv(target_KOs_tsv, sep="\t", index=False)

    return found_KOs, missing_KOs


def copy_over_target_ko_HMMs(found_KOs, full_KO_HMMs_dir, target_KO_profiles_dir):
    for KO in found_KOs:
        source_file = full_KO_HMMs_dir + "/" + KO + ".hmm"
        dest_file = target_KO_profiles_dir + "/" + KO + ".hmm"
        if os.path.exists(source_file):
            shutil.copyfile(source_file, dest_file)
        else:
            print(f"KO HMM file not found for this KO (this shouldn't really happen...): {source_file}")

