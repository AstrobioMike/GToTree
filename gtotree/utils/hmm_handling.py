import os
import sys
from gtotree.utils.messaging import wprint, color_text, report_message, report_early_exit
import pandas as pd
from gtotree.utils.general import download_with_tqdm
import pyhmmer #type: ignore


def check_hmm_file(args):
    if args.hmm == "Universal":
        hmm_arg = "Universal-Hug-et-al"
        args.hmm = "Universal-Hug-et-al"
    else:
        hmm_arg = args.hmm

    # adding .hmm to end if not present
    if not hmm_arg.endswith(".hmm"):
        hmm_file = f"{hmm_arg}.hmm"
    else:
        hmm_file = hmm_arg

    # getting hmm path
    if os.path.isfile(hmm_file): # handles if user-provided full path
        args.hmm_path = hmm_file
    else:
        # getting hmm from prepackaged table
        args.hmm_path = get_hmm_path(hmm_file, args.hmm)

    return args


def get_hmm_path(hmm_file, hmm_arg):
    hmm_data_dir = check_hmm_location_var_is_set()
    if os.path.isfile(os.path.join(hmm_data_dir, hmm_file)):
        return os.path.join(hmm_data_dir, hmm_file)
    else:
        download_prepackaged_hmm(hmm_file, hmm_arg)
        return os.path.join(hmm_data_dir, hmm_file)


def check_hmm_location_var_is_set():
    try:
        hmm_data_dir = os.environ['GToTree_HMM_dir']
    except KeyError:
        wprint(color_text("The environment variable 'GToTree_HMM_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt-data-locations check`.")
        sys.exit(1)
    return hmm_data_dir


def download_prepackaged_hmm(hmm_file, hmm_arg):
    target_hmm_url = get_target_hmm_url(hmm_file, hmm_arg)

    print(color_text(f"    Downloading the prebuilt \"{hmm_arg}\" HMM set (only needs to be done once)...\n", "yellow"))

    try:
        download_with_tqdm(target_hmm_url, f"        {hmm_arg} HMM file", hmm_file)
    except Exception as e:
        report_message(f"Downloading the HMM file failed with the following error:\n{e}", "red")
        report_early_exit()


def read_in_hmm_summary_table():
    hmm_data_dir = check_hmm_location_var_is_set()
    hmm_table_path = os.path.join(hmm_data_dir, "hmm-sources-and-info.tsv")
    df = pd.read_csv(hmm_table_path, sep='\t')
    return(df)


def get_target_hmm_url(hmm_file, hmm_arg):
    df = read_in_hmm_summary_table()
    try:
        target_hmm_url = df.loc[df['file'] == hmm_file, 'link'].values[0]
    except IndexError:
        report_message(f"You specified \"{hmm_arg}\" as the HMM file to use, but that file can't be found.", "red")
        report_message("You can see the available gene-sets packaged with GToTree by running `gtt-hmms`.")
        report_early_exit()
    return target_hmm_url



def get_number_of_targets(hmm_path):
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        number_of_targets = len(list(hmm_file))
    return number_of_targets
