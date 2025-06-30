from gtotree.utils.general import download_and_gunzip, concat_files

def get_additional_pfam_targets(run_data):

    base_link = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/"

    pfam_dict = {} # key: input_id, value: pulled_id

    with open(run_data.target_pfams_file, "r") as pfam_file:
        for line in pfam_file:
            input_id = line.strip()
            pfam_dict[input_id] = None

    found_pfam_targets = []
    failed_pfam_targets = []
    found_pfam_paths = []

    combined_pfam_hmm_path = f"{run_data.pfam_results_dir}/target-pfam-profiles/all-pfam-targets.hmm"
    run_data.all_pfam_targets_hmm_path = combined_pfam_hmm_path

    for input_id in pfam_dict.keys():

        target = input_id.split(".")[0]

        target_url = f"{base_link}{target}?annotation=hmm"
        pfam_hmm_out_path = f"{run_data.pfam_results_dir}/target-pfam-profiles/{target}.hmm"

        success = download_and_gunzip(target_url, pfam_hmm_out_path)

        if not success:
            failed_pfam_targets.append(input_id)
            continue

        found_pfam_targets.append(target)
        found_pfam_paths.append(pfam_hmm_out_path)
        # getting the full pfam id (with version) from the downloaded file
        with open(pfam_hmm_out_path, "r") as pfam_hmm_file:
            for line in pfam_hmm_file:
                if line.startswith("ACC"):
                    pfam_id = line.split()[1]
                    pfam_dict[input_id] = pfam_id
                    break

    run_data.found_pfam_targets = found_pfam_targets
    run_data.failed_pfam_targets = failed_pfam_targets
    run_data.pfam_dict = pfam_dict

    write_out_failed_pfams(run_data)
    write_requested_and_pulled_pfams(run_data, pfam_dict)
    concat_files(found_pfam_paths, combined_pfam_hmm_path)

    return run_data


def write_out_failed_pfams(run_data):
    if len(run_data.failed_pfam_targets) > 0:
        with open(run_data.run_files_dir + "/failed-pfam-targets.txt", "w") as failed_pfams_file:
            for entry in run_data.failed_pfam_targets:
                failed_pfams_file.write(entry + "\n")


def write_requested_and_pulled_pfams(run_data, pfam_dict):
    with open(run_data.pfam_results_dir + "/info/requested-and-pulled-pfams.tsv", "w") as pfam_file:
        pfam_file.write("specified_pfam\tpulled_pfam\n")
        for input_id, pulled_id in pfam_dict.items():
            if pulled_id is None:
                pulled_id = "NA"
            pfam_file.write(f"{input_id}\t{pulled_id}\n")
