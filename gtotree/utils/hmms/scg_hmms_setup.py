import os
import sys
from gtotree.utils.misc.messaging import wprint, color_text, report_message, report_early_exit
import pandas as pd # type: ignore
from gtotree.utils.misc.general import download_with_tqdm, SCGset, decode_pyhmmer_text
from gtotree.utils.hmms.hmm_searching_engine import profiles_missing_gathering_cutoffs
import pyhmmer #type: ignore


# aliases that aren't simply the file stem of a packaged set
HMM_SET_ALIASES = {"universal": "Universal-Hug-et-al"}


def with_hmm_suffix(name):
    return name if name.endswith(".hmm") else f"{name}.hmm"


def find_local_hmm_file(hmm_arg):
    """
    Absolute path if `-H` points at a real file on disk, else None
    """
    for candidate in (with_hmm_suffix(hmm_arg), hmm_arg):
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def packaged_scg_set_names():
    """
    Official names of the pre-packaged SCG-sets, as file stems
    """
    try:
        df = read_in_hmm_summary_table()
    except Exception:
        return []
    names = []
    for entry in df["file"]:
        entry = str(entry)
        names.append(entry[:-4] if entry.lower().endswith(".hmm") else entry)
    return names


def canonicalize_hmm_arg(hmm_arg):
    """
    Resolve a requested SCG-set to the official spelling of its name

    e.g., `-H bacteria`, `-H BACTERIA`, and `-H Bacteria.HMM` all -> `Bacteria.hmm`
    Matching case-insensitively is not sufficient on its own: the value is reused as a
    filename under GToTree_HMM_dir, as the lookup key into hmm-sources-and-info.tsv, and
    in run reporting. Left as typed, `-H bacteria` would download a second copy to
    `bacteria.hmm` alongside the existing `Bacteria.hmm` on a case-sensitive filesystem.

    Returned unchanged when it doesn't name a packaged set, so user paths and genuine
    typos both flow on to the existing handling.
    """
    stem = hmm_arg[:-4] if hmm_arg.lower().endswith(".hmm") else hmm_arg
    key = stem.lower()

    if key in HMM_SET_ALIASES:
        return HMM_SET_ALIASES[key]

    for name in packaged_scg_set_names():
        if name.lower() == key:
            return name

    return hmm_arg


def resolve_hmm_arg(args):
    """
    Settle `-H` on the official spelling of the requested set
    """
    if find_local_hmm_file(args.hmm) is None:
        args.hmm = canonicalize_hmm_arg(args.hmm)
    return args


def resolve_hmm_source(args, run_data):
    """
    Point run_data at the HMM file and validate it. Runs on every invocation.

    Cheap on a resume: the packaged set is already in GToTree_HMM_dir, so `get_hmm_path`
    is a stat rather than a download.
    """
    local_path = find_local_hmm_file(args.hmm)

    if local_path is not None:
        run_data.hmm_path = local_path
    else:
        run_data.hmm_path = get_hmm_path(with_hmm_suffix(args.hmm), args.hmm)

    check_gathering_cutoffs(run_data.hmm_path, args.hmm)

    return run_data


def populate_SCG_targets(run_data):
    """
    Build the initial SCG-target list. Fresh runs only, a resume reads its targets,
    with all their accumulated per-SCG state, back out of run-data.json.
    """
    initial_SCG_targets = get_SCG_hmm_targets(run_data.hmm_path)
    run_data.SCG_targets = [SCGset.from_id(target) for target in initial_SCG_targets]

    return run_data


def check_hmm_file(args, run_data):
    """
    The whole fresh-run HMM setup: resolve the name, resolve the file, list the targets
    """
    args = resolve_hmm_arg(args)
    run_data = resolve_hmm_source(args, run_data)

    return populate_SCG_targets(run_data)


def check_gathering_cutoffs(hmm_path, hmm_arg):
    """
    Fail early if any profile lacks a gathering (GA) threshold

    Searching uses gathering thresholds (the `--cut_ga` equivalent), so a profile
    without one can't be searched
    """
    try:
        missing = profiles_missing_gathering_cutoffs(hmm_path)
    except Exception:
        # unreadable HMM is reported elsewhere; don't mask it with a cutoffs message
        return

    if not missing:
        return

    shown = ", ".join(missing[:5])
    if len(missing) > 5:
        shown += f", ... ({len(missing)} total)"

    report_message(
        f'The HMM file for "{hmm_arg}" has {len(missing)} profile(s) with no gathering '
        f"(GA) threshold: {shown}.", "red")
    report_message(
        "GToTree identifies target genes using each profile's gathering threshold, so "
        "every profile needs a `GA` line. If you built this HMM yourself, you can add "
        "gathering thresholds with `hmmbuild --cut_ga` inputs or by editing the `GA` "
        "lines in directly.")
    report_early_exit(None, copy_log=False)


def get_hmm_path(hmm_file, hmm_arg):
    hmm_data_dir = check_hmm_location_var_is_set()
    dest_path = os.path.join(hmm_data_dir, hmm_file)
    if os.path.isfile(dest_path):
        return dest_path
    download_prepackaged_hmm(dest_path, hmm_arg)
    return dest_path


def check_hmm_location_var_is_set():
    try:
        hmm_data_dir = os.environ['GToTree_HMM_dir']
    except KeyError:
        wprint(color_text("The environment variable 'GToTree_HMM_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt data locations check`.")
        sys.exit(1)
    return hmm_data_dir


def download_prepackaged_hmm(dest_path, hmm_arg):
    hmm_file = os.path.basename(dest_path)
    target_hmm_url = get_target_hmm_url(hmm_file, hmm_arg)

    print(color_text(f"    Downloading the prebuilt \"{hmm_arg}\" HMM set (only needs to be done once)...\n", "yellow"))

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)

    try:
        download_with_tqdm(target_hmm_url, f"        {hmm_arg} HMM file", dest_path)
    except Exception as e:
        report_message(f"Downloading the HMM file failed with the following error:\n{e}", "red")
        report_early_exit(None, copy_log=False)


def read_in_hmm_summary_table():
    hmm_data_dir = check_hmm_location_var_is_set()
    hmm_table_path = os.path.join(hmm_data_dir, "hmm-sources-and-info.tsv")
    df = pd.read_csv(hmm_table_path, sep='\t')
    return(df)


def get_target_hmm_url(hmm_file, hmm_arg):
    df = read_in_hmm_summary_table()
    match = df.loc[df["file"].astype(str).str.lower() == hmm_file.lower(), "link"]
    if len(match) == 0:
        report_message(f"You specified \"{hmm_arg}\" as the HMM file to use, but that file can't be found.", "red")
        report_message("You can see the available gene-sets packaged with GToTree by running `gtt hmms`.")
        report_early_exit(None, copy_log = False)
    return match.values[0]


def get_SCG_hmm_targets(hmm_path):
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        hmms = list(hmm_file)
    SCG_targets = [decode_pyhmmer_text(hmm.name) for hmm in hmms]
    return SCG_targets
