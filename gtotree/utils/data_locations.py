#!/usr/bin/env python

import sys
import os
import argparse
import shutil
from gtotree.utils.messaging import report_message, color_text, wprint


ENV_VARIABLES = [
    "NCBI_ASSEMBLY_DATA_DIR",
    "GToTree_HMM_dir",
    "TAXONKIT_DB",
    "GTDB_DIR",
    "KO_data_dir",
    "Pfam_data_dir",
]


class PathNotWritable(Exception):
    pass

class PathNotAbsolute(Exception):
    pass


def main():
    args = parse_args()

    if args.task == "check":
        check_and_report_env_variables()

    elif args.task == "set":
        paths_dict = set_env_variables()
        modify_conda_activate_startup_script(paths_dict)
        notify_to_reactivate_conda()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Helper program to check and set GToTree database environment variables.",
        epilog="Ex. usage: gtt-data-locations check\n",
    )
    parser.add_argument(
        "task", choices=["check", "set"],
        help="check or set database environment variables",
    )

    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv)


### check ###
def check_and_report_env_variables():

    writable_dict = {}
    paths = {}

    for variable in ENV_VARIABLES:
        path, writable = check_location_var_is_set_and_writable(variable)
        paths[variable] = path
        writable_dict[variable] = writable

    report_message("GToTree database environment variables are set as follows:", "yellow")
    print()

    print("\t  {:<38} {:>10}".format("\033[4mvariable\033[0m", "\033[4mpath\033[0m"))
    for variable in ENV_VARIABLES:
        print("\t  {:<30} {:>10}".format(variable, paths[variable]))
    print()

    for variable, writable in writable_dict.items():
        if not writable:
            report_message(
                "The path set for the '" + variable + "' variable is not writable. "
                "This may cause problems, so you might want to put it somewhere else "
                "(with `gtt-data-locations set`).",
                "red",
            )
            print()


def check_location_var_is_set_and_writable(variable):

    path = os.environ.get(variable, "")

    if not path:
        report_message(
            "The environment variable '" + variable + "' does not seem to be set :(",
            "yellow",
        )
        wprint("Try to set it with `gtt-data-locations set`.")
        print()
        sys.exit(1)

    path_writable = os.access(path, os.W_OK)

    return path, path_writable


### set ###
def get_variable_path(variable):
    return os.environ.get(variable, False)


def set_variable_path(variable, curr_path):

    if curr_path:
        print("\n  " + "-" * 80)
        report_message("The current path set for '" + variable + "' is:", "yellow")
        print("\n    " + str(curr_path) + "\n")

        while True:
            want_to_change = input("  Would you like to change it? (y/n): ")
            if want_to_change in ("y", "n"):
                break
            report_message("Must respond with 'y' or 'n'.", "red")

        if want_to_change == "n":
            print("\n  " + "-" * 80)
            return curr_path

    else:
        print("\n  " + "-" * 80)
        report_message("There is no current path set for '" + variable + "'.", "yellow")

    while True:
        try:
            new_path = input(color_text(
                "\n  Enter the wanted full path for the '" + variable + "' variable: ",
                "yellow",
            ))
            new_path = os.path.join(new_path, "")

            if not os.path.isdir(new_path):
                try:
                    os.makedirs(new_path)
                except Exception:
                    raise PathNotWritable

            if not os.access(new_path, os.W_OK):
                raise PathNotWritable

            if not os.path.isabs(new_path):
                raise PathNotAbsolute

            break

        except PathNotWritable:
            report_message(
                "That location is not writable for you (press `ctrl + c` if wanting to exit).",
                "red",
            )
        except PathNotAbsolute:
            report_message(
                "Please provide an absolute path (press `ctrl + c` if wanting to exit).",
                "red",
            )

    print("\n  " + "-" * 80 + "\n")
    return new_path


def set_env_variables():

    paths_dict = {}
    for variable in ENV_VARIABLES:
        curr_path = get_variable_path(variable)
        paths_dict[variable] = set_variable_path(variable, curr_path)

    return paths_dict


def modify_conda_activate_startup_script(paths_dict):
    """
    Adjust the conda activate startup script so it holds the appropriate
    environment variables and loads them on future activation.

    Some users lack write permission to the conda env area; in that case a
    script is written to the user's ~/.config/gtotree/ location instead (the
    conda activate script sources it if present).
    """

    path_to_startup_script = os.path.join(
        os.environ["CONDA_PREFIX"], "etc/conda/activate.d/gtotree.sh"
    )
    new_one_tmp = path_to_startup_script + ".tmp"

    conda_writable = os.access(path_to_startup_script, os.W_OK)

    if not conda_writable:
        user_config_location = os.path.join(
            os.path.expanduser("~"), ".config/gtotree/", ""
        )
        try:
            os.makedirs(user_config_location, exist_ok=True)
        except Exception:
            report_message(
                "We can't seem to find a place to store startup environment variables :(",
                "yellow",
            )
            wprint("This is likely due to permission restrictions in the conda "
                   "environment location, and in your home location.")
            wprint("If you can't sort this out, please feel free to post an issue here:")
            print("        github.com/AstrobioMike/GToTree/issues\n\n")
            sys.exit(1)

        new_one = os.path.join(user_config_location, "gtotree.sh")

    with open(path_to_startup_script) as initial_file:
        initial_lines = [line.strip() for line in initial_file.readlines()]

    new_lines = []
    for line in initial_lines:

        # if writing to the user's home location, skip the block that sources
        # that same file to avoid a self-sourcing singularity
        if not conda_writable:
            if line.startswith("if") or line.startswith(".") or line.startswith("fi"):
                continue

        try:
            curr_var = line.split(" ")[1].split("=")[0]
        except Exception:
            new_lines.append(line)
            continue

        if curr_var not in paths_dict:
            new_lines.append(line)

    for variable, path in paths_dict.items():
        new_lines.append("export " + variable + "=" + os.path.join(path))

    if conda_writable:
        with open(new_one_tmp, "w") as outfile:
            outfile.write("\n".join(new_lines) + "\n")
        shutil.move(new_one_tmp, path_to_startup_script)
    else:
        with open(new_one, "w") as outfile:
            outfile.write("\n".join(new_lines) + "\n")


def notify_to_reactivate_conda():

    curr_conda_name = os.environ["CONDA_DEFAULT_ENV"]

    print(color_text("\n  " + "-" * 80, "green"))
    print(color_text("  " + "-" * 80, "green"))
    report_message(
        "Environment variables have been updated. But for the changes to take effect, "
        "be sure to reactivate the conda environment, e.g.:",
        "yellow",
    )
    print(f"\n        `conda deactivate && conda activate {curr_conda_name}`\n")
    wprint("Then you can double-check with `gtt-data-locations check`.")
    print(color_text("  " + "-" * 80, "green"))
    print(color_text("  " + "-" * 80 + "\n", "green"))


if __name__ == "__main__":
    main()
