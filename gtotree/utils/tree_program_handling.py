import os
import shutil
import subprocess
from gtotree.utils.messaging import report_early_exit, report_notice, color_text
from gtotree.utils.helper_scripts.gtt_midpoint_root_tree import midpoint_root_tree


def run_tree_building(args, run_data):
    orig_out_tree, final_out_tree = get_out_tree_file(run_data)

    try:
        builder = command_building_dict[args.tree_program]
    except KeyError:
        report_early_exit(run_data, f"Tree building program {args.tree_program!r} not recognized.", suggest_help=True)

    cmd, log_file = builder(args, run_data, orig_out_tree)

    print(f"""
    The tree is being built with {args.tree_program}. You can check the log file for
    progress at:
        {color_text(log_file, "yellow")}
    """)

    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        report_early_exit(run_data, "Tree building failed. Please check the log file noted above.")

    if args.tree_program == "IQTREE":
        shutil.copyfile(os.path.join(run_data.run_files_dir, "iqtree-out", "iqtree.treefile"),
                        orig_out_tree)

    run_data.original_tree_path = orig_out_tree
    midpoint_root_tree(orig_out_tree, final_out_tree)
    run_data.final_tree_path = final_out_tree

    return run_data


def get_out_tree_file(run_data, suffix="aligned-SCGs"):
    output_tree_basename = os.path.basename(run_data.output_dir)
    orig_out_tree = os.path.join(run_data.run_files_dir, f"{output_tree_basename}-non-mid-point-rooted.tre")
    final_out_tree = os.path.join(run_data.output_dir, f"{output_tree_basename}.tre")
    return orig_out_tree, final_out_tree


def build_fasttree_cmd(args, run_data, out_tree):
    log = os.path.join(run_data.logs_dir_rel, "fasttree.log")
    exe = "FastTreeMP" if args.tree_program == "FastTreeMP" else "FastTree"
    if exe == "FastTreeMP":
        exe = f"OMP_NUM_THREADS={args.num_jobs} FastTreeMP"
    nt = "-nt -gtr" if args.nucleotide_mode else ""
    cmd = f"{exe} {nt} {run_data.final_alignment_path} > {out_tree} 2> {log}"
    return cmd, log


def build_veryfasttree_cmd(args, run_data, out_tree):
    log = os.path.join(run_data.logs_dir_rel, "veryfasttree.log")
    threads = f"-threads {args.num_jobs}"
    nt = "-nt -gtr" if args.nucleotide_mode else ""
    cmd = f"VeryFastTree {threads} {nt} {run_data.final_alignment_path} > {out_tree} 2> {log}"
    return cmd, log


def build_iqtree_cmd(args, run_data, out_tree):
    log = os.path.join(run_data.logs_dir_rel, "iqtree.log")
    n = len(run_data.get_all_remaining_input_genomes())
    out_dir = os.path.join(run_data.run_files_dir, "iqtree-out")
    os.makedirs(out_dir, exist_ok=True)

    if n >= 4:
        boot = "-B 1000"
    else:
        report_notice("Bootstrapping not performed with IQTREE with fewer than 4 genomes.")
        boot = ""
    pre  = os.path.join(out_dir, "iqtree")
    cmd  = (
        f"iqtree -s {run_data.final_alignment_path}"
        f" -nt {args.num_jobs} -m MFP {boot}"
        f" -pre {pre} -redo > {log}"
    )
    return cmd, log


command_building_dict = {
    "FastTree":     build_fasttree_cmd,
    "FastTreeMP":   build_fasttree_cmd,
    "VeryFastTree": build_veryfasttree_cmd,
    "IQTREE":       build_iqtree_cmd,
}
