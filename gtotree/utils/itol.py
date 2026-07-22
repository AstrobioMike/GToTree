import os
import pandas as pd  # type: ignore


ITOL_COLOR = "#0000ff"


def _read_hit_counts(counts_path):
    """
    Return (target_ids, {assembly_id: {target_id: count}}).

    The hit-counts tsv has header: assembly_id, total_gene_count, <target ...>.
    Target columns are everything after the first two.
    """
    df = pd.read_csv(counts_path, sep="\t", header=0, dtype={"assembly_id": str})
    target_ids = [c for c in df.columns if c not in ("assembly_id", "total_gene_count")]

    counts_by_genome = {}
    for _, row in df.iterrows():
        assembly_id = row["assembly_id"]
        counts_by_genome[assembly_id] = {t: int(row[t]) for t in target_ids}

    return target_ids, counts_by_genome


def _read_summary(summary_path):
    """
    Return (label_map, in_tree_set) from genomes-summary-info.tsv.

    label_map: assembly_id -> label (label differs from assembly_id only when
               headers were swapped; if identical, no swap effectively happened)
    in_tree_set: assembly_ids whose in_final_tree == "Yes"
    """
    df = pd.read_csv(summary_path, sep="\t", header=0, dtype={"assembly_id": str})

    label_map = dict(zip(df["assembly_id"], df["label"].astype(str)))
    in_tree_set = set(df.loc[df["in_final_tree"] == "Yes", "assembly_id"])

    # a swap occurred if any label differs from its assembly_id
    labels_swapped = any(label_map[a] != a for a in label_map)

    return label_map, in_tree_set, labels_swapped


def _write_itol_file(out_path, target, leaf_ids):
    header = (
        "DATASET_STYLE\n"
        "SEPARATOR SPACE\n"
        f"DATASET_LABEL {target}\n"
        f"COLOR {ITOL_COLOR}\n"
        "DATA\n"
    )
    with open(out_path, "w") as f:
        f.write(header)
        for leaf in leaf_ids:
            f.write(f"{leaf} branch node {ITOL_COLOR} 3 normal\n")


def generate_search_itol_files(counts_path, summary_path, itol_dir):
    """
    Generate one DATASET_STYLE iToL file per target that has >=1 hit among
    genomes retained in the final tree.

    counts_path : path to {pfam,ko}-hit-counts.tsv
    summary_path: path to genomes-summary-info.tsv
    itol_dir    : output directory for the iToL files (created if needed)

    Returns the list of targets for which a file was written.
    """
    if not (os.path.isfile(counts_path) and os.path.isfile(summary_path)):
        return []

    target_ids, counts_by_genome = _read_hit_counts(counts_path)
    label_map, in_tree_set, labels_swapped = _read_summary(summary_path)

    os.makedirs(itol_dir, exist_ok=True)

    written = []
    for target in target_ids:
        # genomes with a hit to this target AND retained in the final tree
        leaf_ids = []
        for assembly_id, counts in counts_by_genome.items():
            if counts.get(target, 0) > 0 and assembly_id in in_tree_set:
                leaf = label_map.get(assembly_id, assembly_id) if labels_swapped else assembly_id
                leaf_ids.append(leaf)

        if not leaf_ids:
            continue

        out_path = os.path.join(itol_dir, f"{target}-iToL.txt")
        _write_itol_file(out_path, target, leaf_ids)
        written.append(target)

    return written


def generate_all_search_itol_files(args, run_data):
    """
    Post-tree orchestration: generate iToL files for whichever additional
    searches were run (Pfam and/or KO). Reads the summary table written by
    generate_primary_summary_table plus each search's hit-counts table.

    Fine to call unconditionally; it no-ops for any search that wasn't run or
    whose counts table is absent (e.g. no targets were found).
    """
    summary_path = f"{run_data.output_dir}/genomes-summary-info.tsv"

    if run_data.target_pfams_file:
        generate_search_itol_files(
            counts_path=f"{run_data.pfam_results_dir}/pfam-hit-counts.tsv",
            summary_path=summary_path,
            itol_dir=f"{run_data.pfam_results_dir}/iToL-files",
        )

    if run_data.target_kos_file:
        generate_search_itol_files(
            counts_path=f"{run_data.ko_results_dir}/ko-hit-counts.tsv",
            summary_path=summary_path,
            itol_dir=f"{run_data.ko_results_dir}/iToL-files",
        )
