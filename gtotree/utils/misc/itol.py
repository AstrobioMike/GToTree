"""
Writers for the iToL annotation files GToTree produces

Two entry points use these:

  - the automatic per-target files written after a run that used `-p`/`-K`
    (`generate_all_search_itol_files`), and
  - the `gtt itol` subcommand, which decorates an arbitrary list of genome IDs.
"""

import os
import pandas as pd  # type: ignore

from gtotree.utils.misc.general import (SOURCE_ACCESSION, SOURCE_AMINO_ACID,
                                        SOURCE_FASTA, SOURCE_GENBANK,
                                        WANTED_REF_TAX_SOURCE_LABEL,
                                        atomic_write_text)


class ItolError(Exception):
    """Anything that went wrong that we can explain in one sentence."""


SEPARATOR_LINE = "SEPARATOR TAB"

ITOL_COLOR = "#0000ff"

COLOR_MAP = {
    "blue": "#434da7",
    "green": "#48a743",
    "red": "#c01820",
    "purple": "#512f9c",
    "black": "#000000",
    "orange": "#e07b28",
}

SHAPE_MAP = {
    "square": "1",
    "circle": "2",
    "star": "3",
    "rtriangle": "4",
    "ltriangle": "5",
    "check": "6",
}

BINARY_PRESENT = "1"

SEARCH_STYLE_LINE_WEIGHT = 3


################################################################################
# input
################################################################################

def read_targets(path):
    """
    Read a single-column file of genome IDs into an ordered, de-duplicated list.
    """
    if not os.path.isfile(path):
        raise ItolError(f'The input file "{path}" (passed to `-g`) was not found.')

    seen = set()
    targets = []

    with open(path, "r") as in_file:
        for line in in_file:
            target = line.strip()

            if not target or target in seen:
                continue

            seen.add(target)
            targets.append(target)

    if not targets:
        raise ItolError(
            f'The input file "{path}" (passed to `-g`) has no genome IDs in it. It '
            "should hold one ID per line, matching the labels in the tree.")

    return targets


def resolve_color(color):
    """Accept either a named color or a literal hex code."""
    if color in COLOR_MAP:
        return COLOR_MAP[color]

    if color.startswith("#") and len(color) in (4, 7):
        return color

    raise ItolError(
        f'"{color}" was passed to `-c`, but that\'s neither one of the named colors '
        f"({', '.join(COLOR_MAP)}) nor a hex code like \"#0000ff\".")


################################################################################
# writers
################################################################################

def write_binary_dataset(out_path, targets, label="data", color=ITOL_COLOR,
                         shape="square", height_factor=1):
    """
    DATASET_BINARY: one filled symbol per target genome.
    """
    if shape not in SHAPE_MAP:
        raise ItolError(
            f'"{shape}" was passed to `-s`, but the available shapes are: '
            f"{', '.join(SHAPE_MAP)}.")

    shape_code = SHAPE_MAP[shape]

    def _write(f):
        f.write(f"DATASET_BINARY\n{SEPARATOR_LINE}\n\n")
        f.write(f"DATASET_LABEL\t{label}\n\n")
        f.write(f"COLOR\t{color}\n\n")
        f.write(f"FIELD_LABELS\t{label}\n\n")
        f.write(f"FIELD_SHAPES\t{shape_code}\n\n")
        f.write(f"FIELD_COLORS\t{color}\n\n")
        f.write(f"HEIGHT_FACTOR\t{height_factor}\n\n")
        f.write("DATA\n")
        for target in targets:
            f.write(f"{target}\t{BINARY_PRESENT}\n")

    atomic_write_text(out_path, _write)


def write_colorstrip(out_path, targets, label="data", color=ITOL_COLOR,
                     strip_width=25, color_branches_too=False):
    """
    DATASET_COLORSTRIP: a colored band alongside each target genome.
    """
    def _write(f):
        f.write(f"DATASET_COLORSTRIP\n{SEPARATOR_LINE}\n\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write(f"COLOR\t{color}\n\n")
        f.write(f"COLOR_BRANCHES\t{'1' if color_branches_too else '0'}\n\n")
        f.write(f"STRIP_WIDTH\t{strip_width}\n\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#999999\n\n")
        f.write("DATA\n")
        for target in targets:
            f.write(f"{target}\t{color}\t{label}\n")

    atomic_write_text(out_path, _write)


def write_style_dataset(out_path, targets, label="data", color=ITOL_COLOR,
                        what_to_color="branches", line_weight=SEARCH_STYLE_LINE_WEIGHT):
    """
    DATASET_STYLE: color branches and/or labels for a set of genomes
    """

    def _write(f):
        f.write(f"DATASET_STYLE\n{SEPARATOR_LINE}\n\n")
        f.write(f"DATASET_LABEL\t{label}\n\n")
        f.write(f"COLOR\t{color}\n\n")
        f.write("DATA\n")
        for target in targets:
            if what_to_color in ("both", "branches"):
                f.write(f"{target}\tbranch\tnode\t{color}\t{line_weight}\tnormal\n")
            if what_to_color in ("both", "labels"):
                f.write(f"{target}\tlabel\tnode\t{color}\t1\tbold\n")

    atomic_write_text(out_path, _write)


def write_text_dataset(out_path, targets, text, label="data", color=ITOL_COLOR):
    """
    DATASET_TEXT: attach a text label to each target genome
    """
    if "\t" in text:
        raise ItolError(
            "The text passed to `-t` contains a tab, which is the field separator "
            "in the file being written. Please use spaces instead.")

    def _write(f):
        f.write(f"DATASET_TEXT\n{SEPARATOR_LINE}\n\n")
        f.write(f"DATASET_LABEL\t{label}\n\n")
        f.write(f"COLOR\t{color}\n\n")
        f.write("DATA\n")
        for target in targets:
            f.write(f"{target}\t{text}\t-1\t{color}\tnormal\t1\t0\n")

    atomic_write_text(out_path, _write)


################################################################################
# resolving a selection against a finished run
################################################################################

SOURCE_CHOICES = {
    SOURCE_ACCESSION:   SOURCE_ACCESSION,
    SOURCE_GENBANK:     SOURCE_GENBANK,
    SOURCE_FASTA:       SOURCE_FASTA,
    SOURCE_AMINO_ACID:  SOURCE_AMINO_ACID,
    "wanted-ref-tax":   WANTED_REF_TAX_SOURCE_LABEL,
}

SUMMARY_FILENAME = "genomes-summary-info.tsv"

# Columns a wanted-genome can be named by, in the order they're consulted. The
# user's list is whatever they happened to keep (accessions, the file paths
# they passed in, or labels copied off the tree) and any of those should work.
LOOKUP_COLUMNS = ("genome_id", "label", "input", "taxid")


class Selection:
    """
    The outcome of resolving a selection: what to annotate, and what fell out
    """

    def __init__(self, labels, dropped_not_in_tree):
        self.labels = labels
        self.dropped_not_in_tree = dropped_not_in_tree


def load_summary(input_dir):
    """Read a finished run's summary table, or explain why we can't."""
    summary_path = os.path.join(input_dir, SUMMARY_FILENAME)

    if not os.path.isdir(input_dir):
        raise ItolError(
            f'The run directory "{input_dir}" (passed to `-i`) was not found.')

    if not os.path.isfile(summary_path):
        raise ItolError(
            f'"{input_dir}" (passed to `-i`) doesn\'t have a {SUMMARY_FILENAME} in '
            "it, so it doesn't look like a completed GToTree run.")

    df = pd.read_csv(summary_path, sep="\t", header=0, dtype=str).fillna("NA")

    missing = [c for c in ("genome_id", "label", "in_final_tree")
               if c not in df.columns]
    if missing:
        raise ItolError(
            f'The {SUMMARY_FILENAME} in "{input_dir}" is missing expected '
            f"column(s): {', '.join(missing)}.")

    return df


def _build_lookup(df):
    """
    Map every usable identifier to its row index, skipping ambiguous ones
    """
    lookup = {}

    for column in LOOKUP_COLUMNS:
        if column not in df.columns:
            continue

        per_column = {}
        for idx, value in df[column].items():
            if not value or value == "NA":
                continue
            per_column.setdefault(value, set()).add(idx)

        for value, idxs in per_column.items():
            if len(idxs) == 1 and value not in lookup:
                lookup[value] = next(iter(idxs))

    return lookup


def _finish_selection(df, idxs):
    """Split matched rows into annotatable labels and genomes with no leaf."""
    labels = []
    dropped = []

    for idx in idxs:
        row = df.loc[idx]
        if row["in_final_tree"] == "Yes":
            labels.append(row["label"])
        else:
            dropped.append(row["genome_id"])

    if not labels:
        raise ItolError(
            "None of the selected genomes made it into the final tree, so there's "
            "nothing to annotate. Check the `in_final_tree` and `reason_removed` "
            f"columns of {SUMMARY_FILENAME}.")

    return Selection(labels, dropped)


def resolve_wanted_genomes(input_dir, wanted):
    """
    Translate a user's list of genome names into the labels used in the tree
    """
    df = load_summary(input_dir)
    lookup = _build_lookup(df)

    idxs = []
    unmatched = []

    for name in wanted:
        if name in lookup:
            idxs.append(lookup[name])
        else:
            unmatched.append(name)

    if unmatched:
        shown = ", ".join(unmatched[:10])
        more = f" (and {len(unmatched) - 10} more)" if len(unmatched) > 10 else ""
        raise ItolError(
            f"{len(unmatched)} of the genomes passed to `-g` weren't found in "
            f'"{input_dir}": {shown}{more}. They should match the genome IDs, '
            f"labels, or original inputs listed in {SUMMARY_FILENAME}.")

    return _finish_selection(df, idxs)


def select_by_source(input_dir, sources):
    """
    Select every genome that came from one or more given input sources
    """
    df = load_summary(input_dir)

    if "source" not in df.columns:
        raise ItolError(
            f'The {SUMMARY_FILENAME} in "{input_dir}" has no `source` column, so '
            "genomes can't be selected with `--source`.")

    wanted_labels = {SOURCE_CHOICES[s] for s in sources}
    idxs = [idx for idx, value in df["source"].items() if value in wanted_labels]

    if not idxs:
        present = sorted({v for v in df["source"] if v and v != "NA"})
        raise ItolError(
            f"No genomes in \"{input_dir}\" came from: {', '.join(sorted(sources))}. "
            f"Sources present in this run: {', '.join(present) or 'none'}.")

    return _finish_selection(df, idxs)


################################################################################
# automatic per-target files after a Pfam/KO search
################################################################################

def _read_hit_counts(counts_path):
    """
    Return (target_ids, {genome_id: {target_id: count}}).

    The hit-counts tsv has header: genome_id, total_gene_count, <target ...>.
    Target columns are everything after the first two.
    """
    df = pd.read_csv(counts_path, sep="\t", header=0, dtype={"genome_id": str})
    target_ids = [c for c in df.columns if c not in ("genome_id", "total_gene_count")]

    counts_by_genome = {}
    for _, row in df.iterrows():
        genome_id = row["genome_id"]
        counts_by_genome[genome_id] = {t: int(row[t]) for t in target_ids}

    return target_ids, counts_by_genome


def _read_summary(summary_path):
    """
    Return (label_map, in_tree_set, labels_swapped) from genomes-summary-info.tsv.

    label_map: genome_id -> label (label differs from genome_id only when
               headers were swapped; if identical, no swap effectively happened)
    in_tree_set: genome_ids whose in_final_tree == "Yes"
    """
    df = pd.read_csv(summary_path, sep="\t", header=0, dtype={"genome_id": str})

    label_map = dict(zip(df["genome_id"], df["label"].astype(str), strict=True))
    in_tree_set = set(df.loc[df["in_final_tree"] == "Yes", "genome_id"])

    # a swap occurred if any label differs from its genome_id
    labels_swapped = any(label_map[a] != a for a in label_map)

    return label_map, in_tree_set, labels_swapped


def generate_search_itol_files(counts_path, summary_path, itol_dir):
    """
    For each target with >=1 hit among genomes retained in the final tree,
    generate a pair of iToL files: a DATASET_STYLE (`-branches-`) coloring the
    branches, and a DATASET_BINARY (`-symbols-`) placing a filled square at each
    hit genome.

    counts_path : path to {pfam,ko}-hit-counts.tsv
    summary_path: path to genomes-summary-info.tsv
    itol_dir    : output directory for the iToL files (created if needed)

    Returns the list of targets for which files were written.
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
        for genome_id, counts in counts_by_genome.items():
            if counts.get(target, 0) > 0 and genome_id in in_tree_set:
                leaf = label_map.get(genome_id, genome_id) if labels_swapped else genome_id
                leaf_ids.append(leaf)

        if not leaf_ids:
            continue

        write_style_dataset(
            os.path.join(itol_dir, f"{target}-branches-itol.txt"), leaf_ids,
            label=target, color=ITOL_COLOR, what_to_color="branches",
            line_weight=SEARCH_STYLE_LINE_WEIGHT)
        write_binary_dataset(
            os.path.join(itol_dir, f"{target}-symbols-itol.txt"), leaf_ids,
            label=target, color=ITOL_COLOR, shape="square")
        written.append(target)

    return written


################################################################################
# automatic per-input-source files
################################################################################

SOURCE_FILE_SPECS = {
    SOURCE_FASTA:                (ITOL_COLOR,          "nucleotide-fasta"),
    SOURCE_AMINO_ACID:           (COLOR_MAP["purple"], "amino-acid-fasta"),
    SOURCE_GENBANK:              (COLOR_MAP["green"],  "genbank-file"),
    SOURCE_ACCESSION:            (COLOR_MAP["red"],    "ncbi-accession"),
    WANTED_REF_TAX_SOURCE_LABEL: (COLOR_MAP["orange"], "wanted-ref-tax"),
}


def generate_input_source_itol_files(summary_path, itol_dir):
    """
    Write one DATASET_STYLE file per input source, when a run drew from more than
    one
    """
    if not os.path.isfile(summary_path):
        return []

    df = pd.read_csv(summary_path, sep="\t", header=0, dtype=str).fillna("NA")

    if "source" not in df.columns:
        return []

    in_tree = df[df["in_final_tree"] == "Yes"]

    # iterate the spec rather than the frame so file order is stable run to run
    present = [src for src in SOURCE_FILE_SPECS
               if (in_tree["source"] == src).any()]

    if len(present) < 2:
        return []

    os.makedirs(itol_dir, exist_ok=True)

    written = []
    for source in present:
        color, slug = SOURCE_FILE_SPECS[source]
        labels = list(in_tree.loc[in_tree["source"] == source, "label"])

        write_style_dataset(os.path.join(itol_dir, f"{slug}-branches-itol.txt"),
                            labels, label=source, color=color,
                            what_to_color="branches",
                            line_weight=SEARCH_STYLE_LINE_WEIGHT)
        write_binary_dataset(os.path.join(itol_dir, f"{slug}-symbols-itol.txt"),
                             labels, label=source, color=color, shape="square")
        written.append(source)

    return written


def generate_all_search_itol_files(args, run_data):
    """
    Post-tree orchestration: generate iToL files for whichever additional
    searches were run (Pfam and/or KO). Reads the summary table written by
    generate_primary_summary_table plus each search's hit-counts table.
    """
    summary_path = f"{run_data.output_dir}/genomes-summary-info.tsv"

    generate_input_source_itol_files(
        summary_path=summary_path,
        itol_dir=f"{run_data.output_dir}/itol-files",
    )

    if run_data.target_pfams_file:
        generate_search_itol_files(
            counts_path=f"{run_data.pfam_results_dir}/pfam-hit-counts.tsv",
            summary_path=summary_path,
            itol_dir=f"{run_data.pfam_results_dir}/itol-files",
        )

    if run_data.target_kos_file:
        generate_search_itol_files(
            counts_path=f"{run_data.ko_results_dir}/ko-hit-counts.tsv",
            summary_path=summary_path,
            itol_dir=f"{run_data.ko_results_dir}/itol-files",
        )
