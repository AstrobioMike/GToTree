"""
Output tables for `gtt search-pfams` / `gtt search-kos`.

The main driver's `generate_primary_summary_table` can't be reused: its columns are
SCG-hit counts and `in_final_tree`, neither of which exists without a tree. What a
search-only run can say about each genome is what it was, whether it made it through
processing, how many genes it had, and whether its search finished -- so that's the
table.

The counts matrix and per-target hit fastas are written by the main driver's own
helpers (`write_pfam_counts_table` / `write_ko_counts_table` and
`combine_all_*_hits`), unmodified, so those files are byte-for-byte what a full GToTree
run with `-p`/`-K` would have produced.
"""

import os

from gtotree.utils.misc.general import (atomic_write_text, genome_source_label,
                                        genome_input_label)
from gtotree.utils.misc.summary_info import search_completed_value


SUMMARY_FILENAME = "genomes-summary-info.tsv"

COLUMNS = ["genome_id", "input_source", "input_path", "num_genes",
           "prodigal_used", "search_completed", "reason_removed"]


def _row(gd, spec, run_data=None):
    completed = search_completed_value(gd, spec.search_done_flag,
                                       spec.search_failed_flag)

    return [
        gd.id,
        genome_source_label(gd),
        genome_input_label(gd, run_data),
        "NA" if gd.num_genes is None else str(gd.num_genes),
        "Yes" if gd.prodigal_used else "No",
        completed,
        gd.reason_removed or "NA",
    ]


def write_genomes_summary(out_dir, run_data, spec):
    """
    Write one row per input genome, in input order.

    Every input genome appears, including ones that dropped out -- a genome missing
    from the counts matrix is only interpretable if there's somewhere that says why it
    isn't there.
    """
    path = os.path.join(out_dir, SUMMARY_FILENAME)

    rows = [_row(gd, spec, run_data) for gd in run_data.all_input_genomes]

    def write(f):
        f.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    atomic_write_text(path, write)

    return path


def summarize_counts(run_data, spec):
    """
    Return (num_searched, num_removed, num_failed_search) for the finish banner.
    """
    searched = 0
    removed = 0
    failed = 0

    failed_attr = spec.search_failed_flag

    for gd in run_data.all_input_genomes:
        if gd.removed:
            removed += 1
            continue
        if getattr(gd, failed_attr, False):
            failed += 1
        elif getattr(gd, spec.search_done_flag, False):
            searched += 1
        else:
            failed += 1

    return searched, removed, failed
