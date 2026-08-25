"""
Output tables for `gtt search-annotations`.

The main driver's `generate_primary_summary_table` can't be reused: its columns are
SCG-hit counts and `in_final_tree`, neither of which exists without a tree.

A combined run writes two kinds of genomes summary:

  * a per-target-type table under each type's subdirectory
    (`pfam/pfam-genomes-summary-info.tsv`, `ko/ko-genomes-summary-info.tsv`), carrying
    that type's per-genome hit counts and whether that type's search completed; and
  * one run-level table at the top (`genomes-summary-info.tsv`) covering only what's
    shared across target types -- what each genome was, whether it made it through
    preprocessing, its gene count, and why it was dropped if it was. The per-type hit
    and search-completion columns are deliberately absent here: they'd be ambiguous
    (which target type does `search_completed` refer to?) and are already answered,
    unambiguously, by the per-type tables.

The counts matrix and per-target hit fastas are written by the main driver's own
helpers (`write_pfam_counts_table` / `write_ko_counts_table` and
`combine_all_*_hits`), so those files are byte-for-byte what a full GToTree run with
`-p`/`-K` would have produced. The counts writer additionally returns the per-genome
hit tallies it already computed, which feed the per-type table's `num_hits` /
`num_targets_hit` columns.
"""

import os

from gtotree.utils.misc.general import (atomic_write_text, genome_source_label,
                                        genome_input_label)
from gtotree.utils.misc.summary_info import search_completed_value


# run-level table: what's true of a genome regardless of target type
ROOT_SUMMARY_FILENAME = "genomes-summary-info.tsv"
ROOT_COLUMNS = ["genome_id", "input", "source", "num_genes", "prodigal_used",
                "reason_removed"]

# per-type table: the above plus this target type's hit counts and search status
SPEC_COLUMNS = ["genome_id", "input", "source", "num_genes", "num_hits",
                "num_targets_hit", "prodigal_used", "search_completed",
                "reason_removed"]


def _root_row(gd, run_data=None):
    return [
        gd.id,
        genome_input_label(gd, run_data),
        genome_source_label(gd),
        "NA" if gd.num_genes is None else str(gd.num_genes),
        "Yes" if gd.prodigal_used else "No",
        gd.reason_removed or "NA",
    ]


def _spec_row(gd, spec, run_data=None, hit_tallies=None):
    completed = search_completed_value(gd, spec.search_done_flag,
                                       spec.search_failed_flag)

    num_hits, num_targets_hit = (hit_tallies or {}).get(gd.id, (None, None))

    return [
        gd.id,
        genome_input_label(gd, run_data),
        genome_source_label(gd),
        "NA" if gd.num_genes is None else str(gd.num_genes),
        "NA" if num_hits is None else str(num_hits),
        "NA" if num_targets_hit is None else str(num_targets_hit),
        "Yes" if gd.prodigal_used else "No",
        completed,
        gd.reason_removed or "NA",
    ]


def _write_table(path, columns, rows):
    def write(f):
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    atomic_write_text(path, write)
    return path


def write_spec_genomes_summary(spec_out_dir, run_data, spec, hit_tallies=None):
    """
    Write this target type's per-genome summary into its subdirectory, in input order.

    Includes the type-specific hit counts and this type's `search_completed`, since the
    whole point of the per-type file is to answer "how did the <type> search go for
    each genome" without the ambiguity a shared table would have.
    """
    path = os.path.join(spec_out_dir, spec.summary_filename)
    rows = [_spec_row(gd, spec, run_data, hit_tallies)
            for gd in run_data.all_input_genomes]
    return _write_table(path, SPEC_COLUMNS, rows)


def write_root_genomes_summary(out_dir, run_data):
    """
    Write the run-level per-genome summary at the top of the output directory.

    Only the target-type-independent columns: what the genome was, whether it made it
    through preprocessing, its gene count, and why it was removed if it was. No hit or
    search-completion columns -- those live in the per-type tables where they have an
    unambiguous meaning.
    """
    path = os.path.join(out_dir, ROOT_SUMMARY_FILENAME)
    rows = [_root_row(gd, run_data) for gd in run_data.all_input_genomes]
    return _write_table(path, ROOT_COLUMNS, rows)


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
