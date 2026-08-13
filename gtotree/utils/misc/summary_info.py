import os

import pandas as pd # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.misc.general import (atomic_write_text, genome_source_label,
                                        genome_input_label)
from gtotree.utils.misc.stages import (GenomeRemovalStage,
                                       GENOME_REMOVAL_STAGE_ORDER)
from gtotree.utils.misc.messaging import (REMOVED_GENOMES_FILENAME,
                                          SCG_INFO_FILENAME)


SCG_INFO_COLUMNS = (
    "target_SCG",
    "num_genomes_with_any_hit",
    "num_genomes_with_hits",
    "perc_genomes_with_hits",
    "num_genomes_after_length_filtering",
    "perc_genomes_after_length_filtering",
    "num_genomes_after_genome_filtering",
    "perc_genomes_after_genome_filtering",
    "retained",
    "stage_removed",
    "reason_removed",
)


REMOVED_GENOMES_COLUMNS = (
    "genome_id",
    "input",
    "source",
    "stage_removed",
    "reason_removed",
)


def write_removed_genomes_report(run_data):
    """
    Record every input genome that has left the run so far, and why
    """
    removed = [gd for gd in run_data.all_input_genomes if gd.removed]
    if not removed:
        return

    removed.sort(key=lambda gd: (GENOME_REMOVAL_STAGE_ORDER.index(gd.removed_at),
                                 gd.id))

    out_path = os.path.join(run_data.run_files_dir, REMOVED_GENOMES_FILENAME)

    def _write(f):
        f.write("\t".join(REMOVED_GENOMES_COLUMNS) + "\n")
        for gd in removed:
            f.write("\t".join([
                gd.id,
                genome_input_label(gd, run_data),
                genome_source_label(gd),
                gd.removed_at,
                gd.reason_removed or "NA",
            ]) + "\n")

    atomic_write_text(out_path, _write)


def _perc(count, total):
    if count is None or not total:
        return "NA"
    return f"{count / total * 100:.1f}"


def _na(value):
    return "NA" if value is None else str(value)


def scg_info_denominators(run_data):
    """
    The two genome pools the representation percentages are against

    `searched` is the pool alive coming out of the HMM search, which is also exactly the
    pool `filter_genes` evaluates `-r` against. So perc_genomes_after_length_filtering
    is directly comparable to the user's `-r` value. `retained` is the pool left after
    `-G`, making perc_genomes_after_genome_filtering the same quantity afterwards. Both
    derive from `removed_at` rather than being stored, which is what keeps them right
    on a resume.
    """
    searched = len(run_data.genomes_alive_through(GenomeRemovalStage.HMM_SEARCH))
    retained = len(run_data.genomes_alive_through(GenomeRemovalStage.SCG_HIT_FILTER))
    return searched, retained


def write_SCG_info_table(run_data):
    """
    Record every target SCG-set, how well represented it was at each filtering stage,
    and whether (and why) it left the run
    """
    if not run_data.SCG_targets:
        return

    searched, retained = scg_info_denominators(run_data)
    out_path = f"{run_data.run_files_dir}/{SCG_INFO_FILENAME}"

    def _write(f):
        f.write("\t".join(SCG_INFO_COLUMNS) + "\n")
        for scg in run_data.SCG_targets:
            after_len = scg.num_genomes_with_hits_after_len_filtering
            after_genome = scg.num_genomes_with_hits_after_genome_filtering
            f.write("\t".join([
                scg.id,
                _na(scg.num_genomes_with_any_hit),
                _na(scg.num_genomes_with_hits),
                _perc(scg.num_genomes_with_hits, searched),
                _na(after_len),
                _perc(after_len, searched),
                _na(after_genome),
                _perc(after_genome, retained),
                "no" if scg.removed else "yes",
                scg.removed_at or "NA",
                scg.reason_removed or "NA",
            ]) + "\n")

    atomic_write_text(out_path, _write)


def SCG_sets_below_representation_cutoff(run_data):
    """
    Retained SCG-sets whose representation fell below `-r` once `-G` had run

    `-r` is applied in `filter_genes` against the genomes alive at that point, and never
    re-checked afterwards, so a set can pass it and then lose most of its genomes to
    `-G`
    """
    cutoff = run_data.gene_representation_cutoff
    if not cutoff:
        return []

    _searched, retained = scg_info_denominators(run_data)
    if not retained:
        return []

    below = []
    for scg in run_data.get_all_SCG_targets_remaining():
        after_genome = scg.num_genomes_with_hits_after_genome_filtering
        if after_genome is None:
            continue
        if after_genome / retained < cutoff:
            below.append(scg)
    return below


def generate_primary_summary_table(args, run_data):

    if run_data.ncbi_sub_table_path:
        ncbi_df     = pd.read_csv(run_data.ncbi_sub_table_path, sep="\t", header=0)
        taxid_map   = dict(zip(ncbi_df["input_accession"], ncbi_df["taxid"], strict=True))
    else:
        taxid_map   = {}

    m = run_data.mapping_dict

    rows = []
    for g in run_data.all_input_genomes:
        label = m.get(g.id) or g.id

        # lookup taxid only if an NCBI sub-table was produced
        taxid = taxid_map.get(g.id, "NA") if run_data.ncbi_sub_table_path else "NA"

        rows.append({
            "assembly_id":                   g.id,
            "label":                         label,
            "source":                        genome_source_label(g),
            "taxid":                         taxid,
            "num_SCG_hits":                  g.num_SCG_hits,
            "num_uniq_SCG_hits":             g.num_unique_SCG_hits,
            "num_SCG_hits_after_filtering":  g.num_SCG_hits_after_filtering,
            "num_total_genes":               g.num_genes,
            "in_final_tree":                 "Yes" if not g.removed else "No",
            "reason_removed":                g.reason_removed,
        })

    df = (pd.DataFrame(rows).convert_dtypes())

    if args.add_gtdb_tax or args.add_ncbi_tax:
        df = add_tax_info(df, run_data, args)

    out = f"{run_data.output_dir}/genomes-summary-info.tsv"
    df.to_csv(out, sep="\t", index=False, na_rep="NA")


def add_tax_info(df, run_data, args):

    # taxonomy was resolved from the hosted Parquet asset during header-updating and
    # stashed on run_data as {input_acc: {rank: value, ...}}; build the summary
    # columns from that rather than re-reading anything from disk
    tax_info = run_data.tax_info_dict or {}

    tax_cols = list(RANKS)

    tax_df = pd.DataFrame(
        [
            {"input_acc": acc, **{rank: ranks.get(rank, "NA") for rank in tax_cols}}
            for acc, ranks in tax_info.items()
        ],
        columns=["input_acc"] + tax_cols,
    )

    # merging
    df = df.merge(
        tax_df,
        how = "left",
        left_on = "assembly_id",
        right_on = "input_acc",
    )

    df = df.drop(columns="input_acc")

    df[tax_cols] = df[tax_cols].fillna("NA")

    return df
