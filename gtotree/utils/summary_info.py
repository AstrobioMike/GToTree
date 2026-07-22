import pandas as pd # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS

def generate_primary_summary_table(args, run_data):

    if run_data.ncbi_sub_table_path:
        ncbi_df     = pd.read_csv(run_data.ncbi_sub_table_path, sep="\t", header=0)
        taxid_map   = dict(zip(ncbi_df["input_accession"], ncbi_df["taxid"]))
    else:
        taxid_map   = {}

    m = run_data.mapping_dict

    rows = []
    for g in run_data.all_input_genomes:
        label = (
            m.get(g.id)
            or m.get(g.full_path)
            or m.get(g.provided_path)
            or g.id
        )

        # lookup taxid only if an NCBI sub-table was produced
        taxid = taxid_map.get(g.id, "NA") if run_data.ncbi_sub_table_path else "NA"

        rows.append({
            "assembly_id":                   g.id,
            "label":                         label,
            "source":                        g.source,
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
