import pandas as pd

def generate_primary_summary_table(args, run_data):
    # building a lookup if needed
    if args.ncbi_accessions:
        ncbi_df     = pd.read_csv(run_data.ncbi_sub_table_path, sep="\t", header=0)
        taxid_map   = dict(zip(ncbi_df["input_accession"], ncbi_df["taxid"]))
    else:
        taxid_map   = {}

    # quick alias for mapping dict
    m = run_data.mapping_dict

    rows = []
    for g in run_data.all_input_genomes:
        label = (
            m.get(g.id)
            or m.get(g.full_path)
            or m.get(g.provided_path)
            or g.id
        )

        # lookup taxid only if requested
        taxid = taxid_map.get(g.id, "NA") if args.ncbi_accessions else "NA"

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
        })

    df = (pd.DataFrame(rows).convert_dtypes())

    if args.add_gtdb_tax or args.add_ncbi_tax:
        df = add_tax_info(df, run_data, args)

    out = f"{run_data.output_dir}/genomes-summary-info.tsv"
    df.to_csv(out, sep="\t", index=False, na_rep="NA")


def add_tax_info(df, run_data, args):

    if args.add_gtdb_tax:
        tax_df_path = run_data.tmp_dir + "/gtdb-subset-tax.tsv"
    else:
        tax_df_path = run_data.tmp_dir + "/ncbi-accession-tax.tsv"

    tax_df = pd.read_csv(tax_df_path, sep="\t", header=0)

    # getting taxonomy columns
    start_col = list(tax_df.columns).index("domain")
    tax_cols = tax_df.columns[start_col:].tolist()

    # merging
    df = df.merge(
        tax_df[["input_acc"] + tax_cols],
        how = "left",
        left_on = "assembly_id",
        right_on = "input_acc",
    )

    df = df.drop(columns="input_acc")

    df[tax_cols] = df[tax_cols].fillna("NA")

    return df
