import sys
import os
import argparse
import pyarrow.parquet as pq # type: ignore

from gtotree.utils.messaging import wprint, color_text
from gtotree.utils.gtdb.get_gtdb_data import (get_gtdb_data, gtdb_data_table_path,
                                              report_gtdb_version_info as _read_gtdb_version_info)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import TaxonNotFound, AmbiguousTaxon
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes


_RANK_COLUMNS = list(RANKS)


################################################################################

def main():
    args = parse_args()
    get_accessions_from_gtdb(args)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(

        description="This is a helper program to facilitate using taxonomy and genomes "
                    "from the Genome Taxonomy Database (gtdb.ecogenomic.org) with GToTree. "
                    "It primarily returns NCBI accessions and GTDB metadata subsets based "
                    "on GTDB-taxonomy searches, with optional filtering to GTDB "
                    "representative species or RefSeq reference genomes, and optional "
                    "dereplication down to one genome per rank.",
        epilog="Ex. usage: gtt-get-accessions-from-GTDB -t Archaea --GTDB-representatives-only\n")

    parser.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        help=("Target taxon (enter 'all' for all). Not needed with `--get-rank-counts`."),
        action="store",
    )

    parser.add_argument(
        "-r",
        "--target-rank",
        metavar="<STR>",
        help=("Target rank (if needed to disambiguate a taxon name that exists at multiple ranks)"),
        action="store",
    )
    parser.add_argument(
        "--derep-rank",
        action="store",
        default="off",
        help="Dereplicate the pulled genomes down to a single best genome "
             "per unique value of this rank (e.g., '--derep-rank family' "
             "keeps one genome per family within the target taxon). "
             "Default: off (all matching genomes are returned). Use 'auto' "
             "for two ranks finer than the target, or pass an explicit rank."
    )

    parser.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help="Add this flag along with a specified taxon to the `-t` parameter "
             "to see how many of that taxon are in the database.",
    )

    parser.add_argument(
        "--get-rank-counts",
        action="store_true",
        help="Provide just this flag alone to see counts of how many "
             "unique taxa there are for each rank.",
    )

    parser.add_argument(
        "-G",
        "--gtdb-representatives-only",
        action="store_true",
        help="Add this flag to only pull accessions for genomes "
             "designated as GTDB species representatives (see, e.g., "
             "https://gtdb.ecogenomic.org/faq#gtdb_species_clusters).",
    )

    parser.add_argument(
        "-R",
        "--refseq-reference-genomes-only",
        action="store_true",
        help="Add this flag to only pull accessions for genomes designated as "
             "RefSeq \"reference\" genomes (these used to be called \"representative\" genomes, see, e.g., "
             "https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes).",
    )

    parser.add_argument(
        "--get-table",
        action="store_true",
        help="Provide just this flag alone to write out a tsv of GToTree's "
             "GTDB metadata table.")

    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args(argv)

################################################################################


def get_accessions_from_gtdb(args):

    # make sure the prepared GTDB Parquet is present, then work against it
    get_gtdb_data()
    gtdb_path = gtdb_data_table_path()

    if args.get_table:
        copy_gtdb_table(gtdb_path)
        sys.exit(0)

    _validate_flag_combo(args)

    _report_gtdb_version(gtdb_path)

    representatives_source = _representatives_source(args)

    if args.get_rank_counts:
        full = _read_rank_columns(gtdb_path)
        rep = _apply_reps_filter(gtdb_path, representatives_source)
        get_unique_taxa_counts_of_all_ranks(full, rep,
                                            representatives_source=representatives_source)
        sys.exit(0)

    if not args.target_taxon:
        return

    if args.get_taxon_counts:
        full = _read_rank_columns(gtdb_path)
        rep = _apply_reps_filter(gtdb_path, representatives_source)
        _report_taxon_counts_or_exit(args.target_taxon, full, rep, representatives_source)
        sys.exit(0)

    _select_and_write(gtdb_path, args, representatives_source)
    sys.exit(0)


def _validate_flag_combo(args):
    if args.get_taxon_counts and not args.target_taxon:
        print("")
        wprint(color_text("A specific taxon needs to also be provided to the `-t` flag "
                          "in order to use `--get-taxon-counts`.", "yellow"))
        print("")
        wprint("  E.g.,: gtt-get-accessions-from-GTDB --get-taxon-counts -t Alteromonas")
        print("")
        sys.exit(1)

    if args.gtdb_representatives_only and args.refseq_reference_genomes_only:
        print("")
        wprint(color_text("Only one of `--GTDB-representatives-only` or "
                          "`--RefSeq-reference-genomes-only` can be provided.", "yellow"))
        print("")
        sys.exit(1)


def _representatives_source(args):
    if args.gtdb_representatives_only:
        return "GTDB"
    if args.refseq_reference_genomes_only:
        return "RefSeq"
    return None


def _select_and_write(gtdb_path, args, representatives_source):
    """
    Pull accessions for the target taxon through the shared taxonomy core, then write
    the accessions list and a metadata-subset TSV. The 'all' taxon is a bulk dump that
    the core's taxon-resolution doesn't cover, so it's handled directly.
    """
    if args.target_taxon.lower() == "all":
        _write_all(gtdb_path, representatives_source)
        return

    # reps-only requested -> pass through to the core; RefSeq maps to the NCBI-style
    # reference-genome filter that SOURCES['gtdb'] understands via reps_only, while GTDB
    # representatives are the core's default representative pool.
    reps_only = representatives_source is not None

    try:
        selection = select_ref_genomes(
            gtdb_path, "gtdb", args.target_taxon,
            target_rank=args.target_rank, derep_rank=args.derep_rank,
            reps_only=reps_only)
    except AmbiguousTaxon:
        wprint(color_text("Since '" + str(args.target_taxon) + "' occurs at more than 1 "
               "rank, we'll need to specify which rank is wanted as well before we pull "
               "the accessions. This can be specified with the `-r` flag.", "yellow"))
        print("")
        sys.exit(0)
    except TaxonNotFound:
        wprint(color_text("Input taxon '" + str(args.target_taxon) + "' doesn't seem to "
               "exist at any rank :(", "yellow"))
        print("")
        sys.exit(0)
    except ValueError as err:
        # a derep rank coarser than the target's rank
        wprint(color_text(str(err), "yellow"))
        print("")
        sys.exit(0)

    if selection.canonical != args.target_taxon:
        wprint(color_text("Matched input '" + str(args.target_taxon) + "' to GTDB taxon '"
               + selection.canonical + "'.", "yellow"))
        print("")

    for warning in selection.warnings:
        wprint(color_text(warning, "yellow"))
        print("")

    if not selection.accessions:
        wprint(color_text("No accessions were found for the given target :(", "yellow"))
        print("")
        sys.exit(0)

    _write_outputs(selection, representatives_source)


def _write_outputs(selection, representatives_source):
    taxon_for_filename = selection.canonical.replace(" ", "-")
    rank = selection.resolved_rank

    if representatives_source:
        suffix = "-" + representatives_source + "-rep"
    else:
        suffix = ""

    acc_out_filename = f"GTDB-{taxon_for_filename}-{rank}{suffix}-accs.txt"
    tab_out_filename = f"GTDB-{taxon_for_filename}-{rank}{suffix}-metadata.tsv"

    with open(acc_out_filename, "w") as out:
        for acc in selection.accessions:
            out.write(acc + "\n")

    _write_metadata_tsv(selection.rows, tab_out_filename)

    print("")
    wprint(f"Wrote {len(selection.accessions):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def _write_all(gtdb_path, representatives_source):
    """Bulk dump of every accession (optionally reps-only), plus a metadata TSV."""
    if representatives_source == "GTDB":
        filt = [("gtdb_representative", "=", "t")]
    elif representatives_source == "RefSeq":
        filt = [("ncbi_refseq_category", "=", "reference genome")]
    else:
        filt = None

    if filt:
        table = pq.read_table(gtdb_path, filters=filt)
    else:
        table = pq.read_table(gtdb_path)

    accessions = table.column("ncbi_genbank_assembly_accession").to_pylist()

    if representatives_source:
        acc_out_filename = "GTDB-arc-and-bac-" + representatives_source + "-rep-accessions.txt"
        tab_out_filename = "GTDB-arc-and-bac-" + representatives_source + "-rep-metadata.tsv"
        table.to_pandas().to_csv(tab_out_filename, sep="\t", index=False)
    else:
        acc_out_filename = "GTDB-arc-and-bac-accessions.txt"
        tab_out_filename = None

    with open(acc_out_filename, "w") as out:
        for acc in accessions:
            out.write(str(acc) + "\n")

    print("")
    wprint(f"Wrote {len(accessions):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    if tab_out_filename:
        wprint("Associated taxonomy and metadata written to:")
        wprint("  " + color_text(tab_out_filename))
        print("")


def _write_metadata_tsv(rows, out_filename):
    """Write selected genome rows to a TSV, columns in the asset's natural order."""
    if not rows:
        open(out_filename, "w").close()
        return
    # preserve a stable, readable column order: accession + ranks first, then the rest
    first = ["ncbi_genbank_assembly_accession"] + list(RANKS)
    seen = set(first)
    header = [c for c in first if c in rows[0]] + [c for c in rows[0] if c not in seen]
    with open(out_filename, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in header) + "\n")


################################################################################
# counts + table helpers (read straight from the Parquet, no selection core needed)
################################################################################

def _read_rank_columns(gtdb_path):
    return pq.read_table(gtdb_path, columns=_RANK_COLUMNS).to_pandas()


def _apply_reps_filter(gtdb_path, representatives_source):
    if not representatives_source:
        return None
    if representatives_source == "GTDB":
        filt = [("gtdb_representative", "=", "t")]
    else:
        filt = [("ncbi_refseq_category", "=", "reference genome")]
    return pq.read_table(gtdb_path, columns=_RANK_COLUMNS, filters=filt).to_pandas()


def find_ranks_for_taxon(taxon, tab):
    """ ranks (of the 7) whose column contains `taxon`, within a DataFrame """
    return [rank for rank in RANKS if taxon in tab[rank].unique()]


def _report_taxon_counts_or_exit(taxon, full_df, rep_df, representatives_source):
    if taxon.lower() == "all":
        count = len(full_df.index)
        print("")
        wprint("  There are " + str(count) + " total genomes in the database.")
        print("")
        if representatives_source:
            wprint(color_text("In considering only " + representatives_source
                              + " representative genomes:", "yellow"))
            print("")
            wprint("  There are " + str(len(rep_df.index))
                   + " total representative genomes in the database.")
            print("")
        return

    ranks_found_in = find_ranks_for_taxon(taxon, full_df)
    if len(ranks_found_in) == 0:
        wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any "
               "rank :(", "yellow"))
        wprint("Be aware, to be on the safe side, our searching is case-sensitive.")
        print("")
        sys.exit(0)

    print("")
    for rank in ranks_found_in:
        count = len(full_df[full_df[rank] == taxon].index)
        wprint("  The rank '" + rank + "' has " + str(count) + " " + taxon + " entries.")
    print("")

    if representatives_source:
        wprint(color_text("In considering only " + representatives_source
                          + " representative genomes:", "yellow"))
        print("")
        ranks_found_in_rep = find_ranks_for_taxon(taxon, rep_df)
        for rank in ranks_found_in_rep:
            count = len(rep_df[rep_df[rank] == taxon].index)
            wprint("  The rank '" + rank + "' has " + str(count) + " " + taxon
                   + " representative genome entries.")
            print("")
        if len(ranks_found_in_rep) == 0:
            wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any "
                   "rank as a representative genome :(", "yellow"))
            print("")
            sys.exit(0)


def get_unique_taxa_counts_of_all_ranks(gtdb_tab, gtdb_rep_tab=None,
                                        representatives_source=None):
    """ counts of unique taxa at each rank, from a DataFrame """
    print("\n    {:<10} {:}".format("Rank", "Num. Unique Taxa"))
    for rank in RANKS:
        print("    {:<10} {:}".format(rank, str(gtdb_tab[rank].nunique())))
    print("")

    if representatives_source == "GTDB":
        wprint(color_text("(The `--GTDB-representatives-only` flag doesn't change these "
                          "counts: every GTDB taxon has a representative genome, so the "
                          "number of unique taxa per rank is the same with or without it.)",
                          "yellow"))
        print("")
    elif representatives_source == "RefSeq":
        wprint(color_text("In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print("    {:<10} {:}".format("Rank", "Num. Unique Ref. Taxa"))
        for rank in RANKS:
            print("    {:<10} {:}".format(rank, str(gtdb_rep_tab[rank].nunique())))
        print("")


def _report_gtdb_version(gtdb_path):
    version, release_date = _read_gtdb_version_info(os.path.dirname(gtdb_path))
    print("\n    Using GTDB " + version + ": " + release_date)


def copy_gtdb_table(gtdb_path):
    """
    Materialize the full GTDB metadata table (all columns in the asset) to the current
    directory as a TSV -- the `--get-table` escape hatch.
    """
    _report_gtdb_version(gtdb_path)
    out_name = "GTDB-arc-and-bac-metadata.tsv"
    pq.read_table(gtdb_path).to_pandas().to_csv(out_name, sep="\t", index=False)
    print("")
    wprint("GTDB table written to:")
    print(color_text("    " + out_name + "\n"))


if __name__ == "__main__":
    main()
