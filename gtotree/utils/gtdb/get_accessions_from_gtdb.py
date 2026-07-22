import sys
import os
import argparse
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.messaging import wprint, color_text, report_message
from gtotree.utils.gtdb.get_gtdb_data import (get_gtdb_data, gtdb_data_table_path,
                                              report_gtdb_version_info as _read_gtdb_version_info)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, AmbiguousTaxon,
                                               find_ranks_for_taxon as _resolve_ranks)
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes


_RANK_COLUMNS = list(RANKS)


################################################################################

def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate using taxonomy and genomes "
            "from the Genome Taxonomy Database (gtdb.ecogenomic.org) with GToTree. "
            "It primarily returns NCBI accessions and GTDB metadata subsets based "
            "on GTDB-taxonomy searches, with optional filtering to GTDB "
            "representative species or RefSeq reference genomes, plus optional "
            "dereplication down to one genome per specified rank.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "get-accs-from-gtdb",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt get-accs-from-gtdb -t Archaea --gtdb-representatives-only`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        help=("Target taxon (enter 'all' for all). Not needed with `--get-rank-counts`."),
        action="store",
    )

    optional.add_argument(
        "-r",
        "--target-rank",
        choices=list(RANKS),
        help=("Target rank (if needed to disambiguate a taxon name that exists at multiple ranks)"),
        action="store",
    )
    optional.add_argument(
        "--derep-rank",
        choices=["auto", "off"] + list(RANKS),
        default="off",
        help=("Dereplicate the pulled genomes down to a single best genome per unique "
              "value of this rank (default: off). E.g., '--derep-rank family' keeps one genome per "
              "family within the target taxon). Use 'auto' for two ranks finer than the target."),
        action="store",
    )

    optional.add_argument(
        "-G",
        "--gtdb-representatives-only",
        action="store_true",
        help=("Pull only genomes designated as GTDB species representatives."),
    )

    optional.add_argument(
        "-R",
        "--refseq-reference-genomes-only",
        action="store_true",
        help=("Pull only genomes designated as RefSeq reference genomes."),
    )

    optional.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help=("Provide this flag along with a specified taxon to `-t` to see how many "
              "genomes match the set parameters (excluding --derep-rank)"),
    )

    optional.add_argument(
        "--get-rank-counts",
        action="store_true",
        help=("Provide just this flag alone to see counts of how many unique taxa there "
              "are for each rank."),
    )

    optional.add_argument(
        "--get-table",
        action="store_true",
        help=("Provide just this flag alone to write out a tsv of GToTree's "
              "GTDB metadata table."),
    )

    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    get_accessions_from_gtdb(args)

################################################################################


def get_accessions_from_gtdb(args):

    preflight_checks(args)

    # make sure the prepared GTDB Parquet is present, then work against it
    get_gtdb_data()
    gtdb_path = gtdb_data_table_path()

    if args.get_table:
        copy_gtdb_table(gtdb_path)
        sys.exit(0)

    _report_gtdb_version(gtdb_path)

    representatives_source = _representatives_source(args)

    if args.get_rank_counts:
        report_unique_taxa_counts_of_all_ranks(
            gtdb_path, representatives_source=representatives_source)
        sys.exit(0)

    if not args.target_taxon:
        return

    if args.get_taxon_counts:
        _report_taxon_counts_or_exit(gtdb_path, args.target_taxon, representatives_source)
        sys.exit(0)

    if args.target_taxon.lower() == "all":
        _write_all(gtdb_path, representatives_source)
        sys.exit(0)

    selection = _select_rows(gtdb_path, args, representatives_source)
    _write_outputs(selection, representatives_source)
    sys.exit(0)


def preflight_checks(args):
    if args.get_taxon_counts and not args.target_taxon:
        report_message("A specific taxon needs to also be provided to the `-t` flag "
                       "in order to use `--get-taxon-counts`.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    if not args.get_rank_counts and not args.get_table and not args.target_taxon:
        report_message("A target taxon needs to be provided to `-t` (or 'all').", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    if args.gtdb_representatives_only and args.refseq_reference_genomes_only:
        print("")
        wprint(color_text("Only one of `--gtdb-representatives-only` or "
                          "`--refseq-reference-genomes-only` can be provided.", "yellow"))
        print("")
        sys.exit(1)


def _representatives_source(args):
    if args.gtdb_representatives_only:
        return "gtdb"
    if args.refseq_reference_genomes_only:
        return "refseq"
    return None


def _select_rows(gtdb_path, args, representatives_source):
    """
    Resolve the target taxon through the shared taxonomy core and return the
    RefGenomeSelection (accessions + metadata rows + resolved rank/canonical). Writing
    is done separately by _write_outputs -- this mirrors the NCBI helper's
    _select_rows / _write_outputs split so the two orchestrators read the same.

    The 'all' taxon is handled by the caller (a bulk dump the core doesn't cover).
    """
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
        report_message(f"Since the input taxon '{args.target_taxon}' occurs at more than 1 rank, "
                        "you'll need to specify which rank is wanted as well before we pull the "
                        "accessions. This can be done with the `-r` parameter.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)
    except TaxonNotFound:
        report_message(f"Input taxon '{args.target_taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)
    except ValueError as err:
        report_message(str(err), "yellow", ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    # can use this if i want to notify about case-insensitive matching (thought i wanted it, but don't feel like it's really needed ATM)
    # if selection.canonical != args.target_taxon:
    #     report_message(f"Matched input '{args.target_taxon}' to GTDB taxon '{selection.canonical}'.",
    #                    "yellow", ii="    ", si="    ", width=100)

    for warning in selection.warnings:
        report_message(warning, "yellow", ii="    ", si="    ", width=100, trailing_newline=True)

    if not selection.accessions:
        report_message(f"No accessions were found for the given target '{selection.canonical}' :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    return selection


def _write_outputs(selection, representatives_source):
    taxon_for_filename = selection.canonical.replace(" ", "-").lower()
    rank = selection.resolved_rank

    if representatives_source:
        if representatives_source == "gtdb":
            suffix = "-" + representatives_source + "-rep"
        else:
            suffix = "-" + representatives_source + "-ref"
    else:
        suffix = ""

    acc_out_filename = f"gtdb-{taxon_for_filename}-{rank}{suffix}-accs.txt"
    tab_out_filename = f"gtdb-{taxon_for_filename}-{rank}{suffix}-metadata.tsv"

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
    if representatives_source == "gtdb":
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
        acc_out_filename = "gtdb-arc-and-bac-" + representatives_source + "-rep-accs.txt"
        tab_out_filename = "gtdb-arc-and-bac-" + representatives_source + "-rep-metadata.tsv"
        table.to_pandas().to_csv(tab_out_filename, sep="\t", index=False)
    else:
        acc_out_filename = "gtdb-arc-and-bac-accs.txt"
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

def _rep_filter_for(representatives_source):
    """The Parquet predicate for a representatives source (or None for no filter)."""
    if representatives_source == "gtdb":
        return ("gtdb_representative", "=", "t")
    if representatives_source == "RefSeq":
        return ("ncbi_refseq_category", "=", "reference genome")
    return None


def _count_at_rank(gtdb_path, rank, taxon, rep_filter=None):
    """Count rows where column `rank` == `taxon` (optionally rep-filtered), via pushdown."""
    filters = [(rank, "=", taxon)]
    if rep_filter:
        filters.append(rep_filter)
    # read a single tiny column; the rank predicate is pushed down to Parquet
    return pq.read_table(gtdb_path, columns=[rank], filters=filters).num_rows


def _count_total(gtdb_path, rep_filter=None):
    filters = [rep_filter] if rep_filter else None
    return pq.read_table(gtdb_path, columns=["ncbi_genbank_assembly_accession"],
                         filters=filters).num_rows


def _report_taxon_counts_or_exit(gtdb_path, taxon, representatives_source):
    """
    Report how many genomes match `taxon` at each rank it occurs at. Resolution goes
    through the shared, case-insensitive resolver (the same one the selection path
    uses), and counts are read straight from the Parquet via predicate pushdown -- so
    this no longer keeps its own taxon-lookup or pandas machinery.
    """
    rep_filter = _rep_filter_for(representatives_source)

    if str(taxon).lower() == "all":
        count = _count_total(gtdb_path)
        print("")
        wprint(f"  There are {count:,} total genomes in the database.")
        print("")
        if representatives_source:
            rep_type = _rep_type_label(representatives_source)
            wprint(color_text(f"  In considering only {rep_type} genomes:", "yellow"))
            print("")
            wprint(f"  There are {_count_total(gtdb_path, rep_filter):,} total "
                   f"{rep_type} genomes in the database.")
            print("")
        return

    # shared resolver: case-insensitive, returns (canonical, [ranks]); raises TaxonNotFound
    try:
        canonical, ranks_found_in = _resolve_ranks(gtdb_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100)
        print("")
        sys.exit(0)

    # can use this if i want to notify about case-insensitive matching (thought i wanted it, but don't feel like it's really needed ATM)
    # if canonical != taxon:
    #     report_message(f"Matched input '{taxon}' to GTDB taxon '{canonical}'.",
    #                    "yellow", ii="    ", si="    ", width=100)

    taxon = canonical

    print("")
    for rank in ranks_found_in:
        count = _count_at_rank(gtdb_path, rank, taxon)
        wprint(f"  The rank '{rank}' has {count:,} {taxon} entries.")
    print("")

    if representatives_source:
        rep_type = _rep_type_label(representatives_source)
        wprint(color_text(f"  In considering only {rep_type} genomes:", "yellow"))
        print("")
        any_rep = False
        for rank in ranks_found_in:
            count = _count_at_rank(gtdb_path, rank, taxon, rep_filter)
            if count:
                any_rep = True
                wprint(f"  The rank '{rank}' has {count:,} {taxon} {rep_type} genome entries.")
                print("")
        if not any_rep:
            wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any "
                   "rank as a representative genome :(", "yellow"))
            print("")
            sys.exit(0)


def _rep_type_label(representatives_source):
    return "refseq reference" if representatives_source == "refseq" else "gtdb representative"


def report_unique_taxa_counts_of_all_ranks(gtdb_path, representatives_source=None):
    """
    Print, for each of the 7 ranks, how many unique taxa exist in the GTDB table.
    Reads the rank columns straight from Parquet and counts distinct values with
    Arrow (matching the NCBI helper's --get-rank-counts), rather than loading pandas.
    """
    ranks = list(RANKS)
    tab = pq.read_table(gtdb_path, columns=ranks)

    print("\n    {:<10} {:}".format("Rank", "Num. Unique Taxa"))
    for rank in ranks:
        n = pc.count_distinct(tab.column(rank)).as_py()
        print("    {:<10} {:}".format(rank, str(n)))
    print("")

    if representatives_source == "gtdb":
        report_message("(The `--gtdb-representatives-only` flag doesn't change these counts: "
                       "every GTDB taxon has a representative genome, so the number of unique "
                       "taxa per rank is the same with or without it.)",
                       "yellow", ii="    ", si="    ", width=100, trailing_newline=True)
    elif representatives_source == "RefSeq":
        rep = pq.read_table(gtdb_path, columns=ranks,
                            filters=[("ncbi_refseq_category", "=", "reference genome")])
        wprint(color_text("  In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print("    {:<10} {:}".format("Rank", "Num. Unique Ref. Taxa"))
        for rank in ranks:
            n = pc.count_distinct(rep.column(rank)).as_py()
            print("    {:<10} {:}".format(rank, str(n)))
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
    out_name = "gtdb-arc-and-bac-metadata.tsv"
    pq.read_table(gtdb_path).to_pandas().to_csv(out_name, sep="\t", index=False)
    print("")
    wprint("GTDB table written to:")
    print(color_text("    " + out_name + "\n"))


if __name__ == "__main__":
    main()
