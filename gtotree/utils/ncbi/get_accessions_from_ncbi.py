import sys
import os
import argparse
from collections import namedtuple
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.messaging import report_message, wprint, color_text, spinner
from gtotree.utils.ncbi.get_ncbi_assembly_data import (get_ncbi_assembly_data,
                                                       ncbi_data_table_path,
                                                       read_date_retrieved)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, AmbiguousTaxon,
                                               find_ranks_for_taxon as _resolve_ranks)
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes


_COLUMNS = [
    "assembly_accession", "organism_name", "taxid", "asm_name", "assembly_level",
    "refseq_category", "checkm_completeness", "checkm_contamination", "genome_size",
] + list(RANKS)

# accepted --assembly-level values mapped to the strings NCBI uses in the table
_ASSEMBLY_LEVELS = {
    "complete": "Complete Genome",
    "chromosome": "Chromosome",
    "scaffold": "Scaffold",
    "contig": "Contig",
}

# result of _select_rows: the metadata rows, a human label for messages, the resolved
# rank, and the canonical taxon string used to build output filenames
_NcbiSelection = namedtuple("_NcbiSelection", ["rows", "label", "rank", "taxon"])


################################################################################

def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate using taxonomy and genomes "
            "from NCBI with GToTree. It primarily returns NCBI accessions and "
            "metadata subsets based on NCBI-taxonomy searches, with optional "
            "filtering by source (RefSeq/GenBank), assembly level, and/or RefSeq 'reference' genomes "
            "only, plus optional dereplication down to one genome per specified rank.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "get-accs-from-ncbi",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt get-accs-from-ncbi -t Nitrospirota --source refseq`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        help=("Target taxon (a name, an NCBI taxid, or 'all'). Not needed with "
              "`--get-rank-counts`."),
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
              "family within the target taxon). Use 'auto' for two ranks finer than the target. "
              "Only applies to a taxon-name search (not a taxid or 'all')."),
        action="store",
    )

    optional.add_argument(
        "-s",
        "--source",
        default="refseq",
        choices=["refseq", "genbank", "both"],
        help=("Specify which section of NCBI to pull from (default: refseq)"),
        action="store",
    )

    optional.add_argument(
        "-a",
        "--assembly-level",
        choices=list(_ASSEMBLY_LEVELS),
        nargs="+",
        help=("Restrict to one or more assembly levels (can be multiple space-separated)"),
        action="store",
    )

    optional.add_argument(
        "-R",
        "--refseq-reference-genomes-only",
        dest="refseq_reference_genomes_only",
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
              "NCBI assembly-summary metadata table."),
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
    get_accessions_from_ncbi(args)

################################################################################


def get_accessions_from_ncbi(args):

    preflight_checks(args)

    # make sure the prepared NCBI Parquet is present, then work against it
    get_ncbi_assembly_data()
    table_path = ncbi_data_table_path()

    _report_ncbi_date(table_path)

    if args.get_table:
        copy_ncbi_table(table_path)
        sys.exit(0)

    if args.get_rank_counts:
        report_unique_taxa_counts_of_all_ranks(
            table_path, source=args.source,
            reps_only=args.refseq_reference_genomes_only)
        sys.exit(0)

    try:
        assembly_levels = parse_assembly_levels(args.assembly_level)
    except ValueError as err:
        wprint(color_text(str(err), "yellow"))
        print("")
        sys.exit(0)

    target = str(args.target_taxon)
    if args.get_taxon_counts and target.lower() != "all" and not target.isdigit():
        _report_taxon_counts_or_exit(table_path, target, args, assembly_levels)
        sys.exit(0)

    selection = _select_rows(table_path, args)
    rows, label = selection.rows, selection.label

    # assembly-level is the only post-filter left; --source scoping already happened
    # up front inside _select_rows (before any dereplication)
    if assembly_levels:
        rows = [r for r in rows if r.get("assembly_level") in assembly_levels]

    if args.get_taxon_counts:
        print("")
        wprint(f"There are {len(rows):,} genome(s) under {label} with any specified "
               "filters.")
        print("")
        sys.exit(0)

    if not rows:
        wprint(color_text(f"No genomes were found under {label} with any specified "
                          "filters.", "yellow"))
        print("")
        sys.exit(0)

    _write_outputs(rows, args, selection)


def preflight_checks(args):

    if args.get_taxon_counts and not args.target_taxon:
        report_message("A specific taxon needs to also be provided to the `-t` flag "
                       "in order to use `--get-taxon-counts`.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    if not args.get_rank_counts and not args.get_table and not args.target_taxon:
        report_message("A target taxon needs to be provided to `-t` (a name, a taxid, or 'all').", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)


def _count_at_rank(table_path, rank, taxon, prefixes=None, reps_only=False,
                   assembly_levels=None):
    """
    Count rows where column `rank` == `taxon`, with the set POOL filters applied
    (source prefix, RefSeq-reference-only, assembly level). --derep-rank is
    intentionally NOT applied: counts report how many genomes MATCH the filters, not
    how many survive dereplication (a pull-time reduction). Reads via pushdown where
    possible, then applies the prefix / level filters in Arrow.
    """
    filters = [(rank, "=", taxon)]
    if reps_only:
        filters.append(("refseq_category", "=", "reference genome"))
    if assembly_levels:
        filters.append(("assembly_level", "in", set(assembly_levels)))

    cols = [rank, "assembly_accession"]
    tab = pq.read_table(table_path, columns=cols, filters=filters)

    if prefixes:
        tab = tab.filter(_prefix_mask(tab.column("assembly_accession"), prefixes))

    return tab.num_rows


def _prefix_mask(acc_col, prefixes):
    mask = None
    for p in prefixes:
        m = pc.starts_with(acc_col, p)
        mask = m if mask is None else pc.or_(mask, m)
    return mask


def _report_taxon_counts_or_exit(table_path, taxon, args, assembly_levels):
    """
    Report how many genomes match `taxon` at each rank it occurs at, matching the GTDB
    helper's format: a primary per-rank block for the base pool (scoped by --source and
    --assembly-level), then if --refseq-reference-genomes-only is set a separate
    "In considering only RefSeq reference genomes:" block, like GTDB's reps block.

    The wording is explicit about WHICH filters each block reflects: the primary block
    reflects --source and --assembly-level (but not the reference-genome filter, which
    is applied only in the second block), so the two numbers aren't confused.
    """
    prefixes = _source_prefixes(args.source)

    # a short human description of the filters folded into the PRIMARY block, so the
    # count line says what it actually reflects rather than a vague "any filters"
    scope_note = _counts_scope_note(args, assembly_levels)

    try:
        canonical, ranks_found_in = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    # can use this if i want to notify about case-insensitive matching (thought i wanted it, but don't feel like it's really needed ATM)
    # if canonical != taxon:
    #     report_message(f"Matched input '{taxon}' to NCBI taxon '{canonical}'.", "yellow",
    #                    ii="    ", si="    ", width=100)
    taxon = canonical

    print("")
    for rank in ranks_found_in:
        count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                               assembly_levels=assembly_levels)
        report_message(f"The rank '{rank}' has {count:,} {taxon} entries{scope_note}.", color=None,
                       ii="    ", si="    ", width=100, newline=False, trailing_newline=True)

    if args.refseq_reference_genomes_only:
        report_message("Of those, in considering only RefSeq reference genomes:", "yellow",
                       ii="    ", si="    ", width=100, newline=False, trailing_newline=True)
        any_rep = False
        for rank in ranks_found_in:
            count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                                   reps_only=True, assembly_levels=assembly_levels)
            if count:
                any_rep = True
                report_message(f"The rank '{rank}' has {count:,} {taxon} RefSeq reference genome entries.", color=None,
                               ii="    ", si="    ", width=100, newline=False, trailing_newline=True)
        if not any_rep:
            report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank as a RefSeq reference genome :(", "yellow",
                           ii="    ", si="    ", width=100, newline=False, trailing_newline=True)
            sys.exit(0)


def _counts_scope_note(args, assembly_levels):
    bits = []
    if args.source == "refseq":
        bits.append("in refseq")
    elif args.source == "genbank":
        bits.append("in genbank")
    if assembly_levels:
        levels = ", ".join(sorted(assembly_levels))
        bits.append(f"at assembly level {levels}")
    if not bits:
        return ""
    return " (" + ", ".join(bits) + ")"


def _source_prefixes(source):
    """Accession prefixes for a --source value (None means no restriction)."""
    if source == "refseq":
        return ("GCF_",)
    if source == "genbank":
        return ("GCA_",)
    return None


def _select_rows(table_path, args):
    """
    Resolve the target and return an _NcbiSelection(rows, label, rank, taxon) where
    `rank` is the resolved rank for a taxon-name search (None for 'all'/taxid, which
    don't resolve to a single rank) and `taxon` is the canonical name used for output
    filenames. Three modes:
      - 'all'          -> every genome (optionally reps-only)
      - a numeric taxid-> lineage-taxid lookup
      - a taxon name   -> the shared select_ref_genomes core (honours --derep-rank)
    Only the taxon-name mode dereplicates; taxid and 'all' don't resolve to a single
    rank, so derep (a one-per-finer-rank operation) doesn't apply there.

    --source scoping (refseq/genbank -> GCF_/GCA_ prefix) is applied up front in every
    mode: inside the core for the taxon path, and at read time for all/taxid.
    """
    target = str(args.target_taxon)
    reps_only = args.refseq_reference_genomes_only
    prefixes = _source_prefixes(args.source)

    if target.lower() == "all":
        filters = [("refseq_category", "=", "reference genome")] if reps_only else None
        tab = pq.read_table(table_path, columns=_COLUMNS, filters=filters)
        rows = _apply_source_prefix(tab.to_pylist(), prefixes)
        return _NcbiSelection(rows, "all genomes", None, "all")

    if target.isdigit():
        rows = _select_by_taxid(table_path, target, reps_only=reps_only)
        rows = _apply_source_prefix(rows, prefixes)
        return _NcbiSelection(rows, f"taxid {target}", None, f"taxid-{target}")

    # taxon name -> shared core (taxon resolution + optional dereplication).
    # --source scopes the candidate pool BY ACCESSION PREFIX up front (inside the
    # core, before dereplication), so a best-per-rank pick is made within the
    # requested pool rather than dropped afterward.
    try:
        selection = select_ref_genomes(
            table_path, "ncbi", target,
            target_rank=args.target_rank, derep_rank=args.derep_rank,
            reps_only=reps_only,
            accession_prefixes=_source_prefixes(args.source))
    except AmbiguousTaxon as e:
        report_message(f"Since the input taxon '{e.taxon}' occurs at more than 1 rank, "
                       "you'll need to specify which rank is wanted as well before we pull the "
                       "accessions. This can be done with the `-r` parameter, or you can try passing "
                       "the NCBI taxid to `-t` instead.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)
    except TaxonNotFound:
        report_message(f"Input taxon '{target}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)
    except ValueError as err:
        report_message(str(err), "yellow", ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    # can use this if i want to notify about case-insensitive matching (thought i wanted it, but don't feel like it's really needed ATM)
    # if selection.canonical != target:
    #     report_message(f"Matched input '{target}' to NCBI taxon '{selection.canonical}'.",
    #                    "yellow", ii="    ", si="    ", width=100)

    for warning in selection.warnings:
        report_message(warning, "yellow", ii="    ", si="    ", width=100, trailing_newline=True)

    if not selection.accessions:
        report_message(f"No accessions were found for the given target '{selection.canonical}' :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    return _NcbiSelection(selection.rows,
                          f"{selection.resolved_rank} '{selection.canonical}'",
                          selection.resolved_rank, selection.canonical)


def _select_by_taxid(table_path, taxid, reps_only=False):
    """
    Select rows whose lineage carries `taxid` at any rank (checking each rank's
    *_taxid column, then the row's own taxid), returning full metadata rows.
    """
    base_filters = [("refseq_category", "=", "reference genome")] if reps_only else []
    for rank in list(RANKS) + [None]:
        col = f"{rank}_taxid" if rank else "taxid"
        tab = pq.read_table(table_path, columns=_COLUMNS,
                            filters=base_filters + [(col, "=", str(taxid))])
        if tab.num_rows:
            return tab.to_pylist()
    return []


def parse_assembly_levels(value):
    if not value:
        return []
    if isinstance(value, str):
        value = value.split(",")
    parts = [v.strip().lower() for v in value if str(v).strip()]
    unknown = [p for p in parts if p not in _ASSEMBLY_LEVELS]
    if unknown:
        raise ValueError(
            f"unrecognised --assembly-level value(s): {', '.join(unknown)}. "
            f"Choose from: {', '.join(_ASSEMBLY_LEVELS)}")
    return [_ASSEMBLY_LEVELS[p] for p in parts]


def _apply_source_prefix(rows, prefixes):
    """Scope rows to a source by accession prefix (None -> unfiltered)."""
    if not prefixes:
        return rows
    prefixes = tuple(prefixes)
    return [r for r in rows
            if str(r.get("assembly_accession") or "").startswith(prefixes)]


def report_unique_taxa_counts_of_all_ranks(table_path, source="refseq", reps_only=False):
    """
    Print, for each of the 7 ranks, how many unique taxa exist in the NCBI table,
    scoped to `source` (refseq -> GCF_ only, genbank -> GCA_ only, both -> no filter).
    If reps_only, also print counts among RefSeq reference genomes.
    """
    ranks = list(RANKS)
    tab = pq.read_table(table_path, columns=ranks + ["assembly_accession"])
    tab = _filter_table_by_source(tab, source)

    label = {"refseq": "refseq", "genbank": "genbank", "both": "all"}.get(source, source)
    print("\n    {:<10} {:}".format("Rank", f"Num. Unique Taxa ({label})"))
    for rank in ranks:
        n = pc.count_distinct(tab.column(rank)).as_py()
        print("    {:<10} {:}".format(rank, str(n)))
    print("")

    if reps_only:
        rep = pq.read_table(table_path, columns=ranks,
                            filters=[("refseq_category", "=", "reference genome")])
        wprint(color_text("In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print("    {:<10} {:}".format("Rank", "Num. Unique Ref. Taxa"))
        for rank in ranks:
            n = pc.count_distinct(rep.column(rank)).as_py()
            print("    {:<10} {:}".format(rank, str(n)))
        print("")


def _filter_table_by_source(tab, source):
    """Arrow-table variant of the source prefix filter (for the counts path)."""
    if source in ("refseq", "genbank"):
        prefix = "GCF_" if source == "refseq" else "GCA_"
        tab = tab.filter(pc.starts_with(tab.column("assembly_accession"), prefix))
    return tab


def _report_ncbi_date(table_path):
    date_str = read_date_retrieved(os.path.dirname(table_path))
    print("\n    Date NCBI assembly-data retrieved: " + date_str)


def copy_ncbi_table(table_path):
    out_name = "ncbi-assembly-summary-metadata.tsv"
    print("")
    with spinner("Writing NCBI table...", "", clear_on_done=True):
        pq.read_table(table_path).to_pandas().to_csv(out_name, sep="\t", index=False)
    wprint("  NCBI table written to:")
    print(color_text("      " + out_name + "\n"))


def _write_outputs(rows, args, selection):

    taxon_for_filename = selection.taxon.replace(" ", "-").replace("/", "-").lower()

    rank_bit = f"-{selection.rank}" if selection.rank else ""

    suffix_bits = []
    if args.refseq_reference_genomes_only:
        suffix_bits.append("refseq-ref")
    elif args.source != "both":
        suffix_bits.append(args.source)
    suffix = ("-" + "-".join(suffix_bits)) if suffix_bits else ""

    acc_out = f"ncbi-{taxon_for_filename}{rank_bit}{suffix}-accs.txt"
    tab_out = f"ncbi-{taxon_for_filename}{rank_bit}{suffix}-metadata.tsv"

    _write_metadata_tsv(rows, tab_out)

    accs = [r.get("assembly_accession") for r in rows if r.get("assembly_accession")]
    with open(acc_out, "w") as out:
        for acc in accs:
            out.write(acc + "\n")

    print("")
    wprint(f"Wrote {len(accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out))
    print("")


def _write_metadata_tsv(rows, out_filename):
    """Write selected genome rows to a TSV, accession + ranks first then the rest."""
    if not rows:
        open(out_filename, "w").close()
        return
    first = ["assembly_accession"] + list(RANKS)
    seen = set(first)
    header = [c for c in first if c in rows[0]] + [c for c in rows[0] if c not in seen]
    with open(out_filename, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in header) + "\n")


if __name__ == "__main__":
    main()
