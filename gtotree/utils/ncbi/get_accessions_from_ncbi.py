import sys
import os
import argparse
from collections import namedtuple
import pyarrow.parquet as pq # type: ignore
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc.messaging import report_message, wprint, color_text, spinner
from gtotree.utils.misc.general import (write_table_tsv, write_accessions,
                                        atomic_write_text)
from gtotree.utils.ncbi.get_ncbi_assembly_data import (get_ncbi_assembly_data,
                                                       ncbi_data_table_path,
                                                       read_date_retrieved)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, AmbiguousTaxon,
                                               find_ranks_for_taxon as _resolve_ranks)
from gtotree.utils.taxonomy.tax_derep import (select_ref_genomes, select_all_domains,
                                              resolve_derep_rank)
from gtotree.utils.taxonomy.tax_counts import (representatives_filter, count_genomes,
                                               derep_size, rank_counts,
                                               render_rank_count_table)


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
            epilog="Ex. usage: `gtt get-accs-from-ncbi -t Alteromonas`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-w",
        "--wanted-ref-tax",
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
              "family within the target taxon). Use 'auto' for two ranks finer than the target."),
        action="store",
    )

    optional.add_argument(
        "--source",
        type=str.lower,
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
        help=("Provide this flag along with a specified taxon to `-w` to see how many "
              "genomes match the set parameters. If `--derep-rank` is also set, the "
              "number of genomes following dereplication is reported too."),
    )

    optional.add_argument(
        "--get-rank-counts",
        action="store_true",
        help=("Provide this flag to see counts of how many unique taxa there are for each rank. "
              "By itself, that'd be the whole database, but it can also be combined with "
              "`-w` and `--derep-rank`."),
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

    try:
        assembly_levels = parse_assembly_levels(args.assembly_level)
    except ValueError as err:
        wprint(color_text(str(err), "yellow"))
        print("")
        sys.exit(0)

    target = str(args.wanted_ref_tax) if args.wanted_ref_tax else ""
    named_taxon = bool(target) and target.lower() != "all" and not target.isdigit()

    if args.get_rank_counts:
        if named_taxon:
            _report_rank_counts_for_taxon_or_exit(table_path, target, args,
                                                  assembly_levels)
        else:
            report_unique_taxa_counts_of_all_ranks(
                table_path, source=args.source,
                reps_only=args.refseq_reference_genomes_only,
                assembly_levels=assembly_levels)
        sys.exit(0)

    if args.get_taxon_counts and named_taxon:
        _report_taxon_counts_or_exit(table_path, target, args, assembly_levels)
        sys.exit(0)

    selection = _select_rows(table_path, args, assembly_levels)
    rows, label = selection.rows, selection.label

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

    if args.get_taxon_counts and not args.wanted_ref_tax:
        report_message("A specific taxon needs to also be provided to the `-w` flag "
                       "in order to use `--get-taxon-counts`.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    if not args.get_rank_counts and not args.get_table and not args.wanted_ref_tax:
        report_message("A target taxon needs to be provided to `-w` (a name, a taxid, or 'all').", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    check_derep_rank_is_applicable(args)


def check_derep_rank_is_applicable(args):
    """
    A taxid's rank is known only after the lookup, and the taxid path doesn't run
    through the selection core, so there is no resolved rank to group beneath. So it's
    not applicable with --derep-rank
    """
    target = str(args.wanted_ref_tax or "")
    derep_rank = getattr(args, "derep_rank", "off")

    if derep_rank in (None, "off", "none", "None"):
        return

    if not target.isdigit():
        return

    report_message(
        f"`--derep-rank` can't be applied with a taxid. Pass the taxon "
        f"name instead if you also want to dereplicate",
        "yellow", ii="    ", si="    ", width=100, trailing_newline=True)
    sys.exit(0)


def _derep_is_on(args):
    """True when --derep-rank asks for actual dereplication."""
    return getattr(args, "derep_rank", "off") not in (None, "off", "none", "None")


def _count_at_rank(table_path, rank, taxon, prefixes=None, reps_only=False,
                   assembly_levels=None):
    """
    Count rows where column `rank` == `taxon`, with the set POOL filters applied
    (source prefix, RefSeq-reference-only, assembly level). This is the number of
    genomes that MATCH. What a --derep-rank pull would return is reported separately
    by _derep_count_at_rank(), since the two are different questions.
    """
    return count_genomes(table_path, "ncbi", rank=rank, taxon=taxon,
                         rep_filter=_rep_filter(reps_only),
                         accession_prefixes=prefixes,
                         assembly_levels=assembly_levels)


def _rep_filter(reps_only):
    """The RefSeq-reference-genome predicate, or None."""
    return representatives_filter("ncbi", "refseq" if reps_only else None)


def _derep_count_at_rank(table_path, rank, taxon, derep_rank, prefixes=None,
                         reps_only=False, assembly_levels=None):
    """How many genomes survive dereplication at `derep_rank`, under the same pool."""
    return derep_size(table_path, "ncbi", rank, taxon, derep_rank,
                      rep_filter=_rep_filter(reps_only),
                      accession_prefixes=prefixes,
                      assembly_levels=assembly_levels)


def _derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                reps_only=False):
    """
    The "...and dereplicated at X, that would be N" line for one rank, or None when
    dereplication is off. Any advisory from resolving 'auto' is returned alongside so
    the caller can print it (e.g. auto on a species target turns derep off).

    Returns (line_or_None, warnings).
    """
    derep_rank = getattr(args, "derep_rank", "off")
    if derep_rank in (None, "off", "none", "None"):
        return None, []

    try:
        effective, warnings = resolve_derep_rank(rank, derep_rank)
    except ValueError as err:
        # an explicit rank coarser than this one -- can't dereplicate to it
        return str(err), []

    if effective is None:
        return None, warnings

    n = _derep_count_at_rank(table_path, rank, taxon, effective, prefixes=prefixes,
                             reps_only=reps_only, assembly_levels=assembly_levels)
    return (f"Dereplicated at '{effective}', that would be {n:,} genome(s).",
            warnings)


def _report_rank_counts_for_taxon_or_exit(table_path, taxon, args, assembly_levels):
    """
    `--get-rank-counts` scoped to a taxon, how many unique taxa sit under it at each
    rank from its own rank down
    """
    prefixes = _source_prefixes(args.source)
    reps_only = args.refseq_reference_genomes_only
    scope_note = _counts_scope_note(args, assembly_levels)

    try:
        canonical, ranks_found_in = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    for rank in ranks_found_in:
        total = _count_at_rank(table_path, rank, canonical, prefixes=prefixes,
                               reps_only=reps_only, assembly_levels=assembly_levels)
        rows = rank_counts(table_path, "ncbi", scope_rank=rank, scope_taxon=canonical,
                           rep_filter=_rep_filter(reps_only),
                           accession_prefixes=prefixes,
                           assembly_levels=assembly_levels)

        print("")
        report_message(f"The rank '{rank}' has {total:,} {canonical} entries{scope_note}.",
                       color=None, ii="    ", si="    ", width=100, newline=False)
        print("")
        print(render_rank_count_table(
            rows, count_header=f"Num. Unique Taxa under '{canonical}'"))

    print("")
    report_message("Each count above is also how many genomes `--derep-rank <rank>` "
                   "would return, since dereplication keeps one genome per unique "
                   "taxon at that rank.", "yellow",
                   ii="    ", si="    ", width=100, newline=False, trailing_newline=True)


def _report_taxon_counts_or_exit(table_path, taxon, args, assembly_levels):
    """
    Report how many genomes match `taxon` at each rank it occurs at

    """
    prefixes = _source_prefixes(args.source)

    scope_note = _counts_scope_note(args, assembly_levels)

    try:
        canonical, ranks_found_in = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    taxon = canonical

    print("")
    for rank in ranks_found_in:
        count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                               assembly_levels=assembly_levels)
        report_message(f"The rank '{rank}' has {count:,} {taxon} entries{scope_note}.", color=None,
                       ii="    ", si="    ", width=100, newline=False, trailing_newline=True)
        _report_derep_note(table_path, rank, taxon, args, prefixes, assembly_levels)

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
                _report_derep_note(table_path, rank, taxon, args, prefixes,
                                   assembly_levels, reps_only=True)
        if not any_rep:
            report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank as a RefSeq reference genome :(", "yellow",
                           ii="    ", si="    ", width=100, newline=False, trailing_newline=True)
            sys.exit(0)


def _report_derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                       reps_only=False):
    """Print the dereplicated-count line (and any 'auto' advisory) for one rank."""
    line, warnings = _derep_note(table_path, rank, taxon, args, prefixes,
                                 assembly_levels, reps_only=reps_only)
    if line:
        report_message(line, color=None, ii="      ", si="      ", width=100,
                       newline=False, trailing_newline=True)
    for warning in warnings:
        report_message(warning, "yellow", ii="      ", si="      ", width=100,
                       newline=False, trailing_newline=True)


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


def _select_rows(table_path, args, assembly_levels=None):
    """
    Resolve the target and return an _NcbiSelection(rows, label, rank, taxon) where
    `rank` is the resolved rank for a taxon-name search
    """
    target = str(args.wanted_ref_tax)
    reps_only = args.refseq_reference_genomes_only
    prefixes = _source_prefixes(args.source)

    if target.lower() == "all":
        if _derep_is_on(args):
            selection = select_all_domains(
                table_path, "ncbi", derep_rank=args.derep_rank, reps_only=reps_only,
                accession_prefixes=prefixes, assembly_levels=assembly_levels)
            report_message(f"Dereplicating within each domain "
                           f"({', '.join(selection.domains)}).", "yellow",
                           ii="    ", si="    ", width=100, trailing_newline=True)
            for warning in selection.warnings:
                report_message(warning, "yellow", ii="    ", si="    ", width=100,
                               trailing_newline=True)
            label = ("all genomes (dereplicated within each domain: "
                     + ", ".join(selection.domains) + ")")
            return _NcbiSelection(selection.rows, label, None, "all")

        filters = [("refseq_category", "=", "reference genome")] if reps_only else []
        if assembly_levels:
            filters.append(("assembly_level", "in", set(assembly_levels)))
        tab = pq.read_table(table_path, columns=_COLUMNS, filters=filters or None)
        rows = _apply_source_prefix(tab.to_pylist(), prefixes)
        return _NcbiSelection(rows, "all genomes", None, "all")

    if target.isdigit():
        rows = _select_by_taxid(table_path, target, reps_only=reps_only,
                                assembly_levels=assembly_levels)
        rows = _apply_source_prefix(rows, prefixes)
        return _NcbiSelection(rows, f"taxid {target}", None, f"taxid-{target}")

    try:
        selection = select_ref_genomes(
            table_path, "ncbi", target,
            target_rank=args.target_rank, derep_rank=args.derep_rank,
            reps_only=reps_only,
            accession_prefixes=_source_prefixes(args.source),
            assembly_levels=assembly_levels)
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

    for warning in selection.warnings:
        report_message(warning, "yellow", ii="    ", si="    ", width=100, trailing_newline=True)

    if not selection.accessions:
        report_message(f"No accessions were found for the given target '{selection.canonical}' :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    return _NcbiSelection(selection.rows,
                          f"{selection.resolved_rank} '{selection.canonical}'",
                          selection.resolved_rank, selection.canonical)


def _select_by_taxid(table_path, taxid, reps_only=False, assembly_levels=None):
    """
    Select rows whose lineage carries `taxid` at any rank (checking each rank's
    *_taxid column, then the row's own taxid), returning full metadata rows.
    """
    base_filters = [("refseq_category", "=", "reference genome")] if reps_only else []
    if assembly_levels:
        base_filters = base_filters + [("assembly_level", "in", set(assembly_levels))]
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


def report_unique_taxa_counts_of_all_ranks(table_path, source="refseq", reps_only=False,
                                           assembly_levels=None):
    """
    Print, for each of the 7 ranks, how many unique taxa exist in the NCBI table,
    scoped to `source` (refseq -> GCF_ only, genbank -> GCA_ only, both -> no filter).
    If reps_only, also print counts among RefSeq reference genomes.
    """
    prefixes = _source_prefixes(source)
    label = {"refseq": "refseq", "genbank": "genbank", "both": "all"}.get(source, source)

    rows = rank_counts(table_path, "ncbi", accession_prefixes=prefixes,
                       assembly_levels=assembly_levels)
    print("")
    print(render_rank_count_table(rows, count_header=f"Num. Unique Taxa ({label})"))
    print("")

    if reps_only:
        rep_rows = rank_counts(table_path, "ncbi", accession_prefixes=prefixes,
                               rep_filter=_rep_filter(True),
                               assembly_levels=assembly_levels)
        wprint(color_text("In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print(render_rank_count_table(rep_rows, count_header="Num. Unique Ref. Taxa"))
        print("")


def _report_ncbi_date(table_path):
    date_str = read_date_retrieved(os.path.dirname(table_path))
    print("\n    Date NCBI assembly-data retrieved: " + date_str)


def copy_ncbi_table(table_path):
    out_name = "ncbi-assembly-summary-metadata.tsv"
    print("")
    with spinner("Writing NCBI table...", "", clear_on_done=True):
        write_table_tsv(pq.read_table(table_path), out_name)
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
    write_accessions(acc_out, accs)

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
        atomic_write_text(out_filename, lambda f: None)
        return
    first = ["assembly_accession"] + list(RANKS)
    seen = set(first)
    header = [c for c in first if c in rows[0]] + [c for c in rows[0] if c not in seen]

    def _write(out):
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in header) + "\n")

    atomic_write_text(out_filename, _write)


if __name__ == "__main__":
    main()
