import sys
import os
import argparse
import pyarrow.parquet as pq # type: ignore

from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc.messaging import wprint, color_text, report_message, spinner
from gtotree.utils.misc.general import (write_table_tsv, write_accessions,
                                        atomic_write_text)
from gtotree.utils.gtdb.get_gtdb_data import (get_gtdb_data, gtdb_data_table_path,
                                              report_gtdb_version_info as _read_gtdb_version_info)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, AmbiguousTaxon,
                                               CrossDomainTaxon,
                                               find_ranks_for_taxon as _resolve_ranks)
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes, select_all_domains
from gtotree.utils.taxonomy.tax_counts import (representatives_filter, count_genomes,
                                               rank_counts,
                                               render_rank_count_table)
from gtotree.utils.taxonomy.tax_targets import (is_all_target,
                                                unassigned_domain_summary)
from gtotree.utils.taxonomy.get_accs_shared import (PoolSpec,
                                                    add_common_get_accs_args,
                                                    all_derep_size,
                                                    derep_note as _shared_derep_note,
                                                    is_derep_on, scoped_counts_note)


_RANK_COLUMNS = list(RANKS)


################################################################################
from gtotree.utils.taxonomy.exclusion_list import (load_exclusion_cores,
                                                   filter_table_by_exclusion,
                                                   exclusion_warning)


def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate retrieving accessions from the "
            "Genome Taxonomy Database (gtdb.ecogenomic.org). It returns NCBI "
            "accessions and GTDB metadata subsets based on GTDB-taxonomy searches, "
            "with optional filtering and dereplication settings.")

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
            epilog="Ex. usage: `gtt get-accs-from-gtdb -w Archaea --gtdb-representatives-only`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    add_common_get_accs_args(
        required, optional, "GTDB")

    optional.add_argument(
        "-G",
        "--gtdb-representatives-only",
        action="store_true",
        help=("Pull only genomes designated as GTDB species representatives."),
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

    # the counts flags preview what a pull would return, so they honor the same
    # exclusion list the pull itself would apply
    exclude_cores = load_exclusion_cores(getattr(args, "exclusion_list", None))

    preflight_checks(args)

    # make sure the prepared GTDB Parquet is present, then work against it
    get_gtdb_data()
    gtdb_path = gtdb_data_table_path()

    if args.get_table:
        copy_gtdb_table(gtdb_path)
        sys.exit(0)

    _report_gtdb_version(gtdb_path)

    representatives_source = _representatives_source(args)

    named_taxon = bool(args.wanted_ref_tax) and not is_all_target(args.wanted_ref_tax)

    if args.get_rank_counts:
        if named_taxon:
            _report_rank_counts_for_taxon_or_exit(
                gtdb_path, args.wanted_ref_tax, representatives_source,
                                                  exclude_cores=exclude_cores)
        else:
            report_unique_taxa_counts_of_all_ranks(
                gtdb_path, representatives_source=representatives_source,
                scoped_to_all=is_all_target(args.wanted_ref_tax),
                                            exclude_cores=exclude_cores)
        sys.exit(0)

    if not args.wanted_ref_tax:
        return

    if args.get_taxon_counts:
        _report_taxon_counts_or_exit(gtdb_path, args.wanted_ref_tax,
                                     representatives_source, args,
                                     exclude_cores=exclude_cores)
        sys.exit(0)

    if is_all_target(args.wanted_ref_tax):
        if _derep_is_on(args):
            _write_all_dereplicated(gtdb_path, args, representatives_source)
        else:
            _write_all(gtdb_path, representatives_source,
                       exclude_cores=load_exclusion_cores(
                           getattr(args, "exclusion_list", None)))
        sys.exit(0)

    selection = _select_rows(gtdb_path, args, representatives_source)
    _write_outputs(selection, representatives_source)
    sys.exit(0)


def preflight_checks(args):
    if args.get_taxon_counts and not args.wanted_ref_tax:
        report_message("A specific taxon needs to also be provided to the `-w` flag "
                       "in order to use `--get-taxon-counts`.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    if not args.get_rank_counts and not args.get_table and not args.wanted_ref_tax:
        report_message("A target taxon needs to be provided to `-w` (or 'all').", "yellow",
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
    RefGenomeSelection (accessions + metadata rows + resolved rank/canonical)
    """
    reps_only = representatives_source is not None

    try:
        selection = select_ref_genomes(
            gtdb_path, "gtdb", args.wanted_ref_tax,
            target_rank=args.target_rank, derep_rank=args.derep_rank,
            reps_only=reps_only,
            exclude_cores=load_exclusion_cores(
                getattr(args, "exclusion_list", None)),
            target_domain=getattr(args, "target_domain", None))
    except AmbiguousTaxon:
        report_message(f"Since the input taxon '{args.wanted_ref_tax}' occurs at more than 1 rank, "
                        "you'll need to specify which rank is wanted as well before we pull the "
                        "accessions. This can be done with the `-r` parameter.", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)
    except CrossDomainTaxon as e:
        report_message(f"The input taxon '{e.taxon}' occurs in more than one domain "
                       f"({', '.join(e.domains_found)}), so pulling on the name alone "
                       "would mix genomes from different domains. Specify which domain "
                       "is wanted with `--target-domain` "
                       f"(e.g. `--target-domain {e.domains_found[0]}`).", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)
    except TaxonNotFound:
        report_message(f"Input taxon '{args.wanted_ref_tax}' doesn't seem to exist at any rank :(", "yellow",
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

    write_accessions(acc_out_filename, selection.accessions)
    _write_metadata_tsv(selection.rows, tab_out_filename)

    print("")
    wprint(f"Wrote {len(selection.accessions):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def _derep_is_on(args):
    """True when --derep-rank asks for actual dereplication."""
    return is_derep_on(getattr(args, "derep_rank", "off"))


def _pool(gtdb_path, representatives_source=None):
    """
    The genomes this invocation is working within, as a PoolSpec
    """
    return PoolSpec(gtdb_path, "gtdb",
                    rep_filter=_rep_filter_for(representatives_source),
                    label="GTDB", taxon_flag="-w")


def _write_all_dereplicated(gtdb_path, args, representatives_source):
    """`-w all` WITH --derep-rank: one selection per domain, merged."""
    reps_only = representatives_source is not None

    try:
        selection = select_all_domains(gtdb_path, "gtdb", derep_rank=args.derep_rank,
                                       reps_only=reps_only,
                                       exclude_cores=load_exclusion_cores(
                                           getattr(args, "exclusion_list", None)))
    except ValueError as err:
        report_message(str(err), "yellow", ii="    ", si="    ", width=100,
                       trailing_newline=True)
        sys.exit(0)

    report_message(f"Dereplicating within each domain "
                   f"({', '.join(selection.domains)}).", "yellow",
                   ii="    ", si="    ", width=100, trailing_newline=False)

    _report_unassigned_domains(getattr(selection, "unassigned", None))

    for warning in selection.warnings:
        report_message(warning, "yellow", ii="    ", si="    ", width=100,
                       trailing_newline=True)

    if not selection.accessions:
        report_message("No accessions were found :(", "yellow",
                       ii="    ", si="    ", width=100, trailing_newline=True)
        sys.exit(0)

    suffix = f"-{representatives_source}-rep" if representatives_source else ""
    acc_out_filename = f"gtdb-arc-and-bac{suffix}-accs.txt"
    tab_out_filename = f"gtdb-arc-and-bac{suffix}-metadata.tsv"

    write_accessions(acc_out_filename, selection.accessions)
    _write_metadata_tsv(selection.rows, tab_out_filename)

    print("")
    wprint(f"Wrote {len(selection.accessions):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def _write_all(gtdb_path, representatives_source, exclude_cores=None):
    """
    Bulk dump of every accession (optionally reps-only), plus a metadata TSV.

    The representatives predicate comes from _rep_filter_for()

    This path reads the asset directly rather than going through the selection core,
    so it applies any `--exclusion-list` itself. Nothing is dereplicated here, so
    the listed genomes are simply absent from the dump.
    """
    rep_filter = _rep_filter_for(representatives_source)
    filt = [(rep_filter[0], "=", rep_filter[1])] if rep_filter else None

    if filt:
        table = pq.read_table(gtdb_path, filters=filt)
    else:
        table = pq.read_table(gtdb_path)

    table, n_excluded = filter_table_by_exclusion(
        table, "ncbi_genbank_assembly_accession", exclude_cores)
    note = exclusion_warning(n_excluded)
    if note:
        report_message(note, "yellow", ii="    ", si="    ", width=100,
                       trailing_newline=True)

    accessions = table.column("ncbi_genbank_assembly_accession").to_pylist()

    if representatives_source:
        acc_out_filename = "gtdb-arc-and-bac-" + representatives_source + "-rep-accs.txt"
        tab_out_filename = "gtdb-arc-and-bac-" + representatives_source + "-rep-metadata.tsv"
        write_table_tsv(table, tab_out_filename)
    else:
        acc_out_filename = "gtdb-arc-and-bac-accs.txt"
        tab_out_filename = None

    write_accessions(acc_out_filename, accessions)

    print("")
    wprint(f"Wrote {len(accessions):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    if tab_out_filename:
        wprint("Associated taxonomy and metadata written to:")
        wprint("  " + color_text(tab_out_filename))
        print("")


def _report_unassigned_domains(summary):
    """
    Say what an 'all' pull leaves behind, if anything
    """
    message = summary.message("GTDB") if summary else None
    if message:
        report_message(message, "yellow", ii="    ", si="    ", width=100,
                       trailing_newline=True)


def _write_metadata_tsv(rows, out_filename):
    """Write selected genome rows to a TSV, columns in the asset's natural order."""
    if not rows:
        atomic_write_text(out_filename, lambda f: None)
        return
    # preserve a stable, readable column order: accession + ranks first, then the rest
    first = ["ncbi_genbank_assembly_accession"] + list(RANKS)
    seen = set(first)
    header = [c for c in first if c in rows[0]] + [c for c in rows[0] if c not in seen]

    def _write(out):
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in header) + "\n")

    atomic_write_text(out_filename, _write)


################################################################################
# counts + table helpers (read straight from the Parquet, no selection core needed)
################################################################################

def _rep_filter_for(representatives_source):
    """
    The Parquet predicate for a representatives source (or None for no filter)
    """
    if representatives_source == "gtdb":
        return representatives_filter("gtdb", "source")
    if representatives_source == "refseq":
        return representatives_filter("gtdb", "refseq")
    return None


def _count_at_rank(gtdb_path, rank, taxon, rep_filter=None, exclude_cores=None):
    """Count rows where column `rank` == `taxon` (optionally rep-filtered), via pushdown."""
    return count_genomes(gtdb_path, "gtdb", rank=rank, taxon=taxon,
                         rep_filter=rep_filter,
                         exclude_cores=exclude_cores)


def _count_total(gtdb_path, rep_filter=None, exclude_cores=None):
    return count_genomes(gtdb_path, "gtdb", rep_filter=rep_filter,
                         exclude_cores=exclude_cores)


def _derep_note(gtdb_path, rank, taxon, args, rep_filter=None, exclude_cores=None):
    """
    The "...dereplicated at X, that would be N" line for one rank, or None when
    dereplication is off. Returns (line_or_None, warnings)
    """
    pool = PoolSpec(gtdb_path, "gtdb", rep_filter=rep_filter, label="GTDB",
                    taxon_flag="-w",
                    exclude_cores=exclude_cores)
    return _shared_derep_note(pool, rank, taxon, getattr(args, "derep_rank", "off"))


def _report_derep_note(gtdb_path, rank, taxon, args, rep_filter=None, exclude_cores=None):
    """Print the dereplicated-count line (and any 'auto' advisory) for one rank."""
    line, warnings = _derep_note(gtdb_path, rank, taxon, args, rep_filter=rep_filter,
                                 exclude_cores=exclude_cores)
    if line:
        wprint("    " + line)
    for warning in warnings:
        report_message(warning, "yellow", ii="      ", si="      ", width=100,
                       newline=False, trailing_newline=True)


def _report_rank_counts_for_taxon_or_exit(gtdb_path, taxon, representatives_source, exclude_cores=None):
    """
    `--get-rank-counts` scoped to a taxon, how many unique taxa sit under it at each
    rank from its own rank down
    """
    rep_filter = _rep_filter_for(representatives_source)

    try:
        canonical, ranks_found_in, _domains_found = _resolve_ranks(gtdb_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100)
        print("")
        sys.exit(0)

    for rank in ranks_found_in:
        total = _count_at_rank(gtdb_path, rank, canonical, rep_filter,
                               exclude_cores=exclude_cores)
        rows = rank_counts(gtdb_path, "gtdb", scope_rank=rank, scope_taxon=canonical,
                           rep_filter=rep_filter,
                           exclude_cores=exclude_cores)
        print("")
        wprint(f"  The rank '{rank}' has {total:,} {canonical} entries.")
        print("")
        print(render_rank_count_table(
            rows, count_header=f"Num. Unique Taxa under '{canonical}'"))

    print("")
    report_message("Each count above is also how many genomes `--derep-rank <rank>` "
                   "would return, since dereplication keeps one genome per unique "
                   "taxon at that rank.", "yellow",
                   ii="    ", si="    ", width=100, newline=False, trailing_newline=True)


def _report_taxon_counts_or_exit(gtdb_path, taxon, representatives_source, args=None, exclude_cores=None):
    """
    Report how many genomes match `taxon` at each rank it occurs at
    """
    rep_filter = _rep_filter_for(representatives_source)

    if str(taxon).lower() == "all":
        count = _count_total(gtdb_path,
                             exclude_cores=exclude_cores)
        print("")
        wprint(f"  There are {count:,} total genomes in the database.")
        if args is not None and _derep_is_on(args):
            derep_n = _all_derep_size(gtdb_path, "gtdb", args,
                                      exclude_cores=exclude_cores)
            wprint(f"    Dereplicated within each domain, that would be "
                   f"{derep_n:,} genome(s).")
        print("")
        if representatives_source:
            rep_type = _rep_type_label(representatives_source)
            wprint(color_text(f"  In considering only {rep_type} genomes:", "yellow"))
            print("")
            rep_n = _count_total(gtdb_path, rep_filter,
                                 exclude_cores=exclude_cores)
            wprint(f"  There are {rep_n:,} total "
                   f"{rep_type} genomes in the database.")
            print("")
        return

    # shared resolver: case-insensitive, returns (canonical, [ranks]); raises TaxonNotFound
    try:
        canonical, ranks_found_in, _domains_found = _resolve_ranks(gtdb_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(", "yellow",
                       ii="    ", si="    ", width=100)
        print("")
        sys.exit(0)

    taxon = canonical

    print("")
    for rank in ranks_found_in:
        count = _count_at_rank(gtdb_path, rank, taxon,
                               exclude_cores=exclude_cores)
        wprint(f"  The rank '{rank}' has {count:,} {taxon} entries.")
        if args is not None:
            _report_derep_note(gtdb_path, rank, taxon, args,
                               exclude_cores=exclude_cores)
    print("")

    if representatives_source:
        rep_type = _rep_type_label(representatives_source)
        wprint(color_text(f"  In considering only {rep_type} genomes:", "yellow"))
        print("")
        any_rep = False
        for rank in ranks_found_in:
            count = _count_at_rank(gtdb_path, rank, taxon, rep_filter,
                                   exclude_cores=exclude_cores)
            if count:
                any_rep = True
                wprint(f"  The rank '{rank}' has {count:,} {taxon} {rep_type} genome entries.")
                if args is not None:
                    _report_derep_note(gtdb_path, rank, taxon, args,
                                       rep_filter=rep_filter,
                           exclude_cores=exclude_cores)
                print("")
        if not any_rep:
            wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any "
                   "rank as a representative genome :(", "yellow"))
            print("")
            sys.exit(0)


def _all_derep_size(path, source, args, rep_filter=None, exclude_cores=None):
    """How many genomes `-w all --derep-rank X` would return."""
    pool = PoolSpec(path, source, rep_filter=rep_filter, label="GTDB",
                    taxon_flag="-w", exclude_cores=exclude_cores)
    return all_derep_size(pool, args.derep_rank)


def _rep_type_label(representatives_source):
    return "refseq reference" if representatives_source == "refseq" else "gtdb representative"


def report_unique_taxa_counts_of_all_ranks(gtdb_path, representatives_source=None,
                                           scoped_to_all=False, exclude_cores=None):
    """
    Print, for each of the 7 ranks, how many unique taxa exist in the GTDB table
    """
    domain_assigned = True if scoped_to_all else None

    print("")
    print(render_rank_count_table(rank_counts(gtdb_path, "gtdb",
                                              domain_assigned=domain_assigned,
                                              exclude_cores=exclude_cores)))
    print("")

    if scoped_to_all:
        report_message(scoped_counts_note("-w"), "yellow", ii="    ", si="    ",
                       width=100, trailing_newline=True)
        _report_unassigned_domains(unassigned_domain_summary(gtdb_path, "gtdb"))

    if representatives_source == "gtdb":
        report_message("(The `--gtdb-representatives-only` flag doesn't change these counts: "
                       "every GTDB taxon has a representative genome, so the number of unique "
                       "taxa per rank is the same with or without it.)",
                       "yellow", ii="    ", si="    ", width=100, newline=False, trailing_newline=True)
    elif representatives_source == "refseq":
        rep_rows = rank_counts(gtdb_path, "gtdb",
                               rep_filter=_rep_filter_for("refseq"),
                                                          exclude_cores=exclude_cores)
        wprint(color_text("  In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print(render_rank_count_table(rep_rows, count_header="Num. Unique Ref. Taxa"))
        print("")


def _report_gtdb_version(gtdb_path):
    version, release_date = _read_gtdb_version_info(os.path.dirname(gtdb_path))
    print("\n    Using GTDB " + version + ": " + release_date)


def copy_gtdb_table(gtdb_path):
    """
    Write the full GTDB metadata table (all columns in the asset) to the current
    directory as a tsv
    """
    _report_gtdb_version(gtdb_path)
    out_name = "gtdb-arc-and-bac-metadata.tsv"
    print("")
    with spinner("Writing GTDB table...", "", clear_on_done=True):
        write_table_tsv(pq.read_table(gtdb_path), out_name)
    wprint("  GTDB table written to:")
    print(color_text("      " + out_name + "\n"))


if __name__ == "__main__":
    main()
