"""
The surface `gtt get-accs-from-gtdb` and `gtt get-accs-from-ncbi` share
"""

from gtotree.utils.taxonomy.tax_counts import (count_distinct_taxa, count_genomes,
                                               derep_size, rank_counts, read_pool,
                                               representatives_filter)
from gtotree.utils.taxonomy.tax_derep import (resolve_derep_rank, select_all_domains,
                                              select_ref_genomes)
from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_targets import (domains_in_asset,
                                                unassigned_domain_summary)
from gtotree.utils.taxonomy.exclusion_list import exclusion_list_help

ASSEMBLY_LEVELS = {
    "complete": "Complete Genome",
    "chromosome": "Chromosome",
    "scaffold": "Scaffold",
    "contig": "Contig",
}

DEREP_OFF = (None, "off", "none", "None")

# --derep-rank's real default, and the placeholder argparse actually stores.
#
# The sentinel exists so a caller can tell "the user typed --derep-rank auto" from
# "the user typed nothing". The taxid path needs that distinction: dereplication
# can't be applied to a taxid, so an inherited default should quietly step aside
# while an explicit request should still be an error.
#
# It is deliberately NOT None: None is in DEREP_OFF, so an unresolved None would
# silently mean "off". An unresolved sentinel instead reaches resolve_derep_rank()
# and raises, which is the failure mode we want.
DEREP_DEFAULT = "auto"
DEREP_UNSET = "__derep_unset__"


def representatives_suffix(representatives_source):
    """
    The filename suffix for a pull restricted to a source's curated genomes.

    Empty string when the pull wasn't restricted.
    """
    if not representatives_source:
        return ""
    word = "rep" if str(representatives_source).strip().lower() == "gtdb" else "ref"
    return f"-{representatives_source}-{word}"


def resolved_derep_rank(args):
    """
    The --derep-rank sitting on `args`, already resolved by apply_derep_default().
    """
    derep_rank = getattr(args, "derep_rank", DEREP_DEFAULT)
    return DEREP_DEFAULT if derep_rank == DEREP_UNSET else derep_rank


def is_derep_on(derep_rank):
    """True when --derep-rank asks for actual dereplication."""
    return derep_rank not in DEREP_OFF


def apply_derep_default(args):
    """
    Swap the --derep-rank sentinel on `args` for the real default, in place.

    Returns True when the user set --derep-rank themselves. Safe to call more than
    once, but only the first call can report that faithfully, so the caller that
    needs the answer (the taxid guard) has to be the first one in.
    """
    if getattr(args, "derep_rank", DEREP_UNSET) == DEREP_UNSET:
        args.derep_rank = DEREP_DEFAULT
        return False
    return True


def parse_assembly_levels(value):
    """--assembly-level values -> the NCBI table's strings. Raises ValueError."""
    if not value:
        return []
    if isinstance(value, str):
        value = value.split(",")
    parts = [str(v).strip().lower() for v in value if str(v).strip()]
    unknown = [p for p in parts if p not in ASSEMBLY_LEVELS]
    if unknown:
        raise ValueError(
            f"unrecognised --assembly-level value(s): {', '.join(unknown)}. "
            f"Choose from: {', '.join(ASSEMBLY_LEVELS)}")
    return [ASSEMBLY_LEVELS[p] for p in parts]


def source_prefixes(source):
    """Accession prefixes for an NCBI --ncbi-section value (None means no restriction)."""
    if source == "refseq":
        return ("GCF_",)
    if source == "genbank":
        return ("GCA_",)
    return None


SOURCE_OVERLAP_NOTE = (
    "There is a lot of overlap between RefSeq and GenBank (most RefSeq entries are "
    "derived from a GenBank one), so `--ncbi-section both` will return many duplicate "
    "assemblies. You most likely only want one of them, but you're the boss!")


def source_overlap_note(source):
    """The `--ncbi-section both` warning, or None for a single-section source."""
    if str(source or "").strip().lower() != "both":
        return None
    return SOURCE_OVERLAP_NOTE


class PoolSpec:
    """
    The set of genomes a helper is working within, and every question asked of it

    Attributes:
        table_path         : the Parquet asset
        source             : "gtdb" | "ncbi" (the tax_select SOURCES key)
        rep_filter         : representatives predicate, or None
        accession_prefixes : e.g. ("GCF_",), or None
        assembly_levels    : NCBI level strings, or None
        label              : "GTDB" / "NCBI", for user-facing messages
        taxon_flag         : the flag this surface calls the target, for messages
        exclude_cores      : accession cores from `--exclusion-list`, or None
    """

    def __init__(self, table_path, source, rep_filter=None, accession_prefixes=None,
                 assembly_levels=None, label=None, taxon_flag="-w",
                 exclude_cores=None):
        self.table_path = table_path
        self.source = source
        self.rep_filter = rep_filter
        self.accession_prefixes = accession_prefixes
        self.assembly_levels = assembly_levels or None
        self.label = label or source.upper()
        self.taxon_flag = taxon_flag
        self.exclude_cores = exclude_cores or None

    # -- the pool filters as kwargs, so every call site passes the same set ---------

    def _filters(self, with_reps=True):
        return {
            "rep_filter": self.rep_filter if with_reps else None,
            "accession_prefixes": self.accession_prefixes,
            "assembly_levels": self.assembly_levels,
            "exclude_cores": self.exclude_cores,
        }

    def without_reps(self):
        """The same pool minus the representatives predicate."""
        return PoolSpec(self.table_path, self.source,
                        accession_prefixes=self.accession_prefixes,
                        assembly_levels=self.assembly_levels, label=self.label,
                        taxon_flag=self.taxon_flag,
                        exclude_cores=self.exclude_cores)

    def with_reps(self, kind):
        """The same pool plus a representatives predicate ('source' | 'refseq')."""
        return PoolSpec(self.table_path, self.source,
                        rep_filter=representatives_filter(self.source, kind),
                        accession_prefixes=self.accession_prefixes,
                        assembly_levels=self.assembly_levels, label=self.label,
                        taxon_flag=self.taxon_flag,
                        exclude_cores=self.exclude_cores)

    # -- counting ------------------------------------------------------------------

    def count_genomes(self, rank=None, taxon=None, domain_assigned=None):
        return count_genomes(self.table_path, self.source, rank=rank, taxon=taxon,
                             domain_assigned=domain_assigned, **self._filters())

    def count_distinct(self, rank_col, scope_rank=None, scope_taxon=None,
                       domain_assigned=None):
        return count_distinct_taxa(self.table_path, self.source, rank_col,
                                   scope_rank=scope_rank, scope_taxon=scope_taxon,
                                   domain_assigned=domain_assigned, **self._filters())

    def rank_counts(self, scope_rank=None, scope_taxon=None, domain_assigned=None):
        return rank_counts(self.table_path, self.source, scope_rank=scope_rank,
                           scope_taxon=scope_taxon, domain_assigned=domain_assigned,
                           **self._filters())

    def derep_size(self, scope_rank, scope_taxon, derep_rank):
        """How many genomes a --derep-rank pull of this taxon would return."""
        return derep_size(self.table_path, self.source, scope_rank, scope_taxon,
                          derep_rank, **self._filters())

    # -- the pool itself -----------------------------------------------------------

    def read(self, columns, rank=None, taxon=None, domain_assigned=None):
        """The pool as an Arrow Table, filters applied before materializing."""
        return read_pool(self.table_path, self.source, columns, rank=rank, taxon=taxon,
                         domain_assigned=domain_assigned, **self._filters())

    def domains(self):
        return domains_in_asset(self.table_path, self.source, **self._filters())

    def unassigned(self, rank=None):
        """What a `<taxon_flag> all` pull leaves behind. See tax_targets."""
        return unassigned_domain_summary(self.table_path, self.source, rank=rank,
                                         **self._filters())

    # -- selection -----------------------------------------------------------------

    def select(self, taxon, target_rank=None, derep_rank="off", reps_only=None,
               **kwargs):
        return select_ref_genomes(self.table_path, self.source, taxon,
                                  target_rank=target_rank, derep_rank=derep_rank,
                                  reps_only=reps_only,
                                  accession_prefixes=self.accession_prefixes,
                                  assembly_levels=self.assembly_levels, **kwargs)

    def select_all_domains(self, derep_rank="auto", reps_only=None, **kwargs):
        return select_all_domains(self.table_path, self.source, derep_rank=derep_rank,
                                  reps_only=reps_only,
                                  accession_prefixes=self.accession_prefixes,
                                  assembly_levels=self.assembly_levels, **kwargs)


def derep_note(pool, rank, taxon, derep_rank):
    """
    The "...dereplicated at X, that would be N" line for one rank, or None when
    dereplication is off

    Returns (line_or_None, warnings).
    """
    if not is_derep_on(derep_rank):
        return None, []

    try:
        effective, warnings = resolve_derep_rank(rank, derep_rank)
    except ValueError as err:
        # an explicit rank coarser than this one -- can't dereplicate to it
        return str(err), []

    if effective is None:
        return None, warnings

    n = pool.derep_size(rank, taxon, effective)
    return f"Dereplicated at '{effective}', that would be {n:,} genome(s).", warnings


FILTERS_APPLIED_NOTE = " (after any specified filters)"
FILTERS_IMPLICATED_NOTE = " with the specified filters"


def _with_note(text, note):
    """Splice `note` in just inside a line's trailing period, if it has one."""
    if text.endswith("."):
        return text[:-1] + note + "."
    return text + note


def with_filters_note(text):
    """
    Append FILTERS_APPLIED_NOTE just inside a line's trailing period, if it has one.

    "The rank 'genus' has 2 X entries." ->
    "The rank 'genus' has 2 X entries (after any specified filters)."
    """
    return _with_note(text, FILTERS_APPLIED_NOTE)


def with_filters_implicated_note(text):
    """
    The stronger sibling of with_filters_note(), for the case where we KNOW a filter
    is what emptied the result rather than only that filters may have applied.

    "No accessions were found for the --wanted-ref-tax target 'Trichodesmium'." ->
    "No accessions were found for the --wanted-ref-tax target 'Trichodesmium' with
     the specified filters."
    """
    return _with_note(text, FILTERS_IMPLICATED_NOTE)


def pull_count_lines(pool, resolved_rank, taxon, effective_derep_rank, kept_n):
    """
    The count block a completed pull prints above its "Wrote N accession(s)" lines
    """
    total = pool.count_genomes(rank=resolved_rank, taxon=taxon)
    header = with_filters_note(
        f"The rank '{resolved_rank}' has {total:,} {taxon} entries.")
    if effective_derep_rank is None:
        return header, None
    derep = (f"Dereplicated at '{effective_derep_rank}', that is "
             f"{kept_n:,} genome(s).")
    return header, derep


def all_derep_size(pool, derep_rank):
    """
    How many genomes a `<taxon_flag> all --derep-rank X` pull would return.

    The effective rank is resolved PER DOMAIN, because `auto` is domain-aware:
    Eukaryota dereplicates one rank coarser than the prokaryotic domains, so a single
    asset-wide rank would over- or under-count it. An explicit --derep-rank resolves
    to the same rank for every domain, so this still matches that case exactly.
    """
    total = 0
    for domain in pool.domains():
        effective, _warnings = resolve_derep_rank(RANKS[0], derep_rank, domain=domain)
        if effective is None:
            total += pool.count_genomes(rank=RANKS[0], taxon=domain,
                                        domain_assigned=True)
        else:
            total += pool.derep_size(RANKS[0], domain, effective)
    return total


def scoped_counts_note(taxon_flag="-w"):
    """Why a rank-counts table scoped to 'all' is the number that matters."""
    return (f"Counts are scoped to what `{taxon_flag} all` pulls, so each is also how "
            f"many accessions `{taxon_flag} all --derep-rank <rank>` returns.")


def add_common_get_accs_args(required, optional, source_label,
                             taxon_flags=("-w", "--wanted-tax"),
                             taxon_help=None):
    """
    Declare the flags both helpers take
    """
    required.add_argument(
        *taxon_flags,
        metavar="<STR>",
        dest="wanted_ref_tax",
        help=taxon_help or ("Wanted tax to get accessions for (a name, or 'all'). Not needed with `--get-rank-counts`."),
        action="store",
    )

    optional.add_argument(
        "-r", "--target-rank",
        choices=list(RANKS),
        help=("Target rank (if needed to disambiguate a taxon name that exists at "
              "multiple ranks)"),
        action="store",
    )

    optional.add_argument(
        "--target-domain",
        dest="target_domain",
        metavar="<STR>",
        default=None,
        help=("Target domain (if needed to disambiguate a taxon name that exists in "
              "multiple domains)"),
        action="store",
    )

    optional.add_argument(
        "--derep-rank",
        choices=["auto", "off"] + list(RANKS),
        default=DEREP_UNSET,
        help=("Dereplicate the pulled genomes down to a single best genome per unique "
              "value of this rank (default: auto). 'auto' is two ranks finer than the "
              "target (one rank finer for eukaryotes). E.g., '--derep-rank family' "
              "keeps one genome per family within the target taxon. 'off' returns "
              "every genome under the target taxon, so use with care :)"),
        action="store",
    )

    optional.add_argument(
        "--exclusion-list",
        metavar="<FILE>",
        dest="exclusion_list",
        default=None,
        help=exclusion_list_help(taxon_flags[0]),
        action="store",
    )

    optional.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help=(f"Provide this flag along with a specified taxon to `{taxon_flags[0]}` "
              f"to see how many genomes match the set parameters."),
    )

    optional.add_argument(
        "--get-rank-counts",
        action="store_true",
        help=(f"Provide this flag to see counts of how many unique taxa there are for "
              f"each rank. By itself, that'd be the whole database, but it can also "
              f"be combined with `{taxon_flags[0]}` and `--derep-rank`."),
    )

    optional.add_argument(
        "--get-table",
        action="store_true",
        help=(f"Provide just this flag alone to write out a tsv of GToTree's {source_label} "
              f"metadata table."),
    )

    optional.add_argument(
        "-R", "--refseq-ref-genomes-only",
        dest="refseq_reference_genomes_only",
        action="store_true",
        help="Pull only genomes designated as RefSeq reference genomes.",
    )

    return optional
