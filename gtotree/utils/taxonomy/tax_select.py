"""
Taxon-based genome/accession selection from the GTDB / NCBI Parquet assets.

The filter/select piece of the taxonomy toolkit: resolve a user-supplied taxon to
a (canonical, rank) pair and pull the rows under it from the hosted Parquet assets
(bit's ncbi-data.parquet / gtdb-data.parquet). Column names are the only thing that
differs between sources; that difference is isolated in SourceSpec / SOURCES.
"""

import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA, REFERENCE_VALUE, accession_core, rank_index


class SourceSpec:
    """Per-source column names. The only thing that differs between GTDB and NCBI."""

    def __init__(self, name, acc_col, rep_filter, quality_cols,
                 ref_col, default_reps_only, level_col=None, contig_col=None,
                 size_col=None, size_fallback_col=None):
        self.name = name
        self.acc_col = acc_col
        self.rep_filter = rep_filter
        self.quality_cols = quality_cols
        self.ref_col = ref_col
        self.default_reps_only = default_reps_only
        self.level_col = level_col
        self.contig_col = contig_col
        self.size_col = size_col
        self.size_fallback_col = size_fallback_col


SOURCES = {
    "gtdb": SourceSpec(
        name="gtdb",
        acc_col="ncbi_genbank_assembly_accession",
        rep_filter=("gtdb_representative", "t"),
        quality_cols=("checkm2_completeness", "checkm2_contamination"),
        ref_col="ncbi_refseq_category",
        default_reps_only=True,
        size_col="genome_size",
        contig_col="contig_count",
    ),
    "ncbi": SourceSpec(
        name="ncbi",
        acc_col="assembly_accession",
        rep_filter=("refseq_category", REFERENCE_VALUE),
        quality_cols=("checkm_completeness", "checkm_contamination"),
        ref_col="refseq_category",
        default_reps_only=False,
        level_col="assembly_level",
        contig_col="contig_count",
        size_col="genome_size_ungapped",
        size_fallback_col="genome_size",
    ),
}


class TaxonNotFound(Exception):
    pass


class AmbiguousTaxon(Exception):
    """The taxon name exists at more than one rank; caller must disambiguate."""

    def __init__(self, taxon, ranks_found):
        self.taxon = taxon
        self.ranks_found = ranks_found
        super().__init__(
            f"'{taxon}' occurs at more than one rank ({', '.join(ranks_found)}); "
            f"specify which rank is wanted")


# ---------------------------------------------------------------------------
# resolving a user-supplied taxon
# ---------------------------------------------------------------------------

def find_ranks_for_taxon(path, taxon):
    """
    Which of the 7 ranks contain `taxon`? Case-insensitive.
    Returns (canonical_name, [ranks]). Reads one column at a time, so this stays
    cheap even on the 4M-row NCBI table.
    """
    target = str(taxon).strip().lower()
    canonical = None
    found = []
    for rank in RANKS:
        col = pq.read_table(path, columns=[rank]).column(rank)
        uniq = pc.unique(col).to_pylist()
        for name in uniq:
            if name is not None and name != NA and name.lower() == target:
                canonical = name
                found.append(rank)
                break
    if canonical is None:
        raise TaxonNotFound(f"'{taxon}' doesn't exist at any rank in this source")
    return canonical, found


def resolve_taxon(path, taxon, rank=None):
    """
    Resolve a user-supplied taxon (+ optional explicit rank) to (canonical, rank).

    Raises AmbiguousTaxon if the name lives at multiple ranks and no rank was
    given. E.g., a name used as both an order and a family. On the NCBI side a user can
    sidestep this entirely by passing a taxid for now (though i might remove this; revisit Mike)
    """
    canonical, found = find_ranks_for_taxon(path, taxon)

    if rank:
        r = str(rank).strip().lower()
        rank_index(r)
        if r not in found:
            raise TaxonNotFound(
                f"'{canonical}' exists at rank(s) {', '.join(found)}, not '{r}'")
        return canonical, r

    if len(found) > 1:
        raise AmbiguousTaxon(canonical, found)

    return canonical, found[0]


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------

def select(path, source, rank, taxon, reps_only=False, columns=None):
    """
    All rows under `taxon` at `rank`. Returns a pyarrow Table.

    The rank predicate is pushed down to Parquet, so on the lineage-sorted assets
    this skips whole row groups rather than scanning.
    """
    spec = SOURCES[source]
    filters = [(rank, "=", taxon)]
    if reps_only:
        filters.append((spec.rep_filter[0], "=", spec.rep_filter[1]))

    cols = columns or [spec.acc_col]
    # always need the accession back
    if spec.acc_col not in cols:
        cols = [spec.acc_col] + list(cols)

    return pq.read_table(path, columns=cols, filters=filters)


def select_accessions(path, source, rank, taxon, reps_only=False):
    """The common case: just the accession list."""
    spec = SOURCES[source]
    tab = select(path, source, rank, taxon, reps_only=reps_only)
    return [a for a in tab.column(spec.acc_col).to_pylist() if a and a != NA]


def select_by_taxid(path, rank, taxid):
    """
    NCBI only: select by a lineage TAXID rather than a name
    """
    col = f"{rank}_taxid"
    tab = pq.read_table(path, columns=["assembly_accession"],
                        filters=[(col, "=", str(taxid))])
    return tab.column("assembly_accession").to_pylist()


# ---------------------------------------------------------------------------
# liveness screening (suppressed / removed / version-drifted assemblies)
# ---------------------------------------------------------------------------

def live_accession_cores(ncbi_table_path):
    """
    The set of assembly core accs currently present in the NCBI assembly summary

    NCBI drops suppressed/removed assemblies from the assembly summary files, so
    ABSENCE == suppressed / removed. GTDB is a snapshot and can still have them
    """
    tab = pq.read_table(ncbi_table_path, columns=["assembly_accession"])
    return {accession_core(a) for a in tab.column("assembly_accession").to_pylist() if a}
