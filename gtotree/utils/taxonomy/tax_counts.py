"""
COUNTING genomes and unique taxa in the GTDB / NCBI Parquet assets
"""

import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA, REFERENCE_VALUE, rank_index
from gtotree.utils.taxonomy.tax_select import SOURCES

# values in a rank column that mean "no taxon assigned here". derep() skips these
# rows entirely (`if not group or group == NA`), so they must not be counted as groups
UNASSIGNED = (NA, "")


################################################################################
# pool filters
################################################################################

def representatives_filter(source, kind):
    """
    The (column, value) predicate for a representatives-only pool, or None.

    `kind` is one of:
        None       -- no restriction
        "source"   -- the source's OWN representative flag (GTDB representatives;
                      for NCBI this is the RefSeq reference-genome flag)
        "refseq"   -- RefSeq 'reference genome', whatever the source

    Both CLIs route through here rather than writing the column names out, so a
    source's column naming lives in exactly one place (SourceSpec).
    """
    if kind in (None, False):
        return None

    spec = SOURCES[source]

    if kind == "source":
        return spec.rep_filter
    if kind == "refseq":
        return (spec.ref_col, REFERENCE_VALUE)

    raise ValueError(
        f"unknown representatives kind {kind!r}; expected None, 'source', or 'refseq'")


################################################################################
# reading
################################################################################

def _read(path, source, columns, rank=None, taxon=None, rep_filter=None,
          accession_prefixes=None, assembly_levels=None):
    """
    Read `columns` with the pool filters applied.

    rank/taxon scope the read to one taxon (both None -> the whole table). The rank
    predicate, the representatives predicate and the assembly-level predicate are all
    pushed down to Parquet; only the accession-prefix filter has to happen in Arrow,
    since prefix matching isn't expressible as a Parquet predicate.
    """
    spec = SOURCES[source]

    filters = []
    if rank and taxon:
        filters.append((rank, "=", taxon))
    if rep_filter:
        filters.append((rep_filter[0], "=", rep_filter[1]))
    if assembly_levels and spec.level_col:
        filters.append((spec.level_col, "in", set(assembly_levels)))

    cols = list(columns)
    if accession_prefixes and spec.acc_col not in cols:
        cols.append(spec.acc_col)

    tab = pq.read_table(path, columns=cols, filters=filters or None)

    if accession_prefixes:
        tab = tab.filter(_prefix_mask(tab.column(spec.acc_col), accession_prefixes))

    return tab


def _prefix_mask(acc_col, prefixes):
    """OR of starts_with over each prefix. Kept flat -- see the note in general.py
    about deeply chained Arrow predicates overflowing the stack."""
    mask = None
    for p in prefixes:
        m = pc.starts_with(acc_col, p)
        mask = m if mask is None else pc.or_(mask, m)
    return mask


################################################################################
# counts
################################################################################

def count_genomes(path, source, rank=None, taxon=None, rep_filter=None,
                  accession_prefixes=None, assembly_levels=None):
    """
    How many genome rows match. rank/taxon both None counts the whole table.
    """
    spec = SOURCES[source]
    tab = _read(path, source, [spec.acc_col], rank=rank, taxon=taxon,
                rep_filter=rep_filter, accession_prefixes=accession_prefixes,
                assembly_levels=assembly_levels)
    return tab.num_rows


def count_distinct_taxa(path, source, rank_col, scope_rank=None, scope_taxon=None,
                        rep_filter=None, accession_prefixes=None,
                        assembly_levels=None):
    """
    How many DISTINCT taxa appear in `rank_col`, optionally scoped to a taxon.

    Unassigned values (null, "NA", "") are excluded, matching what derep() does with
    them. The pre-existing global rank-counts tables didn't exclude them, so a rank
    with any unassigned lineage was reported one too high.
    """
    tab = _read(path, source, [rank_col], rank=scope_rank, taxon=scope_taxon,
                rep_filter=rep_filter, accession_prefixes=accession_prefixes,
                assembly_levels=assembly_levels)
    return _count_distinct_assigned(tab.column(rank_col))


def _count_distinct_assigned(col):
    """Distinct values of a rank column, ignoring nulls and the unassigned markers."""
    if len(col) == 0:
        return 0
    unassigned = pa.array(list(UNASSIGNED), type=col.type)
    keep = pc.and_(pc.is_valid(col), pc.invert(pc.is_in(col, value_set=unassigned)))
    return pc.count_distinct(col.filter(keep)).as_py()


def distinct_taxa(path, source, rank_col, scope_rank=None, scope_taxon=None,
                  rep_filter=None, accession_prefixes=None, assembly_levels=None):
    """
    The assigned taxa present in `rank_col`, sorted. Same exclusions as
    count_distinct_taxa (null / "NA" / ""), so the list and the count agree.

    Used to discover which domains an asset actually holds, rather than hardcoding
    Bacteria + Archaea: the NCBI table also carries eukaryotes (and anything else NCBI
    has assemblies for), and a hardcoded pair would silently drop them.
    """
    tab = _read(path, source, [rank_col], rank=scope_rank, taxon=scope_taxon,
                rep_filter=rep_filter, accession_prefixes=accession_prefixes,
                assembly_levels=assembly_levels)
    col = tab.column(rank_col)
    if len(col) == 0:
        return []
    unassigned = pa.array(list(UNASSIGNED), type=col.type)
    keep = pc.and_(pc.is_valid(col), pc.invert(pc.is_in(col, value_set=unassigned)))
    return sorted(pc.unique(col.filter(keep)).to_pylist())


def derep_size(path, source, scope_rank, scope_taxon, derep_rank, rep_filter=None,
               accession_prefixes=None, assembly_levels=None):
    """
    How many genomes a --derep-rank pull of this taxon would return.

    One genome per distinct value of `derep_rank`, so this is count_distinct_taxa()
    under a name that says what the number MEANS at the call site. The equivalence
    holds only because the same pool filters are applied here as in the pull; see the
    module docstring.
    """
    return count_distinct_taxa(path, source, derep_rank,
                               scope_rank=scope_rank, scope_taxon=scope_taxon,
                               rep_filter=rep_filter,
                               accession_prefixes=accession_prefixes,
                               assembly_levels=assembly_levels)


def rank_counts(path, source, scope_rank=None, scope_taxon=None, rep_filter=None,
                accession_prefixes=None, assembly_levels=None):
    """
    [(rank, num_unique_taxa), ...] for the ranks worth reporting.

    Unscoped, that's all seven ranks (the global summary). Scoped to a taxon it's the
    taxon's OWN rank and everything finer
    """
    ranks = RANKS if scope_rank is None else RANKS[rank_index(scope_rank):]
    return [(rank,
             count_distinct_taxa(path, source, rank,
                                 scope_rank=scope_rank, scope_taxon=scope_taxon,
                                 rep_filter=rep_filter,
                                 accession_prefixes=accession_prefixes,
                                 assembly_levels=assembly_levels))
            for rank in ranks]


################################################################################
# rendering
################################################################################

def render_rank_count_table(rows, count_header="Num. Unique Taxa", indent="    "):
    """
    The rank-counts table as a string, so the GTDB and NCBI helpers can't drift into
    printing it two different ways. Pure formatting -- the CLIs still decide what to
    say around it and when to say it.
    """
    lines = ["{}{:<10} {:}".format(indent, "Rank", count_header)]
    lines += ["{}{:<10} {:,}".format(indent, rank, n) for rank, n in rows]
    return "\n".join(lines)
