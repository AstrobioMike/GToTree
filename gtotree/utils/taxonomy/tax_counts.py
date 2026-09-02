"""
COUNTING genomes and unique taxa in the GTDB / NCBI Parquet assets
"""

import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, REFERENCE_VALUE, rank_index
from gtotree.utils.taxonomy.tax_select import (SOURCES, assigned_mask,
                                               prefix_mask)
from gtotree.utils.taxonomy.exclusion_list import filter_table_by_exclusion


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

def read_pool(path, source, columns, rank=None, taxon=None, rep_filter=None,
          accession_prefixes=None, assembly_levels=None, domain_assigned=None,
          exclude_cores=None):
    """
    Read `columns` with the pool filters applied, as an Arrow Table

    domain_assigned:
        None  : no restriction (the whole table)
        True  : only rows with an assigned domain, i.e. the pool a `-w all` pull can
                actually reach
        False : only rows WITHOUT one (the viral / metagenome / unclassified stuff)
    exclude_cores:
        Optional accession cores from `--exclusion-list`, dropped from the pool the
        same way they're dropped from a real selection. Every count in this module
        funnels through here, so honoring it in this one place keeps the previewed
        numbers matching what a pull would actually return
    """
    spec = SOURCES[source]
    domain_col = RANKS[0]

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
    if domain_assigned is not None and domain_col not in cols:
        cols.append(domain_col)

    tab = pq.read_table(path, columns=cols, filters=filters or None)

    if accession_prefixes:
        tab = tab.filter(prefix_mask(tab.column(spec.acc_col), accession_prefixes))

    if domain_assigned is not None:
        tab = tab.filter(assigned_mask(tab.column(domain_col),
                                       assigned=domain_assigned))

    if exclude_cores:
        if spec.acc_col not in tab.column_names:
            tab = read_pool(path, source, cols + [spec.acc_col], rank=rank,
                            taxon=taxon, rep_filter=rep_filter,
                            accession_prefixes=accession_prefixes,
                            assembly_levels=assembly_levels,
                            domain_assigned=domain_assigned)
        tab, _n = filter_table_by_exclusion(tab, spec.acc_col, exclude_cores)

    return tab


################################################################################
# counts
################################################################################

def count_genomes(path, source, rank=None, taxon=None, rep_filter=None,
                  accession_prefixes=None, assembly_levels=None,
                  domain_assigned=None,
                  exclude_cores=None):
    """
    How many genome rows match. rank/taxon both None counts the whole table.
    """
    spec = SOURCES[source]
    tab = read_pool(path, source, [spec.acc_col], rank=rank, taxon=taxon,
                rep_filter=rep_filter, accession_prefixes=accession_prefixes,
                assembly_levels=assembly_levels, domain_assigned=domain_assigned,
                exclude_cores=exclude_cores)
    return tab.num_rows


def count_distinct_taxa(path, source, rank_col, scope_rank=None, scope_taxon=None,
                        rep_filter=None, accession_prefixes=None,
                        assembly_levels=None, domain_assigned=None,
                        exclude_cores=None):
    """
    How many DISTINCT taxa appear in `rank_col`, optionally scoped to a taxon
    """
    tab = read_pool(path, source, [rank_col], rank=scope_rank, taxon=scope_taxon,
                rep_filter=rep_filter, accession_prefixes=accession_prefixes,
                assembly_levels=assembly_levels, domain_assigned=domain_assigned,
                exclude_cores=exclude_cores)
    return _count_distinct_assigned(tab.column(rank_col))


def _count_distinct_assigned(col):
    """
    Distinct values of a rank column, ignoring nulls and the unassigned markers
    """
    if len(col) == 0:
        return 0
    return pc.count_distinct(col.filter(assigned_mask(col))).as_py()


def distinct_taxa(path, source, rank_col, scope_rank=None, scope_taxon=None,
                  rep_filter=None, accession_prefixes=None, assembly_levels=None,
                  domain_assigned=None,
                  exclude_cores=None):
    """
    The assigned taxa present in `rank_col`, sorted
    """
    tab = read_pool(path, source, [rank_col], rank=scope_rank, taxon=scope_taxon,
                rep_filter=rep_filter, accession_prefixes=accession_prefixes,
                assembly_levels=assembly_levels, domain_assigned=domain_assigned,
                exclude_cores=exclude_cores)
    col = tab.column(rank_col)
    if len(col) == 0:
        return []
    return sorted(pc.unique(col.filter(assigned_mask(col))).to_pylist())


def derep_size(path, source, scope_rank, scope_taxon, derep_rank, rep_filter=None,
               accession_prefixes=None, assembly_levels=None,
               exclude_cores=None):
    """
    How many genomes a --derep-rank pull of this taxon would return
    """
    return count_distinct_taxa(path, source, derep_rank,
                               scope_rank=scope_rank, scope_taxon=scope_taxon,
                               rep_filter=rep_filter,
                               accession_prefixes=accession_prefixes,
                               assembly_levels=assembly_levels,
                               exclude_cores=exclude_cores)


def rank_counts(path, source, scope_rank=None, scope_taxon=None, rep_filter=None,
                accession_prefixes=None, assembly_levels=None, domain_assigned=None,
                exclude_cores=None):
    """
    [(rank, num_unique_taxa), ...] for the ranks worth reporting

    Unscoped, that's all seven ranks. Scoped to a taxon it's the
    taxon's OWN rank and everything finer
    """
    ranks = RANKS if scope_rank is None else RANKS[rank_index(scope_rank):]
    return [(rank,
             count_distinct_taxa(path, source, rank,
                                 scope_rank=scope_rank, scope_taxon=scope_taxon,
                                 rep_filter=rep_filter,
                                 accession_prefixes=accession_prefixes,
                                 assembly_levels=assembly_levels,
                                 domain_assigned=domain_assigned,
                                 exclude_cores=exclude_cores))
            for rank in ranks]


################################################################################
# rendering
################################################################################

def render_rank_count_table(rows, count_header="Num. Unique Taxa", indent="    "):
    lines = ["{}{:<10} {:}".format(indent, "Rank", count_header)]
    lines += ["{}{:<10} {:,}".format(indent, rank, n) for rank, n in rows]
    return "\n".join(lines)
