"""
The `-w` TARGET vocabulary shared by every surface that takes a taxon

`-w` accepts a taxon name (see tax_ranks.TAXON_ALIASES for the spellings that get
redirected), an NCBI taxid (the NCBI helper only), or the keyword 'all'. This module
owns 'all': what it expands to, and what it can't reach.

'all' expands to each domain present in the current asset (gtdb/ncbi)

it does not take things beyond the primary domains, though those can be grabbed
explicitly by themselves (e.g., `-w Uroviricota` works fine)
"""

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_counts import (count_genomes, count_distinct_taxa,
                                               distinct_taxa)

ALL_KEYWORD = "all"

DOMAIN_RANK = RANKS[0]


def is_all_target(value):
    """True for the 'all' keyword, in any casing/spacing."""
    return str(value or "").strip().lower() == ALL_KEYWORD


def domains_in_asset(path, source, rep_filter=None, accession_prefixes=None,
                     assembly_levels=None):
    """
    The domains present in the asset, under the same pool filters the pull will use.

    Sorted, and with unassigned ('NA' / '' / null) excluded
    """
    return distinct_taxa(path, source, DOMAIN_RANK, rep_filter=rep_filter,
                         accession_prefixes=accession_prefixes,
                         assembly_levels=assembly_levels)


def expand_all_targets(path, source, taxa, rep_filter=None, accession_prefixes=None,
                       assembly_levels=None):
    """
    Replace any 'all' in `taxa` with the asset's domains, preserving order and
    dropping duplicates (so `-w all -w bacteria` doesn't resolve Bacteria twice).

    Returns (expanded_taxa, domains_used). `domains_used` is empty when 'all' wasn't
    among the targets, which is also how a caller knows whether to say anything.
    """
    taxa = [t for t in (taxa or []) if str(t).strip()]

    if not any(is_all_target(t) for t in taxa):
        return list(taxa), []

    domains = domains_in_asset(path, source, rep_filter=rep_filter,
                               accession_prefixes=accession_prefixes,
                               assembly_levels=assembly_levels)

    expanded = []
    seen = set()
    for taxon in taxa:
        for name in (domains if is_all_target(taxon) else [taxon]):
            key = str(name).strip().lower()
            if key not in seen:
                seen.add(key)
                expanded.append(name)

    return expanded, domains


class UnassignedDomainSummary:
    """
    What an 'all' pull leaves behind, if anything
    """

    def __init__(self, n_genomes, rank=None, n_taxa_at_rank=0):
        self.n_genomes = n_genomes
        self.rank = rank
        self.n_taxa_at_rank = n_taxa_at_rank

    def __bool__(self):
        return bool(self.n_genomes)

    def message(self, source_label="NCBI"):
        """One human-facing line, or None when there's nothing to report."""
        if not self.n_genomes:
            return None

        line = (f"{self.n_genomes:,} genome(s) in the {source_label} table have no "
                f"assigned domain (likely viral and/or uncultured entries) and so "
                f"aren't included in 'all'.")

        if self.rank and self.n_taxa_at_rank:
            line += (f" They hold {self.n_taxa_at_rank:,} '{self.rank}' value(s) found "
                     f"nowhere else, which is why an unscoped `--get-rank-counts` "
                     f"reports more than 'all' can pull.")

        line += (" They can still be pulled by name explicitly if you're looking "
                 "for something like that (e.g., `-w Uroviricota`).")

        return line


def unassigned_domain_summary(path, source, rank=None, rep_filter=None,
                              accession_prefixes=None, assembly_levels=None):
    """
    Describe the domain-less slice of the asset under the given pool filters.

    `rank` is the rank the caller is about to dereplicate at (or report on); when
    given, the summary also carries how many distinct values of that rank live only
    in the domain-less rows.
    """
    n_genomes = count_genomes(path, source, rep_filter=rep_filter,
                              accession_prefixes=accession_prefixes,
                              assembly_levels=assembly_levels,
                              domain_assigned=False)

    if not n_genomes or not rank:
        return UnassignedDomainSummary(n_genomes, rank=rank)

    n_taxa = count_distinct_taxa(path, source, rank, rep_filter=rep_filter,
                                 accession_prefixes=accession_prefixes,
                                 assembly_levels=assembly_levels,
                                 domain_assigned=False)

    return UnassignedDomainSummary(n_genomes, rank=rank, n_taxa_at_rank=n_taxa)
