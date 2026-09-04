"""
Driver-side --wanted-ref-tax (-w) resolution
"""

from gtotree.utils.taxonomy.tax_derep import select_ref_genomes, size_advice
from gtotree.utils.taxonomy.empty_selection import empty_selection_message
from gtotree.utils.taxonomy.tax_targets import expand_all_targets, is_all_target
from gtotree.utils.gtdb.get_gtdb_data import (gtdb_data_table_path,
                                              check_gtdb_location_var_is_set,
                                              report_gtdb_version_info)
from gtotree.utils.ncbi.get_ncbi_assembly_data import (ncbi_data_table_path,
                                                       check_ncbi_assembly_info_location_var_is_set,
                                                       read_date_retrieved)


_SOURCE_ASSETS = {
    "gtdb": ("gtdb", gtdb_data_table_path),
    "ncbi": ("ncbi", ncbi_data_table_path),
}


NCBI_SECTION_PREFIXES = {
    "refseq": ("GCF_",),
    "genbank": ("GCA_",),
    "both": None,
}


def section_prefixes(source, ncbi_section):
    """
    Accession prefixes for an --ncbi-section value, or None for no restriction.
    """
    if str(source).strip().lower() == "gtdb":
        return None
    return NCBI_SECTION_PREFIXES.get(str(ncbi_section or "refseq").strip().lower())


class WantedRefTaxError(Exception):
    """A --wanted-ref-tax request that resolved to nothing usable."""


def _table_path_for_source(source):
    """The asset path for a driver `--source` value. Raises WantedRefTaxError."""
    key = str(source).strip().lower()
    if key not in _SOURCE_ASSETS:
        raise WantedRefTaxError(
            f"'{source}' is not a recognized --source for --wanted-ref-tax "
            f"(expected one of: {', '.join(s.upper() for s in _SOURCE_ASSETS)}).")
    core_source, table_path_fn = _SOURCE_ASSETS[key]
    return core_source, table_path_fn()


def expand_wanted_ref_tax(source, taxa):
    taxa = list(taxa or [])
    if not any(is_all_target(t) for t in taxa):
        return taxa, []

    core_source, table_path = _table_path_for_source(source)

    return expand_all_targets(table_path, core_source, taxa)


def describe_all_expansion(source, domains):
    """The one-line 'this is what all meant' note, or None."""
    if not domains:
        return None
    return (f"`-w all` was expanded to the domains in the {str(source).upper()} "
            f"table: {', '.join(domains)}.")


def describe_source_version(source):
    try:
        if source == "gtdb":
            location = check_gtdb_location_var_is_set()
            version, release_date = report_gtdb_version_info(location)
            if release_date:
                return f"GTDB {version} (released {release_date})"
            return f"GTDB {version}"
        if source == "ncbi":
            location = check_ncbi_assembly_info_location_var_is_set()
            accessed = read_date_retrieved(location)
            return f"NCBI (accessed {accessed})"
    except (OSError, ValueError, IndexError, SystemExit):
        return None
    return None


def resolve_wanted_ref_tax_accessions(source, taxon, target_rank=None,
                                      derep_rank="auto", min_completeness=None,
                                      max_contamination=None, target_domain=None,
                                      ncbi_section="refseq", include_rows=True,
                                      building_tree=False, exclude_cores=None,
                                      reps_only=None, assembly_levels=None):
    """
    Resolve `-w <taxon>` to a list of assembly accessions plus the RefGenomeSelection
    it came from (for warnings / provenance the caller may want to surface).

    Parameters
    ----------
    source : str
        The driver's --source value ('gtdb' or 'ncbi'; case-insensitive here).
    taxon : str
        The -w taxon name.
    target_rank : str or None
        --target-rank; disambiguates a name that lives at multiple ranks.
    derep_rank : str
        --derep-rank ('auto' | 'off' | a rank name).
    min_completeness : float or None
        --min-completeness; drops candidates below this checkm completeness before
        selection. None (the default) means no floor.
    max_contamination : float or None
        --max-contamination; drops candidates above this checkm contamination before
        selection. None (the default) means no ceiling.
    ncbi_section : str
        --ncbi-section ('refseq' | 'genbank' | 'both')
    target_domain : str or None
        --target-domain; disambiguates a name shared across domains (e.g. `Bacillus`).
        Without it such a name raises CrossDomainTaxon. Also scopes selection to that
        domain and drives the domain-aware `auto` derep step.
    include_rows : bool
        Whether to carry the selected genomes' metadata rows back on the returned
        selection. Default True, and the main driver needs them -- scg_hmms_setup
        reads the pulled rows' lineage to auto-select an HMM set. Pass False only from
        a caller that wants accessions and counts alone (`gtt dl-ncbi-assemblies`); it
        skips a second filtered read of the asset per taxon without changing which
        genomes are selected.
    reps_only : bool or None
        Restrict the candidate pool to the source's representative/reference
        genomes. None (the default) defers to the source's own default, which is
        representatives-only for GTDB and all genomes for NCBI. Callers that want
        the same meaning for both sources pass an explicit bool.
    assembly_levels : list or None
        Optional assembly_level restriction, applied to the candidate pool before
        selection.
    exclude_cores : set or None
        Accession cores from `--exclusion-list`. Handed to the selection core rather
        than applied to its result, so the listed genomes leave the candidate pool
        before dereplication: excluding a group's chosen representative promotes the
        next-best genome in that group instead of dropping the group entirely.
    building_tree : bool
        Whether the selected genomes become a tree. Only the main GToTree driver sets
        this; `gtt gen-scg-hmms` and `gtt search-annotations` share this
        entry point but use the genomes as a search set, so tree-framed advice ("that
        is a very large tree", "consider a coarser --derep-rank") is wrong for them.

    Returns
    -------
    (accessions, selection) : (list[str], RefGenomeSelection)

    Raises
    ------
    WantedRefTaxError
        Unknown source, or the selection produced no accessions.
    TaxonNotFound, AmbiguousTaxon, CrossDomainTaxon, ValueError
        Propagated from the taxonomy core for the CLI layer to translate.
    """

    core_source, table_path = _table_path_for_source(source)

    # liveness screening: ensuring what we pick from gtdb is still present in ncbi (e.g., not suppressed)
    screen_against = ncbi_data_table_path() if core_source == "gtdb" else None

    selection = select_ref_genomes(
        table_path, core_source, taxon,
        target_rank=target_rank, derep_rank=derep_rank,
        screen_against=screen_against,
        accession_prefixes=section_prefixes(source, ncbi_section),
        reps_only=reps_only,
        assembly_levels=assembly_levels,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        target_domain=target_domain,
        exclude_cores=exclude_cores,
        include_rows=include_rows,
        diagnose_empty=True)

    if not selection.accessions:
        raise WantedRefTaxError(empty_selection_message(
            selection, taxon_flag="--wanted-ref-tax",
            assembly_levels=assembly_levels,
            ncbi_section=(None if str(source).strip().lower() == "gtdb"
                          else ncbi_section),
            reps_only_requested=bool(reps_only),
            min_completeness=min_completeness,
            max_contamination=max_contamination))

    if building_tree:
        selection.warnings.extend(
            size_advice(len(selection.accessions),
                        selection.resolved_rank,
                        selection.effective_derep_rank))

    return selection.accessions, selection
