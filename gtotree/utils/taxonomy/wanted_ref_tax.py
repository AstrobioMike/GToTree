"""
Driver-side --wanted-ref-tax (-w) resolution
"""

from gtotree.utils.taxonomy.tax_derep import select_ref_genomes, size_advice
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
                                      max_contamination=None, building_tree=False):
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
    TaxonNotFound, AmbiguousTaxon, ValueError
        Propagated from the taxonomy core for the CLI layer to translate.
    """

    core_source, table_path = _table_path_for_source(source)

    # liveness screening: ensuring what we pick from gtdb is still present in ncbi (e.g., not suppressed)
    screen_against = ncbi_data_table_path() if core_source == "gtdb" else None

    selection = select_ref_genomes(
        table_path, core_source, taxon,
        target_rank=target_rank, derep_rank=derep_rank,
        screen_against=screen_against,
        min_completeness=min_completeness,
        max_contamination=max_contamination)

    if not selection.accessions:
        detail = ""
        if min_completeness is not None or max_contamination is not None:
            detail = (" No genomes cleared the requested quality floor, so you may "
                      "want to relax `--min-completeness` / `--max-contamination`.")
        raise WantedRefTaxError(
            f"No accessions were found for the --wanted-ref-tax target "
            f"'{selection.canonical}'.{detail}")

    if building_tree:
        selection.warnings.extend(
            size_advice(len(selection.accessions),
                        selection.resolved_rank,
                        selection.effective_derep_rank))

    return selection.accessions, selection
