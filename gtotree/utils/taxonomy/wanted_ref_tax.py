"""
Driver-side --wanted-ref-tax (-w) resolution.

This is the GToTree driver's counterpart to the standalone get-accs-from-*
helpers: it turns a `-w <taxon>` request (scoped by `-S {gtdb,ncbi}`, and honouring
`--target-rank` / `--derep-rank`) into a list of assembly accessions to fold into the
user's other input genomes

Library/CLI seam: everything here RAISES (WantedRefTaxError, or the taxonomy core's
TaxonNotFound / AmbiguousTaxon / ValueError). The CLI translation -- friendly message
+ early exit -- lives in preflight_checks.resolve_wanted_ref_tax(), so this module can
stay thin and be exercised without a process exit.
"""

from gtotree.utils.taxonomy.tax_derep import select_ref_genomes, size_advice
from gtotree.utils.gtdb.get_gtdb_data import (gtdb_data_table_path,
                                              check_gtdb_location_var_is_set,
                                              report_gtdb_version_info)
from gtotree.utils.ncbi.get_ncbi_assembly_data import (ncbi_data_table_path,
                                                       check_ncbi_assembly_info_location_var_is_set,
                                                       read_date_retrieved)


# source-of-taxonomy (-S) -> (core source name, table-path resolver). The driver's
# -S is GTDB/NCBI (which asset supplies the reference genomes), distinct from the NCBI
# helper's --source (refseq/genbank accession-prefix scoping).
_SOURCE_ASSETS = {
    "gtdb": ("gtdb", gtdb_data_table_path),
    "ncbi": ("ncbi", ncbi_data_table_path),
}


class WantedRefTaxError(Exception):
    """A --wanted-ref-tax request that resolved to nothing usable."""


def describe_source_version(source):
    """
    Best-effort one-line description of the taxonomy source actually in use, for the
    Phase 1 banner (the counterpart to Phase 3's "Pfam version being used: ...").

    GTDB -> "GTDB <release> (released <date>)"; NCBI -> "NCBI (accessed <date>)".

    Returns None if the source is unknown or its version/date sidecar can't be read.
    This is a display label, so a missing sidecar shouldn't abort the run -- the caller
    just skips the line.
    """
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
                                      max_contamination=None):
    """
    Resolve `-w <taxon>` to a list of assembly accessions plus the RefGenomeSelection
    it came from (for warnings / provenance the caller may want to surface).

    Parameters
    ----------
    source : str
        The driver's -S value ('GTDB' or 'NCBI'; case-insensitive here).
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
    key = str(source).strip().lower()
    if key not in _SOURCE_ASSETS:
        raise WantedRefTaxError(
            f"'{source}' is not a recognized --source for --wanted-ref-tax "
            f"(expected one of: {', '.join(s.upper() for s in _SOURCE_ASSETS)}).")

    core_source, table_path_fn = _SOURCE_ASSETS[key]
    table_path = table_path_fn()

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

    selection.warnings.extend(
        size_advice(len(selection.accessions),
                    selection.resolved_rank,
                    selection.effective_derep_rank))

    return selection.accessions, selection
