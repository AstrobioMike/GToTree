"""
`--exclusion-list` handling, shared by every surface that pulls genomes by taxonomy.

An exclusion list names assembly accessions that must not be used. It is applied to
the CANDIDATE POOL before dereplication and before any best-per-group pick

Matching is on the accession core alone

The list only ever constrains genomes selected BY TAXONOMY. Accessions a user names
directly (`-a`) are always used as provided.
"""

from gtotree.utils.taxonomy.tax_ranks import accession_core


def exclusion_list_help(taxon_flag="-w"):
    """
    The `--exclusion-list` help string, reused in multiple places
    """
    return (f"single-column file of assembly accessions to exclude from what "
            f"`{taxon_flag}` pulls")


def read_exclusion_list(path):
    entries = []
    seen = set()
    with open(path) as f:
        for line in f:
            entry = line.strip()
            if not entry or entry.startswith("#") or entry in seen:
                continue
            seen.add(entry)
            entries.append(entry)
    return entries


def load_exclusion_cores(path):
    """
    An exclusion-list file -> the set of accession cores it names.
    """
    if not path:
        return set()

    cores = set()
    for entry in read_exclusion_list(path):
        core = accession_core(entry)
        if core:
            cores.add(core)
    return cores


def filter_rows_by_exclusion(rows, acc_col, exclude_cores):
    """
    Drop excluded genomes from a list of candidate metadata rows.

    Returns (kept_rows, num_excluded), preserving order. This is the form the
    selection core works in; see filter_accessions_by_exclusion for the plain
    accession-list form.
    """
    rows = list(rows or [])
    if not exclude_cores:
        return rows, 0

    kept = [r for r in rows
            if accession_core(r.get(acc_col) or "") not in exclude_cores]
    return kept, len(rows) - len(kept)


def filter_accessions_by_exclusion(accessions, exclude_cores):
    """
    Drop excluded accessions from a plain accession list.

    Returns (kept, num_excluded), preserving order. For the surfaces that hold
    accessions rather than metadata rows at the point the exclusion applies.
    """
    accessions = list(accessions or [])
    if not exclude_cores:
        return accessions, 0

    kept = [acc for acc in accessions
            if accession_core(acc) not in exclude_cores]
    return kept, len(accessions) - len(kept)


def filter_table_by_exclusion(table, acc_col, exclude_cores):
    """
    Drop excluded genomes from an Arrow table of candidates.

    Returns (kept_table, num_excluded). For the bulk `all` paths that read the asset
    directly rather than going through the selection core; the core's own paths use
    filter_rows_by_exclusion instead. Kept here so all three forms normalise
    accessions the same way.
    """
    if not exclude_cores:
        return table, 0

    accs = table.column(acc_col).to_pylist()
    keep = [i for i, a in enumerate(accs)
            if accession_core(a or "") not in exclude_cores]
    num_excluded = table.num_rows - len(keep)
    if num_excluded:
        table = table.take(keep)
    return table, num_excluded


def exclusion_warning(num_excluded):
    """
    The advisory line for a pool the exclusion list shrank, or None if it removed none.
    """
    if not num_excluded:
        return None
    word = "genome" if num_excluded == 1 else "genomes"
    verb = "was" if num_excluded == 1 else "were"
    return (f"{num_excluded:,} candidate {word} {verb} removed by the "
            f"--exclusion-list before selection.")
