"""
Optional lineage annotation for the `dl-ncbi-assemblies` info table

This is mirrored in bit at bit/modules/taxonomy/lineage_lookup.py, mike,
apply any changes there too :)
"""

import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.dataset as ds # type: ignore
from gtotree.utils.taxonomy.tax_ranks import RANKS, NA, accession_core


# the GTDB asset's NCBI-accession column (GenBank/GCA_ side only)
GTDB_ACCESSION_COLUMN = "ncbi_genbank_assembly_accession"

NCBI_LINEAGE_COLUMNS = [f"ncbi_{rank}" for rank in RANKS]
GTDB_LINEAGE_COLUMNS = [f"gtdb_{rank}" for rank in RANKS]

# what a row gets when the taxonomy in question has no entry for it
NO_LINEAGE = tuple(NA for _ in RANKS)

_ACCESSION_CORE_REGEX = r"^(?:RS_|GB_)?GC[AF]_([0-9]+).*$"


def lineage_columns(add_ncbi_tax=False, add_gtdb_tax=False):
    """The lineage column names to append to the info table, in output order."""
    columns = []
    if add_ncbi_tax:
        columns.extend(NCBI_LINEAGE_COLUMNS)
    if add_gtdb_tax:
        columns.extend(GTDB_LINEAGE_COLUMNS)
    return columns


def clean_rank_value(value):
    if value is None:
        return NA
    value = str(value).strip()
    return value if value else NA


def lineage_from_row(row):
    """
    The seven rank values off a row of either asset.

    Both Parquet assets are built on the same lineage schema (bare names, no `d__`
    prefixes), so one accessor covers both.
    """
    return tuple(clean_rank_value(row.get(rank)) for rank in RANKS)


def _accession_cores(column):
    return pc.replace_substring_regex(column.cast(pa.string()),
                                      _ACCESSION_CORE_REGEX, r"\1")


def gtdb_lineage_map(gtdb_table_path, accessions):
    """
    {accession_core -> (domain, ..., species)} from the GTDB asset, for `accessions`.

    Only the wanted rows are ever materialized: the cores are computed and matched
    inside Arrow, so the full table never becomes Python objects.

    GTDB holds bacteria and archaea only, so anything eukaryotic or viral -- and any
    prokaryote not in the current release, or that didn't clear its QC -- simply won't
    be in the returned map. Callers fill those rows with NO_LINEAGE rather than
    treating a miss as an error; a partly-NA column is the honest answer there.
    """
    wanted_cores = {accession_core(acc) for acc in (accessions or [])}
    wanted_cores.discard("")
    if not wanted_cores:
        return {}

    table = ds.dataset(str(gtdb_table_path), format="parquet").to_table(
        columns=[GTDB_ACCESSION_COLUMN] + list(RANKS))

    cores = _accession_cores(table.column(GTDB_ACCESSION_COLUMN))
    mask = pc.is_in(cores,
                    value_set=pa.array(sorted(wanted_cores), type=pa.string()))

    table = table.filter(mask)
    cores = pc.filter(cores, mask).to_pylist()

    return {core: lineage_from_row(row)
            for core, row in zip(cores, table.to_pylist())}
