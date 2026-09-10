"""
Optional lineage annotation for the dl-ncbi-assemblies info table.

The interesting seam is the GTDB join: GTDB keys assemblies on their GenBank (GCA_)
accession, so a RefSeq (GCF_) accession list has to match through accession_core()
rather than on the raw string. Mirrored as
bit/tests/taxonomy/test_lineage_lookup.py.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA, accession_core
from gtotree.utils.taxonomy.lineage_lookup import (GTDB_ACCESSION_COLUMN,
                                                 GTDB_LINEAGE_COLUMNS,
                                                 NCBI_LINEAGE_COLUMNS,
                                                 NO_LINEAGE,
                                                 gtdb_lineage_map,
                                                 lineage_columns,
                                                 lineage_from_row)


_ECOLI = ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
          "Enterobacteriaceae", "Escherichia", "Escherichia coli")


def _write_gtdb(tmp_path, rows):
    """rows: list of (accession, 7 rank names)."""
    cols = {GTDB_ACCESSION_COLUMN: [acc for acc, _ in rows]}
    for i, rank in enumerate(RANKS):
        cols[rank] = [lineage[i] for _, lineage in rows]
    path = tmp_path / "gtdb-data.parquet"
    pq.write_table(pa.table(cols), str(path))
    return path


# --- column naming --------------------------------------------------------

def test_lineage_columns_are_prefixed_and_ordered():
    assert lineage_columns(add_ncbi_tax=True) == NCBI_LINEAGE_COLUMNS
    assert lineage_columns(add_gtdb_tax=True) == GTDB_LINEAGE_COLUMNS
    assert lineage_columns(True, True) == NCBI_LINEAGE_COLUMNS + GTDB_LINEAGE_COLUMNS
    assert lineage_columns() == []


def test_column_names_never_collide():
    """
    Both taxonomies can be on at once, so a bare `phylum` would be a column whose
    values don't all mean the same thing.
    """
    assert not set(NCBI_LINEAGE_COLUMNS) & set(GTDB_LINEAGE_COLUMNS)
    assert NCBI_LINEAGE_COLUMNS[0] == "ncbi_domain"
    assert GTDB_LINEAGE_COLUMNS[-1] == "gtdb_species"


# --- row accessor ---------------------------------------------------------

def test_lineage_from_row_reads_all_ranks():
    row = dict(zip(RANKS, _ECOLI))
    assert lineage_from_row(row) == _ECOLI


def test_lineage_from_row_fills_missing_and_blank_with_NA():
    row = dict(zip(RANKS, _ECOLI))
    row["genus"] = ""
    del row["species"]
    got = lineage_from_row(row)
    assert got[RANKS.index("genus")] == NA
    assert got[RANKS.index("species")] == NA


# --- the GTDB join --------------------------------------------------------

def test_gtdb_lineage_map_basic_hit(tmp_path):
    path = _write_gtdb(tmp_path, [("GCA_000005845.2", _ECOLI)])
    mapping = gtdb_lineage_map(path, ["GCA_000005845.2"])
    assert mapping[accession_core("GCA_000005845.2")] == _ECOLI


def test_gcf_accession_matches_gtdb_gca_row(tmp_path):
    """The whole reason the join goes through accession_core()."""
    path = _write_gtdb(tmp_path, [("GCA_000005845.2", _ECOLI)])
    mapping = gtdb_lineage_map(path, ["GCF_000005845.3"])
    assert mapping[accession_core("GCF_000005845.3")] == _ECOLI


def test_version_differences_do_not_block_the_match(tmp_path):
    path = _write_gtdb(tmp_path, [("GCA_000005845.9", _ECOLI)])
    mapping = gtdb_lineage_map(path, ["GCA_000005845.1"])
    assert mapping[accession_core("GCA_000005845.1")] == _ECOLI


def test_accession_not_in_gtdb_is_simply_absent(tmp_path):
    """
    GTDB is bacteria and archaea only, so a euk/viral/unreleased genome legitimately
    has no row. Callers fill the gap with NO_LINEAGE rather than erroring.
    """
    path = _write_gtdb(tmp_path, [("GCA_000005845.2", _ECOLI)])
    mapping = gtdb_lineage_map(path, ["GCA_999999999.1"])
    assert mapping == {}
    assert mapping.get(accession_core("GCA_999999999.1"), NO_LINEAGE) == NO_LINEAGE
    assert len(NO_LINEAGE) == len(RANKS)


def test_only_wanted_rows_come_back(tmp_path):
    path = _write_gtdb(tmp_path, [
        ("GCA_000005845.2", _ECOLI),
        ("GCA_000001405.1", ("Bacteria",) + ("Other",) * 6),
        ("GCA_000002345.1", ("Archaea",) + ("Other",) * 6),
    ])
    mapping = gtdb_lineage_map(path, ["GCA_000005845.2", "GCA_000002345.1"])
    assert set(mapping) == {accession_core("GCA_000005845.2"),
                            accession_core("GCA_000002345.1")}


def test_junk_and_empty_accessions_are_tolerated(tmp_path):
    path = _write_gtdb(tmp_path, [("GCA_000005845.2", _ECOLI), ("na", (NA,) * 7)])
    assert gtdb_lineage_map(path, []) == {}
    assert gtdb_lineage_map(path, ["", "not-an-accession"]) == {}


def test_over_match_on_shared_digit_prefix_is_rejected(tmp_path):
    """'000000001' must not also match the longer '0000000019'."""
    path = _write_gtdb(tmp_path, [
        ("GCA_000000001.1", ("Bacteria",) + ("A",) * 6),
        ("GCA_0000000019.1", ("Bacteria",) + ("B",) * 6),
    ])
    mapping = gtdb_lineage_map(path, ["GCA_000000001.1"])
    assert set(mapping) == {"000000001"}
    assert mapping["000000001"][1] == "A"
