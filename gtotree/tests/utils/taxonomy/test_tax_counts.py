"""
Unit tests for gtotree/utils/taxonomy/tax_counts.py.

The load-bearing property specified here is the one the counts surfaces rely on:
derep_size() must equal the number of accessions an equivalent select_ref_genomes()
pull actually returns. If those two ever drift, `--get-taxon-counts` starts quoting a
number the pull won't produce, which is worse than not reporting it at all.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA
from gtotree.utils.taxonomy.tax_counts import (count_distinct_taxa, count_genomes,
                                               derep_size, rank_counts,
                                               render_rank_count_table,
                                               representatives_filter)
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes

_EXTRA_COLS = ["gtdb_representative", "ncbi_refseq_category",
               "checkm2_completeness", "checkm2_contamination",
               "genome_size", "contig_count"]


def _rec(acc, lineage, rep="t", refseq="", comp="99.0"):
    d = {"ncbi_genbank_assembly_accession": acc,
         "gtdb_representative": rep, "ncbi_refseq_category": refseq,
         "checkm2_completeness": comp, "checkm2_contamination": "0.5",
         "genome_size": "4000000", "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


# Two classes under Testophyla. ClassB's genus is unassigned ("NA") on one genome,
# which derep() skips -- so the distinct-genus count must skip it too.
_RECORDS = [
    _rec("GCA_000000001.1", ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", "GenA sp1"),
         refseq="reference genome"),
    _rec("GCA_000000002.1", ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", "GenA sp2")),
    _rec("GCA_000000003.1", ("Bacteria", "Testophyla", "ClassB", "OrdB", "FamB", "GenB", "GenB sp1")),
    _rec("GCA_000000004.1", ("Bacteria", "Testophyla", "ClassB", "OrdB", "FamB", NA, NA),
         rep="f"),
]


@pytest.fixture
def gtdb_path(tmp_path):
    path = tmp_path / "gtdb-data.parquet"
    keys = ["ncbi_genbank_assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in _RECORDS] for k in keys}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}), str(path))
    return str(path)


class TestCountGenomes:

    def test_counts_the_whole_table_when_unscoped(self, gtdb_path):
        assert count_genomes(gtdb_path, "gtdb") == 4

    def test_scopes_to_a_taxon(self, gtdb_path):
        assert count_genomes(gtdb_path, "gtdb", rank="class", taxon="ClassA") == 2

    def test_applies_a_representatives_filter(self, gtdb_path):
        rep_filter = representatives_filter("gtdb", "source")
        # GCA_000000004.1 has gtdb_representative="f"
        assert count_genomes(gtdb_path, "gtdb", rep_filter=rep_filter) == 3

    def test_applies_an_accession_prefix_filter(self, gtdb_path):
        assert count_genomes(gtdb_path, "gtdb", accession_prefixes=("GCF_",)) == 0
        assert count_genomes(gtdb_path, "gtdb", accession_prefixes=("GCA_",)) == 4


class TestRepresentativesFilter:

    def test_none_means_no_filter(self):
        assert representatives_filter("gtdb", None) is None

    def test_source_kind_is_the_sources_own_flag(self):
        assert representatives_filter("gtdb", "source") == ("gtdb_representative", "t")

    def test_refseq_kind_is_the_refseq_reference_column(self):
        """The GTDB helper compared against a literal "RefSeq" that never matched."""
        assert representatives_filter("gtdb", "refseq") == ("ncbi_refseq_category",
                                                            "reference genome")

    def test_an_unknown_kind_is_rejected(self):
        with pytest.raises(ValueError):
            representatives_filter("gtdb", "gtdb")


class TestCountDistinctTaxa:

    def test_counts_unique_values_of_a_rank(self, gtdb_path):
        assert count_distinct_taxa(gtdb_path, "gtdb", "class") == 2

    def test_unassigned_values_are_not_counted_as_a_taxon(self, gtdb_path):
        """One genome has genus "NA"; that is not a genus."""
        assert count_distinct_taxa(gtdb_path, "gtdb", "genus") == 2

    def test_scoping_restricts_the_pool(self, gtdb_path):
        assert count_distinct_taxa(gtdb_path, "gtdb", "genus",
                                   scope_rank="class", scope_taxon="ClassA") == 1


class TestDerepSize:
    """The count must equal what an equivalent pull actually returns."""

    @pytest.mark.parametrize("derep_rank", ["class", "order", "family", "genus"])
    def test_matches_the_number_of_accessions_a_pull_returns(self, gtdb_path, derep_rank):
        predicted = derep_size(gtdb_path, "gtdb", "phylum", "Testophyla", derep_rank)
        selection = select_ref_genomes(gtdb_path, "gtdb", "Testophyla",
                                       derep_rank=derep_rank, reps_only=False)
        assert predicted == len(selection.accessions)

    def test_matches_a_pull_restricted_to_representatives(self, gtdb_path):
        rep_filter = representatives_filter("gtdb", "source")
        predicted = derep_size(gtdb_path, "gtdb", "phylum", "Testophyla", "genus",
                               rep_filter=rep_filter)
        selection = select_ref_genomes(gtdb_path, "gtdb", "Testophyla",
                                       derep_rank="genus", reps_only=True)
        assert predicted == len(selection.accessions)


class TestRankCounts:

    def test_unscoped_covers_every_rank(self, gtdb_path):
        assert [r for r, _ in rank_counts(gtdb_path, "gtdb")] == list(RANKS)

    def test_scoped_starts_at_the_taxons_own_rank(self, gtdb_path):
        rows = rank_counts(gtdb_path, "gtdb", scope_rank="class", scope_taxon="ClassA")
        assert [r for r, _ in rows] == ["class", "order", "family", "genus", "species"]

    def test_the_taxons_own_rank_reports_itself(self, gtdb_path):
        rows = dict(rank_counts(gtdb_path, "gtdb",
                                scope_rank="class", scope_taxon="ClassA"))
        assert rows["class"] == 1


def test_render_rank_count_table_lays_out_rank_and_count():
    text = render_rank_count_table([("genus", 1234)], count_header="Num. Unique Taxa")
    assert "Rank" in text
    assert "Num. Unique Taxa" in text
    assert "genus" in text
    assert "1,234" in text
