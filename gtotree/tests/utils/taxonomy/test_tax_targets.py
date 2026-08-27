"""
Unit tests for gtotree/utils/taxonomy/tax_targets.py.

The property specified here is the one that broke: a rank-counts table and a `-w all`
pull have to describe the SAME pool. The NCBI asset carries rows with no assigned
domain (viral / metagenome entries); walking the domains can't reach them, so counting
them makes the table quote a number the pull will never produce.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA
from gtotree.utils.taxonomy.tax_counts import count_distinct_taxa, rank_counts
from gtotree.utils.taxonomy.tax_derep import select_all_domains
from gtotree.utils.taxonomy.tax_targets import (ALL_KEYWORD, domains_in_asset,
                                                expand_all_targets, is_all_target,
                                                unassigned_domain_summary)

_EXTRA_COLS = ["assembly_level", "refseq_category", "checkm_completeness",
               "checkm_contamination", "genome_size", "genome_size_ungapped",
               "contig_count"]


def _rec(acc, lineage, comp="99.0"):
    d = {"assembly_accession": acc, "assembly_level": "Complete Genome",
         "refseq_category": "", "checkm_completeness": comp,
         "checkm_contamination": "0.5", "genome_size": "4000000",
         "genome_size_ungapped": "4000000", "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


# Three domains plus two domain-less (viral) rows. The viral rows carry a class
# ("Caudoviricetes") that exists nowhere else, which is the shape of the real gap:
# 4 classes in the table, but only 3 reachable by walking the domains.
_RECORDS = [
    _rec("GCF_000000001.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam", "Bgen", "Bgen sp1")),
    _rec("GCF_000000002.1", ("Archaea", "Aphy", "Aclass", "Aord", "Afam", "Agen", "Agen sp1")),
    _rec("GCF_000000003.1", ("Eukaryota", "Ephy", "Eclass", "Eord", "Efam", "Egen", "Egen sp1")),
    # 'Xgen' is a genus in BOTH Bacteria and Eukaryota (as `Bacillus` is in the real
    # asset), so it belongs to no single domain
    _rec("GCF_000000006.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam", "Xgen", "Xgen sp1")),
    _rec("GCF_000000007.1", ("Eukaryota", "Ephy", "Eclass", "Eord", "Efam", "Xgen", "Xgen sp2")),
    _rec("GCF_000000004.1", (NA, "Uroviricota", "Caudoviricetes", NA, NA, NA, "phage sp1")),
    _rec("GCF_000000005.1", (NA, "Uroviricota", "Caudoviricetes", NA, NA, NA, "phage sp2")),
]


@pytest.fixture
def ncbi_path(tmp_path):
    path = tmp_path / "ncbi-data.parquet"
    keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in _RECORDS] for k in keys}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}), str(path))
    return str(path)


class TestIsAllTarget:

    @pytest.mark.parametrize("value", ["all", "ALL", " All "])
    def test_recognizes_the_keyword_in_any_casing(self, value):
        assert is_all_target(value)

    @pytest.mark.parametrize("value", ["Bacteria", "", None, "allomonas"])
    def test_leaves_everything_else_alone(self, value):
        assert not is_all_target(value)


class TestDomainsInAsset:

    def test_reads_the_domains_from_the_asset(self, ncbi_path):
        assert domains_in_asset(ncbi_path, "ncbi") == ["Archaea", "Bacteria", "Eukaryota"]

    def test_excludes_unassigned_domains(self, ncbi_path):
        """The domain-less viral rows are not a domain, so they mustn't become one."""
        assert NA not in domains_in_asset(ncbi_path, "ncbi")


class TestExpandAllTargets:

    def test_all_becomes_one_target_per_domain(self, ncbi_path):
        expanded, domains = expand_all_targets(ncbi_path, "ncbi", [ALL_KEYWORD])
        assert expanded == ["Archaea", "Bacteria", "Eukaryota"]
        assert domains == expanded

    def test_named_taxa_pass_through_untouched(self, ncbi_path):
        expanded, domains = expand_all_targets(ncbi_path, "ncbi", ["Bacteria", "Aphy"])
        assert expanded == ["Bacteria", "Aphy"]
        assert domains == []

    def test_all_alongside_a_named_domain_does_not_duplicate_it(self, ncbi_path):
        """`-w all -w bacteria` must resolve Bacteria once, not twice."""
        expanded, _domains = expand_all_targets(ncbi_path, "ncbi",
                                                [ALL_KEYWORD, "bacteria"])
        assert expanded == ["Archaea", "Bacteria", "Eukaryota"]

    def test_empty_input_is_not_an_expansion(self, ncbi_path):
        assert expand_all_targets(ncbi_path, "ncbi", []) == ([], [])


class TestUnassignedDomainSummary:

    def test_counts_the_genomes_all_cannot_reach(self, ncbi_path):
        summary = unassigned_domain_summary(ncbi_path, "ncbi")
        assert summary.n_genomes == 2
        assert bool(summary) is True

    def test_counts_the_taxa_that_exist_only_there(self, ncbi_path):
        """
        The number that reconciles the two views: 4 classes in the table, 3 reachable
        by walking the domains, 1 (Caudoviricetes) only in the domain-less rows.
        """
        summary = unassigned_domain_summary(ncbi_path, "ncbi", rank="class")
        assert summary.n_taxa_at_rank == 1

        unscoped = dict(rank_counts(ncbi_path, "ncbi"))["class"]
        scoped = dict(rank_counts(ncbi_path, "ncbi", domain_assigned=True))["class"]
        assert unscoped - scoped == summary.n_taxa_at_rank

    def test_message_names_the_numbers(self, ncbi_path):
        message = unassigned_domain_summary(ncbi_path, "ncbi", rank="class").message()
        assert "2 genome(s)" in message
        assert "1 'class' value(s)" in message

    def test_no_message_when_every_genome_has_a_domain(self, tmp_path):
        path = tmp_path / "clean.parquet"
        keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
        records = _RECORDS[:3]
        cols = {k: [rec[k] for rec in records] for k in keys}
        pq.write_table(
            pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
            str(path))

        summary = unassigned_domain_summary(str(path), "ncbi", rank="class")
        assert not summary
        assert summary.message() is None


class TestCountsReconcileWithThePull:

    def test_scoped_class_count_equals_what_all_actually_pulls(self, ncbi_path):
        """
        The regression this whole module exists for: `--get-rank-counts` scoped to
        'all' has to equal the number of accessions `-w all --derep-rank class`
        writes. Against the unscoped count it was 4 vs 3.
        """
        selection = select_all_domains(ncbi_path, "ncbi", derep_rank="class")
        scoped = count_distinct_taxa(ncbi_path, "ncbi", "class", domain_assigned=True)

        assert len(selection.accessions) == scoped

    def test_the_pull_reports_what_it_left_behind(self, ncbi_path):
        selection = select_all_domains(ncbi_path, "ncbi", derep_rank="class")
        assert selection.unassigned.n_genomes == 2
        assert selection.unassigned.rank == "class"
