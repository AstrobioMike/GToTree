"""
`-w all` on the driver side (the main GToTree run, gen-scg-hmms, search-annotations).

'all' is not a taxon and has no rank of its own, so it can't go through resolve_taxon()
-- it used to reach it anyway and die with "doesn't exist at any rank in the gtdb
taxonomy". It is expanded to one target per DOMAIN in the source's asset instead, which
keeps every downstream piece (per-taxon warnings, SCG-set auto-selection, per-source
iTOL files) working with no 'all'-shaped special case in it.
"""

import types

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA
import gtotree.utils.taxonomy.wanted_ref_tax as wrt
from gtotree.utils.taxonomy.wanted_ref_tax import (WantedRefTaxError,
                                                   describe_all_expansion,
                                                   expand_wanted_ref_tax,
                                                   resolve_wanted_ref_tax_accessions)

_NCBI_EXTRA = ["assembly_level", "refseq_category", "checkm_completeness",
               "checkm_contamination", "genome_size", "genome_size_ungapped",
               "contig_count"]
_GTDB_EXTRA = ["gtdb_representative", "ncbi_refseq_category", "checkm2_completeness",
               "checkm2_contamination", "genome_size", "contig_count"]


def _ncbi_rec(acc, lineage):
    d = {"assembly_accession": acc, "assembly_level": "Complete Genome",
         "refseq_category": "", "checkm_completeness": "99.0",
         "checkm_contamination": "0.5", "genome_size": "4000000",
         "genome_size_ungapped": "4000000", "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


def _gtdb_rec(acc, lineage):
    d = {"ncbi_genbank_assembly_accession": acc, "gtdb_representative": "t",
         "ncbi_refseq_category": "", "checkm2_completeness": "99.0",
         "checkm2_contamination": "0.5", "genome_size": "4000000",
         "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


def _write(path, records, acc_col, extra_cols):
    keys = [acc_col] + list(RANKS) + extra_cols
    cols = {k: [rec[k] for rec in records] for k in keys}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}), str(path))
    return str(path)


@pytest.fixture
def assets(tmp_path, monkeypatch):
    """Point both source resolvers at small mock tables."""
    ncbi = _write(tmp_path / "ncbi-data.parquet", [
        _ncbi_rec("GCF_000000001.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam",
                                      "Bgen", "Bgen sp1")),
        _ncbi_rec("GCF_000000002.1", ("Archaea", "Aphy", "Aclass", "Aord", "Afam",
                                      "Agen", "Agen sp1")),
        _ncbi_rec("GCF_000000003.1", ("Eukaryota", "Ephy", "Eclass", "Eord", "Efam",
                                      "Egen", "Egen sp1")),
        _ncbi_rec("GCF_000000004.1", (NA, "Uroviricota", "Caudoviricetes", NA, NA, NA,
                                      "Phage sp1")),
    ], "assembly_accession", _NCBI_EXTRA)

    gtdb = _write(tmp_path / "gtdb-data.parquet", [
        _gtdb_rec("GCA_000000001.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam",
                                      "Bgen", "Bgen sp1")),
        _gtdb_rec("GCA_000000002.1", ("Archaea", "Aphy", "Aclass", "Aord", "Afam",
                                      "Agen", "Agen sp1")),
    ], "ncbi_genbank_assembly_accession", _GTDB_EXTRA)

    monkeypatch.setitem(wrt._SOURCE_ASSETS, "ncbi", ("ncbi", lambda: ncbi))
    monkeypatch.setitem(wrt._SOURCE_ASSETS, "gtdb", ("gtdb", lambda: gtdb))
    # the GTDB path screens its picks against the NCBI table for liveness, and reaches
    # for that path directly rather than through _SOURCE_ASSETS
    monkeypatch.setattr(wrt, "ncbi_data_table_path", lambda: ncbi)
    return types.SimpleNamespace(ncbi=ncbi, gtdb=gtdb)


class TestExpansion:

    def test_gtdb_all_is_the_two_prokaryotic_domains(self, assets):
        taxa, domains = expand_wanted_ref_tax("gtdb", ["all"])
        assert taxa == ["Archaea", "Bacteria"]
        assert domains == taxa

    def test_ncbi_all_includes_eukaryota(self, assets):
        taxa, _domains = expand_wanted_ref_tax("ncbi", ["all"])
        assert taxa == ["Archaea", "Bacteria", "Eukaryota"]

    def test_the_domain_list_comes_from_the_asset_not_a_hardcode(self, assets):
        """
        The two sources differ only because their tables differ -- neither pair is
        written down anywhere, so an asset that gains a domain is picked up.
        """
        assert expand_wanted_ref_tax("gtdb", ["all"])[0] != \
            expand_wanted_ref_tax("ncbi", ["all"])[0]

    def test_domainless_rows_do_not_become_a_domain(self, assets):
        """The viral row has no domain, so 'all' can't pull it (see tax_targets)."""
        assert NA not in expand_wanted_ref_tax("ncbi", ["all"])[0]

    def test_repeated_w_is_preserved_and_deduped(self, assets):
        taxa, _domains = expand_wanted_ref_tax("gtdb", ["all", "Bacteria"])
        assert taxa == ["Archaea", "Bacteria"]

    def test_no_expansion_without_all(self, assets):
        taxa, domains = expand_wanted_ref_tax("gtdb", ["Bphy"])
        assert taxa == ["Bphy"]
        assert domains == []

    def test_unknown_source_is_a_wanted_ref_tax_error(self, assets):
        with pytest.raises(WantedRefTaxError):
            expand_wanted_ref_tax("nonsense", ["all"])


class TestExpansionNote:

    def test_note_names_the_source_and_the_domains(self, assets):
        _taxa, domains = expand_wanted_ref_tax("ncbi", ["all"])
        note = describe_all_expansion("ncbi", domains)
        assert "NCBI" in note
        assert "Archaea, Bacteria, Eukaryota" in note

    def test_no_note_when_nothing_was_expanded(self):
        assert describe_all_expansion("gtdb", []) is None


class TestSelectionAfterExpansion:

    def test_each_expanded_domain_resolves_like_any_other_target(self, assets):
        taxa, _domains = expand_wanted_ref_tax("gtdb", ["all"])
        for taxon in taxa:
            accessions, selection = resolve_wanted_ref_tax_accessions(
                "gtdb", taxon, derep_rank="class")
            assert accessions
            assert selection.resolved_rank == "domain"
            assert selection.effective_derep_rank == "class"
