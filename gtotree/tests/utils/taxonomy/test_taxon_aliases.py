"""
Unit tests for the taxon-name aliases (tax_ranks.TAXON_ALIASES).

These are applied in tax_select.find_ranks_for_taxon(), which is the single funnel
every surface resolves a `-w` name through -- so specifying them here is specifying
them for the main driver, both get-accs helpers, gen-scg-hmms, search-annotations at once.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, normalize_taxon_name
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, find_ranks_for_taxon,
                                               resolve_taxon)

_EXTRA_COLS = ["assembly_level", "refseq_category", "checkm_completeness",
               "checkm_contamination", "genome_size", "genome_size_ungapped",
               "contig_count"]


def _rec(acc, lineage):
    d = {"assembly_accession": acc, "assembly_level": "Complete Genome",
         "refseq_category": "", "checkm_completeness": "99.0",
         "checkm_contamination": "0.5", "genome_size": "4000000",
         "genome_size_ungapped": "4000000", "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


_RECORDS = [
    _rec("GCF_000000001.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam", "Bgen", "Bgen sp1")),
    _rec("GCF_000000002.1", ("Archaea", "Aphy", "Aclass", "Aord", "Afam", "Agen", "Agen sp1")),
    _rec("GCF_000000003.1", ("Eukaryota", "Ephy", "Eclass", "Eord", "Efam", "Egen", "Egen sp1")),
]


@pytest.fixture
def ncbi_path(tmp_path):
    path = tmp_path / "ncbi-data.parquet"
    keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in _RECORDS] for k in keys}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}), str(path))
    return str(path)


class TestNormalizeTaxonName:

    @pytest.mark.parametrize("spelling", ["eukarya", "Eukarya", "eukaryotes"])
    def test_eukaryote_spellings_map_to_the_ncbi_name(self, spelling):
        assert normalize_taxon_name(spelling) == "Eukaryota"

    def test_unknown_names_pass_through(self):
        assert normalize_taxon_name("Nitrospirota") == "Nitrospirota"
        assert normalize_taxon_name("  Nitrospirota  ") == "Nitrospirota"


class TestAliasesResolveAgainstTheAsset:

    @pytest.mark.parametrize("spelling", ["eukarya", "eukaryotes"])
    def test_aliases_resolve_to_the_canonical_asset_name(self, ncbi_path, spelling):
        assert resolve_taxon(ncbi_path, spelling) == ("Eukaryota", "domain")

    def test_find_ranks_returns_the_assets_spelling_not_the_users(self, ncbi_path):
        canonical, ranks = find_ranks_for_taxon(ncbi_path, "eukarya")
        assert canonical == "Eukaryota"
        assert ranks == ["domain"]

    def test_an_alias_for_something_absent_still_fails_cleanly(self, tmp_path):
        """
        GTDB has no eukaryotes. `-w eukarya` there must still be a plain
        "doesn't exist at any rank", not a crash.
        """
        path = tmp_path / "prokaryotes-only.parquet"
        keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
        records = _RECORDS[:2]
        cols = {k: [rec[k] for rec in records] for k in keys}
        pq.write_table(
            pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
            str(path))

        with pytest.raises(TaxonNotFound):
            resolve_taxon(str(path), "eukarya")
