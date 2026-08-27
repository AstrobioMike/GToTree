"""
Cross-domain taxon guard, and --target-domain to disambiguate.

Taxonomic names are not unique across domains in NCBI: `Bacillus`, for instance, is
both a bacterial genus and a stick-insect (eukaryotic) genus. Selecting on the bare
name pulled genomes from BOTH domains into one alignment/tree -- a silent, chimeric
result nobody asked for.

resolve_taxon() now refuses such a name (CrossDomainTaxon) unless --target-domain says
which domain is wanted. When it is given, resolution scopes to that domain, selection
pulls only that domain's genomes, and (because the pulled rows carry the domain) the
existing HMM-set auto-selection routes euk targets to the Universal set for free.

These tests use a mock asset where the genus 'Shared' occurs in both Bacteria and
Eukaryota, mirroring the real `Bacillus` collision.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, NA
from gtotree.utils.taxonomy.tax_select import (TaxonNotFound, CrossDomainTaxon,
                                               resolve_taxon)
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes

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
    # 'Shared' is a genus in BOTH Bacteria and Eukaryota (like `Bacillus`)
    _rec("GCF_000000004.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam", "Shared", "Shared sp1")),
    _rec("GCF_000000005.1", ("Eukaryota", "Ephy", "Eclass", "Eord", "Efam", "Shared", "Shared sp2")),
    # a domain-less viral row
    _rec("GCF_000000006.1", (NA, "Uroviricota", "Caudoviricetes", NA, NA, NA, "phage sp1")),
]


@pytest.fixture
def ncbi_path(tmp_path):
    path = tmp_path / "ncbi-data.parquet"
    keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in _RECORDS] for k in keys}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}), str(path))
    return str(path)


class TestResolveTaxonReturnsDomain:

    def test_a_single_domain_taxon_returns_its_domain(self, ncbi_path):
        assert resolve_taxon(ncbi_path, "Bgen") == ("Bgen", "genus", "Bacteria")

    def test_a_domain_target_returns_itself_as_domain(self, ncbi_path):
        assert resolve_taxon(ncbi_path, "Eukaryota") == ("Eukaryota", "domain", "Eukaryota")

    def test_a_euk_lineage_target_returns_eukaryota(self, ncbi_path):
        assert resolve_taxon(ncbi_path, "Ephy") == ("Ephy", "phylum", "Eukaryota")

    def test_a_domainless_taxon_returns_none_domain(self, ncbi_path):
        # the viral row carries no assigned domain
        canonical, rank, domain = resolve_taxon(ncbi_path, "Caudoviricetes")
        assert canonical == "Caudoviricetes"
        assert domain is None


class TestCrossDomainGuard:

    def test_a_name_in_two_domains_is_refused(self, ncbi_path):
        with pytest.raises(CrossDomainTaxon) as exc:
            resolve_taxon(ncbi_path, "Shared")
        assert exc.value.domains_found == ["Bacteria", "Eukaryota"]

    def test_the_message_names_the_domains_and_the_flag(self, ncbi_path):
        with pytest.raises(CrossDomainTaxon) as exc:
            resolve_taxon(ncbi_path, "Shared")
        msg = str(exc.value)
        assert "Bacteria" in msg and "Eukaryota" in msg
        assert "--target-domain" in msg

    def test_the_guard_also_fires_through_select_ref_genomes(self, ncbi_path):
        with pytest.raises(CrossDomainTaxon):
            select_ref_genomes(ncbi_path, "ncbi", "Shared", derep_rank="off")


class TestTargetDomainDisambiguates:

    def test_target_domain_picks_one_side(self, ncbi_path):
        assert resolve_taxon(ncbi_path, "Shared", None, "Bacteria") == \
            ("Shared", "genus", "Bacteria")
        assert resolve_taxon(ncbi_path, "Shared", None, "Eukaryota") == \
            ("Shared", "genus", "Eukaryota")

    def test_target_domain_accepts_aliases(self, ncbi_path):
        # 'eukarya' -> 'Eukaryota' via the same alias table -w uses
        assert resolve_taxon(ncbi_path, "Shared", None, "eukarya") == \
            ("Shared", "genus", "Eukaryota")

    def test_target_domain_is_case_insensitive(self, ncbi_path):
        assert resolve_taxon(ncbi_path, "Shared", None, "bacteria") == \
            ("Shared", "genus", "Bacteria")

    def test_a_domain_the_taxon_isnt_in_is_rejected(self, ncbi_path):
        with pytest.raises(TaxonNotFound) as exc:
            resolve_taxon(ncbi_path, "Shared", None, "Archaea")
        assert "Archaea" in str(exc.value)

    def test_target_domain_on_a_single_domain_taxon_is_fine(self, ncbi_path):
        # redundant but harmless: naming the domain a taxon is already unambiguously in
        assert resolve_taxon(ncbi_path, "Bgen", None, "Bacteria") == \
            ("Bgen", "genus", "Bacteria")

    def test_target_domain_that_contradicts_a_single_domain_taxon_is_rejected(self, ncbi_path):
        with pytest.raises(TaxonNotFound):
            resolve_taxon(ncbi_path, "Bgen", None, "Eukaryota")


class TestTargetDomainScopesSelection:
    """
    The disambiguated selection must pull ONLY the chosen domain's genomes -- otherwise
    the tree would still be chimeric. Checked on the pulled rows' domains.
    """

    def test_bacteria_side_pulls_only_bacterial_rows(self, ncbi_path):
        sel = select_ref_genomes(ncbi_path, "ncbi", "Shared",
                                 target_domain="Bacteria", derep_rank="off")
        assert {r.get("domain") for r in sel.rows} == {"Bacteria"}
        assert len(sel.accessions) == 1

    def test_eukaryota_side_pulls_only_eukaryotic_rows(self, ncbi_path):
        sel = select_ref_genomes(ncbi_path, "ncbi", "Shared",
                                 target_domain="eukarya", derep_rank="off")
        assert {r.get("domain") for r in sel.rows} == {"Eukaryota"}
        assert len(sel.accessions) == 1

    def test_combined_with_target_rank(self, ncbi_path):
        # both disambiguators together
        sel = select_ref_genomes(ncbi_path, "ncbi", "Shared",
                                 target_rank="genus", target_domain="Bacteria",
                                 derep_rank="off")
        assert {r.get("domain") for r in sel.rows} == {"Bacteria"}


################################################################################
# CLI wiring: --target-domain is exposed on every -w surface and folded into each
# resume fingerprint so a resume invalidates when it changes.
################################################################################

class TestTargetDomainOnParsers:

    def test_main_gtotree_parser(self):
        from gtotree.cli.parser import parser
        args = parser().parse_args(["-w", "Bacillus", "--target-domain", "Bacteria",
                                    "-o", "out"])
        assert args.target_domain == "Bacteria"

    def test_main_gtotree_default_is_none(self):
        from gtotree.cli.parser import parser
        args = parser().parse_args(["-w", "Bacillus", "-o", "out"])
        assert args.target_domain is None

    def test_gen_scg_hmms_parser(self):
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import build_parser
        args = build_parser().parse_args(["-w", "Bacillus", "--target-domain", "eukarya"])
        assert args.target_domain == "eukarya"

    def test_search_annotations_parser(self):
        from gtotree.utils.target_search.target_search_cli import build_parser
        args = build_parser().parse_args(
            ["-p", "pfams.txt", "-w", "Bacillus", "--target-domain", "Bacteria"])
        assert args.target_domain == "Bacteria"

    def test_get_accs_shared_parser(self):
        import argparse
        from gtotree.utils.taxonomy.get_accs_shared import add_common_get_accs_args
        p = argparse.ArgumentParser()
        req = p.add_argument_group("req")
        opt = p.add_argument_group("opt")
        add_common_get_accs_args(req, opt, "ncbi")
        args = p.parse_args(["-w", "Bacillus", "--target-domain", "Bacteria"])
        assert args.target_domain == "Bacteria"


class TestTargetDomainInFingerprints:

    def test_main_gtotree_field_labels(self):
        from gtotree.utils.misc.preflight_checks import RESUME
        assert "target_domain" in RESUME.field_labels

    def test_gen_scg_hmms_field_labels(self):
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import RESUME
        assert "target_domain" in RESUME.field_labels

    def test_search_annotations_field_labels(self):
        from gtotree.utils.target_search.target_search_cli import RESUME
        assert "target_domain" in RESUME.field_labels
