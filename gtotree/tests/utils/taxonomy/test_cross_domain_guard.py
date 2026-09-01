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
from unittest.mock import patch

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


################################################################################
# --ncbi-section: which part of NCBI the `-w` path draws from.
#
# resolve_wanted_ref_tax_accessions() previously never passed accession_prefixes to
# the selection core, so there was no way to restrict a `-w` pull to one section.
#
# The DEFAULT stays "both" (no filter), which is what the `-w` path has always done.
# Restricting to "refseq" is a genuine narrowing rather than mere de-duplication:
# GenBank-only assemblies vanish, taking their whole group with them. Near-duplicate
# GCF_/GCA_ pairs are already collapsed by dereplication when it is on.
################################################################################

class TestNcbiSectionScopesTheWantedRefTaxPath:

    @pytest.fixture
    def both_sections(self, tmp_path):
        """One genus holding a RefSeq and a GenBank copy of the same assembly."""
        recs = [
            _rec("GCF_000000010.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam",
                                     "Dup", "Dup sp1")),
            _rec("GCA_000000010.1", ("Bacteria", "Bphy", "Bclass", "Bord", "Bfam",
                                     "Dup", "Dup sp1")),
        ]
        path = tmp_path / "ncbi-data.parquet"
        keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
        cols = {k: [r[k] for r in recs] for k in keys}
        pq.write_table(
            pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
            str(path))
        return str(path)

    def _resolve(self, path, section):
        from gtotree.utils.taxonomy import wanted_ref_tax as wrt
        with patch.object(wrt, "_table_path_for_source",
                          return_value=("ncbi", path)), \
             patch.object(wrt, "ncbi_data_table_path", return_value=path):
            accs, _sel = wrt.resolve_wanted_ref_tax_accessions(
                "ncbi", "Dup", derep_rank="off", ncbi_section=section)
        return accs

    def test_refseq_takes_only_the_gcf_copy(self, both_sections):
        assert self._resolve(both_sections, "refseq") == ["GCF_000000010.1"]

    def test_genbank_takes_only_the_gca_copy(self, both_sections):
        assert self._resolve(both_sections, "genbank") == ["GCA_000000010.1"]

    def test_both_is_the_default_and_returns_the_pair_when_derep_is_off(self, both_sections):
        assert sorted(self._resolve(both_sections, "both")) == \
            ["GCA_000000010.1", "GCF_000000010.1"]

    def test_dereplication_collapses_the_pair_without_any_section_filter(self, both_sections):
        """
        The GCF_/GCA_ records of one assembly share a lineage, so with derep on they
        land in the same group and only one survives. This is why the default doesn't
        need to be "refseq" -- de-duplication is derep's job, not this flag's.
        """
        from gtotree.utils.taxonomy import wanted_ref_tax as wrt
        with patch.object(wrt, "_table_path_for_source",
                          return_value=("ncbi", both_sections)), \
             patch.object(wrt, "ncbi_data_table_path", return_value=both_sections):
            accs, _sel = wrt.resolve_wanted_ref_tax_accessions(
                "ncbi", "Dup", derep_rank="species", ncbi_section="both")
        assert len(accs) == 1

    def test_refseq_drops_genbank_only_genomes_and_their_whole_group(self, tmp_path):
        """
        The reason "refseq" is NOT the default. A GenBank-only genome (a MAG never
        promoted to RefSeq) is the sole member of its family; filtering to RefSeq
        before grouping removes that family from the selection entirely.
        """
        recs = [
            _rec("GCF_000000010.1", ("Bacteria", "Bphy", "Bcls", "Bord", "Fam1",
                                     "G1", "G1 sp1")),
            _rec("GCA_000000020.1", ("Bacteria", "Bphy", "Bcls", "Bord", "Fam2",
                                     "G2", "G2 sp1")),
        ]
        path = tmp_path / "ncbi-data.parquet"
        keys = ["assembly_accession"] + list(RANKS) + _EXTRA_COLS
        cols = {k: [r[k] for r in recs] for k in keys}
        pq.write_table(
            pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
            str(path))

        from gtotree.utils.taxonomy import wanted_ref_tax as wrt
        def pull(section):
            with patch.object(wrt, "_table_path_for_source",
                              return_value=("ncbi", str(path))), \
                 patch.object(wrt, "ncbi_data_table_path", return_value=str(path)):
                accs, _ = wrt.resolve_wanted_ref_tax_accessions(
                    "ncbi", "Bphy", derep_rank="family", ncbi_section=section)
            return sorted(accs)

        assert pull("both") == ["GCA_000000020.1", "GCF_000000010.1"]
        # Fam2 is gone entirely, not just de-duplicated
        assert pull("refseq") == ["GCF_000000010.1"]

    def test_gtdb_ignores_the_section_rather_than_filtering_everything_out(self):
        """
        GTDB's accession column is GenBank-based, so a 'refseq' prefix filter would
        empty the pool. section_prefixes() must return None for GTDB.
        """
        from gtotree.utils.taxonomy.wanted_ref_tax import section_prefixes
        assert section_prefixes("gtdb", "refseq") is None
        assert section_prefixes("gtdb", "both") is None
        assert section_prefixes("ncbi", "refseq") == ("GCF_",)


class TestNcbiSectionIsWiredIntoTheDriver:

    def test_main_gtotree_parser_accepts_it(self):
        from gtotree.cli.parser import parser
        args = parser().parse_args(["-w", "Bacteria", "--ncbi-section", "genbank",
                                    "-o", "out"])
        assert args.ncbi_section == "genbank"

    def test_main_gtotree_default_is_both_preserving_prior_behavior(self):
        from gtotree.cli.parser import parser
        assert parser().parse_args(["-w", "Bacteria", "-o", "out"]).ncbi_section == "both"

    def test_it_is_in_the_resume_fingerprint(self):
        """Changing which section genomes come from must invalidate a resume."""
        from gtotree.utils.misc.preflight_checks import RESUME
        assert "ncbi_section" in RESUME.field_labels


class TestNcbiSectionOnEveryWantedRefTaxSurface:
    """
    Every surface that takes `-w` routes through resolve_wanted_ref_tax_accessions,
    so every one of them needs the flag -- otherwise `--source ncbi` silently pulls
    GCF_/GCA_ near-duplicates on that surface.
    """

    def test_gen_scg_hmms_parser(self):
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import build_parser
        args = build_parser().parse_args(["-w", "Bacillus", "--ncbi-section", "genbank"])
        assert args.ncbi_section == "genbank"

    def test_gen_scg_hmms_default_is_genbank(self):
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import build_parser
        assert build_parser().parse_args(["-w", "Bacillus"]).ncbi_section == "genbank"

    def test_gen_scg_hmms_source_flag_still_works(self):
        """The --source arg sits next to the new one; make sure it survived."""
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import build_parser
        args = build_parser().parse_args(["-w", "Bacillus", "--source", "ncbi"])
        assert args.source == "ncbi"

    def test_search_annotations_parser(self):
        from gtotree.utils.target_search.target_search_cli import build_parser
        args = build_parser().parse_args(
            ["-p", "pfams.txt", "-w", "Bacillus", "--ncbi-section", "both"])
        assert args.ncbi_section == "both"

    def test_search_annotations_default_is_both(self):
        from gtotree.utils.target_search.target_search_cli import build_parser
        args = build_parser().parse_args(["-p", "pfams.txt", "-w", "Bacillus"])
        assert args.ncbi_section == "both"

    def test_get_accs_from_ncbi_keeps_its_long_standing_refseq_default(self):
        """
        The one surface that should NOT change: it has always defaulted to refseq and
        defaults --derep-rank to off, so its section filter is load-bearing.
        """
        from gtotree.utils.ncbi.get_accessions_from_ncbi import build_parser
        args = build_parser().parse_args(["-w", "Bacillus"])
        assert args.ncbi_section == "refseq"
        assert args.derep_rank == "off"

    def test_search_annotations_source_flag_still_works(self):
        from gtotree.utils.target_search.target_search_cli import build_parser
        args = build_parser().parse_args(
            ["-p", "pfams.txt", "-w", "Bacillus", "--source", "ncbi"])
        assert args.source == "ncbi"

    def test_dl_ncbi_assemblies_parser(self):
        from gtotree.utils.ncbi.dl_ncbi_assemblies import build_parser
        assert build_parser().parse_args(
            ["-w", "Bacillus", "--ncbi-section", "genbank"]).ncbi_section == "genbank"

    def test_it_is_in_both_subcommand_resume_fingerprints(self):
        from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_cli import RESUME as SCG
        from gtotree.utils.target_search.target_search_cli import RESUME as SEARCH
        assert "ncbi_section" in SCG.field_labels
        assert "ncbi_section" in SEARCH.field_labels

    def test_the_shared_resolver_receives_it(self):
        """
        gen-scg-hmms and search-annotations both select genomes through
        general.resolve_input_genomes, so the flag has to reach the resolver there
        rather than at each call site.
        """
        import inspect
        from gtotree.utils.misc.general import resolve_input_genomes
        src = inspect.getsource(resolve_input_genomes)
        assert "ncbi_section=" in src
