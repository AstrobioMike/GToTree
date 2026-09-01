"""
Domain-aware `auto` dereplication.

`-w all --source ncbi` includes Eukaryota, which is far more genome-heavy per rank
than the prokaryotic domains: dereplicated at the standard two-ranks-finer default it
contributes ~250 genomes, most of the selection. So any EUKARYOTIC-LINEAGE target's
`auto` default is one rank coarser -- one step finer than the target instead of two --
while Bacteria/Archaea keep the two-step default. This follows the target's domain, so
`-w Eukaryota`, a eukaryotic phylum, and a eukaryotic class are all coarsened.

The step count is chosen in resolve_derep_rank via a `domain` argument; select_ref_genomes
supplies it from resolve_taxon (which now returns the target's sole domain, or raises
CrossDomainTaxon for a name shared across domains). These tests pin the step logic
directly. The invariant that an EXPLICIT --derep-rank is never coarsened is covered here;
the guard that a name spanning >1 domain is refused (not silently scoped) is covered in
the resolution tests (test_cross_domain_guard.py).
"""

import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_derep import (
    DEREP_STEPS,
    DEREP_STEPS_BY_DOMAIN,
    default_derep_rank,
    default_derep_steps,
    resolve_derep_rank,
)


class TestDefaultDerepSteps:

    def test_unlisted_domain_uses_the_standard_step_count(self):
        assert default_derep_steps("Bacteria") == DEREP_STEPS
        assert default_derep_steps("Archaea") == DEREP_STEPS

    def test_none_domain_uses_the_standard_step_count(self):
        # a non-domain target (the domain isn't known / doesn't apply)
        assert default_derep_steps(None) == DEREP_STEPS

    def test_eukaryota_uses_one_fewer_step(self):
        assert default_derep_steps("Eukaryota") == 1
        assert DEREP_STEPS_BY_DOMAIN["Eukaryota"] == 1


class TestDefaultDerepRankIsDomainAware:

    def test_eukaryota_domain_target_defaults_to_phylum(self):
        assert default_derep_rank("domain", domain="Eukaryota") == "phylum"

    def test_prokaryote_domain_target_defaults_to_class(self):
        assert default_derep_rank("domain", domain="Bacteria") == "class"
        assert default_derep_rank("domain", domain="Archaea") == "class"

    def test_no_domain_keeps_the_standard_class_default(self):
        # unchanged behavior: a domain target with no domain hint is class
        assert default_derep_rank("domain") == "class"

    @pytest.mark.parametrize("target,expected", [
        ("domain", "phylum"),   # was class (two finer); euk -> one finer
        ("phylum", "class"),    # was order
        ("class", "order"),     # was family
        ("order", "family"),    # was genus
        ("family", "genus"),    # was species
        ("genus", "species"),   # clamped either way
    ])
    def test_eukaryotic_lineage_is_one_rank_finer_at_every_level(self, target, expected):
        assert default_derep_rank(target, domain="Eukaryota") == expected

    @pytest.mark.parametrize("target,expected", [
        ("domain", "class"),
        ("phylum", "order"),
        ("class", "family"),
    ])
    def test_prokaryotic_lineage_keeps_two_finer_at_every_level(self, target, expected):
        assert default_derep_rank(target, domain="Bacteria") == expected


class TestResolveDerepRankDomainAware:

    def test_auto_eukaryota_domain_resolves_to_phylum(self):
        assert resolve_derep_rank("domain", "auto", domain="Eukaryota") == ("phylum", [])

    def test_auto_prokaryote_domain_resolves_to_class(self):
        assert resolve_derep_rank("domain", "auto", domain="Bacteria") == ("class", [])
        assert resolve_derep_rank("domain", "auto", domain="Archaea") == ("class", [])

    def test_auto_without_domain_is_unchanged(self):
        # the original behavior for a domain target with no hint
        assert resolve_derep_rank("domain", "auto") == ("class", [])

    def test_explicit_rank_is_never_coarsened_for_eukaryota(self):
        # --derep-rank is domain-blind: a user who asks for class gets class
        assert resolve_derep_rank("domain", "class", domain="Eukaryota") == ("class", [])
        assert resolve_derep_rank("domain", "order", domain="Eukaryota") == ("order", [])

    def test_euk_lineage_below_domain_is_coarsened_one_step(self):
        """
        A eukaryotic phylum/class target is coarsened too (one step finer), because the
        step count follows the target's DOMAIN, not just a domain-level target.
        """
        assert resolve_derep_rank("phylum", "auto", domain="Eukaryota") == ("class", [])
        assert resolve_derep_rank("class", "auto", domain="Eukaryota") == ("order", [])

    def test_no_domain_hint_keeps_the_two_step_default(self):
        # when the caller can't resolve a sole domain, the ordinary default stands
        assert resolve_derep_rank("phylum", "auto", domain=None) == ("order", [])

    def test_off_is_unaffected_by_domain(self):
        assert resolve_derep_rank("domain", "off", domain="Eukaryota") == (None, [])


################################################################################
# include_rows: skipping metadata must not change WHICH genomes are selected, and
# the default must stay True because HMM auto-selection reads the pulled rows.
#
# Also pins why counts come from the selection path rather than the cheaper
# tax_counts.derep_size(): that counts distinct groups present in the table, so it
# sees neither liveness screening nor the quality floor, and over-reports.
################################################################################

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS
from gtotree.utils.taxonomy.tax_derep import select_ref_genomes
from gtotree.utils.taxonomy.tax_counts import derep_size


_NCBI_EXTRA = ["assembly_level", "refseq_category", "checkm_completeness",
               "checkm_contamination", "genome_size", "genome_size_ungapped",
               "contig_count"]


def _ncbi_rec(acc, lineage, completeness="99.0"):
    d = {"assembly_accession": acc, "assembly_level": "Complete Genome",
         "refseq_category": "", "checkm_completeness": completeness,
         "checkm_contamination": "0.5", "genome_size": "4000000",
         "genome_size_ungapped": "4000000", "contig_count": "1"}
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


@pytest.fixture
def three_families(tmp_path):
    """
    One phylum, three families of one genome each:
      Fam1 -- fine
      Fam2 -- its only genome is absent from the liveness-screen asset (suppressed)
      Fam3 -- its only genome is low completeness
    Returns (table_path, screen_path).
    """
    recs = [
        _ncbi_rec("GCF_000000010.1", ("Bacteria", "Bphy", "Bcls", "Bord", "Fam1",
                                      "G1", "G1 sp1")),
        _ncbi_rec("GCF_000000020.1", ("Bacteria", "Bphy", "Bcls", "Bord", "Fam2",
                                      "G2", "G2 sp1")),
        _ncbi_rec("GCF_000000030.1", ("Bacteria", "Bphy", "Bcls", "Bord", "Fam3",
                                      "G3", "G3 sp1"), completeness="40.0"),
    ]
    table = tmp_path / "ncbi-data.parquet"
    keys = ["assembly_accession"] + list(RANKS) + _NCBI_EXTRA
    cols = {k: [r[k] for r in recs] for k in keys}
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
        str(table))

    screen = tmp_path / "ncbi-screen.parquet"
    pq.write_table(
        pa.table({"assembly_accession": pa.array(
            ["GCF_000000010.1", "GCF_000000030.1"], type=pa.string())}), str(screen))

    return str(table), str(screen)


class TestIncludeRowsDoesNotChangeSelection:

    @pytest.mark.parametrize("derep", ["off", "family"])
    def test_accessions_are_identical_with_and_without_rows(self, three_families, derep):
        table, _screen = three_families
        with_rows = select_ref_genomes(table, "ncbi", "Bphy", derep_rank=derep,
                                       include_rows=True)
        without = select_ref_genomes(table, "ncbi", "Bphy", derep_rank=derep,
                                     include_rows=False)
        assert with_rows.accessions == without.accessions

    def test_rows_are_carried_only_when_asked_for(self, three_families):
        table, _screen = three_families
        assert select_ref_genomes(table, "ncbi", "Bphy", derep_rank="family",
                                  include_rows=True).rows
        assert select_ref_genomes(table, "ncbi", "Bphy", derep_rank="family",
                                  include_rows=False).rows == []

    def test_the_default_keeps_rows_because_hmm_autopick_reads_them(self, three_families):
        """
        scg_hmms_setup reads selection.rows for the pulled genomes' lineage to choose
        an SCG set (a euk target routing to Universal depends on it). If the default
        ever flipped to False, that would silently stop working.
        """
        table, _screen = three_families
        sel = select_ref_genomes(table, "ncbi", "Bphy", derep_rank="family")
        assert sel.rows, "default must carry rows"
        assert all("domain" in row for row in sel.rows)

    def test_warnings_and_resolution_survive_without_rows(self, three_families):
        table, screen = three_families
        sel = select_ref_genomes(table, "ncbi", "Bphy", derep_rank="family",
                                 screen_against=screen, include_rows=False)
        assert sel.canonical == "Bphy"
        assert sel.resolved_rank == "phylum"
        assert sel.effective_derep_rank == "family"
        assert any("no longer available at NCBI" in w for w in sel.warnings)


class TestCheapCountPathWouldOverreport:
    """
    Why `--dry-run` resolves the selection instead of calling derep_size().

    derep_size() counts distinct assigned values of the derep-rank column in the
    table, never applying `screen_against` or the quality floor. A dry run built on
    it would promise more genomes than the real run downloads.
    """

    def test_derep_size_counts_every_group_in_the_table(self, three_families):
        table, _screen = three_families
        assert derep_size(table, "ncbi", "phylum", "Bphy", "family") == 3

    def test_liveness_screening_makes_the_real_selection_smaller(self, three_families):
        table, screen = three_families
        sel = select_ref_genomes(table, "ncbi", "Bphy", derep_rank="family",
                                 screen_against=screen, include_rows=False)
        assert len(sel.accessions) == 2
        assert derep_size(table, "ncbi", "phylum", "Bphy", "family") == 3

    def test_a_quality_floor_makes_it_smaller_still(self, three_families):
        table, screen = three_families
        sel = select_ref_genomes(table, "ncbi", "Bphy", derep_rank="family",
                                 screen_against=screen, min_completeness=90.0,
                                 include_rows=False)
        assert len(sel.accessions) == 1
        assert derep_size(table, "ncbi", "phylum", "Bphy", "family") == 3
