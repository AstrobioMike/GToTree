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
