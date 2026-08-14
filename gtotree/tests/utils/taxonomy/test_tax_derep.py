"""
Unit tests for gtotree/utils/taxonomy/tax_derep.py.

These are the pure decision functions: how the derep rank is chosen, and which single
genome is picked to represent a group. Both determine what ends up in a reference tree,
so they're specified here rather than only exercised through a live Parquet query.
"""

import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_select import SOURCES
from gtotree.utils.taxonomy.tax_derep import (DEREP_STEPS,
                                              _num,
                                              default_derep_rank,
                                              first_key,
                                              mean_contig_length,
                                              quality_key,
                                              resolve_derep_rank,
                                              resolve_picker,
                                              size_advice)

NCBI = SOURCES["ncbi"]
GTDB = SOURCES["gtdb"]


class TestDefaultDerepRank:
    """Default is DEREP_STEPS ranks finer than the target, clamped at species."""

    @pytest.mark.parametrize("target,expected", [
        ("domain", "class"),
        ("phylum", "order"),
        ("class", "family"),
        ("order", "genus"),
        ("family", "species"),
        ("genus", "species"),     # clamped
    ])
    def test_is_two_ranks_finer_clamped_at_species(self, target, expected):
        assert default_derep_rank(target) == expected
        assert DEREP_STEPS == 2

    def test_species_target_has_no_default(self):
        """Nothing is finer than species, so derep would return a single genome."""
        assert default_derep_rank("species") is None


class TestResolveDerepRank:

    def test_auto_uses_the_default_for_the_target_rank(self):
        assert resolve_derep_rank("phylum", "auto") == ("order", [])

    @pytest.mark.parametrize("value", [None, "off", "none", "None"])
    def test_off_disables_dereplication(self, value):
        assert resolve_derep_rank("phylum", value) == (None, [])

    def test_an_explicit_finer_rank_is_honored(self):
        rank, warnings = resolve_derep_rank("phylum", "species")
        assert rank == "species"
        assert warnings == []

    def test_a_rank_coarser_than_the_target_is_rejected(self):
        """Dereplicating above the target would collapse the whole selection."""
        with pytest.raises(ValueError):
            resolve_derep_rank("family", "phylum")

    def test_auto_on_a_species_target_warns_and_turns_derep_off(self):
        rank, warnings = resolve_derep_rank("species", "auto")
        assert rank is None
        assert len(warnings) == 1
        assert "--derep-rank species" in warnings[0]


class TestSizeAdvice:

    def test_no_advice_for_a_sensible_selection(self):
        assert size_advice(500, "phylum", "order") == []

    def test_suggests_a_coarser_rank_for_a_huge_selection(self):
        advice = size_advice(200_000, "domain", "species")
        assert len(advice) == 1
        assert "coarser" in advice[0]

    def test_suggests_a_finer_rank_for_a_tiny_selection(self):
        advice = size_advice(1, "domain", "phylum")
        assert len(advice) == 1
        assert "finer" in advice[0]

    def test_huge_selection_with_derep_off_suggests_turning_it_on(self):
        advice = size_advice(200_000, "phylum", None)
        assert len(advice) == 1
        assert "--derep-rank" in advice[0] or "--derep-rank off" in advice[0]

    def test_species_target_with_derep_off_explains_there_is_nothing_finer(self):
        advice = size_advice(200_000, "species", None)
        assert len(advice) == 1
        assert "no rank finer" in advice[0]


class TestNum:

    @pytest.mark.parametrize("value", [None, "", "na", "NA"])
    def test_missing_values_become_the_default(self, value):
        assert _num(value, default=7) == 7

    def test_parses_numeric_strings(self):
        assert _num("98.5") == 98.5

    def test_unparseable_values_become_the_default(self):
        assert _num("not a number", default=0.0) == 0.0


class TestMeanContigLength:
    """Assembly-quality proxy for genomes with no checkm values."""

    def test_divides_genome_size_by_contig_count(self):
        row = {"genome_size_ungapped": "4000000", "contig_count": "100"}
        assert mean_contig_length(row, NCBI) == 40000.0

    def test_falls_back_to_the_secondary_size_column(self):
        row = {"genome_size": "4000000", "contig_count": "100"}
        assert mean_contig_length(row, NCBI) == 40000.0

    @pytest.mark.parametrize("row", [
        {"contig_count": "0", "genome_size_ungapped": "4000000"},
        {"contig_count": "", "genome_size_ungapped": "4000000"},
        {"contig_count": "100"},
    ])
    def test_returns_zero_when_it_cannot_be_computed(self, row):
        assert mean_contig_length(row, NCBI) == 0.0

    def test_returns_zero_for_a_source_without_the_columns(self):
        spec = SOURCES["gtdb"]
        spec_no_cols = type(spec)(name="x", acc_col="a", rep_filter=("r", "t"),
                                  quality_cols=("c", "d"), ref_col="r",
                                  default_reps_only=False)
        assert mean_contig_length({"contig_count": "10"}, spec_no_cols) == 0.0


class TestQualityKey:
    """
    Lowest key wins. Ordering: reference genomes, then checkm'd genomes, then highest
    completeness / lowest contamination, then assembly level and mean contig length for
    genomes without checkm, then accession as a stable tiebreak.
    """

    def _row(self, acc="GCF_1", ref="", comp=None, cont=None, level=None,
             contigs=None, size=None):
        return {"assembly_accession": acc, "refseq_category": ref,
                "checkm_completeness": comp, "checkm_contamination": cont,
                "assembly_level": level, "contig_count": contigs,
                "genome_size_ungapped": size}

    def _best(self, rows):
        return min(rows, key=lambda r: quality_key(r, NCBI))["assembly_accession"]

    def test_reference_genomes_win(self):
        rows = [self._row("GCF_plain", comp="99", cont="0"),
                self._row("GCF_ref", ref="reference genome", comp="90", cont="2")]
        assert self._best(rows) == "GCF_ref"

    def test_a_poor_quality_reference_loses_its_privilege(self):
        """
        The reference flag is gated on not being junk: below 85% complete or above 10%
        contaminated and it competes on quality like anything else.
        """
        rows = [self._row("GCF_good", comp="99", cont="1"),
                self._row("GCF_badref", ref="reference genome", comp="50", cont="1")]
        assert self._best(rows) == "GCF_good"

    def test_genomes_with_checkm_beat_genomes_without(self):
        rows = [self._row("GCF_nocheckm", level="complete genome"),
                self._row("GCF_checkm", comp="80", cont="5")]
        assert self._best(rows) == "GCF_checkm"

    def test_higher_completeness_wins(self):
        rows = [self._row("GCF_low", comp="90", cont="1"),
                self._row("GCF_high", comp="99", cont="1")]
        assert self._best(rows) == "GCF_high"

    def test_lower_contamination_breaks_a_completeness_tie(self):
        rows = [self._row("GCF_dirty", comp="99", cont="8"),
                self._row("GCF_clean", comp="99", cont="0.5")]
        assert self._best(rows) == "GCF_clean"

    def test_unknown_contamination_sorts_worst_among_checkm_genomes(self):
        rows = [self._row("GCF_known", comp="99", cont="9"),
                self._row("GCF_unknown", comp="99", cont=None)]
        assert self._best(rows) == "GCF_known"

    def test_without_checkm_better_assembly_level_wins(self):
        rows = [self._row("GCF_contig", level="contig"),
                self._row("GCF_complete", level="complete genome")]
        assert self._best(rows) == "GCF_complete"

    def test_without_checkm_longer_mean_contig_breaks_a_level_tie(self):
        rows = [self._row("GCF_frag", level="contig", contigs="1000", size="4000000"),
                self._row("GCF_solid", level="contig", contigs="10", size="4000000")]
        assert self._best(rows) == "GCF_solid"

    def test_accession_is_the_final_tiebreak_so_results_are_reproducible(self):
        rows = [self._row("GCF_2", comp="99", cont="1"),
                self._row("GCF_1", comp="99", cont="1")]
        assert self._best(rows) == "GCF_1"
        # and the ordering is total: same input, same answer, regardless of input order
        assert self._best(list(reversed(rows))) == "GCF_1"


class TestResolvePicker:

    def test_named_pickers_resolve(self):
        assert resolve_picker("quality") is quality_key
        assert resolve_picker("first") is first_key

    def test_a_callable_is_passed_through(self):
        f = lambda row, spec: (0,)
        assert resolve_picker(f) is f

    def test_an_unknown_name_is_rejected(self):
        with pytest.raises(ValueError, match="unknown pick policy"):
            resolve_picker("best-vibes")


def test_first_key_is_the_accession():
    assert first_key({"assembly_accession": "GCF_9"}, NCBI) == ("GCF_9",)
