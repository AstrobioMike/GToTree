"""
Unit tests for gtotree/utils/taxonomy/empty_selection.py.

This module is what tells a user WHY a taxonomy pull came back with nothing, and which
flag to loosen to get genomes back. It is pure string assembly over a PoolAttrition
record, so the failure mode isn't a crash -- it's confidently naming the wrong cause,
which sends someone off loosening a filter that was never the problem.

The two things worth holding onto while reading:

  * `filtered` is not "were any filters set?", it's "did one of the USER's filters do
    this?". A source's default representatives-only pool, NCBI suppressions and an
    all-unclassified derep rank are all causes no flag produced, and saying "with the
    specified filters" for those is a lie.

  * the emoticon must be applied AFTER the filters note, because the note splices
    itself in just inside a trailing period. Swapping the period out first strands the
    note on the far end of the sentence.

Built on the real PoolAttrition/RefGenomeSelection rather than stubs, so a field rename in
tax_derep.py breaks these instead of silently passing against a hand-rolled shape.
"""

import pytest  # type: ignore

from gtotree.utils.taxonomy.tax_derep import PoolAttrition, RefGenomeSelection
from gtotree.utils.taxonomy.empty_selection import (_genomes,
                                                    _oxford,
                                                    _restricting_section,
                                                    _with_emoticon,
                                                    empty_pull_message,
                                                    empty_selection_message,
                                                    explain_empty_selection)


def _selection(canonical="Testophyla", resolved_rank="phylum", **attrition_fields):
    """
    An empty RefSelection whose attrition record says what the caller wants it to.

    Stage defaults are all 0, i.e. "emptied at stage 1", so each test only sets the
    counts that carry it past the stages it isn't about.
    """
    derep_rank = attrition_fields.pop("derep_rank", None)
    reps_only = attrition_fields.pop("reps_only", False)

    attrition = PoolAttrition(resolved_rank=resolved_rank, derep_rank=derep_rank,
                              reps_only=reps_only)
    for key, value in attrition_fields.items():
        assert hasattr(attrition, key), f"PoolAttrition has no field {key!r}"
        setattr(attrition, key, value)

    return RefGenomeSelection(accessions=[], rows=[], canonical=canonical,
                        resolved_rank=resolved_rank,
                        effective_derep_rank=derep_rank, warnings=[],
                        attrition=attrition)


# ---------------------------------------------------------------------------
# the small string helpers
# ---------------------------------------------------------------------------

class TestOxford:

    @pytest.mark.parametrize("items,expected", [
        ([], ""),
        (["a"], "a"),
        (["a", "b"], "a and b"),
        (["a", "b", "c"], "a, b, and c"),
    ])
    def test_joins_with_a_serial_comma(self, items, expected):
        assert _oxford(items) == expected

    def test_conjunction_is_configurable(self):
        assert _oxford(["a", "b"], "or") == "a or b"
        assert _oxford(["a", "b", "c"], "or") == "a, b, or c"

    def test_empty_items_are_dropped(self):
        # a clause builder returning "" must not leave a dangling ", and"
        assert _oxford(["a", "", "c"]) == "a and c"


class TestGenomes:

    def test_one_genome_is_singular(self):
        assert _genomes(1) == "1 genome"

    def test_more_than_one_is_plural(self):
        assert _genomes(2) == "2 genomes"

    def test_zero_is_plural(self):
        assert _genomes(0) == "0 genomes"

    def test_large_counts_are_comma_grouped(self):
        assert _genomes(1234567) == "1,234,567 genomes"


class TestWithEmoticon:

    def test_replaces_a_trailing_period(self):
        assert _with_emoticon("Nothing found.", ":(") == "Nothing found :("

    def test_appends_when_there_is_no_period(self):
        assert _with_emoticon("Nothing found", ":(") == "Nothing found :("

    def test_is_a_noop_without_an_emoticon(self):
        assert _with_emoticon("Nothing found.", None) == "Nothing found."


class TestRestrictingSection:
    """`--ncbi-section both` narrows nothing, so it can never be what emptied a pool."""

    @pytest.mark.parametrize("value", ["refseq", "genbank", "RefSeq", "  GENBANK "])
    def test_a_real_restriction_is_reported_lowercased(self, value):
        assert _restricting_section(value) == value.strip().lower()

    @pytest.mark.parametrize("value", ["both", "", None, "anything-else"])
    def test_a_non_restriction_is_none(self, value):
        assert _restricting_section(value) is None


# ---------------------------------------------------------------------------
# explain_empty_selection -- stage by stage
# ---------------------------------------------------------------------------

class TestNothingToSay:

    def test_no_attrition_record_explains_nothing(self):
        class Bare:
            attrition = None
            resolved_rank = "phylum"
        assert explain_empty_selection(Bare()) == ("", False)

    def test_an_undiagnosed_empty_pool_explains_nothing(self):
        # diagnose_empty was off, so there is no n_unfiltered to reason from
        selection = _selection(diagnosed=False, n_unfiltered=0)
        assert explain_empty_selection(selection) == ("", False)

    def test_a_genuinely_empty_taxon_blames_no_filter(self):
        # diagnosed, but there was nothing under the taxon even unfiltered
        selection = _selection(diagnosed=True, n_unfiltered=0)
        assert explain_empty_selection(selection) == ("", False)

    def test_a_selection_that_survived_every_stage_explains_nothing(self):
        # empty for some reason none of the stages account for
        selection = _selection(n_candidates=5, n_after_liveness=5,
                               n_after_exclusion=5, n_after_floor=5)
        assert explain_empty_selection(selection) == ("", False)


class TestPoolFilterStage:

    def test_assembly_level_is_named_with_what_the_taxon_actually_has(self):
        selection = _selection(diagnosed=True, n_unfiltered=12,
                               pool_culprits=["assembly_levels"],
                               present_levels=["Contig", "Scaffold"])

        detail, filtered = explain_empty_selection(
            selection, assembly_levels=["Complete Genome"])

        assert "12 genomes at 'phylum' rank" in detail
        assert "`--assembly-level complete`" in detail
        # the NCBI strings must come back as the words a user would type
        assert "it has: contig and scaffold" in detail
        assert filtered is True

    def test_several_requested_levels_are_joined_with_or(self):
        selection = _selection(diagnosed=True, n_unfiltered=3,
                               pool_culprits=["assembly_levels"],
                               present_levels=["Contig"])

        detail, _ = explain_empty_selection(
            selection, assembly_levels=["Complete Genome", "Chromosome"])

        assert "`--assembly-level complete or chromosome`" in detail

    def test_ncbi_section_is_named_when_it_actually_restricts(self):
        selection = _selection(diagnosed=True, n_unfiltered=4,
                               pool_culprits=["accession_prefixes"])

        detail, filtered = explain_empty_selection(selection, ncbi_section="refseq")

        assert "`--ncbi-section refseq`" in detail
        assert filtered is True

    def test_section_blame_stays_generic_when_the_section_restricts_nothing(self):
        # 'both' imposes no prefix restriction, so there is no value to quote back
        selection = _selection(diagnosed=True, n_unfiltered=4,
                               pool_culprits=["accession_prefixes"])

        detail, filtered = explain_empty_selection(selection, ncbi_section="both")

        assert "requested part of NCBI" in detail
        assert "`--ncbi-section" not in detail
        assert filtered is True

    def test_a_requested_representatives_only_filter_is_blamed_on_the_flag(self):
        selection = _selection(diagnosed=True, n_unfiltered=9,
                               pool_culprits=["reps_only"])

        detail, filtered = explain_empty_selection(selection,
                                                   reps_only_requested=True)

        assert "`--representatives-only`" in detail
        assert filtered is True

    def test_a_sources_default_representatives_pool_is_not_the_users_filter(self):
        """
        GTDB pulls representatives by default. Blaming "the specified filters" when the
        user specified nothing would send them looking for a flag they never passed.
        """
        selection = _selection(diagnosed=True, n_unfiltered=9,
                               pool_culprits=["reps_only"])

        detail, filtered = explain_empty_selection(selection,
                                                   reps_only_requested=False)

        assert "default" in detail
        assert "`--representatives-only`" not in detail
        assert filtered is False

    def test_a_flag_alongside_the_default_reps_pool_still_implicates_the_flag(self):
        # mixed culprits: one is the user's, one is the source's -> the user's wins
        selection = _selection(diagnosed=True, n_unfiltered=9,
                               pool_culprits=["reps_only", "accession_prefixes"])

        detail, filtered = explain_empty_selection(selection, ncbi_section="genbank",
                                                   reps_only_requested=False)

        assert filtered is True
        assert "`--ncbi-section genbank`" in detail

    def test_filters_that_only_empty_it_together_are_named_together(self):
        # no single filter is the culprit, so none of them gets accused alone
        selection = _selection(diagnosed=True, n_unfiltered=20, pool_culprits=[])

        detail, filtered = explain_empty_selection(
            selection, assembly_levels=["Contig"], ncbi_section="refseq",
            reps_only_requested=True)

        assert "none survive" in detail
        assert "`--assembly-level contig`" in detail
        assert "`--ncbi-section refseq`" in detail
        assert "`--representatives-only`" in detail
        assert "Loosening any one of them" in detail
        assert filtered is True

    def test_no_culprit_and_no_active_flags_explains_nothing(self):
        selection = _selection(diagnosed=True, n_unfiltered=20, pool_culprits=[])
        assert explain_empty_selection(selection) == ("", False)


class TestLivenessStage:

    def test_all_suppressed_at_ncbi_is_not_the_users_filter(self):
        selection = _selection(n_candidates=7, n_after_liveness=0)

        detail, filtered = explain_empty_selection(selection)

        assert "All 7 genomes under it are suppressed or removed at NCBI" in detail
        assert filtered is False


class TestExclusionListStage:

    def test_everything_named_in_the_exclusion_list_is_a_user_filter(self):
        selection = _selection(n_candidates=7, n_after_liveness=5,
                               n_after_exclusion=0)

        detail, filtered = explain_empty_selection(selection)

        assert "Every one of its 5 genomes" in detail
        assert "`--exclusion-list`" in detail
        assert filtered is True


class TestQualityFloorStage:

    def _floored(self, **extra):
        return _selection(n_candidates=7, n_after_liveness=7, n_after_exclusion=6,
                          n_after_floor=0, **extra)

    def test_both_thresholds_are_quoted_back(self):
        detail, filtered = explain_empty_selection(
            self._floored(n_below_floor=6),
            min_completeness=90, max_contamination=5)

        assert "None of its 6 genomes cleared the quality floor" in detail
        assert "completeness >= 90" in detail
        assert "contamination <= 5" in detail
        assert "`--min-completeness`" in detail
        assert filtered is True

    def test_thresholds_render_without_trailing_zeros(self):
        # %g, so 90.0 must not come out as "90.000000"
        detail, _ = explain_empty_selection(self._floored(n_below_floor=6),
                                            min_completeness=90.0,
                                            max_contamination=2.5)
        assert "completeness >= 90" in detail
        assert "contamination <= 2.5" in detail

    def test_missing_checkm_values_are_reported_separately(self):
        detail, _ = explain_empty_selection(self._floored(n_missing_quality=6))

        assert "6 have no checkm values recorded" in detail
        assert "fell outside" not in detail

    def test_both_causes_are_reported_together(self):
        detail, _ = explain_empty_selection(
            self._floored(n_below_floor=4, n_missing_quality=2),
            min_completeness=90)

        assert "4 fell outside" in detail
        assert "2 have no checkm values recorded" in detail

    def test_a_floor_with_no_thresholds_passed_still_reads_sensibly(self):
        detail, _ = explain_empty_selection(self._floored(n_below_floor=6))
        assert "fell outside the floor" in detail


class TestDereplicationStage:

    def test_all_unclassified_at_the_derep_rank_is_not_a_filter(self):
        selection = _selection(derep_rank="genus", n_candidates=7,
                               n_after_liveness=7, n_after_exclusion=7,
                               n_after_floor=3, n_unassigned_group=3)

        detail, filtered = explain_empty_selection(selection)

        assert "unclassified at 'genus'" in detail
        assert "`--derep-rank off`" in detail
        assert filtered is False

    def test_only_some_unclassified_is_not_this_cause(self):
        selection = _selection(derep_rank="genus", n_candidates=7,
                               n_after_liveness=7, n_after_exclusion=7,
                               n_after_floor=3, n_unassigned_group=1)
        assert explain_empty_selection(selection) == ("", False)

    def test_unclassified_survivors_without_derep_on_is_not_this_cause(self):
        selection = _selection(derep_rank=None, n_candidates=7, n_after_liveness=7,
                               n_after_exclusion=7, n_after_floor=3,
                               n_unassigned_group=3)
        assert explain_empty_selection(selection) == ("", False)


# ---------------------------------------------------------------------------
# the two message wrappers
# ---------------------------------------------------------------------------

class TestEmptySelectionMessage:

    def test_headline_names_the_taxon_and_the_surfaces_own_flag(self):
        selection = _selection(canonical="Trichodesmium", diagnosed=False)

        assert empty_selection_message(selection) == (
            "No accessions were found for the --wanted-ref-tax target "
            "'Trichodesmium'.")

    def test_the_taxon_flag_is_configurable_for_bit(self):
        selection = _selection(canonical="Trichodesmium", diagnosed=False)
        assert "-t target 'Trichodesmium'" in empty_selection_message(
            selection, taxon_flag="-t")

    def test_the_filters_note_lands_before_the_emoticon(self):
        """
        The note splices in just inside a trailing period, so it has to go on before
        the emoticon eats that period -- otherwise it ends up stranded after it.
        """
        selection = _selection(n_candidates=7, n_after_liveness=5,
                               n_after_exclusion=0)

        message = empty_selection_message(selection, emoticon=":(")

        # note, then emoticon, then the detail sentence -- never ".) :(" ordering
        assert "with the specified filters :( Every one of its" in message
        assert not message.startswith("No accessions were found for the "
                                      "--wanted-ref-tax target 'Testophyla' :(")

    def test_an_unimplicated_cause_gets_no_filters_note(self):
        selection = _selection(n_candidates=7, n_after_liveness=0)

        message = empty_selection_message(selection)

        assert "with the specified filters" not in message
        assert "suppressed or removed at NCBI" in message

    def test_detail_is_appended_to_the_headline(self):
        selection = _selection(n_candidates=7, n_after_liveness=0)
        message = empty_selection_message(selection)
        assert message.startswith("No accessions were found for the")
        assert message.rstrip().endswith("so none can be downloaded.")


class TestEmptyPullMessage:
    """
    The get-accs surfaces word their own headline, so this only adds the explanation.
    """

    HEADLINE = "No genomes were found under phylum 'Testophyla'."

    def test_no_selection_keeps_the_honest_hedge(self):
        # an 'all' or taxid pull: nothing in scope knows what was filtered
        assert empty_pull_message(self.HEADLINE) == (
            "No genomes were found under phylum 'Testophyla' "
            "(after any specified filters).")

    def test_no_selection_still_takes_an_emoticon(self):
        message = empty_pull_message(self.HEADLINE, emoticon=":(")
        assert message.endswith("(after any specified filters) :(")

    def test_a_known_filter_cause_upgrades_the_hedge_to_an_accusation(self):
        selection = _selection(n_candidates=7, n_after_liveness=5,
                               n_after_exclusion=0)

        message = empty_pull_message(self.HEADLINE, selection)

        assert "with the specified filters" in message
        assert "(after any specified filters)" not in message
        assert "`--exclusion-list`" in message

    def test_a_non_filter_cause_gets_the_detail_but_no_filter_language(self):
        selection = _selection(n_candidates=7, n_after_liveness=0)

        message = empty_pull_message(self.HEADLINE, selection)

        assert "suppressed or removed at NCBI" in message
        assert "specified filters" not in message

    def test_an_unexplainable_selection_falls_back_to_the_hedge(self):
        selection = _selection(diagnosed=False)

        message = empty_pull_message(self.HEADLINE, selection)

        assert message.endswith("(after any specified filters).")

    def test_the_headline_is_never_restated(self):
        # the surfaces already printed the rank/canonical/derep label themselves
        selection = _selection(n_candidates=7, n_after_liveness=0)
        message = empty_pull_message(self.HEADLINE, selection)
        assert message.count("No genomes were found") == 1
