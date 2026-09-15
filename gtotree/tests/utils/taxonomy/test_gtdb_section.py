"""
Unit tests for `--gtdb-section`, the GTDB analogue of `--ncbi-section`.

The thing this flag exists to fix: the `-w` drivers used to pass `reps_only=None`,
which fell through to `SourceSpec.default_reps_only` -- True for GTDB. So every
`-w` selection on the main driver was drawn from the species-representative pool,
with no flag anywhere to widen it and nothing in the output saying so.

Two properties are easy to break and expensive to notice:

  * `--gtdb-section` must stay INERT under `--source ncbi`. `reps_only_for` returns
    None there, which defers to NCBI's own default (all genomes). Returning False
    would look equivalent but isn't -- it would pin NCBI's pool to a value this flag
    has no business setting.

  * "the user typed `--gtdb-section reps`" and "the user got the default" must stay
    distinguishable, which is why the argparse default is None rather than "reps".
    Collapse them and every empty GTDB pull starts accusing the user of a flag they
    never passed. That's what `resolved_gtdb_section` is for: it fills the default in
    at the point of USE, so the raw args value stays honest about what was typed.
"""

import pytest  # type: ignore

from gtotree.utils.taxonomy.wanted_ref_tax import (DEFAULT_GTDB_SECTION,
                                                   GTDB_SECTIONS,
                                                   describe_gtdb_section,
                                                   reps_only_for,
                                                   resolved_gtdb_section)
from gtotree.utils.taxonomy.tax_select import SOURCES


class TestResolvedSection:

    def test_an_unset_flag_resolves_to_the_default(self):
        assert resolved_gtdb_section(None) == DEFAULT_GTDB_SECTION

    def test_the_default_is_the_narrow_slice(self):
        # the whole point: a defaulted GToTree run pulls representatives
        assert GTDB_SECTIONS[DEFAULT_GTDB_SECTION] is True

    @pytest.mark.parametrize("value", ["reps", "all"])
    def test_a_set_flag_resolves_to_itself(self, value):
        assert resolved_gtdb_section(value) == value

    @pytest.mark.parametrize("value", ["REPS", "  All  "])
    def test_case_and_whitespace_are_tolerated(self, value):
        assert resolved_gtdb_section(value) == value.strip().lower()

    def test_an_unrecognized_value_falls_back_rather_than_raising(self):
        # argparse `choices` is the real gate; this is just belt-and-braces for the
        # getattr() paths, and it must not widen the pool by accident
        assert resolved_gtdb_section("nonsense") == DEFAULT_GTDB_SECTION


class TestRepsOnlyFor:

    def test_gtdb_reps_asks_for_the_representatives_pool(self):
        assert reps_only_for("gtdb", "reps") is True

    def test_gtdb_all_asks_for_every_genome(self):
        assert reps_only_for("gtdb", "all") is False

    def test_an_unset_flag_under_gtdb_still_means_reps(self):
        assert reps_only_for("gtdb", None) is True

    def test_source_is_case_insensitive(self):
        # the main driver lowercases via type=str.lower, but gen-scg-hmms and
        # search-annotations hand their own args through, and RunData round-trips
        # `--source` as an uppercased string
        assert reps_only_for("GTDB", "all") is False

    def test_it_is_inert_under_ncbi(self):
        """
        None, not False: None defers to the source's own default, which for NCBI is
        already all-genomes. Pinning it to False would give the same answer today and
        silently override NCBI if that default ever changed.
        """
        assert reps_only_for("ncbi", "all") is None
        assert reps_only_for("ncbi", "reps") is None
        assert reps_only_for("ncbi", None) is None

    def test_none_really_does_defer_to_the_ncbi_default(self):
        # guards the reasoning above against a SourceSpec change
        assert SOURCES["ncbi"].default_reps_only is False


class TestSectionDescription:

    def test_gtdb_names_the_pool_for_the_source_header(self):
        assert describe_gtdb_section("gtdb", None) == "species representatives"
        assert describe_gtdb_section("gtdb", "all") == "all genomes"

    def test_ncbi_gets_no_gtdb_pool_label(self):
        assert describe_gtdb_section("ncbi", "all") is None
