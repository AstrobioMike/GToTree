"""
Each surface must name ITS OWN representatives-only flag in empty-result messages.

`--representatives-only` is `gtt dl-ncbi-assemblies`' spelling and nothing else's.
`gtt get-accs-from-gtdb` narrows with `-G/--gtdb-representatives-only` or
`-R/--refseq-ref-genomes-only`, `gtt get-accs-from-ncbi` with `-R`, and the `-w`
drivers with `--gtdb-section reps`. All four used to emit `--representatives-only`,
so a user who hit an empty pull was told to loosen a flag their subcommand does not
have -- and on get-accs-from-gtdb the `-G` and `-R` cases produced byte-identical
text despite being different filters, which is why that surface resolves the name at
runtime instead of holding a constant.
"""

import pytest  # type: ignore

from gtotree.utils.gtdb.get_accessions_from_gtdb import _reps_flag_for
from gtotree.utils.taxonomy.empty_selection import (DEFAULT_REPS_FLAG,
                                                    empty_pull_message)
from gtotree.utils.taxonomy.tax_derep import PoolAttrition, RefGenomeSelection


def _emptied_by_reps_only():
    """A selection whose pool was emptied solely by the representatives filter."""
    attrition = PoolAttrition(resolved_rank="class", derep_rank=None, reps_only=True)
    attrition.diagnosed = True
    attrition.n_unfiltered = 412
    attrition.n_candidates = 0
    attrition.pool_culprits = ["reps_only"]
    attrition.present_levels = []
    return RefGenomeSelection(accessions=[], rows=[], canonical="Nitrososphaeria",
                              resolved_rank="class", effective_derep_rank=None,
                              warnings=[], attrition=attrition)


class TestGtdbHelperFlagNames:

    def test_gtdb_representatives_are_named_as_dash_g(self):
        assert _reps_flag_for("gtdb") == "-G/--gtdb-representatives-only"

    def test_refseq_references_are_named_as_dash_r(self):
        assert _reps_flag_for("refseq") == "-R/--refseq-ref-genomes-only"

    def test_the_two_sources_are_not_named_identically(self):
        # the actual pre-existing bug: both rendered as `--representatives-only`
        assert _reps_flag_for("gtdb") != _reps_flag_for("refseq")

    @pytest.mark.parametrize("source", ["gtdb", "refseq"])
    def test_neither_names_dl_ncbi_assemblies_flag(self, source):
        assert _reps_flag_for(source) != DEFAULT_REPS_FLAG


class TestMessageRendering:

    @pytest.mark.parametrize("source,expected", [
        ("gtdb", "-G/--gtdb-representatives-only"),
        ("refseq", "-R/--refseq-ref-genomes-only"),
    ])
    def test_the_message_carries_the_surfaces_own_flag(self, source, expected):
        message = empty_pull_message("No genomes were found.", _emptied_by_reps_only(),
                                     reps_only_requested=True,
                                     reps_flag=_reps_flag_for(source))

        assert f"`{expected}`" in message
        assert DEFAULT_REPS_FLAG not in message

    def test_an_unpassed_filter_is_still_never_mentioned(self):
        """
        These helpers default to all genomes, so with no `-G`/`-R` the representatives
        filter never enters the culprit list and must not be named at all.
        """
        selection = _emptied_by_reps_only()
        selection.attrition.reps_only = False
        selection.attrition.pool_culprits = []

        message = empty_pull_message("No genomes were found.", selection,
                                     reps_only_requested=False)

        assert "representative" not in message.lower()
