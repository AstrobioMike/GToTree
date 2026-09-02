"""
`--derep-rank` now defaults to 'auto' on both get-accs helpers, matching the main
GToTree program and dl-ncbi-assemblies.

argparse stores a sentinel rather than 'auto' directly, so preflight can tell an
inherited default from a value the user actually typed. The taxid path is the reason
that distinction has to exist, and these cover both sides of it.
"""

import argparse

import pytest

from gtotree.utils.gtdb.get_accessions_from_gtdb import (
    build_parser as gtdb_parser, preflight_checks as gtdb_preflight)
from gtotree.utils.ncbi.get_accessions_from_ncbi import (
    build_parser as ncbi_parser, preflight_checks as ncbi_preflight)
from gtotree.utils.taxonomy.get_accs_shared import (DEREP_DEFAULT, DEREP_OFF,
                                                    DEREP_UNSET,
                                                    apply_derep_default, is_derep_on)


def _ncbi(argv):
    args = ncbi_parser().parse_args(argv)
    ncbi_preflight(args)
    return args


# -- the sentinel itself -------------------------------------------------------


def test_sentinel_is_not_read_as_off():
    """
    The reason the sentinel isn't None: None is in DEREP_OFF, so an unresolved None
    would quietly mean 'off' instead of failing.
    """
    assert DEREP_UNSET not in DEREP_OFF
    assert is_derep_on(DEREP_UNSET)


def test_apply_derep_default_fills_in_the_default():
    args = argparse.Namespace(derep_rank=DEREP_UNSET)
    assert apply_derep_default(args) is False
    assert args.derep_rank == DEREP_DEFAULT


def test_apply_derep_default_leaves_an_explicit_value_alone():
    args = argparse.Namespace(derep_rank="family")
    assert apply_derep_default(args) is True
    assert args.derep_rank == "family"


def test_explicitly_typed_auto_still_counts_as_explicit():
    """
    The whole point of the sentinel: `--derep-rank auto` and typing nothing both end
    up as 'auto', but only one of them is a request.
    """
    assert apply_derep_default(argparse.Namespace(derep_rank="auto")) is True


# -- what the helpers actually default to --------------------------------------


def test_ncbi_helper_defaults_to_auto():
    assert _ncbi(["-w", "Alteromonas"]).derep_rank == "auto"


def test_gtdb_helper_defaults_to_auto():
    args = gtdb_parser().parse_args(["-w", "Archaea"])
    gtdb_preflight(args)
    assert args.derep_rank == "auto"


def test_ncbi_section_defaults_to_both():
    """
    The refseq default only existed to keep an underep'd pull from returning
    everything. --derep-rank auto does that job now, so the section filter goes back
    to matching every other wanted-ref-tax surface.
    """
    assert _ncbi(["-w", "Alteromonas"]).ncbi_section == "both"


# -- the taxid path ------------------------------------------------------------


def test_taxid_quietly_drops_the_inherited_default():
    """A plain taxid pull has to keep working now that the default dereplicates."""
    assert _ncbi(["-w", "28108"]).derep_rank == "off"


def test_taxid_with_an_explicit_derep_rank_still_errors():
    with pytest.raises(SystemExit):
        _ncbi(["-w", "28108", "--derep-rank", "genus"])


def test_taxid_with_explicit_auto_errors_too():
    """
    Not distinguishable from the default by value, only by provenance -- so this is
    the case that would regress if the sentinel were dropped.
    """
    with pytest.raises(SystemExit):
        _ncbi(["-w", "28108", "--derep-rank", "auto"])


def test_taxid_with_explicit_off_is_fine():
    assert _ncbi(["-w", "28108", "--derep-rank", "off"]).derep_rank == "off"


def test_a_named_target_keeps_an_explicit_derep_rank():
    assert _ncbi(["-w", "Alteromonas", "--derep-rank", "family"]).derep_rank == "family"
