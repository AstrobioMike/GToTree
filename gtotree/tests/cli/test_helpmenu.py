"""
The help menu is hand-written, so it can drift from the parser. These tests are
what make that safe: they fail if a flag exists in one place and not the other.

This is the failure mode the v1 menu actually had -- 20 long flags were defined in
the parser and documented nowhere, and `-r`'s stated default had fallen out of sync
with the real one.
"""

import argparse
import re

import pytest

from gtotree.cli import helpmenu as hm
from gtotree.cli.parser import parser


ANSI = re.compile(r"\033\[[0-9;]*m")
FLAG = re.compile(r"-{1,2}[A-Za-z][\w-]*")

# Flags that are intentionally only in the parser, or only in the menu.
# Keep this list short and justified -- it's the escape hatch that lets drift back in.
NOT_IN_MENU = set()


def parser_flags():
    """Every option string the parser accepts."""
    flags = set()
    for action in parser()._actions:
        flags.update(action.option_strings)
    return flags


def menu_flags(detailed=True):
    """Every flag named in an entry's flag field (not prose references)."""
    flags = set()
    for _, items in hm.ENTRIES:
        for entry_flags, _, _, tier in items:
            if not detailed and tier is not hm.BRIEF:
                continue
            flags.update(FLAG.findall(entry_flags))
    for entry_flags, _ in hm.REQUIRED_INPUTS:
        flags.update(FLAG.findall(hm.pick(entry_flags, detailed)))
    flags.update(FLAG.findall(hm.pick(hm.HMM_INPUT[0], detailed)))
    flags.update({"-h", "--help", "-S", "--show-detailed-help", "-v", "--version"})
    return flags


def test_every_documented_flag_is_defined():
    extra = menu_flags() - parser_flags()
    assert not extra, (
        f"documented in the help menu but not defined in cli/parser.py: {sorted(extra)}"
    )


def test_condensed_flags_are_a_subset_of_detailed():
    assert menu_flags(detailed=False) <= menu_flags(detailed=True)


@pytest.mark.parametrize("detailed", [False, True])
def test_menu_renders_and_respects_width(detailed):
    text = ANSI.sub("", hm.build_helpmenu(detailed=detailed))
    assert text.strip()
    too_wide = [ln for ln in text.splitlines() if len(ln) > hm.WIDTH]
    assert not too_wide, f"lines exceeding {hm.WIDTH} columns: {too_wide[:3]}"


def test_detailed_menu_is_a_superset_of_condensed():
    condensed = ANSI.sub("", hm.build_helpmenu(detailed=False))
    detailed = ANSI.sub("", hm.build_helpmenu(detailed=True))
    assert len(detailed) > len(condensed)
    # the "run -S for more" pointer belongs only on the condensed menu
    assert "condensed help menu" in condensed
    assert "This is the condensed help menu" not in detailed


def test_stated_defaults_match_the_parser():
    """
    Catches the class of bug the v1 menu had: a summary line saying "default: 0.2"
    for a parameter the parser defaults to 0.1.
    """
    defaults = {}
    for action in parser()._actions:
        if action.default in (None, argparse.SUPPRESS, False):
            continue
        for opt in action.option_strings:
            defaults[opt] = action.default

    mismatches = []
    for _, items in hm.ENTRIES:
        for entry_flags, brief, _, _ in items:
            stated = re.search(r"default:\s*([^\s;)]+)", hm.pick(brief, True))
            if not stated:
                continue
            for opt in FLAG.findall(entry_flags):
                if opt in defaults:
                    if str(defaults[opt]) != stated.group(1):
                        mismatches.append((opt, stated.group(1), defaults[opt]))
                    break

    assert not mismatches, (
        "help menu states a default the parser doesn't use "
        f"(flag, menu says, parser says): {mismatches}"
    )
