"""
Guards for the `gtt` subcommand dispatcher

`build_parser` imports every module named in SUBCOMMAND_MAP, so a single stale entry
takes down the whole CLI, not just the affected subcommand. That's a cheap failure to
introduce during a refactor (moving a module and missing the map) and, without these
tests, an expensive one to notice: nothing else in the suite imports through the map.
"""

import argparse
import importlib
import pytest  # type: ignore
from gtotree.cli.gtt import SUBCOMMAND_MAP, PROGRAM_GROUPS, build_parser


@pytest.mark.parametrize("name,module_path", sorted(SUBCOMMAND_MAP.items()))
def test_every_subcommand_module_is_importable(name, module_path):
    """
    Each mapped module must actually exist. Parametrized so a broken entry names the
    offending subcommand instead of failing the whole map at once.
    """
    importlib.import_module(module_path)


@pytest.mark.parametrize("name,module_path", sorted(SUBCOMMAND_MAP.items()))
def test_every_subcommand_module_exposes_the_dispatch_api(name, module_path):
    """
    Dispatch calls `module.build_parser(parent_subparsers=...)` and `module.main()`, so
    a module that imports cleanly but lacks either would still break at runtime.
    """
    module = importlib.import_module(module_path)
    assert callable(getattr(module, "build_parser", None)), \
        f"{module_path} is missing build_parser()"
    assert callable(getattr(module, "main", None)), \
        f"{module_path} is missing main()"


def test_build_parser_succeeds():
    """
    The end-to-end guard: builds the full parser tree, which imports every mapped
    module. This is the test that fails if any subcommand path goes stale.
    """
    parser = build_parser()
    assert parser is not None


def test_build_parser_registers_every_subcommand():
    """Every mapped name should end up as an actual subparser choice."""
    import argparse

    parser = build_parser()
    choices = set()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            choices.update(action.choices)

    missing = set(SUBCOMMAND_MAP) - choices
    assert not missing, f"subcommands absent from the parser: {sorted(missing)}"


def test_overview_listing_matches_the_map():
    """
    The help overview is maintained by hand alongside SUBCOMMAND_MAP; drift between the
    two means a subcommand that works but is undiscoverable (or is advertised and
    doesn't dispatch).
    """
    listed = {
        prog["name"]
        for group in PROGRAM_GROUPS
        for prog in group["programs"]
    }
    mapped = set(SUBCOMMAND_MAP)

    assert not (listed - mapped), \
        f"advertised but not dispatchable: {sorted(listed - mapped)}"
    # the direction that was previously unchecked, and the likelier drift: adding to
    # SUBCOMMAND_MAP and forgetting PROGRAM_GROUPS gives a subcommand that works but is
    # invisible in `gtt` help
    assert not (mapped - listed), \
        f"dispatchable but not advertised in the overview: {sorted(mapped - listed)}"


def _subcommand_layers_smuggled_in_as_positionals(parser, path):
    """
    Walk the tree and collect any parser that offers its subcommands as a positional
    with `choices` instead of as a real subparser level.
    """
    subparser_actions = [
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]

    if subparser_actions:
        offenders = []
        for name, child in subparser_actions[0].choices.items():
            offenders += _subcommand_layers_smuggled_in_as_positionals(
                child, f"{path} {name}")
        return offenders

    return [
        f"{path} (<{action.dest}>)"
        for action in parser._actions
        if not action.option_strings and getattr(action, "choices", None)
    ]


def test_no_parser_uses_a_positional_with_choices_as_a_subcommand_layer():
    """
    Tab completion is driven by argcomplete walking this tree, and
    _suppress_help_version_on_group_parsers() recognizes a "group" only via
    argparse._SubParsersAction. A parser that instead offers its subcommands as a
    positional with `choices` parses identically and reads fine under `-h`, but TAB
    after it pads the choices with -h/-v and whatever else that parser defines.

    This is exactly how `gtt data get` drifted away from `bit data get`: it took a
    `source` positional with `choices`, so TAB after it offered -f/-h/-v mixed in
    among the four database names.

    Tree-wide on purpose rather than scoped to the parsers that had the problem: the
    shape is an easy one to reach for, and the symptom only ever shows up at a
    terminal, so it goes unnoticed until someone happens to hit TAB there.
    """
    offenders = _subcommand_layers_smuggled_in_as_positionals(build_parser(), "gtt")
    assert not offenders, (
        "these parsers use a positional-with-choices as a subcommand layer, so "
        f"argcomplete will pad their choices with flags: {sorted(offenders)}"
    )
