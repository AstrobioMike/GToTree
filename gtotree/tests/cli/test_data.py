"""
Shape tests for `gtt data`.

Mirrored by bit/tests/test_data_cli.py in the bit repo -- if something changes here,
check there too.

The theme is the failure these were written after: `gtt data get` used to take a
`source` positional with `choices` rather than a real subparser level. It parsed the
same and its help menu read fine, but gtt.py's
_suppress_help_version_on_group_parsers() only recognizes a group via
argparse._SubParsersAction, so TAB after `gtt data get` offered -f/-h/-v mixed in
among the four database names while `bit data get` offered just the names.
"""

import argparse
import io
import contextlib

import pytest  # type: ignore

from gtotree.cli import data as data_cli
from gtotree.cli.gtt import _suppress_help_version_on_group_parsers


HELP_AND_VERSION = {"-h", "--help", "-v", "--version"}


def _subparsers_action(parser):
    """The parser's subparsers action, or None if it isn't a group node."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action
    return None


def _node(*path):
    """Walk down from the `data` parser by subcommand name."""
    node = data_cli.build_parser()
    for name in path:
        action = _subparsers_action(node)
        assert action is not None, f"'{name}' has no subcommand layer to descend into"
        node = action.choices[name]
    return node


def _flags(parser):
    flags = set()
    for action in parser._actions:
        flags.update(action.option_strings)
    return flags


def _help_text(parser):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        parser.print_help()
    return buf.getvalue()


def test_data_get_is_a_real_subparser_level():
    """
    The load-bearing assertion. A `source` positional with `choices` would satisfy any
    test about *parsing*, so this checks the thing tab completion actually keys off.
    """
    action = _subparsers_action(_node("get"))
    assert action is not None
    assert set(action.choices) == set(data_cli.DATA_SOURCES)


def test_group_nodes_hide_help_and_version_from_completion():
    """
    What the shape buys us: the suppressor marks -h/-v SUPPRESS on group parsers, and
    argcomplete skips SUPPRESSed options, so TAB offers only the database names.
    """
    parser = data_cli.build_parser()
    _suppress_help_version_on_group_parsers(parser)

    get_parser = _subparsers_action(parser).choices["get"]
    suppressed = {
        option
        for action in get_parser._actions
        for option in action.option_strings
        if action.help is argparse.SUPPRESS
    }
    assert HELP_AND_VERSION <= suppressed


def test_leaf_nodes_keep_their_flags_visible_to_completion():
    """The flip side -- suppression must stop at the group level, not cascade."""
    parser = data_cli.build_parser()
    _suppress_help_version_on_group_parsers(parser)

    leaf = _subparsers_action(_subparsers_action(parser).choices["get"]).choices["gtdb-data"]
    assert all(action.help is not argparse.SUPPRESS for action in leaf._actions)


@pytest.mark.parametrize("source_name", sorted(data_cli.DATA_SOURCES))
def test_every_source_takes_the_same_flags(source_name):
    """
    No `-q/--quiet` here, unlike bit's mirror of this menu -- GToTree's getter
    functions don't take a quiet argument, so advertising the flag would be a lie.
    """
    assert _flags(_node("get", source_name)) == HELP_AND_VERSION | {"-f", "--force-update"}


@pytest.mark.parametrize("source_name", sorted(data_cli.DATA_SOURCES))
def test_every_source_renders_one_optional_parameters_section(source_name):
    """
    -h and -v have to land in the same explicit argument group. Adding one to the
    group and the other to the parser also "works", but argparse renders the parser's
    default group first and the help menu comes out with two identically-titled
    sections.
    """
    assert _help_text(_node("get", source_name)).count("Optional Parameters:") == 1


@pytest.mark.parametrize("source_name", sorted(data_cli.DATA_SOURCES))
def test_every_source_names_an_importable_worker(source_name):
    """
    DATA_SOURCES holds module paths and function names as strings so that building
    the parser (which gtt.py does on every TAB keypress) imports nothing. That means
    a typo in either string is invisible until someone runs the command.
    """
    import importlib

    source_info = data_cli.DATA_SOURCES[source_name]
    module = importlib.import_module(source_info["module"])
    assert callable(getattr(module, source_info["function"]))


@pytest.mark.parametrize("force", [False, True])
def test_get_dispatches_to_the_selected_source(monkeypatch, force):
    called = {}

    def fake_import(module_path):
        assert module_path == data_cli.DATA_SOURCES["pfam-data"]["module"]

        class Stub:
            pass

        stub = Stub()
        setattr(stub, "get_pfam_data",
                lambda **kwargs: called.update(kwargs))
        return stub

    monkeypatch.setattr(data_cli.importlib, "import_module", fake_import)

    argv = ["get", "pfam-data"] + (["-f"] if force else [])
    args = data_cli.build_parser().parse_args(argv)
    assert args.get_action == "pfam-data"

    data_cli._run_get(args)
    assert called == {"force_update": force}
