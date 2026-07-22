import sys
import shlex
import re
import shutil
import textwrap
import argparse
from rich_argparse import RichHelpFormatter  # type: ignore

# this is basically a port of what i came up with for bit

class CustomRichHelpFormatter(RichHelpFormatter):

    # this was added so help text appears on same line as subcommand for programs with subcommands
    def __init__(self, prog, indent_increment=2, max_help_position=24, width=None, **kwargs):
        super().__init__(prog, indent_increment=indent_increment,
                         max_help_position=max(max_help_position, 36),
                         width=width, **kwargs)

    def add_argument(self, action):
        super().add_argument(action)
        if hasattr(action, '_get_subactions'):
            self._action_max_length += 2

    def start_section(self, heading):
        if heading == "positional arguments":
            heading = "Available Subcommands"
        elif heading == "options":
            heading = "Optional Parameters"
        super().start_section(heading)
    group_name_formatter = lambda name: "Usage" if name.lower() == "usage" else name


def add_help(group):
    group.add_argument(
        "-h",
        "--help",
        action="help",
        help=wrap_help("Show this help message and exit")
    )


class VersionAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from gtotree.utils.messaging import get_version
        print(f"GToTree v{get_version()}")
        parser.exit()


def add_version_arg(group):
    group.add_argument(
        "-v",
        "--version",
        action=VersionAction,
        help=wrap_help("Show GToTree version and exit")
    )


def add_force(group):
    group.add_argument(
        "-F",
        "--force-overwrite",
        help=wrap_help("Force overwrite of output files if they already exist"),
        action="store_true",
    )


def reconstruct_invocation(parser, args):
    """
    This is a helper for capturing the executed command-line (including default args) for logging.
    """
    cmd = [sys.executable, sys.argv[0]]

    # handling if there is a subcommand
    sub = getattr(args, "subcommand", None)
    if sub:
        cmd.append(sub)

    # gathering all actions
    all_actions = list(parser._actions)
    # including subparsers if they exist
    for a in parser._actions:
        if isinstance(a, argparse._SubParsersAction):
            subpar_action = a
            break
    else:
        subpar_action = None

    if sub and subpar_action:
        all_actions += subpar_action.choices[sub]._actions

    # walking through and getting settings for each item
    for action in all_actions:
        # skip help, the subparsers action itself, and version (an exiting flag,
        # never a real run parameter -- avoids a noisy "--version None" in logs)
        if isinstance(action, (argparse._HelpAction,
                               argparse._SubParsersAction,
                               argparse._VersionAction,
                               VersionAction)):
            continue

        # named parameters handled here
        if action.option_strings:
            val = getattr(args, action.dest, None)
            if isinstance(action, argparse._StoreTrueAction):
                if val:
                    cmd.append(action.option_strings[-1])
                continue
            if isinstance(action, argparse._StoreFalseAction):
                if not val:
                    cmd.append(action.option_strings[-1])
                continue

            flag = next((o for o in action.option_strings if o.startswith("--")),
                        action.option_strings[0])

            # this little bit deals with cases where the value is a list or tuple
            if isinstance(val, (list, tuple)):
                for v in val:
                    cmd.append(flag)
                    cmd.append(str(v))
            else:
                cmd.append(flag)
                cmd.append(str(val))

        # positional parameters handled here
        else:
            val = getattr(args, action.dest, None)
            if val is None:
                continue
            if isinstance(val, (list, tuple)):
                cmd.extend(str(v) for v in val)
            else:
                cmd.append(str(val))

    return shlex.join(cmd)


def wrap_help(text, margin=30):
    term_width = shutil.get_terminal_size((80, 24)).columns
    help_width = max(10, term_width - margin)
    cleaned_text = " ".join(text.split())

    return textwrap.fill(cleaned_text, width=help_width)


def wrap_multiline_help(text, margin=30):
    term_width = shutil.get_terminal_size((80, 24)).columns
    help_width = max(10, term_width - margin)

    wrapped_lines = []

    for line in text.splitlines():
        # Separate the leading spaces from the actual text using regex
        match = re.match(r"^(\s*)(.*)$", line)
        lead = match.group(1)
        body = match.group(2)

        # If the line is empty, just keep the whitespace and move on
        if not body:
            wrapped_lines.append(lead)
            continue

        # Wrap the text, applying the extracted leading spaces to every new line it creates
        wrapped = textwrap.fill(
            body,
            width=help_width,
            initial_indent=lead,
            subsequent_indent=lead,
            break_on_hyphens=False
        )
        wrapped_lines.append(wrapped)

    return "\n".join(wrapped_lines)
