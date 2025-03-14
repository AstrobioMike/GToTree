import textwrap
import sys
from importlib.metadata import version

tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


def gtotree_header():
    header = f"""

                                   GToTree v{version('GToTree')}
                         (github.com/AstrobioMike/GToTree)
    """
    return header


def get_version():
    return version('GToTree')


def color_text(text, color = 'green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    print(textwrap.fill(text, width = 80, initial_indent = "  ",
          subsequent_indent = "  ", break_on_hyphens = False))


def report_message(message, color = "yellow"):
    print("")
    wprint(color_text(message, color))


def report_missing_input_genomes_file(path, flag):
    report_message(f'You specified "{path}" as a source of input genomes to use (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_early_exit(message = None, color = "red", suggest_help = False):
    if message:
        print("")
        wprint(color_text(message, color))
    if suggest_help:
        print("\n  See `GToTree -h` for more info.")
    print("\nExiting for now :(\n")
    sys.exit(1)
