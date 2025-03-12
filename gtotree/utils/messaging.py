import textwrap
import sys

tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


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


def report_failure(message, color = "red"):
    print("")
    wprint(color_text(message, color))
    print("\nExiting for now :(\n")
    sys.exit(1)
