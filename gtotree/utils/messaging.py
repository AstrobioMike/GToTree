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


def wprint(text, width = 80):
    # print(textwrap.fill(text, width = 80, initial_indent = "  ",
        #   subsequent_indent = "  ", break_on_hyphens = False))
    wrapper = textwrap.TextWrapper(width=width,
                                   initial_indent="  ",
                                   subsequent_indent="  ",
                                   break_on_hyphens=False)
    # Split the text into lines (or paragraphs)
    paragraphs = text.splitlines()
    # Wrap each paragraph individually
    wrapped_paragraphs = [wrapper.fill(par) if par.strip() else par for par in paragraphs]
    # Join them back preserving newlines
    print("\n".join(wrapped_paragraphs))


def report_message(message, color = "yellow"):
    print("")
    wprint(color_text(message, color))


def report_missing_input_genomes_file(path, flag):
    report_message(f'You specified "{path}" as a source of input genomes to use (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_missing_pfam_targets_file(path, flag):
    report_message(f'You specified "{path}" as a source of Pfam targets to search each genome for (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_missing_ko_targets_file(path, flag):
    report_message(f'You specified "{path}" as a source of KO targets to search each genome for (passed to `{flag}`), but that file can\'t be found.')
    report_early_exit()


def report_early_exit(message = None, color = "red", suggest_help = False):
    if message:
        print("")
        wprint(color_text(message, color))
    if suggest_help:
        print("\n  See `GToTree -h` for more info.")
    print("\nExiting for now :(\n")
    sys.exit(1)


def report_notice(message, color = "yellow"):
    print("")
    print(color_text("  ********************************** NOTICE **********************************  ", "yellow"))
    print(message)
    print(color_text("  ****************************************************************************  ", "yellow"))


### specific notices
def many_genomes_notice(total_input_genomes):
    return (
    f"""    We seem to be aiming to work with {total_input_genomes} genomes. This is quite a bit, and
    the time individual gene alignments can take can quickly become prohibitive
    with many genomes like this. By default, GToTree is going to use the 'super5'
    muscle algorithm to help speed up the alignments for this run. If you don't
    want this to happen, you should cancel this run now with `ctrl + c` now, and
    then add the `-X` flag to the GToTree command and re-run it.

    More info can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider#working-with-many-genomes

    And while we're chatting, you may also want to consider using \"prepresentative\" genomes
    if you're not already. More info on that can be found here:
      github.com/AstrobioMike/GToTree/wiki/things-to-consider#consider-using-representative-genomes

    We will wait 30 seconds before continuing with our regularly scheduled program :)"""
)
