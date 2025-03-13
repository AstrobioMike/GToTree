import os
import sys
from dataclasses import dataclass
from tqdm import tqdm # type: ignore
import urllib.request
from gtotree.utils.messaging import report_missing_input_genomes_file, report_message, report_early_exit
import time

@dataclass
class ToolsUsed:
    parallel_used: bool = False
    prodigal_used: bool = False
    taxonkit_used: bool = False
    gtdb_used: bool = False
    fasttree_used: bool = False
    veryfasttree_used: bool = False
    iqtree_used: bool = False
    universal_SCGs_used: bool = False
    pfam_db_used: bool = False
    kofamscan_used: bool = False


def download_with_tqdm(url, filename, target):
    with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=target, ncols = 90) as t:
        def reporthook(block_num, block_size, total_size):
            if total_size > 0:
                t.total = total_size
            t.update(block_size)
        urllib.request.urlretrieve(url, filename, reporthook=reporthook)
    sys.stdout.write("")


def check_path(path, flag):
    if not os.path.isfile(path):
        report_missing_input_genomes_file(path, flag)


def check_expected_single_column_input(path, flag):

    check_for_whitespace(path, flag)
    print(path)
    path = check_line_endings(path, flag)
    print(path)
    path = check_for_duplicates(path, flag)
    print(path)

    return path


def check_for_whitespace(path, flag):

    # there should be no tabs or spaces in these input files
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    for line in lines:
        if " " in line:
            report_early_exit(f'The specified input genomes file "{path}" (passed to `{flag}`) contains spaces in one or more entries. '
                              f'That is not expected and will break things. Please double-check things and remove any whitespace from the file.', suggest_help=True)
        if "\t" in line:
            report_early_exit(f'The specified input genomes file "{path}" (passed to `{flag}`) contains tabs in one or more entries. '
                              f'That is not expected and will break things. Please double-check things and remove any whitespace from the file.', suggest_help=True)


def check_line_endings(path, flag):

    # checking for any CLRF line endings and creating a new file if so
    with open(path, 'rb') as f:
        content = f.read()

    if b'\r\n' in content:
        new_filename = f"{path}-unix"
        new_content = content.replace(b'\r\n', b'\n')
        with open(new_filename, 'wb') as f:
            f.write(new_content)

        report_message(f'Input file "{path}" (passed to `{flag}`) had Windows-formatting that would have caused problems. '
                       f'A modified version was created, "{new_filename}", which will be used.')

        path = new_filename
        time.sleep(2)

    return path


def check_for_duplicates(path, flag):

    # checking for duplicates in the input file
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) != len(set(lines)):
        new_filename = f"{path}-unique"
        new_lines = list(set(lines))
        with open(new_filename, 'w') as f:
            f.write("\n".join(new_lines))

        report_message(f'Input file "{path}" (passed to `{flag}`) had duplicate entries that would have caused problems. '
                       f'A modified version was created, "{new_filename}", which will be used.')

        path = new_filename
        time.sleep(2)

    return path
