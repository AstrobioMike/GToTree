#!/usr/bin/env python

import os
import sys
import math
import shutil
import argparse
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc.messaging import color_text, report_message
from gtotree.utils.hmms.scg_hmms_setup import (check_hmm_location_var_is_set,
                                              read_in_hmm_summary_table)


LAYOUT_COLUMNS = ("rank", "domain", "parent")

CROSS_DOMAIN_RANKS = ("universal", "multi-domain")

EMPTY = {"", "na", "none", "nan", "-"}

MIN_COLUMN_GAP = 3
SECTION_GAP = 5
LEFT_MARGIN = 4

# Ceiling on how many columns one domain can spread across, and the floor on how full
# each of those columns has to be. Without the floor, a domain with 14 sets would split
# into three columns of five, which is wider and no shorter than two of seven.
MAX_COLUMNS = 3
MIN_CELLS_PER_COLUMN = 6


def build_parser(parent_subparsers=None):

    desc = ("This program lists the pre-packaged HMM SCG-sets available in GToTree. "
            "See github.com/AstrobioMike/GToTree/wiki/SCG-sets for more info.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "hmms",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `gtt hmms`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    optional = parser.add_argument_group("Optional Parameters")
    optional.add_argument(
        "--flat",
        action="store_true",
        help="List every set in one alphabetical column, without grouping them by "
             "taxonomy",
    )
    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()
    args = parser.parse_args()

    print(color_text("\n                   GToTree pre-packaged HMM SCG-sets", "yellow"))
    print("   See github.com/AstrobioMike/GToTree/wiki/SCG-sets for more info\n")

    hmm_data_dir = get_writable_hmm_dir()
    df = read_in_hmm_summary_table()

    if args.flat or not has_layout_info(df):
        report_available_scg_sets(hmm_data_dir, df)
    else:
        report_available_scg_sets_grouped(hmm_data_dir, df)


def get_writable_hmm_dir():
    hmm_data_dir = check_hmm_location_var_is_set()

    if not os.path.isdir(hmm_data_dir):
        try:
            os.makedirs(hmm_data_dir, exist_ok=True)
        except OSError:
            report_message(
                "The 'GToTree_HMM_dir' location does not exist and can't be "
                "created :( Use `gtt data locations check` and `gtt data locations set` to check and configure.",
                "yellow",
            )
            sys.exit(1)

    if not os.access(hmm_data_dir, os.W_OK):
        report_message(
            "The 'GToTree_HMM_dir' location is not writable for you :( "
            "Use `gtt data locations check` and `gtt data locations set` to check and configure.",
            "yellow",
        )
        sys.exit(1)

    return hmm_data_dir


def report_available_scg_sets(hmm_data_dir, df):
    """
    The original flat listing: one line per set, in whatever order the table is in
    """
    print(f"   The environment variable {color_text('GToTree_HMM_dir', 'green')} is set to:")
    print(f"     {hmm_data_dir}\n")

    print(f"   The {len(df)} available pre-packaged HMM SCG-sets include:\n")

    for _, row in df.iterrows():
        gene_set = str(row["file"]).replace(".hmm", "")
        num_genes = row["num_genes"]
        print("\t   {:<30} {:>14}".format(gene_set, f"({num_genes} genes)"))

    print_footer(hmm_data_dir)


################################################################################
# reading the table
################################################################################

def _clean(value):
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in EMPTY else text


def has_layout_info(df):
    if not all(column in df.columns for column in LAYOUT_COLUMNS):
        return False
    return any(_clean(value) for value in df["rank"])


def parse_sets(df):
    sets = []
    for _, row in df.iterrows():
        name = _clean(row.get("file"))
        if not name:
            continue
        if name.lower().endswith(".hmm"):
            name = name[:-4]
        sets.append({
            "name": name,
            "num_genes": _clean(row.get("num_genes")) or "?",
            "rank": _clean(row.get("rank")).lower(),
            "domain": _clean(row.get("domain")),
            "parent": _clean(row.get("parent")),
            "gtdb_release": _clean(row.get("gtdb_release")),
        })
    return sets


def group_sets(sets):
    """
    Sort the sets into the shape they're displayed in
    """
    cross_domain = [s for s in sets if not s["domain"]]
    cross_domain.sort(key=lambda s: (_cross_domain_order(s["rank"]), s["name"].lower()))

    by_domain = {}
    for entry in sets:
        if entry["domain"]:
            by_domain.setdefault(entry["domain"], []).append(entry)

    domains = []
    for domain_name in sorted(by_domain):
        members = by_domain[domain_name]
        domain_set = next((s for s in members if s["rank"] == "domain"), None)
        rest = [s for s in members if s is not domain_set]

        # A set whose parent is also published nests under it; anything else stands on
        # its own line at the top level of the domain, so nothing can go missing just
        # because its parent didn't get built.
        names = {s["name"].lower() for s in rest}
        children = {}
        for entry in rest:
            key = entry["parent"].lower()
            if key in names and key != entry["name"].lower():
                children.setdefault(key, []).append(entry)

        nested = {id(e) for kids in children.values() for e in kids}
        blocks = []
        for entry in sorted((s for s in rest if id(s) not in nested),
                            key=lambda s: s["name"].lower()):
            kids = sorted(children.get(entry["name"].lower(), []),
                          key=lambda s: s["name"].lower())
            blocks.append((entry, kids))

        domains.append((domain_name, domain_set, blocks))

    domains.sort(key=lambda d: (sum(1 + len(k) for _, k in d[2]), d[0]))
    return cross_domain, domains


def _cross_domain_order(rank):
    try:
        return CROSS_DOMAIN_RANKS.index(rank)
    except ValueError:
        return len(CROSS_DOMAIN_RANKS)


################################################################################
# laying it out
################################################################################

def _cell(indent, entry):
    return (indent, entry["name"], f"({entry['num_genes']})")


def _format_cell(cell, label_width, count_width):
    indent, name, count = cell
    label = " " * indent + name
    return f"{label:<{label_width}}{count:>{count_width}}"


def _widths(cells):
    label_width = max(len(" " * i + name) for i, name, _ in cells) + 2
    count_width = max(len(count) for _, _, count in cells) + 1
    return label_width, count_width


def pack_blocks(blocks, n_columns):
    """
    Deal blocks into `n_columns` lists of cells, filling each column top to bottom
    """
    if n_columns <= 1:
        return [[cell for block in blocks for cell in block]]

    total = sum(len(block) for block in blocks)
    per_column = math.ceil(total / n_columns)

    columns, current = [], []
    for block in blocks:
        if current and len(current) + len(block) > per_column \
                and len(columns) < n_columns - 1:
            columns.append(current)
            current = []
        current.extend(block)
    columns.append(current)

    while len(columns) < n_columns:
        columns.append([])
    return columns


def _format_row(columns, row_index, label_width, count_width):
    parts = []
    for column in columns:
        if row_index < len(column):
            parts.append(_format_cell(column[row_index], label_width, count_width))
        else:
            parts.append(" " * (label_width + count_width))
    return (" " * MIN_COLUMN_GAP).join(parts)


################################################################################
# fitting the domains side by side
################################################################################
#
# Every domain gets its own vertical strip, laid out left to right, because the whole
# point of this listing is to be short: a listing taller than the terminal loses its
# own header, and what you're left staring at is the bottom of the alphabet with no
# indication of what it belongs to.
#
# Archaea running out well before Bacteria does leave dead space on the right, but
# dead space costs nothing, whereas the stacked version costs height. Where there's
# room, Bacteria spreads across more than one column of its own strip, which is where
# most of the height goes.

def build_section(domain_name, domain_set, blocks):
    block_cells = [[_cell(2, parent)] + [_cell(4, kid) for kid in kids]
                   for parent, kids in blocks]
    header = _cell(0, domain_set) if domain_set else (0, domain_name, "")
    flat = [cell for block in block_cells for cell in block]
    label_width, count_width = _widths(flat + [header])
    return {
        "header": header,
        "blocks": block_cells,
        "n_cells": len(flat),
        "label_width": label_width,
        "count_width": count_width,
        "cell_width": label_width + count_width,
        "columns": 1,
    }


def _strip_width(section):
    return (section["columns"] * section["cell_width"]
            + (section["columns"] - 1) * MIN_COLUMN_GAP)


def _total_width(sections):
    return (LEFT_MARGIN + sum(_strip_width(s) for s in sections)
            + SECTION_GAP * (len(sections) - 1))


def _strip_height(section):
    """Rows of blocks in this strip at its current column count, header excluded."""
    return max((len(c) for c in pack_blocks(section["blocks"], section["columns"])),
               default=0)


def allocate_columns(sections, terminal_width):
    """
    Widen the tallest strip, one column at a time, for as long as it fits and helps.

    Greedy on height rather than solving it properly: with two or three domains there
    isn't enough of a search space for anything cleverer to find a different answer.
    Returns False if the domains can't sit side by side at one column each, which is
    the caller's cue to stack them instead.
    """
    for section in sections:
        section["columns"] = 1
    if _total_width(sections) > terminal_width:
        return False

    while True:
        # tallest first, so the domain driving the overall height gets the next column
        for section in sorted(sections, key=_strip_height, reverse=True):
            if section["columns"] >= MAX_COLUMNS:
                continue
            wider = section["columns"] + 1
            if math.ceil(section["n_cells"] / wider) < MIN_CELLS_PER_COLUMN:
                continue

            before = _strip_height(section)
            section["columns"] = wider
            if _total_width(sections) > terminal_width or _strip_height(section) >= before:
                section["columns"] = wider - 1
                continue
            break
        else:
            return True


def render_section(section):
    """A domain's strip as a list of equal-width lines, header first."""
    label_width, count_width = section["label_width"], section["count_width"]
    width = _strip_width(section)

    header = _format_cell(section["header"], label_width, count_width)
    lines = [color_text(f"{header:<{width}}", "green")]

    columns = pack_blocks(section["blocks"], section["columns"])
    for row_index in range(max((len(c) for c in columns), default=0)):
        lines.append(_format_row(columns, row_index, label_width, count_width))
    return lines, width


def _side_by_side(sections, terminal_width):
    """Domains as strips laid left to right, or None if they don't fit at all."""
    if not allocate_columns(sections, terminal_width):
        return None

    strips = [render_section(section) for section in sections]
    height = max(len(strip) for strip, _ in strips)

    lines = [""]
    for row_index in range(height):
        parts = [strip[row_index] if row_index < len(strip) else " " * width
                 for strip, width in strips]
        lines.append((" " * LEFT_MARGIN + (" " * SECTION_GAP).join(parts)).rstrip())
    return lines


def _stacked(sections, terminal_width):
    """Domains one under the other, each spread across as many columns as fit."""
    lines = []
    for section in sections:
        section["columns"] = 1
        while (section["columns"] < MAX_COLUMNS
               and math.ceil(section["n_cells"] / (section["columns"] + 1))
                   >= MIN_CELLS_PER_COLUMN
               and _strip_width(dict(section, columns=section["columns"] + 1))
                   + LEFT_MARGIN <= terminal_width):
            section["columns"] += 1
        lines.append("")
        strip, _ = render_section(section)
        lines.extend((" " * LEFT_MARGIN + line).rstrip() for line in strip)
    return lines


def layout_lines(cross_domain, domains, terminal_width=80):
    """
    The whole listing as a list of lines. Pure -- no printing, no terminal poking --
    so the layout can be tested without a terminal.
    """
    lines = []

    if cross_domain:
        cells = [_cell(0, entry) for entry in cross_domain]
        label_width, count_width = _widths(cells)
        for row_index in range(len(cells)):
            lines.append(" " * LEFT_MARGIN
                         + _format_row([cells], row_index, label_width, count_width))

    if not domains:
        return lines

    # Both arrangements, and whichever is shorter wins
    sections = [build_section(*domain) for domain in domains]
    side_by_side = _side_by_side(sections, terminal_width)
    stacked = _stacked([build_section(*domain) for domain in domains], terminal_width)

    if side_by_side is not None and len(side_by_side) <= len(stacked):
        return lines + side_by_side
    return lines + stacked


################################################################################
# printing
################################################################################

def report_available_scg_sets_grouped(hmm_data_dir, df):
    print(f"   The environment variable {color_text('GToTree_HMM_dir', 'green')} is set to:")
    print(f"     {hmm_data_dir}\n")

    sets = parse_sets(df)
    cross_domain, domains = group_sets(sets)

    # release = next((s["gtdb_release"] for s in sets if s["gtdb_release"]), "")
    # built_from = f", built from GTDB {release}" if release else ""

    # print(f"   The {len(sets)} available pre-packaged HMM SCG-sets "
    #       f"(# of genes{built_from}):\n")

    print(f"   The {len(sets)} available pre-packaged HMM SCG-sets include (# of genes):\n")

    terminal_width = shutil.get_terminal_size((80, 24)).columns
    for line in layout_lines(cross_domain, domains, terminal_width):
        print(line)

    print_footer(hmm_data_dir)


def print_footer(hmm_data_dir):
    table_path = os.path.join(hmm_data_dir, "hmm-sources-and-info.tsv")
    print(f"\n   Details can be found in: \n     {table_path}\n\n")


if __name__ == "__main__":
    main()
