#!/usr/bin/env python

import os
import sys
import argparse
from gtotree.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from gtotree.utils.misc.messaging import color_text, report_message
from gtotree.utils.hmms.scg_hmms_setup import (check_hmm_location_var_is_set,
                                              read_in_hmm_summary_table)


LAYOUT_COLUMNS = ("rank", "domain", "parent")

CROSS_DOMAIN_RANKS = ("universal", "multi-domain")

DOMAIN_ORDER = ("Bacteria", "Archaea")

EMPTY = {"", "na", "none", "nan", "-"}

SECTION_GAP = 5
LEFT_MARGIN = 8


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
        help="List every set in one alphabetical column rather than taxonomic groupings",
    )
    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()
    args = parser.parse_args()

    print(color_text("\n                   GToTree pre-packaged HMM SCG-sets", "yellow"))
    print("    See github.com/AstrobioMike/GToTree/wiki/SCG-sets for more info\n")

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
    print(f"    The environment variable {color_text('GToTree_HMM_dir', 'green')} is set to:")
    print(f"      {hmm_data_dir}\n")

    print(f"    The {len(df)} available pre-packaged HMM SCG-sets include:\n")

    for _, row in df.iterrows():
        gene_set = str(row["file"]).replace(".hmm", "")
        num_genes = row["num_genes"]
        print("\t   {:<30} {:>14}".format(gene_set, f"({num_genes} genes)"))

    print_footer(hmm_data_dir)

    print(f"    And see `gtt gen-scg-hmms -h` if you'd like to make your own!\n")



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

    domains.sort(key=lambda d: _domain_order(d[0], sum(1 + len(k) for _, k in d[2])))
    return cross_domain, domains


def _domain_order(domain_name, n_sets):
    """
    Sort key putting the domains in DOMAIN_ORDER first, then anything else biggest
    first
    """
    try:
        return (DOMAIN_ORDER.index(domain_name), 0, "")
    except ValueError:
        return (len(DOMAIN_ORDER), -n_sets, domain_name.lower())


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


def section_cells(domain_name, domain_set, blocks):
    """
    One domain's cells, top to bottom: the domain-level set as the heading, then each
    phylum, then any class-level sets indented under the phylum they sit inside.
    """
    header = _cell(0, domain_set) if domain_set else (0, domain_name, "")
    return [header] + [cell for parent, kids in blocks
                       for cell in [_cell(2, parent)] + [_cell(4, k) for k in kids]]


def build_section(domain_name, domain_set, blocks, label_width, count_width):
    """
    One domain's strip, as a list of equal-width lines
    """
    return [_format_cell(cell, label_width, count_width)
            for cell in section_cells(domain_name, domain_set, blocks)]


def layout_lines(cross_domain, domains):
    if not cross_domain and not domains:
        return []

    cross_cells = [_cell(0, entry) for entry in cross_domain]
    label_width, count_width = _widths(
        cross_cells + [cell for domain in domains for cell in section_cells(*domain)])
    cell_width = label_width + count_width

    lines = [(" " * LEFT_MARGIN + _format_cell(cell, label_width, count_width)).rstrip()
             for cell in cross_cells]

    if not domains:
        return lines

    strips = [build_section(*domain, label_width, count_width) for domain in domains]

    lines.append("")
    for row_index in range(max(len(strip) for strip in strips)):
        parts = [strip[row_index] if row_index < len(strip) else " " * cell_width
                 for strip in strips]
        lines.append((" " * LEFT_MARGIN + (" " * SECTION_GAP).join(parts)).rstrip())

    return lines


################################################################################
# printing
################################################################################

def report_available_scg_sets_grouped(hmm_data_dir, df):
    # print(f"    The environment variable {color_text('GToTree_HMM_dir', 'green')} is set to:")
    # print(f"      {hmm_data_dir}\n")

    sets = parse_sets(df)
    cross_domain, domains = group_sets(sets)

    # release = next((s["gtdb_release"] for s in sets if s["gtdb_release"]), "")
    # built_from = f", built from GTDB {release}" if release else ""

    # print(f"   The {len(sets)} available pre-packaged HMM SCG-sets "
    #       f"(# of genes{built_from}):\n")

    print(f"\n    The {len(sets)} available pre-packaged HMM SCG-sets include (# of genes):\n")

    for line in layout_lines(cross_domain, domains):
        print(line)

    print(f"\n\n    Any of these can be passed to a main gtotree run, e.g., `gtotree -H Bacteria ...`")

    print_footer(hmm_data_dir)

    print(f"    And see `gtt gen-scg-hmms -h` if you'd like to make your own!\n")


def print_footer(hmm_data_dir):
    table_path = os.path.join(hmm_data_dir, "hmm-sources-and-info.tsv")
    print(f"\n    Details of each can be found in: \n      {table_path}\n")


if __name__ == "__main__":
    main()
