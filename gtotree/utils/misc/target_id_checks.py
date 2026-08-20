"""
Shape checks for the target-ID lists passed to `-p`/`--target-pfams` and
`-K`/`--target-kos`
"""

import re
from dataclasses import dataclass


class TargetIDFormatError(Exception):
    """
    A target-ID file whose entries aren't the shape the flag they were given to expects.
    """


@dataclass(frozen=True)
class TargetIDFormat:
    key: str              # "pfam" / "ko"; matches TargetSearchSpec's get_spec() names
    label: str            # "Pfam"
    label_plural: str     # "Pfams"
    flag: str             # "-p"
    flag_long: str        # "--target-pfams"
    example: str          # "PF00789"
    pattern: re.Pattern   # the full canonical form of one ID
    shape: str            # prose describing `pattern`, for the error message


PFAM_FORMAT = TargetIDFormat(
    key="pfam",
    label="Pfam",
    label_plural="Pfams",
    flag="-p",
    flag_long="--target-pfams",
    example="PF00789",
    pattern=re.compile(r"PF\d{4,6}(?:\.\d+)?"),
    shape=('start with "PF" followed by five digits, optionally with a version '
           "suffix (e.g., 'PF00789' or 'PF00789.19')"),
)

KO_FORMAT = TargetIDFormat(
    key="ko",
    label="KO",
    label_plural="KOs",
    flag="-K",
    flag_long="--target-kos",
    example="K01601",
    pattern=re.compile(r"K\d{4,6}"),
    shape='start with "K" followed by five digits (e.g., \'K01601\')',
)

FORMATS = {fmt.key: fmt for fmt in (PFAM_FORMAT, KO_FORMAT)}

# verdicts that aren't "this is an ID of some other target type"
WRONG_CASE = "case"
MALFORMED = "malformed"

# how many bad entries to spell out before summarizing the rest, matching the
# genome-ID-collision reporting's approach
MAX_REPORTED_BAD_IDS = 10


def read_target_ids(path):
    """
    [(1-based line number, stripped entry)] for every non-blank line
    """
    entries = []
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            entry = line.strip()
            if entry:
                entries.append((lineno, entry))
    return entries


def classify_entry(entry, fmt):
    """
    None if `entry` is a well-formed ID for `fmt`, else why it isn't

    A non-None verdict is either another format's key (the entry is a perfectly good ID
    of the wrong type), WRONG_CASE, or MALFORMED. The other-type check comes first
    because it's the diagnosis that leads to an actionable message.
    """
    if fmt.pattern.fullmatch(entry):
        return None

    for other in FORMATS.values():
        if other.key != fmt.key and other.pattern.fullmatch(entry):
            return other.key

    if fmt.pattern.fullmatch(entry.upper()):
        return WRONG_CASE

    return MALFORMED


def find_bad_entries(entries, fmt):
    """
    [(lineno, entry, verdict)] for the entries that aren't well-formed IDs for `fmt`.
    """
    return [(lineno, entry, verdict)
            for lineno, entry in entries
            if (verdict := classify_entry(entry, fmt)) is not None]


def check_target_id_file(path, key):
    """
    Validate one target-ID file, raising TargetIDFormatError if it doesn't hold IDs of
    the expected type.

    `key` is "pfam" or "ko". Returns the parsed entries on success so a caller that
    wants them doesn't have to read the file again.
    """
    fmt = FORMATS[key]
    entries = read_target_ids(path)

    if not entries:
        raise TargetIDFormatError(
            f'The {fmt.label} targets file "{path}" (passed to `{fmt.flag}`) is empty, '
            f"so there'd be nothing to search for.\n\n"
            f"It should hold one {fmt.label} ID per line -- they {fmt.shape}.")

    bad = find_bad_entries(entries, fmt)
    if bad:
        raise TargetIDFormatError(build_message(path, fmt, bad, len(entries)))

    return entries


def wholesale_swap(bad, total):
    """
    The TargetIDFormat this file appears to be entirely, or None.

    All-or-nothing on purpose: a file where every single entry is a well-formed ID of
    one other type was handed to the wrong flag, and saying so plainly beats listing
    its contents back. A file with a couple of stray KOs mixed into a Pfam list is a
    different mistake and gets the itemized message instead.
    """
    if len(bad) != total:
        return None

    verdicts = {verdict for _lineno, _entry, verdict in bad}
    if len(verdicts) != 1:
        return None

    return FORMATS.get(verdicts.pop())


def _sample_block(bad, with_line_numbers):
    shown = bad[:MAX_REPORTED_BAD_IDS]
    if with_line_numbers:
        lines = [f"    line {lineno}: {entry}" for lineno, entry, _v in shown]
    else:
        lines = [f"    {entry}" for _lineno, entry, _v in shown]

    remaining = len(bad) - len(shown)
    if remaining:
        lines.append(f"    ...and {remaining} more.")

    return "\n".join(lines)


def build_message(path, fmt, bad, total):
    """
    The user-facing explanation for a file that failed the check.
    """
    swapped = wholesale_swap(bad, total)
    if swapped is not None:
        return (
            f'The {fmt.label} targets file "{path}" (passed to `{fmt.flag}`) doesn\'t '
            f"hold {fmt.label} IDs.\n\n"
            f"All {total:,} of its entries look like {swapped.label} IDs instead:\n\n"
            f"{_sample_block(bad, with_line_numbers=False)}\n\n"
            f"{fmt.label} IDs {fmt.shape}. If you meant to search for "
            f"{swapped.label_plural}, pass this file to "
            f"`{swapped.flag}`/`{swapped.flag_long}` instead.")

    count = len(bad)
    entry_word = "entry" if count == 1 else "entries"
    verb = "doesn't" if count == 1 else "don't"

    message = (
        f'The {fmt.label} targets file "{path}" (passed to `{fmt.flag}`) has {count:,} '
        f"{entry_word} that {verb} look like {fmt.label} IDs.\n\n"
        f"Problematic ones include:\n\n"
        f"{_sample_block(bad, with_line_numbers=True)}\n\n"
        f"{fmt.label} IDs {fmt.shape}.")

    verdicts = {verdict for _lineno, _entry, verdict in bad}

    if WRONG_CASE in verdicts:
        message += (" They're also matched case-sensitively, so any lowercase ones "
                    "need to be uppercased.")

    for other_key in sorted(verdicts & set(FORMATS)):
        other = FORMATS[other_key]
        message += (f" Some of them look like {other.label} IDs. Those go to "
                    f"`{other.flag}`/`{other.flag_long}`.")

    return message
