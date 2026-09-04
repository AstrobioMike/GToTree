"""
Explaining an EMPTY taxon selection, for every surface that pulls genomes by taxonomy.

The stages, and what an empty result at each one means:

    n_candidates == 0       a POOL filter, or the taxon is genuinely empty here.
                            Which filter comes from diagnose_empty_pool(), which
                            peels them back one at a time.
    n_after_liveness == 0   every candidate is suppressed/removed at NCBI
    n_after_exclusion == 0  every candidate was named in the --exclusion-list
    n_after_floor == 0      no candidate cleared the checkm quality floor
    n_unassigned_group      every survivor is unclassified at the derep rank, so
      == n_after_floor      there was no group for any of them to represent

"""

from gtotree.utils.taxonomy.get_accs_shared import (ASSEMBLY_LEVELS,
                                                  with_filters_implicated_note,
                                                  with_filters_note)

#: the --assembly-level word a user typed, from the NCBI string it maps to
LEVEL_FLAG_VALUES = {v: k for k, v in ASSEMBLY_LEVELS.items()}

#: how each pool filter is named back to the user when it's the one to loosen
POOL_FILTER_FLAGS = {
    "assembly_levels": "--assembly-level",
    "accession_prefixes": "--ncbi-section",
    "reps_only": "--representatives-only",
}


def _oxford(items, conjunction="and"):
    """'a', 'a and b', 'a, b, and c'."""
    items = [str(i) for i in items if i]
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return f"{items[0]} {conjunction} {items[1]}"
    return f"{', '.join(items[:-1])}, {conjunction} {items[-1]}"


def _level_words(levels):
    """NCBI assembly_level strings -> the --assembly-level words that ask for them."""
    return [LEVEL_FLAG_VALUES.get(lvl, str(lvl).lower()) for lvl in (levels or [])]


def _with_emoticon(headline, emoticon):
    """Swap a headline's trailing period for `emoticon` (no-op when None)."""
    if not emoticon:
        return headline
    if headline.endswith("."):
        return f"{headline[:-1]} {emoticon}"
    return f"{headline} {emoticon}"


def _genomes(n):
    return f"{n:,} genome" + ("" if n == 1 else "s")


def _restricting_section(ncbi_section):
    """The --ncbi-section value only when it actually narrows anything."""
    section = str(ncbi_section or "").strip().lower()
    # 'both' imposes no accession-prefix restriction, and GTDB ignores the flag
    # outright, so neither can be what emptied a pool
    return section if section in ("refseq", "genbank") else None


def _pool_filter_clause(name, *, assembly_levels, ncbi_section, present_levels,
                        reps_only_requested):
    """The 'but none are ...' half of a sentence blaming one pool filter."""
    if name == "assembly_levels":
        asked = _oxford(_level_words(assembly_levels), "or")
        clause = f"none are at `--assembly-level {asked}`"
        present = _oxford(_level_words(present_levels), "and")
        if present:
            clause += f" (it has: {present})"
        return clause

    if name == "accession_prefixes":
        section = _restricting_section(ncbi_section)
        if section:
            return f"none are in the `--ncbi-section {section}` part of NCBI"
        return "none are in the requested part of NCBI"

    if name == "reps_only":
        if reps_only_requested:
            return ("none are flagged as representative genomes, so "
                    "`--representatives-only` leaves nothing")
        # the source default did this, not a flag the user typed
        return ("none of them are representative genomes, which is the pool this "
                "source uses by default")

    return f"none survive `{POOL_FILTER_FLAGS.get(name, name)}`"


def _active_pool_flags(*, assembly_levels, ncbi_section, reps_only_requested):
    """The pool-narrowing flags actually in play, named as the user set them."""
    flags = []
    if assembly_levels:
        flags.append(f"`--assembly-level {_oxford(_level_words(assembly_levels), 'or')}`")
    section = _restricting_section(ncbi_section)
    if section:
        flags.append(f"`--ncbi-section {section}`")
    if reps_only_requested:
        flags.append("`--representatives-only`")
    return flags


def explain_empty_selection(selection, *, assembly_levels=None, ncbi_section=None,
                            reps_only_requested=False,
                            min_completeness=None, max_contamination=None):
    """
    Why did this selection come back with nothing?

    Returns (detail, filtered):

        detail   -- one or two sentences naming the stage that emptied it and what to
                    loosen, or "" when nothing more can be said than the headline
        filtered -- True when a user-specified FILTER is what did it, which is the
                    caller's cue to add "with the specified filters" to the headline.
                    False for the causes no flag produced (an empty taxon, everything
                    suppressed at NCBI, a source's default representatives-only pool,
                    everything unclassified at the derep rank), where blaming the
                    user's filters would be wrong.

    `reps_only_requested` is whether the USER asked for representatives-only, as
    opposed to inheriting it from the source (GTDB defaults to it). The effective
    value is read off the attrition record; this argument only decides whether the
    message may point at a flag.

    Deliberately returns text rather than printing: the CLIs wrap and colour their own
    output, and both repos' wrappers raise this as an exception message.
    """
    a = getattr(selection, "attrition", None)
    if a is None:
        return "", False

    rank = a.resolved_rank or selection.resolved_rank

    # ---- stage 1: the pool filters ----------------------------------------
    if not a.n_candidates:
        if not a.diagnosed or not a.n_unfiltered:
            # nothing under the taxon even unfiltered; no filter to blame
            return "", False

        held = f"It has {_genomes(a.n_unfiltered)} at '{rank}' rank, but "

        if a.pool_culprits:
            clauses = [_pool_filter_clause(name,
                                           assembly_levels=assembly_levels,
                                           ncbi_section=ncbi_section,
                                           present_levels=a.present_levels,
                                           reps_only_requested=reps_only_requested)
                       for name in a.pool_culprits]
            # a source's default representatives pool is not one of "the user's
            # filters", so don't let it alone trigger the headline's filter note
            blames_a_flag = any(name != "reps_only" or reps_only_requested
                                for name in a.pool_culprits)
            return held + _oxford(clauses, "and") + ".", blames_a_flag

        # no single filter explains it, so don't accuse one
        flags = _active_pool_flags(assembly_levels=assembly_levels,
                                   ncbi_section=ncbi_section,
                                   reps_only_requested=reps_only_requested)
        if flags:
            return (held + f"none survive {_oxford(flags, 'plus')} together. "
                    f"Loosening any one of them should bring genomes back."), True
        return "", False

    # ---- stage 2: liveness screening --------------------------------------
    if not a.n_after_liveness:
        return (f"All {_genomes(a.n_candidates)} under it are suppressed or removed "
                f"at NCBI, so none can be downloaded."), False

    # ---- stage 3: the exclusion list --------------------------------------
    if not a.n_after_exclusion:
        return (f"Every one of its {_genomes(a.n_after_liveness)} was named in the "
                f"`--exclusion-list`."), True

    # ---- stage 4: the quality floor ---------------------------------------
    if not a.n_after_floor:
        bits = []
        if a.n_below_floor:
            asked = []
            if min_completeness is not None:
                asked.append(f"completeness >= {min_completeness:g}")
            if max_contamination is not None:
                asked.append(f"contamination <= {max_contamination:g}")
            bits.append(f"{a.n_below_floor:,} fell outside "
                        f"{_oxford(asked, 'and') or 'the floor'}")
        if a.n_missing_quality:
            bits.append(f"{a.n_missing_quality:,} have no checkm values recorded")
        detail = (f"None of its {_genomes(a.n_after_exclusion)} cleared the quality "
                  f"floor")
        if bits:
            detail += f" ({_oxford(bits, 'and')})"
        return (detail + ". You may want to relax `--min-completeness` / "
                "`--max-contamination`."), True

    # ---- stage 5: dereplication had nothing to group by -------------------
    if a.derep_rank and a.n_unassigned_group >= a.n_after_floor:
        return (f"All {_genomes(a.n_after_floor)} under it are unclassified at "
                f"'{a.derep_rank}', the rank dereplication groups on, so none of them "
                f"can represent a group. Try `--derep-rank off` to take them all, or a "
                f"coarser `--derep-rank`."), False

    return "", False


def empty_selection_message(selection, *, taxon_flag="--wanted-ref-tax", assembly_levels=None,
                            ncbi_section=None, reps_only_requested=False,
                            min_completeness=None, max_contamination=None,
                            emoticon=None):
    """
    The whole user-facing sentence for an empty selection, headline + detail.

    `taxon_flag` is how the calling surface names its taxon option, since the repos
    differ there (`-t` in bit, `--wanted-ref-tax` in GToTree).

    `emoticon` (e.g. ":(") replaces the headline's full stop, for the surfaces that
    have always signed off that way. It is applied AFTER any filters note, because
    the note splices itself in just inside a trailing period -- swapping the period
    out first would leave the note stranded on the end, after the emoticon.
    """
    detail, filtered = explain_empty_selection(
        selection, assembly_levels=assembly_levels, ncbi_section=ncbi_section,
        reps_only_requested=reps_only_requested,
        min_completeness=min_completeness, max_contamination=max_contamination)

    headline = (f"No accessions were found for the {taxon_flag} target "
                f"'{selection.canonical}'.")
    if filtered:
        headline = with_filters_implicated_note(headline)
    headline = _with_emoticon(headline, emoticon)

    return f"{headline} {detail}".strip() if detail else headline


def empty_pull_message(headline, selection=None, *, assembly_levels=None,
                       ncbi_section=None, reps_only_requested=False,
                       min_completeness=None, max_contamination=None,
                       emoticon=None):
    """
    Attach an explanation to a surface's OWN empty-result headline.

    empty_selection_message() builds the whole sentence for the taxon-flag wrappers,
    which all word the headline identically. The get-accs helpers don't: they say
    "No genomes were found under <label>", where <label> already carries the rank,
    the canonical name and the derep rank. So rather than restate their headline
    here, this takes it as given and adds only what they cannot work out locally.

    `selection` is None on the paths that have no selection behind them (an `all`
    pull, a taxid pull). Those keep with_filters_note()'s hedge, which is the honest
    thing to say when nothing in scope knows what was filtered.
    """
    if selection is None:
        return _with_emoticon(with_filters_note(headline), emoticon)

    detail, filtered = explain_empty_selection(
        selection, assembly_levels=assembly_levels, ncbi_section=ncbi_section,
        reps_only_requested=reps_only_requested,
        min_completeness=min_completeness, max_contamination=max_contamination)

    if filtered:
        # we know a filter did it, so say so rather than hedging
        head = _with_emoticon(with_filters_implicated_note(headline), emoticon)
        return f"{head} {detail}".strip()
    if detail:
        # a cause no filter produced; hedging about filters would be wrong
        return f"{_with_emoticon(headline, emoticon)} {detail}".strip()
    return _with_emoticon(with_filters_note(headline), emoticon)
