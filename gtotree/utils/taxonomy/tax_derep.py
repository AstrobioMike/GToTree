"""
Genome/accession DEREPLICATION and quality-picking from the GTDB / NCBI Parquet
assets

The dereplication piece of the taxonomy toolkit. We need "one best genome per
rank" for tree building, which the filter piece (tax_select.py) doesn't do on its own
"""

from gtotree.utils.taxonomy.tax_ranks import (
    RANKS, NA, REFERENCE_VALUE, accession_core, rank_index, validate_derep_rank)
from gtotree.utils.taxonomy.tax_select import (
    SOURCES, select, resolve_taxon, live_accession_cores)

# Thresholds for the size_advice() nudges at the tails of a selection. A tree far
# above SANE_HIGH is unwieldy; far below SANE_LOW is probably too sparse to be useful.
# Advisory only -- they gate warning messages, never the selection itself.
SANE_LOW = 20
SANE_HIGH = 5000


ASSEMBLY_LEVEL_ORDER = {
    "complete genome": 4,
    "chromosome": 3,
    "scaffold": 2,
    "contig": 1,
}

# A RefSeq "reference genome" is preferred over a merely higher-checkm genome,
# but only if it is not actually bad. These gates are deliberately loose, here to
# catch junk, not to replace all "reference" genomes
REF_MIN_COMPLETENESS = 85.0

REF_MAX_CONTAMINATION = 10.0

MISSING_CONTAMINATION = float("inf")

# How many ranks FINER than the target the default derep-rank
DEREP_STEPS = 2

def _derive_default_derep(target_rank):
    """Two ranks finer, clamped at species; None if that would be the target's own rank."""
    i = rank_index(target_rank)
    j = min(i + DEREP_STEPS, len(RANKS) - 1)
    return None if j == i else RANKS[j]

DEFAULT_DEREP_BY_TARGET_RANK = {r: _derive_default_derep(r) for r in RANKS}

def default_derep_rank(wanted_rank):
    """The default derep rank for a target at `wanted_rank`. None == off."""
    return DEFAULT_DEREP_BY_TARGET_RANK[RANKS[rank_index(wanted_rank)]]

def resolve_derep_rank(wanted_rank, derep_rank="auto"):
    """
    Work out the effective derep rank. Returns (rank_or_None, warnings);
    None means no dereplication.

    derep_rank:
      "auto"            -> DEFAULT_DEREP_BY_TARGET_RANK (two ranks finer, clamped)
      "off"/None/"none" -> no dereplication
      an explicit rank  -> honoured, but must not be COARSER than the target
    """
    warnings = []
    w_rank = RANKS[rank_index(wanted_rank)]

    if derep_rank in (None, "off", "none", "None"):
        return None, warnings

    if derep_rank == "auto":
        eff = DEFAULT_DEREP_BY_TARGET_RANK[w_rank]
        if eff is None:
            warnings.append(
                f"Dereplication is off by default for a target at '{w_rank}' rank: the "
                f"default derep rank would be '{w_rank}' itself, which returns a single "
                f"genome. For this run, all genomes under the taxon will be used. If you DID want a "
                f"single best representative (e.g., as an outgroup), pass "
                f"--derep-rank {w_rank} explicitly.")
        return eff, warnings

    problem = validate_derep_rank(wanted_rank, derep_rank)
    if problem:
        raise ValueError(problem)

    return derep_rank, warnings

def size_advice(n_selected, wanted_rank, derep_rank):
    """
    Nudge at the tails.

    The 2-rank-lower default derep-rank keeps most cases in a sane band
      - small phyla may have very few orders
      - an explicit fine derep on a broad taxon can explode (e.g., Bacteria + species -> ~190k)
      - and derep OFF can explode, e.g., `--wanted-ref-tax "Escherichia coli"`
        defaults to derep off and can return tens of thousands of genomes
    """
    if derep_rank is None:
        if n_selected > SANE_HIGH:
            w = rank_index(wanted_rank)
            if w == len(RANKS) - 1:
                return [f"{n_selected:,} reference genomes selected with dereplication "
                        f"off. That is a very large tree. There is no rank finer than "
                        f"'{wanted_rank}' to group by, so the options are all of them like this, or "
                        f"a single best representative (--derep-rank {wanted_rank})."]
            suggestion = DEFAULT_DEREP_BY_TARGET_RANK[RANKS[w]]
            return [f"{n_selected:,} reference genomes selected with dereplication off. "
                    f"That is a very large tree. Consider --derep-rank "
                    f"'{suggestion}'."]
        return []

    d = rank_index(derep_rank)
    if n_selected > SANE_HIGH and d > 0:
        return [f"{n_selected:,} reference genomes selected. That is a very large "
                f"tree. Consider a coarser --derep-rank (e.g., '{RANKS[d - 1]}')."]
    if n_selected < SANE_LOW and d < len(RANKS) - 1:
        return [f"Only {n_selected:,} reference genome(s) selected. Consider a finer "
                f"--derep-rank (e.g., '{RANKS[d + 1]}') for more genomes to be included."]
    return []

def _num(value, default=None):
    if value in (None, "", NA, "na"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default

def mean_contig_length(row, spec):
    """
    genome_size_ungapped / contig_count for a proxy of quality for those with no checkm values
    """
    if not spec.contig_col or not spec.size_col:
        return 0.0

    contigs = _num(row.get(spec.contig_col))
    if not contigs or contigs <= 0:
        return 0.0

    size = _num(row.get(spec.size_col))
    if size is None and spec.size_fallback_col:
        size = _num(row.get(spec.size_fallback_col))
    if size is None or size <= 0:
        return 0.0

    return size / contigs

def first_key(row, spec):
    """Cheapest possible picker: first row in lineage-sorted order. Deterministic."""
    return (row.get(spec.acc_col) or "",)

#: pickers addressable by name from a CLI flag if i want
PICKERS = {
    "quality": None,
    "first": first_key,
}

def resolve_picker(pick):
    """
    `pick` may be a NAME ("quality" / "first") or a CALLABLE (row, spec) -> sort key
    """
    if callable(pick):
        return pick
    if isinstance(pick, str) and pick in PICKERS:
        return PICKERS[pick]
    raise ValueError(
        f"unknown pick policy {pick!r}. Pass one of {sorted(PICKERS)}, or a callable "
        f"(row, spec) -> sort key (lowest wins; see quality_key).")

def quality_key(row, spec):
    """
    Deterministic "best genome in this group" ordering. Lower tuple sorts first.

    Ordering rationale:

    1. RefSeq REFERENCE genomes first (gated on not being too crappy)

    2. Genomes WITH checkm before genomes without. Prokaryotes have checkm and euks don't.
       In practice a group is shouldn't span those, so this should never matter, but a mixed group
       shouldn't rank a euk above a prokaryote just because its NA sorted low

    3. checkm: highest completeness, then lowest contamination

    4. NO checkm (euks): best assembly_level (Complete > Chromosome >
       Scaffold > Contig), then LONGEST MEAN CONTIG

    5. accession, so the result is stable/reproducible
    """
    comp_col, cont_col = spec.quality_cols
    comp = _num(row.get(comp_col))
    cont = _num(row.get(cont_col))
    has_checkm = comp is not None

    is_ref = str(row.get(spec.ref_col) or "").strip().lower() == REFERENCE_VALUE
    if is_ref and has_checkm:
        is_ref = (comp >= REF_MIN_COMPLETENESS
                  and (cont is None or cont <= REF_MAX_CONTAMINATION))

    lvl = ASSEMBLY_LEVEL_ORDER.get(
        str(row.get(spec.level_col) or "").strip().lower(), 0) if spec.level_col else 0

    return (
        0 if is_ref else 1,                                      # 1. reference genomes first
        0 if has_checkm else 1,                                  # 2. checkm'd before not
        -(comp if has_checkm else 0.0),                          # 3. completeness (higher bette)
        (cont if cont is not None                                #    contamination (lower better)
         else (MISSING_CONTAMINATION if has_checkm else 0.0)),   #    unknown sorts worst)
        -lvl,                                                    # 4. euk fallback, assembly level
        -mean_contig_length(row, spec),                          #    then longest mean contig
        row.get(spec.acc_col) or "",                             # 5. accession tiebreak
    )

PICKERS["quality"] = quality_key

def derep(path, source, wanted_rank, wanted_taxon, derep_rank,
          reps_only=None, pick="quality", screen_against=None):
    """
    One genome per unique value of `derep_rank`, within `wanted_taxon`

    E.g. --wanted-ref-tax bacteria --derep-rank class
         -> one genome for each bacterial class

    reps_only:
        None (default) -> use the SOURCE's default: True for GTDB, False for NCBI.

        This must not be a plain True. GTDB representatives are comprehensive, but
        RefSeq "reference genomes" are sparse

    pick:
        A picker NAME or a CALLABLE (row, spec) -> sort key (lowest wins).

        "quality" (default) -- best-per-group
        "first"             -- first in lineage-sorted order
        callable            -- (row, spec) -> sort key (lowest wins), for a custom
                               selection policy


    screen_against:
        Path to the NCBI Parquet asset. When given, the candidate pool is
        pre-filtered to assemblies that still exist at NCBI BEFORE grouping.

    Returns (accessions, groups_seen, warnings).
    """
    warnings = []

    problem = validate_derep_rank(wanted_rank, derep_rank)
    if problem:
        raise ValueError(problem)

    if rank_index(derep_rank) == rank_index(wanted_rank):
        # Sometimes may be wanted, as in "the single best genome for taxon X" is
        # how you could pick an outgroup. But it is also an easy mis-type for
        # someone who wanted a tree spanning the taxon, and the difference is 1 genome
        # vs hundreds. So mentioning either way
        finer = RANKS[rank_index(derep_rank) + 1] \
            if rank_index(derep_rank) < len(RANKS) - 1 else None
        msg = (f"--derep-rank '{derep_rank}' is the SAME rank as the target taxon, so "
               f"exactly ONE genome will be selected (the single best genome for "
               f"'{wanted_taxon}', which maybe you want for an outgroup or something.")
        if finer:
            msg += (f" But if you wanted a tree spanning '{wanted_taxon}', use a finer "
                    f"--derep-rank (e.g. '{finer}').")
        warnings.append(msg)

    spec = SOURCES[source]
    if reps_only is None:
        reps_only = spec.default_reps_only

    cols = [spec.acc_col, derep_rank, spec.ref_col] + list(spec.quality_cols)
    for extra in (spec.level_col, spec.contig_col, spec.size_col,
                  spec.size_fallback_col):
        if extra and extra not in cols:
            cols.append(extra)

    tab = select(path, source, wanted_rank, wanted_taxon,
                 reps_only=reps_only, columns=cols)

    if screen_against:
        live = live_accession_cores(screen_against)
        before = tab.num_rows
        rows = [r for r in tab.to_pylist()
                if r.get(spec.acc_col) and accession_core(r[spec.acc_col]) in live]
        n_dead = before - len(rows)
        if n_dead:
            warnings.append(
                f"{n_dead:,} candidate genome(s) are no longer available at NCBI "
                f"(suppressed/removed) and were excluded before selection.")
        tab = None
    else:
        rows = tab.to_pylist()

    if not rows:
        warnings.append(f"No genomes found under {wanted_rank} '{wanted_taxon}'"
                        + (" in the representatives pool." if reps_only else "."))
        return [], 0, warnings

    keyfn = resolve_picker(pick)

    best = {}
    n_na_group = 0
    for row in rows:
        group = row.get(derep_rank)
        if not group or group == NA:
            n_na_group += 1
            continue
        k = keyfn(row, spec)
        if group not in best or k < keyfn(best[group], spec):
            best[group] = row

    if n_na_group:
        warnings.append(
            f"{n_na_group:,} genome(s) have no assigned '{derep_rank}' and were "
            f"skipped (unnamed/unclassified at that rank).")

    accs = [best[g][spec.acc_col] for g in sorted(best)]
    return accs, len(best), warnings


class RefGenomeSelection:
    """
    Result of select_ref_genomes(): the kept accessions plus the metadata rows they
    came from, and the resolution details each caller needs to report.

    Attributes:
        accessions   -- list of ncbi_genbank_assembly_accession strings (order stable)
        rows         -- list of dict rows (the selected genomes' metadata) for the
                        accessions above; same order. Suitable for a metadata TSV.
        canonical    -- the taxon name as it appears in the asset (may differ in case
                        from what the user typed)
        resolved_rank-- the rank the taxon was resolved to
        effective_derep_rank -- the rank dereplication collapsed to, or None if derep
                        was off (all genomes under the taxon kept)
        warnings     -- human-facing advisory strings (empty list if none)
    """

    def __init__(self, accessions, rows, canonical, resolved_rank,
                 effective_derep_rank, warnings):
        self.accessions = accessions
        self.rows = rows
        self.canonical = canonical
        self.resolved_rank = resolved_rank
        self.effective_derep_rank = effective_derep_rank
        self.warnings = warnings


def select_ref_genomes(path, source, taxon, target_rank=None, derep_rank="auto",
                       reps_only=None, pick="quality", screen_against=None):
    """
    The one selection entry point shared by every surface that adds reference genomes
    by taxonomy: the standalone get-accessions helpers and the main GToTree driver's
    --wanted-ref-tax path. It runs the full sequence in one place so the surfaces stay
    thin and can't drift:

        1. resolve_taxon()      -- taxon (+ optional target_rank) -> (canonical, rank)
        2. resolve_derep_rank() -- turn "auto"/"off"/<rank> into a concrete rank or None
        3. derep() OR select()  -- dereplicate to one-best-per-rank, or (derep off) take
                                   every genome under the taxon

    This is library-layer code: it RAISES on user-correctable problems
    (TaxonNotFound, AmbiguousTaxon, ValueError for a bad derep rank) and never prints
    or exits. Each CLI surface catches these and translates them to a friendly message,
    per the library-vs-CLI contract.

    Returns a RefGenomeSelection.
    """
    spec = SOURCES[source]

    # 1. resolve the taxon -> (canonical, rank). Raises TaxonNotFound / AmbiguousTaxon.
    canonical, resolved_rank = resolve_taxon(path, taxon, target_rank)

    # 2. work out the effective derep rank. Raises ValueError if the rank is coarser
    #    than the taxon's own rank; may return None (derep off) with a warning.
    effective_derep_rank, warnings = resolve_derep_rank(resolved_rank, derep_rank)

    if reps_only is None:
        reps_only = spec.default_reps_only

    if effective_derep_rank is None:
        # derep off: take every genome under the taxon (optionally NCBI-liveness
        # screened), preserving full metadata rows for the caller's TSV
        cols = _selection_columns(spec, extra_rank=None)
        tab = select(path, source, resolved_rank, canonical,
                     reps_only=reps_only, columns=cols)
        rows = tab.to_pylist()

        if screen_against:
            live = live_accession_cores(screen_against)
            before = len(rows)
            rows = [r for r in rows
                    if r.get(spec.acc_col) and accession_core(r[spec.acc_col]) in live]
            n_dead = before - len(rows)
            if n_dead:
                warnings.append(
                    f"{n_dead:,} candidate genome(s) are no longer available at NCBI "
                    f"(suppressed/removed) and were excluded.")

        accessions = [r.get(spec.acc_col) for r in rows if r.get(spec.acc_col)]
        return RefGenomeSelection(accessions, rows, canonical, resolved_rank,
                                  None, warnings)

    # 3. dereplicate to one best genome per unique value of effective_derep_rank.
    #    derep() returns accessions only, so re-slice to recover metadata rows for the
    #    kept set (keeps a single source of truth for WHICH genomes are picked).
    accessions, _groups, derep_warnings = derep(
        path, source, resolved_rank, canonical, effective_derep_rank,
        reps_only=reps_only, pick=pick, screen_against=screen_against)
    warnings.extend(derep_warnings)

    rows = _rows_for_accessions(path, source, resolved_rank, canonical,
                                effective_derep_rank, reps_only, set(accessions))

    return RefGenomeSelection(accessions, rows, canonical, resolved_rank,
                              effective_derep_rank, warnings)


def _selection_columns(spec, extra_rank=None):
    """The metadata columns worth carrying for a selection's output TSV."""
    cols = [spec.acc_col, spec.ref_col] + list(spec.quality_cols) + list(RANKS)
    for extra in (spec.level_col, spec.contig_col, spec.size_col,
                  spec.size_fallback_col, extra_rank):
        if extra and extra not in cols:
            cols.append(extra)
    return cols


def _rows_for_accessions(path, source, wanted_rank, wanted_taxon, derep_rank,
                         reps_only, wanted_accessions):
    """
    Re-read the taxon slice and return only the rows whose accession is in
    `wanted_accessions` (the derep-kept set), in sorted-accession order.
    """
    spec = SOURCES[source]
    cols = _selection_columns(spec, extra_rank=derep_rank)
    tab = select(path, source, wanted_rank, wanted_taxon,
                 reps_only=reps_only, columns=cols)
    kept = [r for r in tab.to_pylist() if r.get(spec.acc_col) in wanted_accessions]
    kept.sort(key=lambda r: r.get(spec.acc_col) or "")
    return kept
