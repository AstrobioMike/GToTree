import os
import sys
from collections import Counter
from dataclasses import dataclass
from typing import Optional
from gtotree.utils.misc.messaging import wprint, color_text, report_message, report_early_exit
import pandas as pd # type: ignore
from gtotree.utils.misc.general import download_with_tqdm, SCGset, decode_pyhmmer_text
from gtotree.utils.hmms.hmm_searching_engine import profiles_missing_gathering_cutoffs
from gtotree.utils.taxonomy.tax_ranks import RANKS, NA
import pyhmmer #type: ignore


# aliases that aren't simply the file stem of a packaged set
HMM_SET_ALIASES = {"universal": "Universal-Hug-et-al"}

# this is what auto-selection picks when the target reaches outside the Bacteria/Archaea into
# euks, but NOT viruses
UNIVERSAL_SCG_SET = "Universal-Hug-et-al"

# what auto-selection picks when the target taxa span bacteria and archaea
CROSS_DOMAIN_SCG_SET = "Bacteria-and-Archaea"

PREBUILT_DOMAINS = {"bacteria", "archaea"}

# domain names that mean the target is viral, eventually matched as substrings (these get no auto-selection and
# can't use any of the pre-built ones)
VIRAL_DOMAIN_MARKERS = ("virus", "viroid", "viral")

# how much of a linked reference set has to agree on a taxon before that taxon is
# treated as the target group (this is relevant only for --source ncbi -w <taxon> situations where
# NO -H was provided and we are trying to auto-select the best SCG-set to use from the prepackaged gtdb ones)
CONSENSUS_THRESHOLD = 0.9


def with_hmm_suffix(name):
    return name if name.endswith(".hmm") else f"{name}.hmm"


def find_local_hmm_file(hmm_arg):
    """
    Absolute path if `-H` points at a real file on disk, else None
    """
    for candidate in (with_hmm_suffix(hmm_arg), hmm_arg):
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def packaged_scg_set_names():
    """
    Official names of the pre-packaged SCG-sets, as file stems
    """
    try:
        df = read_in_hmm_summary_table()
    except Exception:
        return []
    names = []
    for entry in df["file"]:
        entry = str(entry)
        names.append(entry[:-4] if entry.lower().endswith(".hmm") else entry)
    return names


def canonicalize_hmm_arg(hmm_arg):
    """
    Resolve a requested SCG-set to the official spelling of its name

    e.g., `-H bacteria`, `-H BACTERIA`, and `-H Bacteria.HMM` all -> `Bacteria.hmm`
    Matching case-insensitively is not sufficient on its own since the value is reused as a
    filename under GToTree_HMM_dir, as the lookup key into hmm-sources-and-info.tsv, and
    in run reporting.

    Returned unchanged when it doesn't name a packaged set, so user paths and genuine
    typos both flow on to the existing handling.
    """
    stem = hmm_arg[:-4] if hmm_arg.lower().endswith(".hmm") else hmm_arg
    key = stem.lower()

    if key in HMM_SET_ALIASES:
        return HMM_SET_ALIASES[key]

    for name in packaged_scg_set_names():
        if name.lower() == key:
            return name

    return hmm_arg


def resolve_hmm_arg(args):
    """
    Settle `-H` on the official spelling of the requested set
    """
    if find_local_hmm_file(args.hmm) is None:
        args.hmm = canonicalize_hmm_arg(args.hmm)
    return args


def resolve_hmm_source(args, run_data):
    """
    Point run_data at the HMM file and validate it. Runs on every invocation.

    Cheap on a resume: the packaged set is already in GToTree_HMM_dir, so `get_hmm_path`
    is a stat rather than a download.
    """
    local_path = find_local_hmm_file(args.hmm)

    if local_path is not None:
        run_data.hmm_path = local_path
    else:
        run_data.hmm_path = get_hmm_path(with_hmm_suffix(args.hmm), args.hmm)

    check_gathering_cutoffs(run_data.hmm_path, args.hmm)

    return run_data


def populate_SCG_targets(run_data):
    """
    Build the initial SCG-target list. Fresh runs only, a resume reads its targets,
    with all their accumulated per-SCG state, back out of run-data.json.
    """
    initial_SCG_targets = get_SCG_hmm_targets(run_data.hmm_path)
    run_data.SCG_targets = [SCGset.from_id(target) for target in initial_SCG_targets]

    return run_data


def check_hmm_file(args, run_data):
    """
    The whole fresh-run HMM setup: resolve the name, resolve the file, list the targets
    """
    args = resolve_hmm_arg(args)
    run_data = resolve_hmm_source(args, run_data)

    return populate_SCG_targets(run_data)


################################################################################
# auto-selecting an HMM set from --wanted-ref-tax
################################################################################
#
# I'm adding this for convenience. If `-w` is given without `-H`, we can use the
# wanted target taxon to auto-select the best prebuilt SCG-set we have.
# Only annoying part is the sets are built from GTDB, so if the user also did
# --source ncbi, we have to do some dancing to sort it out.


@dataclass
class AutoPickedSCGSet:
    """
    An auto-selected SCG-set, plus some reason printed out to the user
    """
    name: str
    matched_rank: Optional[str] = None
    matched_taxon: Optional[str] = None
    reason: str = ""


def packaged_sets_by_taxon():
    """
    {(rank, lowercased taxon): official set name} for the packaged sets that sit at one
    of the seven taxonomic ranks

    Sets whose rank is 'universal' or 'multi-domain' are left out on purpose: they don't
    correspond to a node in any lineage, so nothing should reach them by walking one.
    They're only ever chosen explicitly, as a fallback.

    An unreadable summary table yields {} rather than raising -- the caller falls back
    to the universal set, and the table's real problems get reported by the code that
    goes on to look for the HMM file itself.
    """
    try:
        df = read_in_hmm_summary_table()
    except Exception:
        return {}

    index = {}
    for _, row in df.iterrows():
        rank = str(row.get("rank", "")).strip().lower()
        taxon = str(row.get("target_taxa", "")).strip()
        name = str(row.get("file", "")).strip()

        if rank not in RANKS or not taxon or not name:
            continue

        if name.lower().endswith(".hmm"):
            name = name[:-4]
        index[(rank, taxon.lower())] = name

    return index


def pick_set_for_lineage(lineage, resolved_rank, index=None):
    """
    The finest packaged set covering `lineage`, walking upward from `resolved_rank`

    Returns (set_name, matched_rank, matched_taxon), or (None, None, None).
    """
    if index is None:
        index = packaged_sets_by_taxon()

    if resolved_rank not in RANKS:
        return (None, None, None)

    for rank in reversed(RANKS[:RANKS.index(resolved_rank) + 1]):
        taxon = str(lineage.get(rank) or "").strip()
        if not taxon or taxon == NA:
            continue
        match = index.get((rank, taxon.lower()))
        if match:
            return (match, rank, taxon)

    return (None, None, None)


def consensus_lineage(lineages, threshold=CONSENSUS_THRESHOLD):
    """
    Collapse a set of per-genome lineages into the one lineage they share.

    Returns ({rank: taxon}, deepest_agreed_rank); ({}, None) if `lineages` is empty or
    the genomes don't even agree on a domain.

    Walks coarse to fine and stops at the first rank where agreement falls below
    `threshold`, dropping everything below it
    """
    if not lineages:
        return {}, None

    total = len(lineages)
    agreed = {}
    deepest = None

    for rank in RANKS:
        counts = Counter(lineage.get(rank) for lineage in lineages
                         if lineage.get(rank) and lineage.get(rank) != NA)
        if not counts:
            break

        taxon, count = counts.most_common(1)[0]
        if count / total < threshold:
            break

        agreed[rank] = taxon
        deepest = rank

    return agreed, deepest


def domains_outside_prebuilt_scope(lineages):
    """
    The domain names in `lineages` that no pre-built set covers, lowercased.

    On the GTDB path this is always empty: GTDB only classifies Bacteria and Archaea.
    It earns its keep on the NCBI path, where the selection's own rows carry NCBI's
    superkingdom, so e.g., `-w Ascomycota --source ncbi` reports 'eukaryota' here
    """
    found = set()
    for lineage in lineages:
        domain = str(lineage.get("domain") or "").strip()
        if not domain or domain == NA:
            continue
        if domain.lower() not in PREBUILT_DOMAINS:
            found.add(domain.lower())
    return found


def intersect_lineages(resolved):
    """
    Collapse several INDEPENDENTLY resolved taxa into the lineage they share

    `resolved` is [(lineage, deepest_rank), ...], one entry per `-w` taxon. Returns
    (agreed_lineage, deepest_agreed_rank), or ({}, None) when they don't even share a
    domain.

    This is deliberately NOT consensus_lineage(). That tolerates a minority
    disagreeing, which is right within a single provided taxon. But when multiple are
    provided to -w, we want agreement to be unanimous (e.g., pooling `-w Escherichia
    -w Nanoarchaeum` and taking a 90% vote would let the thousands of Escherichia genomes
    outvote the archaea entirely and pick a bacterial set to score them with. Now, with this
    approach, that would select the bacteria-and-archaea set

    """
    if not resolved:
        return {}, None

    limit = min(RANKS.index(rank) for _lineage, rank in resolved
                if rank in RANKS) if all(rank in RANKS for _l, rank in resolved) else -1
    if limit < 0:
        return {}, None

    agreed = {}
    deepest = None

    for rank in RANKS[:limit + 1]:
        values = {str(lineage.get(rank) or "").strip() for lineage, _rank in resolved}
        if len(values) != 1:
            break
        value = values.pop()
        if not value or value == NA:
            break
        agreed[rank] = value
        deepest = rank

    return agreed, deepest


def autopick_scg_set(source, selections):
    """
    Choose an SCG-set for `-w` input that ALSO has no `-H` specified

    `selections` is one RefGenomeSelection or a list of them

    GTDB source: the selection's rows already carry the GTDB lineage, so the walk needs
    no further lookup, and the answer is exact

    NCBI source: since NCBI's names don't line up with GTDB's, the selected accessions
    are linked back to GTDB, where possible, and their lineages collapsed to a consensus.

    Reaching outside Bacteria/Archaea is checked FIRST and takes the universal set if
    it involved eukarya, and blocks viruses from being auto-selected at all. Spanning
    both prokaryote domains but going no further takes 'Bacteria-and-Archaea'

    Returns an AutoPickedSCGSet. Never raises: no answer is a fallback, not an error.
    """
    selections = _as_selection_list(selections)

    if not selections:
        return AutoPickedSCGSet(UNIVERSAL_SCG_SET,
                                reason="there was no taxon to select one from")

    label = _describe_taxa(selections)
    via_genomes = str(source).strip().lower() != "gtdb"
    subject = f"the genomes selected for {label}" if via_genomes else label
    plural = via_genomes or len(selections) > 1

    outside = set()
    for selection in selections:
        outside |= domains_outside_prebuilt_scope(selection.rows or [])
    if outside:
        return AutoPickedSCGSet(
            UNIVERSAL_SCG_SET,
            reason=(f""))

    if via_genomes:
        resolved, unplaced = _resolve_ncbi_selections(selections)
    else:
        resolved, unplaced = _resolve_gtdb_selections(selections)

    if unplaced:
        return AutoPickedSCGSet(
            UNIVERSAL_SCG_SET,
            reason=(f"{_join_names(unplaced)} {'has' if len(unplaced) == 1 else 'have'} "
                    "no counterpart in GTDB, which the pre-built sets are built from"))

    lineage, deepest = intersect_lineages(resolved)

    if deepest is None:
        # they agree on nothing, so they span both domains -- which is precisely what
        # the cross-domain set is for. Eukaryotes were already ruled out above
        return AutoPickedSCGSet(
            CROSS_DOMAIN_SCG_SET,
            reason=(f"{subject} {'span' if plural else 'spans'} both domains "
                    "(Bacteria and Archaea)"))

    name, rank, taxon = pick_set_for_lineage(lineage, deepest)

    if name is None:
        return AutoPickedSCGSet(
            UNIVERSAL_SCG_SET,
            reason=f"no pre-built set was found covering {label}")

    return AutoPickedSCGSet(name, rank, taxon,
                            _pick_reason(selections, rank, taxon, via_genomes))


################################################################################
# viral targets: bring your own HMMs
################################################################################


class ViralTaxonNeedsOwnHMMs(Exception):
    """
    `-w` reached into viruses, and `-H` isn't pointing at the user's own HMM file.
    """


def _selection_is_viral(selection):
    """
    Whether a selection's domain places it in viruses
    """
    for lineage in (getattr(selection, "rows", None) or []):
        domain = str(lineage.get("domain") or "").strip()
        if not domain or domain == NA:
            continue
        if any(marker in domain.lower() for marker in VIRAL_DOMAIN_MARKERS):
            return True
    return False


def viral_selections(selections):
    """
    The subset of `selections` that resolved into viruses (empty list if none did).
    """
    return [s for s in _as_selection_list(selections) if _selection_is_viral(s)]


def check_viruses_have_their_own_hmms(args, selections):
    """
    Viral `-w` targets require the user's own `-H` file
    """
    viral = viral_selections(selections)
    if not viral:
        return

    if args.hmm and find_local_hmm_file(args.hmm) is not None:
        return

    label = _join_names([getattr(s, "canonical", str(s)) for s in viral])
    target, is_are = ("target", "is") if len(viral) == 1 else ("targets", "are")

    message = (f"The `-w` {target} {label} {is_are} viral, and none of GToTree's "
               "pre-built single-copy gene HMM sets are built for viral genomes.")

    if args.hmm:
        message += (f" '{args.hmm}' is one of those pre-built sets, so it can't be "
                    "used for this either.")

    message += (" To build a viral tree, point `-H` at your own HMM file of "
                "single-copy genes suitable for the group (e.g., "
                "`-H my-viral-SCGs.hmm`).")

    raise ViralTaxonNeedsOwnHMMs(message)


def _as_selection_list(selections):
    """One selection, a list of them, or None -- always out as a list."""
    if selections is None:
        return []
    if isinstance(selections, (list, tuple)):
        return [s for s in selections if s is not None]
    return [selections]


def _join_names(names):
    names = list(names)
    if len(names) == 1:
        return f"'{names[0]}'"
    return ", ".join(f"'{n}'" for n in names[:-1]) + f" and '{names[-1]}'"


def _describe_taxa(selections):
    return _join_names([s.canonical for s in selections])


def _pick_reason(selections, rank, taxon, via_genomes):
    """
    Why a set was chosen, phrased for one taxon or several

    The NCBI phrasing says the genomes got there, not the name: NCBI's names don't line
    up with GTDB's, so 'Nitrososphaerota' landing on the Thermoproteota set only makes
    sense to a user if it's clear the match ran through the genomes' GTDB placement.

    An empty string means "no reason worth printing", like when the set carries the taxon's
    own name
    """
    label = _describe_taxa(selections)

    if via_genomes:
        return f"the genomes selected for {label} fall in {rank} {taxon} in GTDB"

    if len(selections) == 1:
        selection = selections[0]
        if (rank == selection.resolved_rank
                and str(taxon).lower() == str(selection.canonical).lower()):
            return ""
        return f"'{selection.canonical}' sits within {rank} {taxon} in GTDB"

    return f"{label} sit within {rank} {taxon} in GTDB"


def _resolve_gtdb_selections(selections):
    """
    ([(lineage, deepest_rank), ...], [taxa with no usable lineage]) for GTDB input.

    The rows already carry the GTDB lineage, and every genome under the taxon shares it
    down to the rank the taxon resolved to, so the first row speaks for the selection.
    """
    resolved, unplaced = [], []

    for selection in selections:
        lineage = selection.rows[0] if selection.rows else {}
        if not lineage or selection.resolved_rank not in RANKS:
            unplaced.append(selection.canonical)
            continue
        resolved.append((lineage, selection.resolved_rank))

    return resolved, unplaced


def _resolve_ncbi_selections(selections):
    """
    ([(lineage, deepest_rank), ...], [taxa GTDB has never heard of]) for NCBI input.

    Each taxon is collapsed on its own, so a taxon's internal noise stays internal and
    can't leak into the across-taxa agreement.
    """
    from gtotree.utils.gtdb.handle_gtdb_tax_info import gtdb_lineages_for_accessions

    resolved, unplaced = [], []

    for selection in selections:
        lineages = list(gtdb_lineages_for_accessions(selection.accessions).values())
        if not lineages:
            unplaced.append(selection.canonical)
            continue

        lineage, deepest = consensus_lineage(lineages)
        if deepest is None:
            # this one taxon's own genomes straddle the domains: it constrains nothing
            # below domain, which intersect_lineages() reads as a cross-domain span
            resolved.append(({}, RANKS[0]))
            continue

        resolved.append((lineage, deepest))

    return resolved, unplaced


def check_gathering_cutoffs(hmm_path, hmm_arg):
    """
    Fail early if any profile lacks a gathering (GA) threshold

    Searching uses gathering thresholds (the `--cut_ga` equivalent), so a profile
    without one can't be searched
    """
    try:
        missing = profiles_missing_gathering_cutoffs(hmm_path)
    except Exception:
        # unreadable HMM is reported elsewhere; don't mask it with a cutoffs message
        return

    if not missing:
        return

    shown = ", ".join(missing[:5])
    if len(missing) > 5:
        shown += f", ... ({len(missing)} total)"

    report_message(
        f'The HMM file for "{hmm_arg}" has {len(missing)} profile(s) with no gathering '
        f"(GA) threshold: {shown}.", "red")
    report_message(
        "GToTree identifies target genes using each profile's gathering threshold, so "
        "every profile needs a `GA` line. If you built this HMM yourself, you can add "
        "gathering thresholds with `hmmbuild --cut_ga` inputs or by editing the `GA` "
        "lines in directly in the HMM file. You may need to do some testing to figure out "
        "what those cutoffs should be!")
    report_early_exit(None, copy_log=False)


def get_hmm_path(hmm_file, hmm_arg):
    hmm_data_dir = check_hmm_location_var_is_set()
    dest_path = os.path.join(hmm_data_dir, hmm_file)
    if os.path.isfile(dest_path):
        return dest_path
    download_prepackaged_hmm(dest_path, hmm_arg)
    return dest_path


def check_hmm_location_var_is_set():
    try:
        hmm_data_dir = os.environ['GToTree_HMM_dir']
    except KeyError:
        wprint(color_text("The environment variable 'GToTree_HMM_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt data locations check`.")
        sys.exit(1)
    return hmm_data_dir


def download_prepackaged_hmm(dest_path, hmm_arg):
    hmm_file = os.path.basename(dest_path)
    target_hmm_url = get_target_hmm_url(hmm_file, hmm_arg)

    print(color_text(f"    Downloading the prebuilt \"{hmm_arg}\" HMM set (only needs to be done once)...\n", "yellow"))

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)

    try:
        download_with_tqdm(target_hmm_url, f"        {hmm_arg} HMM file", dest_path)
    except Exception as e:
        report_message(f"Downloading the HMM file failed with the following error:\n{e}", "red")
        report_early_exit(None, copy_log=False)


def read_in_hmm_summary_table():
    hmm_data_dir = check_hmm_location_var_is_set()
    hmm_table_path = os.path.join(hmm_data_dir, "hmm-sources-and-info.tsv")
    df = pd.read_csv(hmm_table_path, sep='\t')
    return(df)


def get_target_hmm_url(hmm_file, hmm_arg):
    df = read_in_hmm_summary_table()
    match = df.loc[df["file"].astype(str).str.lower() == hmm_file.lower(), "link"]
    if len(match) == 0:
        report_message(f"You specified \"{hmm_arg}\" as the HMM file to use, but that file can't be found.", "red")
        report_message("You can see the available gene-sets packaged with GToTree by running `gtt hmms`.")
        report_early_exit(None, copy_log = False)
    return match.values[0]


def get_SCG_hmm_targets(hmm_path):
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        hmms = list(hmm_file)
    SCG_targets = [decode_pyhmmer_text(hmm.name) for hmm in hmms]
    return SCG_targets
