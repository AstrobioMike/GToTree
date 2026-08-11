#!/usr/bin/env python3
"""
Plan which SCG-HMM sets to pre-build for GToTree from a GTDB release.

This is a maintainer tool, not part of the GToTree CLI. It's expected to be run rarely,
at most once a year, so it favours being readable and auditable over
being fast.

What it does
------------
Reads `scg-sets.toml` and the GTDB parquet, works out which taxa should get a set and
how many genomes each would be built from, and writes out the concrete
`gtt gen-scg-hmms` commands plus the tables backing those decisions.

What it deliberately does NOT do
--------------------------------
It doesn't pick genomes. `gen-scg-hmms -w <taxon> --derep-rank <rank>` already resolves
a taxon to a genome set -> representatives only, best-per-group by quality, NCBI
liveness screened, deterministic. Emitting accession lists here would fork that logic
and guarantee the two drift apart. So the plan names TARGETS, and the build resolves
them. The per-genome record you'd want afterwards is already a `gen-scg-hmms` output
(`target-genomes.tsv`).

The counts here are therefore predictions of what the build will select, computed with
the same rules. They're checked against the real selector in the test suite.

Usage
-----
    python plan_scg_sets.py                      # uses ./scg-sets.toml, writes ./plan/
    python plan_scg_sets.py -c other.toml -o out/
    python plan_scg_sets.py --gtdb-parquet /path/to/gtdb-data.parquet

Requires the GTDB parquet. It's found via, in order: --gtdb-parquet, the GTDB_DIR
environment variable, or GToTree's own downloader (which will fetch it if missing).
"""

import argparse
import csv
import os
import sys
import tomllib # type: ignore
from collections import defaultdict
from datetime import date

import pyarrow.parquet as pq  # type: ignore


RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# GTDB stores everything as strings, including the numerics, and uses these for absent
# values. Matching gtotree.utils.taxonomy.tax_ranks.NA handling.
EMPTY = {None, "", "NA", "na", "none", "None"}

PARQUET_FILENAME = "gtdb-data.parquet"

# `rank` values in the published table that mean "not built from GTDB by this
# toolchain". Currently just the Hug et al. universal set, which is carried forward
# from a published external asset (see publish-scg-sets.py's carried-forward table).
EXTERNAL_RANKS = {"universal"}


################################################################################
# locating the asset
################################################################################

def find_gtdb_parquet(explicit=None):
    """
    Locate the GTDB parquet, preferring an explicit path, then GTDB_DIR, then asking
    GToTree to resolve (and download) it.

    The GToTree fallback is last because it can trigger a multi-minute download, and
    someone running this repeatedly while tuning thresholds shouldn't have to wonder
    whether it's about to.
    """
    if explicit:
        if not os.path.isfile(explicit):
            sys.exit(f"No GTDB parquet at: {explicit}")
        return explicit

    env_dir = os.environ.get("GTDB_DIR")
    if env_dir:
        candidate = os.path.join(env_dir, PARQUET_FILENAME)
        if os.path.isfile(candidate):
            return candidate

    try:
        from gtotree.utils.gtdb.get_gtdb_data import gtdb_data_table_path
        return gtdb_data_table_path()
    except Exception as e:  # noqa: BLE001 - this is a maintainer script, say what broke
        sys.exit(
            f"Couldn't locate the GTDB parquet ({e}).\n"
            "Pass --gtdb-parquet, set GTDB_DIR, or run this inside a GToTree "
            "environment where `gtt data get-gtdb` has been run.")


def gtdb_version(parquet_path):
    """Best-effort release label from the VERSION.txt sidecar next to the parquet."""
    version_file = os.path.join(os.path.dirname(parquet_path), "VERSION.txt")
    try:
        with open(version_file) as f:
            return f.read().strip().splitlines()[0]
    except (OSError, IndexError):
        return "unknown"


################################################################################
# counting
################################################################################

class TaxonStats:
    """Everything the plan needs to know about one candidate taxon."""

    def __init__(self, rank, name, domain, parent="NA"):
        self.rank = rank
        self.name = name
        self.domain = domain
        # the next-coarser taxon (class -> its phylum, phylum -> its domain). Carried
        # through to the published info table so `gtt hmms` can nest a class set under
        # the phylum set it sits inside without re-reading GTDB at display time.
        self.parent = parent
        self.genomes = 0            # all genomes, no filtering
        self.reps = 0               # GTDB species representatives
        self.hq_reps = 0            # representatives clearing the quality floor
        self.groups = defaultdict(set)   # finer rank -> distinct values, floored pool

    def n(self, rank):
        """Distinct values at `rank` within the floored representatives pool. This is
        what a `--derep-rank <rank>` build would produce, one genome per group."""
        return len(self.groups[rank])


def _num(value):
    if value in EMPTY:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def gather_stats(parquet_path, min_completeness, max_contamination, ranks=("phylum", "class")):
    """
    One pass over the parquet, accumulating stats for every taxon at each rank in
    `ranks`.

    The floored pool mirrors what `gen-scg-hmms` selects from: GTDB representatives
    only (the source default), then the quality floor. Groups with no name at a rank
    are skipped, matching derep()'s handling of unclassified values, otherwise the
    plan would promise genomes the build then refuses to group.
    """
    columns = ["gtdb_representative", "checkm2_completeness", "checkm2_contamination"] + RANKS
    table = pq.read_table(parquet_path, columns=columns)
    cols = {c: table.column(c).to_pylist() for c in columns}
    n_rows = len(cols["domain"])

    stats = {}
    for i in range(n_rows):
        domain = cols["domain"][i]
        is_rep = cols["gtdb_representative"][i] == "t"

        completeness = _num(cols["checkm2_completeness"][i])
        contamination = _num(cols["checkm2_contamination"][i])
        clears_floor = is_rep
        if clears_floor and min_completeness is not None:
            clears_floor = completeness is not None and completeness >= min_completeness
        if clears_floor and max_contamination is not None:
            clears_floor = contamination is not None and contamination <= max_contamination

        for rank in ranks:
            name = cols[rank][i]
            if name in EMPTY:
                continue
            key = (rank, name)
            if key not in stats:
                rank_index = RANKS.index(rank)
                parent = cols[RANKS[rank_index - 1]][i] if rank_index > 0 else "NA"
                stats[key] = TaxonStats(rank, name, domain,
                                        "NA" if parent in EMPTY else parent)
            entry = stats[key]
            entry.genomes += 1
            if is_rep:
                entry.reps += 1
            if not clears_floor:
                continue
            entry.hq_reps += 1
            for finer in RANKS[RANKS.index(rank) + 1:]:
                value = cols[finer][i]
                if value not in EMPTY:
                    entry.groups[finer].add(value)

    return stats


################################################################################
# planning
################################################################################

class PlannedSet:
    def __init__(self, rank, name, stats, reason):
        self.rank = rank
        self.name = name
        self.stats = stats
        self.reason = reason          # why it's in the plan
        self.derep_rank = None
        self.n_input = 0
        self.flags = []
        self.taxa = [name]            # what -w values the build uses; 1+ for domain sets
        self.percent_single_copy = None   # None -> build script omits -p, inherits default


def plan_domain_sets(config, parquet_path, quality, default_psc):
    """
    Plan the domain-spanning sets from the [[domain_sets]] config blocks.

    These are handled apart from the phylum machinery because they're a different shape:
    the target is one or more whole domains, the derep rank is fixed by config rather
    than stepped to hit an input floor, and a combined set's input is the sum across its
    domains. Counting is a direct query per (domain, derep_rank) rather than reusing the
    phylum stats, which are keyed on phylum.

    Returns a list of PlannedSet.
    """
    entries = config.get("domain_sets", [])
    if not entries:
        return []

    # one pass, counting distinct values of every rank per domain in the floored pool
    columns = (["gtdb_representative", "checkm2_completeness", "checkm2_contamination"]
               + RANKS)
    table = pq.read_table(parquet_path, columns=columns)
    cols = {c: table.column(c).to_pylist() for c in columns}
    n_rows = len(cols["domain"])

    min_completeness = quality.get("min_completeness")
    max_contamination = quality.get("max_contamination")

    # domain -> rank -> set of distinct values (floored, reps-only)
    per_domain = defaultdict(lambda: defaultdict(set))
    for i in range(n_rows):
        if cols["gtdb_representative"][i] != "t":
            continue
        completeness = _num(cols["checkm2_completeness"][i])
        contamination = _num(cols["checkm2_contamination"][i])
        if min_completeness is not None and (completeness is None
                                             or completeness < min_completeness):
            continue
        if max_contamination is not None and (contamination is None
                                              or contamination > max_contamination):
            continue
        domain = cols["domain"][i]
        for rank in RANKS:
            value = cols[rank][i]
            if value not in EMPTY:
                per_domain[domain][rank].add(value)

    planned = []
    for entry in entries:
        taxa = entry.get("taxa") or [entry["name"]]
        derep_rank = entry.get("derep_rank", "family")
        # input size = distinct derep-rank groups summed over the named domains, since
        # each domain is resolved and dereplicated independently before merging
        n_input = sum(len(per_domain.get(taxon, {}).get(derep_rank, set()))
                      for taxon in taxa)

        planned_set = PlannedSet(
            "domain" if len(taxa) == 1 else "multi-domain",
            entry["name"], None, "domain-set")
        planned_set.taxa = list(taxa)
        planned_set.derep_rank = derep_rank
        planned_set.n_input = n_input
        planned_set.percent_single_copy = entry.get("percent_single_copy")
        planned.append(planned_set)

    return planned


def choose_derep_rank(stats, preferred, min_input, wants_review_above):
    """
    Pick the dereplication rank for one target.

    Start at the preferred rank (genus) and step FINER while the input would be too
    small. Genus is the right default, even breadth, no clinical-isolate pileups,
    but for a narrow clade it can leave 28 genomes, and a 90% single-copy call across
    28 genomes turns on two or three of them. Stepping to species recovers genomes from
    the same floored pool without reintroducing redundancy, since GTDB species clusters
    are already dereplicated.

    Deliberately never steps COARSER on its own. Going genus -> family to trim a large
    input throws away real breadth, so that gets flagged for a human instead.

    Returns (rank, n_genomes, flags).
    """
    flags = []
    candidates = RANKS[RANKS.index(stats.rank) + 1:]
    if preferred not in candidates:
        # e.g. a class target with preferred_rank coarser than class; fall back to the
        # coarsest rank actually finer than the target
        preferred = candidates[0]

    start = candidates.index(preferred)
    chosen, count = preferred, stats.n(preferred)

    for rank in candidates[start:]:
        chosen, count = rank, stats.n(rank)
        if count >= min_input:
            break

    if chosen != preferred:
        flags.append(f"stepped {preferred} -> {chosen} for enough genomes")
    if count < min_input:
        flags.append(f"only {count:,} genomes even at {chosen}; single-copy call "
                     f"will be weak")
    if wants_review_above and count > wants_review_above:
        flags.append(f"{count:,} genomes is large; consider a coarser derep by hand")

    return chosen, count, flags


def build_plan(config, stats):
    """Apply [criteria], [[include]] and [[exclude]] to produce the ordered plan."""
    criteria = config.get("criteria", {})
    auto_rank = criteria.get("rank", "phylum")
    min_genera = criteria.get("min_genera")
    min_genomes = criteria.get("min_genomes")

    excluded = {(e["rank"], e["name"]) for e in config.get("exclude", [])}
    included = {(i["rank"], i["name"]) for i in config.get("include", [])}

    planned = {}

    # automatic, at one rank only
    for (rank, name), entry in stats.items():
        if rank != auto_rank or (rank, name) in excluded:
            continue
        by_genera = min_genera is not None and entry.n("genus") >= min_genera
        by_genomes = min_genomes is not None and entry.genomes >= min_genomes
        if not (by_genera or by_genomes):
            continue
        reasons = []
        if by_genera:
            reasons.append(f"genera>={min_genera}")
        if by_genomes:
            reasons.append(f"genomes>={min_genomes}")
        planned[(rank, name)] = PlannedSet(rank, name, entry, "+".join(reasons))

    # explicit includes, which may also re-flag something already auto-selected
    missing = []
    for (rank, name) in included:
        if (rank, name) in excluded:
            continue
        entry = stats.get((rank, name))
        if entry is None:
            missing.append(f"{name} ({rank})")
            continue
        if (rank, name) in planned:
            existing = planned[(rank, name)]
            existing.reason += "+include"
        else:
            planned[(rank, name)] = PlannedSet(rank, name, entry, "include")

    derep = config.get("derep", {})
    for entry in planned.values():
        entry.derep_rank, entry.n_input, entry.flags = choose_derep_rank(
            entry.stats,
            derep.get("preferred_rank", "genus"),
            derep.get("min_input_genomes", 50),
            derep.get("above_this_wants_review"))

    ordered = sorted(planned.values(),
                     key=lambda p: (p.stats.domain, RANKS.index(p.rank), -p.n_input))
    return ordered, missing


################################################################################
# comparing against what's published
################################################################################

def read_published(path):
    """
    (taxon_lower, rank) -> row, from the hosted-sets info table.

    Rows at an EXTERNAL_RANK are skipped: they're carried-forward external assets that
    this planner never plans, so including them would report them as "published but not
    in this plan" on every single run, forever.
    """
    if not os.path.isfile(path):
        return {}
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            taxa = (row.get("target_taxa") or "").strip().lower()
            rank = (row.get("rank") or "NA").strip()
            if taxa and rank not in EXTERNAL_RANKS:
                out[(taxa, rank)] = row
    return out


def diff_against_published(plan, published):
    """
    What changed relative to the hosted sets.

    Advisory only. Older releases stay hosted, and whether to rebuild at all is a
    separate decision.
    """
    planned_keys = {(p.name.lower(), p.rank) for p in plan}
    published_keys = set(published)

    new = [p for p in plan if (p.name.lower(), p.rank) not in published_keys]
    gone = [published[k] for k in sorted(published_keys - planned_keys)]
    return new, gone


################################################################################
# outputs
################################################################################

def write_plan_table(path, plan):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["set_name", "rank", "taxon", "domain", "parent", "selected_by",
                    "derep_rank", "n_input_genomes", "total_genomes",
                    "species_clusters", "hq_reps", "n_genera_hq", "n_species_hq",
                    "flags"])
        for p in plan:
            if p.stats is None:
                # domain-spanning set: no single-taxon stats row behind it. A
                # single-domain set still names its domain, so the published table can
                # file it under that domain; a set spanning several has no one domain.
                domain = p.taxa[0] if len(p.taxa) == 1 else "NA"
                w.writerow([p.name, p.rank, "+".join(p.taxa), domain, "NA", p.reason,
                            p.derep_rank, p.n_input, "NA", "NA", "NA", "NA", "NA",
                            "; ".join(p.flags)])
            else:
                w.writerow([p.name, p.rank, p.name, p.stats.domain, p.stats.parent,
                            p.reason, p.derep_rank, p.n_input, p.stats.genomes,
                            p.stats.reps, p.stats.hq_reps, p.stats.n("genus"),
                            p.stats.n("species"), "; ".join(p.flags)])


def write_counts_table(path, stats):
    """
    Every taxon at the ranks examined, selected or not.

    The unselected rows are the point: next release you want to see what's approaching
    a threshold, which you can't do from the plan alone.
    """
    rows = sorted(stats.values(),
                  key=lambda s: (s.domain, RANKS.index(s.rank), -s.n("genus")))
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["domain", "rank", "taxon", "total_genomes", "species_clusters",
                    "hq_reps", "n_class_hq", "n_order_hq", "n_family_hq",
                    "n_genus_hq", "n_species_hq"])
        for s in rows:
            w.writerow([s.domain, s.rank, s.name, s.genomes, s.reps, s.hq_reps,
                        s.n("class"), s.n("order"), s.n("family"),
                        s.n("genus"), s.n("species")])


def write_build_script(path, plan, quality, out_root, num_jobs=20, num_threads=16):
    """
    The build stage, as a plain shell script rather than a second program.

    Each set is one independent `gen-scg-hmms` invocation, so a loop is genuinely all
    this is, and a script you can read, edit, and re-run a single line of beats a
    wrapper when a build fails partway.
    """
    floor = []
    if quality.get("min_completeness") is not None:
        floor.append(f"--min-completeness {quality['min_completeness']:g}")
    if quality.get("max_contamination") is not None:
        floor.append(f"--max-contamination {quality['max_contamination']:g}")
    floor_str = (" " + " ".join(floor)) if floor else ""

    lines = [
        "#!/usr/bin/env bash",
        "# Generated by plan_scg_sets.py -> build the planned SCG-HMM sets.",
        "#",
        "# Each set is independent, so a failure only costs that one set. `--resume`",
        "# is safe to add if a run is interrupted; it refuses if parameters changed.",
        "#",
        "# The -o basename becomes the .hmm filename, so these land as <Taxon>.hmm.",
        "set -euo pipefail",
        "",
        f'OUT_ROOT="${{1:-{out_root}}}"',
        'mkdir -p "$OUT_ROOT"',
        "",
    ]
    for p in plan:
        lines.append(f"# {p.name} ({p.rank}) -- {p.n_input:,} genomes, "
                     f"selected by {p.reason}")
        # one -w per taxon (a domain set names two); --target-rank only makes sense for
        # a single-taxon set, so it's omitted when pooling several
        w_flags = " ".join(f'-w "{taxon}"' for taxon in p.taxa)
        target_rank = f" --target-rank {p.rank}" if len(p.taxa) == 1 and p.stats else ""
        psc = (f" -p {p.percent_single_copy}"
               if getattr(p, "percent_single_copy", None) is not None else "")
        lines.append(
            f'gtt gen-scg-hmms {w_flags}{target_rank} '
            f'--derep-rank {p.derep_rank}{floor_str}{psc} '
            f'-j {num_jobs} -t {num_threads} '
            f'-o "$OUT_ROOT/{p.name}"')
        lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))
    os.chmod(path, 0o755)


def write_report(path, plan, stats, new, gone, missing, meta):
    """The human-facing review document, the thing to read before building."""
    L = []
    add = L.append

    add("SCG-HMM set plan")
    add("=" * 76)
    add(f"generated:     {meta['date']}")
    add(f"GTDB release:  {meta['gtdb_version']}")
    add(f"parquet:       {meta['parquet']}")
    add(f"config:        {meta['config']}")
    add(f"quality floor: completeness >= {meta['min_completeness']}, "
        f"contamination <= {meta['max_contamination']}")
    add("")
    add(f"{len(plan)} sets planned, {sum(p.n_input for p in plan):,} genome selections "
        f"in total.")
    add("")

    add("PLANNED SETS")
    add("-" * 76)
    add(f"{'taxon':<26}{'rank':<8}{'derep':<9}{'genomes':>9}  selected by")
    for p in plan:
        add(f"{p.name[:25]:<26}{p.rank:<8}{p.derep_rank:<9}{p.n_input:>9,}  {p.reason}")
    add("")

    flagged = [p for p in plan if p.flags]
    if flagged:
        add("NEEDS A LOOK")
        add("-" * 76)
        for p in flagged:
            add(f"  {p.name} ({p.rank}):")
            for flag in p.flags:
                add(f"      - {flag}")
        add("")

    add("NOT IN THE PUBLISHED TABLE YET")
    add("-" * 76)
    if new:
        for p in new:
            add(f"  {p.name} ({p.rank}) -- {p.n_input:,} genomes")
    else:
        add("  (none)")
    add("")

    add("PUBLISHED BUT NOT IN THIS PLAN")
    add("-" * 76)
    add("  Advisory only. Older releases stay hosted anyway.")
    if gone:
        for row in gone:
            add(f"  {row.get('file', '?')} (target_taxa={row.get('target_taxa')}, "
                f"rank={row.get('rank')})")
    else:
        add("  (none)")
    add("")

    if missing:
        add("CONFIGURED BUT NOT FOUND IN GTDB")
        add("-" * 76)
        add("  These are named in [[include]] but don't exist at that rank in this")
        add("  release.")
        for name in missing:
            add(f"  {name}")
        add("")

    add("APPROACHING THE THRESHOLDS")
    add("-" * 76)
    add("  Unselected phyla with the most genera, for next release's review.")
    planned_names = {(p.rank, p.name) for p in plan}
    near = sorted((s for (rank, name), s in stats.items()
                   if rank == "phylum" and (rank, name) not in planned_names),
                  key=lambda s: -s.n("genus"))[:15]
    for s in near:
        add(f"  {s.name[:25]:<26}{s.domain[:3]:<5}genera={s.n('genus'):>5,}  "
            f"genomes={s.genomes:>8,}")

    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


################################################################################
# driver
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Plan which SCG-HMM sets to pre-build for GToTree from GTDB.")
    here = os.path.dirname(os.path.abspath(__file__))
    parser.add_argument("-j", "--num-jobs", type=int, default=20,
                        help="number of jobs passed to `gtt gen-scg-hmms` (default: 20)")
    parser.add_argument("-t", "--num-threads", type=int, default=16,
                        help="number of threads passed to `gtt gen-scg-hmms` (default: 16)")
    parser.add_argument("-c", "--config", default=os.path.join(here, "scg-sets.toml"),
                        help="TOML config (default: scg-sets.toml beside this script)")
    parser.add_argument("-o", "--output-dir", default=os.path.join(here, "plan"),
                        help="where to write the plan (default: ./plan)")
    parser.add_argument("--gtdb-parquet", default=None,
                        help="path to gtdb-data.parquet (default: GTDB_DIR, then "
                             "GToTree's own resolver)")
    parser.add_argument("--repo-root", default=None,
                        help="repo root, for resolving the published info table "
                             "(default: two levels above this script)")
    args = parser.parse_args()

    with open(args.config, "rb") as f:
        config = tomllib.load(f)

    quality = config.get("quality", {})
    min_completeness = quality.get("min_completeness")
    max_contamination = quality.get("max_contamination")
    num_jobs = args.num_jobs
    num_threads = args.num_threads

    parquet = find_gtdb_parquet(args.gtdb_parquet)
    version = gtdb_version(parquet)
    print(f"GTDB parquet: {parquet}")
    print(f"GTDB release: {version}")
    print(f"quality floor: completeness >= {min_completeness}, "
          f"contamination <= {max_contamination}")

    auto_rank = config.get("criteria", {}).get("rank", "phylum")
    ranks_needed = sorted({auto_rank} | {i["rank"] for i in config.get("include", [])},
                          key=RANKS.index)
    print(f"counting at: {', '.join(ranks_needed)} ...")
    stats = gather_stats(parquet, min_completeness, max_contamination,
                         ranks=tuple(ranks_needed))
    print(f"  {len(stats):,} taxa counted")

    plan, missing = build_plan(config, stats)

    # domain-spanning sets are planned separately (different shape, own count) and
    # appended, ordered first so they head the plan, as the broadest sets
    default_psc = config.get("percent_single_copy")
    domain_sets = plan_domain_sets(config, parquet, quality, default_psc)
    plan = domain_sets + plan

    repo_root = args.repo_root or os.path.dirname(os.path.dirname(here))
    published_path = os.path.join(
        repo_root, config.get("published", {}).get(
            "info_table", "hmm_sets/hmm-sources-and-info.tsv"))
    published = read_published(published_path)
    new, gone = diff_against_published(plan, published)

    os.makedirs(args.output_dir, exist_ok=True)
    plan_tsv = os.path.join(args.output_dir, "scg-set-plan.tsv")
    counts_tsv = os.path.join(args.output_dir, "gtdb-taxon-counts.tsv")
    build_sh = os.path.join(args.output_dir, "build-scg-sets.sh")
    report_txt = os.path.join(args.output_dir, "plan-report.txt")

    write_plan_table(plan_tsv, plan)
    write_counts_table(counts_tsv, stats)
    write_build_script(build_sh, plan, quality, "scg-set-builds", num_jobs=num_jobs, num_threads=num_threads)
    write_report(report_txt, plan, stats, new, gone, missing, {
        "date": date.today().isoformat(),
        "gtdb_version": version,
        "parquet": parquet,
        "config": os.path.abspath(args.config),
        "min_completeness": min_completeness,
        "max_contamination": max_contamination,
    })

    print(f"\n{len(plan)} sets planned "
          f"({sum(p.n_input for p in plan):,} genome selections total)")
    if missing:
        print(f"  {len(missing)} configured taxa NOT found in GTDB: "
              f"{', '.join(missing)}")
    flagged = [p for p in plan if p.flags]
    if flagged:
        print(f"  {len(flagged)} set(s) flagged for review")
    print("\nwritten:")
    for path in (report_txt, plan_tsv, counts_tsv, build_sh):
        print(f"  {path}")
    print(f"\nRead {os.path.basename(report_txt)} first, then run "
          f"{os.path.basename(build_sh)}.")


if __name__ == "__main__":
    main()
