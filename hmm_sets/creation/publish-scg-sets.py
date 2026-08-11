#!/usr/bin/env python3
"""
Publish stage for the pre-built SCG-HMM sets.

Takes a directory of finished `gtt gen-scg-hmms` runs and:

  1. reads each run's outputs to build its row for `hmm_sets/hmm-sources-and-info.tsv`
  2. OVERWRITES that table with the current generation's rows
  3. writes a shell script of `gh release upload` commands

Kept separate from the build because the two fail differently. A build failure costs
hours and is worth resuming; an upload failure costs seconds and shouldn't require
rebuilding anything. Splitting them also leaves a natural point to eyeball the sets --
a set that came back with 12 retained genes is a signal something went wrong upstream,
and you want to notice that before it's public.

This script never uploads anything itself. It writes the commands and stops, so the
irreversible step stays an explicit action.

Usage
-----
    python publish-scg-sets.py --build-root scg-set-builds --release-tag gtdb-scg-sets-r220

Overwriting, not merging
------------------------
v2 is a clean swap: the info table is rewritten to describe only the current
generation. Older releases stay reachable by their own release tags, so nothing is lost
by not carrying their rows forward.
"""

import argparse
import csv
import hashlib
import os
import sys
from datetime import date


# Existing columns, in their published order, plus the provenance columns a generated
# set needs. The v1 rows came from a paper, so a DOI in `source` was the whole story;
# a v2 set comes from a reproducible command, and its provenance is which GTDB release
# and Pfam version produced it. Without those the table can't answer the main question
# anyone will have a year later.
BASE_COLUMNS = ["file", "source", "target_taxa", "rank", "num_genes", "link",
                "pfam_names", "pfam_accessions", "md5"]
# `domain` and `parent` are display columns: they're what lets `gtt hmms` lay the sets
# out as a taxonomy instead of one flat alphabetical list, without shipping GTDB itself
# or hardcoding "Gammaproteobacteria goes under Pseudomonadota" in the CLI. Both come
# from the plan table, which got them from the parquet.
ADDED_COLUMNS = ["domain", "parent", "gtdb_release", "pfam_version", "num_genomes",
                 "derep_rank", "date_built"]
ALL_COLUMNS = BASE_COLUMNS + ADDED_COLUMNS

# Sets that aren't built from GTDB and so have no build directory to scan. Their rows
# are copied verbatim into the published table. Right now that's just the Hug et al.
# universal set, a published external asset carried forward unchanged; keeping it as a
# file rather than a literal means updating it doesn't mean editing this script.
CARRIED_FORWARD = "carried-forward-sets.tsv"

DEFAULT_SOURCE = "https://doi.org/10.1093/bioinformatics/btz188"

# A set this small usually means something went wrong -- too few input genomes, or a
# `--percent-single-copy` that was too strict for a heterogeneous clade. Advisory.
FEW_GENES_THRESHOLD = 30

TARGETS_INFO = "SCG-targets-info.tsv"
TARGET_GENOMES = "target-genomes.tsv"
PFAM_VERSION = "pfam-version-used.txt"


def md5sum(path, chunk_size=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def read_plan(path):
    """
    set_name_lower -> plan row, for the metadata the build doesn't record.

    Keyed on `set_name` rather than `taxon` because that's what the build directories
    are named (`-o "$OUT_ROOT/<set_name>"`). They're the same string for a phylum or
    class set, but a domain-spanning set has taxon `Bacteria+Archaea` and set_name
    `Bacteria-and-Archaea`, so keying on taxon silently missed those rows and published
    them with rank/derep_rank of NA.
    """
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as f:
        return {r["set_name"].lower(): r for r in csv.DictReader(f, delimiter="\t")}


def read_carried_forward(path):
    """
    Rows for sets that have no build directory, copied through as-is.

    Returned as full info-table rows, not build outputs, because there's nothing to
    derive: the .hmm already exists on an old release and its md5 and gene list are
    whatever they were. Unknown columns are dropped and missing ones default to NA, so
    the file only has to carry what's actually known about the set.
    """
    if not path or not os.path.isfile(path):
        return []
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if not (row.get("file") or "").strip():
                continue
            rows.append({column: (row.get(column) or "NA").strip()
                         for column in ALL_COLUMNS})
    return rows


def write_table(path, rows):
    """
    Write a fresh info table.

    v2 is a clean swap: the v1 NCBI-taxonomy sets move to their own hosted release and
    the live table describes only the current generation, so this OVERWRITES rather than
    merging. There's deliberately no read-existing step -- an existing file is treated as
    blank. Older releases stay reachable by their own release tags regardless.

    Atomic, because a half-written info table is the kind of thing that gets committed
    before anyone notices.
    """
    tmp_path = path + ".part"
    try:
        with open(tmp_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=ALL_COLUMNS, delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            for row in sorted(rows, key=lambda r: r["file"].lower()):
                w.writerow({column: row.get(column, "") for column in ALL_COLUMNS})
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise


def _read_first_line(path):
    try:
        with open(path) as f:
            return f.read().strip().splitlines()[0]
    except (OSError, IndexError):
        return "NA"


def scan_build_dir(build_dir):
    """
    Read one finished `gen-scg-hmms` run.

    Returns a dict of everything the info table needs, or None with a reason if the
    run doesn't look complete -- an interrupted build leaves a directory behind, and
    publishing a partial set would be worse than skipping it.
    """
    name = os.path.basename(os.path.normpath(build_dir))

    targets_path = os.path.join(build_dir, TARGETS_INFO)
    if not os.path.isfile(targets_path):
        return None, f"no {TARGETS_INFO} (run didn't finish)"

    hmm_candidates = [f for f in sorted(os.listdir(build_dir)) if f.endswith(".hmm")]
    if not hmm_candidates:
        return None, "no .hmm file"
    if len(hmm_candidates) > 1:
        return None, f"multiple .hmm files ({', '.join(hmm_candidates)})"
    hmm_file = hmm_candidates[0]
    hmm_path = os.path.join(build_dir, hmm_file)

    accessions, names = [], []
    with open(targets_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            accessions.append(row["pfam_id"])
            names.append(row["name"])

    if not accessions:
        return None, f"{TARGETS_INFO} is empty"

    num_genomes = "NA"
    genomes_path = os.path.join(build_dir, TARGET_GENOMES)
    if os.path.isfile(genomes_path):
        with open(genomes_path) as f:
            num_genomes = str(sum(1 for _ in f) - 1)

    return {
        "name": name,
        "hmm_file": hmm_file,
        "hmm_path": hmm_path,
        "num_genes": len(accessions),
        "pfam_names": ",".join(names),
        "pfam_accessions": ",".join(accessions),
        "md5": md5sum(hmm_path),
        "num_genomes": num_genomes,
        "pfam_version": _read_first_line(os.path.join(build_dir, PFAM_VERSION)),
    }, None


def build_row(info, plan_row, release_url, gtdb_release, source, built_on):
    def from_plan(column):
        value = (plan_row or {}).get(column) or "NA"
        return value.strip() or "NA"

    return {
        "file": info["hmm_file"],
        "source": source,
        "target_taxa": info["name"].lower(),
        "rank": from_plan("rank"),
        "domain": from_plan("domain"),
        "parent": from_plan("parent"),
        "num_genes": str(info["num_genes"]),
        "link": f"{release_url}/{info['hmm_file']}",
        "pfam_names": info["pfam_names"],
        "pfam_accessions": info["pfam_accessions"],
        "md5": info["md5"],
        "gtdb_release": gtdb_release,
        "pfam_version": info["pfam_version"],
        "num_genomes": info["num_genomes"],
        "derep_rank": from_plan("derep_rank"),
        "date_built": built_on,
    }


def write_upload_script(path, infos, release_tag, repo):
    lines = [
        "#!/usr/bin/env bash",
        "# Generated by publish_scg_sets.py -- upload the built SCG-HMM sets.",
        "#",
        "# `--clobber` replaces an asset of the same name, so this is safe to re-run",
        "# after a partial upload. Creating the release is left commented out: making",
        "# a release is the one step worth doing deliberately.",
        "set -euo pipefail",
        "",
        f'RELEASE_TAG="{release_tag}"',
        f'REPO="{repo}"',
        "",
        '# gh release create "$RELEASE_TAG" --repo "$REPO" \\',
        '#     --title "$RELEASE_TAG" --notes "SCG-HMM sets"',
        "",
    ]
    for info in infos:
        lines.append(f"# {info['name']}: {info['num_genes']} genes from "
                     f"{info['num_genomes']} genomes")
        lines.append(f'gh release upload "$RELEASE_TAG" "{info["hmm_path"]}" '
                     f'--repo "$REPO" --clobber')
        lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))
    os.chmod(path, 0o755)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="Regenerate the SCG-HMM info table from finished builds and emit "
                    "upload commands.")
    parser.add_argument("--build-root", required=True,
                        help="directory holding one subdirectory per finished build")
    parser.add_argument("--release-tag", required=True,
                        help="GitHub release tag the assets will live under")
    parser.add_argument("--repo", default="AstrobioMike/GToTree",
                        help="owner/repo for the release (default: %(default)s)")
    parser.add_argument("--plan", default=os.path.join(here, "plan", "scg-set-plan.tsv"),
                        help="the plan table, for rank/domain/parent/derep_rank")
    parser.add_argument("--carried-forward",
                        default=os.path.join(here, CARRIED_FORWARD),
                        help="rows for sets with no build directory, copied verbatim "
                             "into the table (default: %(default)s)")
    parser.add_argument("--info-table", default=None,
                        help="info table to write (default: "
                             "<repo-root>/hmm_sets/hmm-sources-and-info.tsv). This is "
                             "OVERWRITTEN, not merged -- see write_table.")
    parser.add_argument("--gtdb-release", default="NA",
                        help="GTDB release label to record for these sets")
    parser.add_argument("--source", default=DEFAULT_SOURCE,
                        help="value for the `source` column (default: %(default)s)")
    parser.add_argument("-o", "--output-dir", default=os.path.join(here, "publish"),
                        help="where to write the upload script (default: ./publish)")
    args = parser.parse_args()

    if not os.path.isdir(args.build_root):
        sys.exit(f"No such build root: {args.build_root}")

    info_table = args.info_table or os.path.join(
        os.path.dirname(os.path.dirname(here)), "hmm_sets", "hmm-sources-and-info.tsv")

    plan = read_plan(args.plan)
    release_url = f"https://github.com/{args.repo}/releases/download/{args.release_tag}"
    built_on = date.today().isoformat()

    subdirs = sorted(
        os.path.join(args.build_root, d) for d in os.listdir(args.build_root)
        if os.path.isdir(os.path.join(args.build_root, d)))
    if not subdirs:
        sys.exit(f"No build directories found under {args.build_root}")

    infos, new_rows, skipped, thin = [], [], [], []
    for subdir in subdirs:
        info, problem = scan_build_dir(subdir)
        if info is None:
            skipped.append((os.path.basename(subdir), problem))
            continue
        plan_row = plan.get(info["name"].lower())
        infos.append(info)
        new_rows.append(build_row(info, plan_row, release_url, args.gtdb_release,
                                  args.source, built_on))
        if info["num_genes"] < FEW_GENES_THRESHOLD:
            thin.append((info["name"], info["num_genes"]))

    if not new_rows:
        sys.exit("No complete builds found -- nothing to publish.")

    # A partial build root would produce a table describing only the sets that finished.
    # For a clean-swap publish that silently shrinks the published set, so compare the
    # count against the plan (if present) and make the operator confirm the shortfall is
    # intended rather than a half-finished build.
    n_planned = len(plan)
    if n_planned and len(new_rows) < n_planned:
        print(f"WARNING: {len(new_rows)} built set(s), but the plan lists {n_planned}.")
        print("The info table is overwritten, so publishing now describes only these")
        print(f"{len(new_rows)}. If the build isn't finished, stop here. Continuing writes")
        print("the table anyway; the upload script only touches the built .hmm files.\n")

    # Rows for sets with no build directory. Appended after the built ones so a name
    # collision resolves in favour of the freshly built set rather than a stale copy.
    carried = [row for row in read_carried_forward(args.carried_forward)
               if row["file"].lower() not in {r["file"].lower() for r in new_rows}]

    os.makedirs(args.output_dir, exist_ok=True)
    write_table(info_table, new_rows + carried)
    upload_script = os.path.join(args.output_dir, "upload-scg-sets.sh")
    write_upload_script(upload_script, infos, args.release_tag, args.repo)

    print(f"scanned {len(subdirs)} build director(ies)")
    print(f"  wrote {len(new_rows)} set(s) to the info table (overwritten fresh)")
    if carried:
        print(f"  plus {len(carried)} carried-forward set(s): "
              f"{', '.join(r['file'] for r in carried)}")
    if skipped:
        print(f"\n  SKIPPED {len(skipped)} incomplete build(s):")
        for name, problem in skipped:
            print(f"    {name}: {problem}")
    if thin:
        print(f"\n  WORTH A LOOK -- fewer than {FEW_GENES_THRESHOLD} retained genes:")
        for name, count in thin:
            print(f"    {name}: {count} genes")
        print("    Usually means too few input genomes or too strict a "
              "--percent-single-copy.")

    print(f"\nwritten:\n  {info_table}\n  {upload_script}")
    print("\nNothing has been uploaded. Review the table, then run the upload script.")


if __name__ == "__main__":
    main()
