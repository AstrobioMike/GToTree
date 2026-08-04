"""
Resume support for `gtt gen-scg-hmms`

Resume is keyed on a fingerprint of everything that affects the result:

  * the resolved target accessions (the actual genome set, not how it was requested)
  * the selection parameters that produced them
  * `--percent-single-copy` and `--min-pfam-coverage`
  * the Pfam version in use

If the fingerprint matches, completed stages are reused. If it doesn't, resume is
refused with an explanation of what changed

Stage state lives in `<work_dir>/run-state.json`, written atomically after each stage
"""

import os
import json
import hashlib
import tempfile


STATE_FILENAME = "run-state.json"
STATE_VERSION = 1

# stage names, in pipeline order
STAGE_GENOMES = "genomes"
STAGE_PFAMS = "pfams"
STAGE_SEARCH = "search"

STAGE_ORDER = [STAGE_GENOMES, STAGE_PFAMS, STAGE_SEARCH]


def _hash_accessions(accessions):
    """
    Order-independent hash of the target accession set
    """
    h = hashlib.sha256()
    for acc in sorted(set(accessions)):
        h.update(acc.encode())
        h.update(b"\n")
    return h.hexdigest()


def _hash_local_genomes(local_genomes):
    """
    Hash the local genome inputs by id, resolved path, size, and mtime

    Size and mtime are included because, unlike an NCBI accession, a local file's
    *contents* can change while its path stays the same. Resuming against an edited
    fasta would silently mix old results with new input, so an edit invalidates the
    resume the same way adding a genome does
    """
    if not local_genomes:
        return None

    h = hashlib.sha256()
    entries = []
    for gd in local_genomes:
        try:
            stat = os.stat(gd.full_path)
            entries.append(f"{gd.id}\t{gd.source}\t{gd.full_path}\t{stat.st_size}\t{int(stat.st_mtime)}")
        except OSError:
            entries.append(f"{gd.id}\t{gd.source}\t{gd.full_path}\tmissing")

    for entry in sorted(entries):
        h.update(entry.encode())
        h.update(b"\n")

    return h.hexdigest()


def build_fingerprint(accessions, args, pfam_version=None, local_genomes=None):
    """
    Build the dict describing everything that affects the final SCG set

    Deliberately does NOT include: --num-cpus, --num-jobs, --keep-working-dir, or the
    output directory name. Those change how the run executes, not what it produces, so
    changing them shouldn't invalidate a resume.
    """
    return {
        "state_version": STATE_VERSION,
        "accessions_sha256": _hash_accessions(accessions),
        "num_accessions": len(set(accessions)),
        "local_genomes_sha256": _hash_local_genomes(local_genomes),
        "num_local_genomes": len(local_genomes or []),
        "percent_single_copy": args.percent_single_copy,
        "min_pfam_coverage": args.min_pfam_coverage,
        "source": (args.source or "").upper(),
        "wanted_ref_tax": args.wanted_ref_tax,
        "target_rank": args.target_rank,
        "derep_rank": args.derep_rank,
        "pfam_version": pfam_version,
    }


# human-readable labels for the fields we compare, for the refusal message
_FIELD_LABELS = {
    "accessions_sha256": "the set of target genomes",
    "local_genomes_sha256": "the local genome files (contents, paths, or set)",
    "percent_single_copy": "--percent-single-copy",
    "min_pfam_coverage": "--min-pfam-coverage",
    "source": "--source",
    "wanted_ref_tax": "--wanted-ref-tax",
    "target_rank": "--target-rank",
    "derep_rank": "--derep-rank",
    "pfam_version": "the Pfam version",
    "state_version": "the run-state format",
}


def compare_fingerprints(old, new):
    """
    Return a list of human-readable descriptions of what differs

    `pfam_version` is only compared when both sides know it -- the Pfam version isn't
    resolved until the Pfam stage runs, so a run interrupted before then legitimately
    has None stored and shouldn't be refused on that basis.
    """
    if not old:
        return ["no previous run state was found"]

    differences = []
    for key, label in _FIELD_LABELS.items():
        old_val = old.get(key)
        new_val = new.get(key)

        if key == "pfam_version" and (old_val is None or new_val is None):
            continue

        if old_val != new_val:
            if key == "accessions_sha256":
                differences.append(
                    f"{label} changed "
                    f"({old.get('num_accessions')} -> {new.get('num_accessions')})")
            elif key == "local_genomes_sha256":
                differences.append(
                    f"{label} changed "
                    f"({old.get('num_local_genomes')} -> {new.get('num_local_genomes')} file(s))")
            else:
                differences.append(f"{label} changed ({old_val!r} -> {new_val!r})")

    return differences


def state_path(work_dir):
    return os.path.join(work_dir, STATE_FILENAME)


def load_state(work_dir):
    """ Read the stored run state, returning None if absent or unreadable. """
    path = state_path(work_dir)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            state = json.load(f)
    except (OSError, ValueError):
        # a corrupt state file means we can't trust anything about the prior run;
        # treat it as no state rather than guessing
        return None

    if not isinstance(state, dict):
        return None

    return state


def save_state(work_dir, state):
    """ Write the run state atomically. """
    path = state_path(work_dir)
    os.makedirs(work_dir, exist_ok=True)

    fd, tmp_path = tempfile.mkstemp(dir=work_dir, suffix=".part")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(state, f, indent=2, sort_keys=True)
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise

    return path


def new_state(fingerprint):
    return {"fingerprint": fingerprint, "completed": {}}


def mark_stage_complete(state, stage, artifacts=None, work_dir=None):
    """
    Record a stage as done, along with the artifacts it produced

    Artifacts are stored RELATIVE to `work_dir`, with their sizes. Relative because the
    whole output directory can legitimately be renamed or moved between runs; absolute
    paths would then all miss, and every stage would silently re-run from scratch even
    though the files are sitting right there. Sizes are stored so `stage_is_reusable`
    can still tell a genuinely finished file from one truncated by a kill -9.
    """
    entry = {"artifacts": {}}
    for path in (artifacts or []):
        key = _relative_artifact(path, work_dir)
        try:
            entry["artifacts"][key] = os.path.getsize(path)
        except OSError:
            entry["artifacts"][key] = None
    state.setdefault("completed", {})[stage] = entry
    return state


def _relative_artifact(path, work_dir):
    """
    Express `path` relative to `work_dir` when it sits inside it
    """
    if not work_dir:
        return path
    try:
        rel = os.path.relpath(path, work_dir)
    except ValueError:
        # different drives on Windows
        return path
    if rel.startswith(os.pardir):
        return path
    return rel


def _resolve_artifact(key, work_dir):
    """ Turn a stored artifact key back into a usable path. """
    if os.path.isabs(key) or not work_dir:
        return key
    return os.path.join(work_dir, key)


def stage_is_reusable(state, stage, work_dir=None):
    """
    True only if the stage was marked complete AND every artifact it recorded is still
    present at the recorded size

    The size check matters: the pipeline writes atomically, so a file at its final path
    should be complete, but a working dir can also be touched by anything else on the
    system between runs, and silently reusing a truncated HMM would corrupt the result
    in a way that's very hard to trace back.
    """
    entry = (state or {}).get("completed", {}).get(stage)
    if not entry:
        return False

    for key, size in (entry.get("artifacts") or {}).items():
        path = _resolve_artifact(key, work_dir)
        if not os.path.exists(path):
            return False
        if size is not None and os.path.getsize(path) != size:
            return False

    return True


def invalidate_from(state, stage):
    """
    Drop the given stage and everything downstream of it

    Used when a stage has to be re-run, anything computed from its output is no longer
    trustworthy, even if it was marked complete
    """
    if stage not in STAGE_ORDER:
        return state

    start = STAGE_ORDER.index(stage)
    completed = state.setdefault("completed", {})
    for name in STAGE_ORDER[start:]:
        completed.pop(name, None)
    return state
