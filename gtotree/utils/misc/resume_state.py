"""
Generic resume-state machinery.

Every GToTree program that supports resuming uses this. The machinery is here; the
configuration lives with each program, declared as a `ResumeProfile`:

  * `GToTree`          -- gtotree/utils/preflight_checks.py
  * `gtt gen-scg-hmms` -- gtotree/utils/hmms/gen_scg_hmms/gen_scg_hmms_cli.py
  * `gtt search-pfams`
    `gtt search-kos`   -- gtotree/utils/target_search/target_search_cli.py

The split is deliberate. What belongs here is program-agnostic: hashing, the state
file, artifact validation, and turning a fingerprint mismatch into a readable
explanation. What belongs with each program is which of *its* arguments affect the
result and what its stages are -- because that's where those arguments are defined, and
where someone adding a new flag will already be looking. Centralizing that config would
mean editing this file every time any program grows a flag.

The contract in every case is the same: if the fingerprint matches, completed work is
reused; if it doesn't, resume is refused with an explanation of what changed.

Not every program needs stage tracking. The main GToTree driver already records its
progress inside run-data.json (per-genome flags, `SCG_hits_filtered`,
`all_SCG_sets_aligned`, and so on), so it uses a profile with no stages -- just the
fingerprint half. Programs that do track stages get the `run-state.json` machinery too.

Stage state is a JSON file written atomically after each stage, so a run killed
mid-write never leaves a half-parsed state behind.
"""

import os
import json
import hashlib
import tempfile


STATE_FILENAME = "run-state.json"
STATE_VERSION = 1


def hash_strings(values):
    """
    Order-independent hash of a set of strings.

    Order-independent because the same genome set requested in a different order is
    the same genome set, and shouldn't invalidate a resume.
    """
    h = hashlib.sha256()
    for value in sorted(set(values)):
        h.update(str(value).encode())
        h.update(b"\n")
    return h.hexdigest()


def hash_file_contents(path):
    """
    Hash a file's contents, or None if it can't be read.

    Used for target-id lists (`-p` / `-K`): unlike an accession set, the meaningful
    unit is what's inside the file, and the path itself is irrelevant.
    """
    if not path:
        return None
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def hash_local_genomes(local_genomes):
    """
    Hash local genome inputs by id, source, resolved path, size, and mtime.

    Size and mtime are included because, unlike an NCBI accession, a local file's
    *contents* can change while its path stays the same. Resuming against an edited
    fasta would silently mix old results with new input, so an edit invalidates the
    resume the same way adding a genome does.
    """
    if not local_genomes:
        return None

    h = hashlib.sha256()
    entries = []
    for gd in local_genomes:
        try:
            stat = os.stat(gd.full_path)
            entries.append(f"{gd.id}\t{gd.source}\t{gd.full_path}\t"
                           f"{stat.st_size}\t{int(stat.st_mtime)}")
        except OSError:
            entries.append(f"{gd.id}\t{gd.source}\t{gd.full_path}\tmissing")

    for entry in sorted(entries):
        h.update(entry.encode())
        h.update(b"\n")

    return h.hexdigest()


def compare_fingerprints(old, new, field_labels, deferred_fields=()):
    """
    Return a list of human-readable descriptions of what differs between two
    fingerprints.

    `field_labels` maps fingerprint key -> the phrase used in the refusal message, and
    also decides which keys are compared at all: anything not in it is ignored.

    `deferred_fields` names keys that aren't known until partway through a run (the
    Pfam version, say). Those are only compared when *both* sides know them, so a run
    interrupted before that point legitimately has None stored and isn't refused on
    that basis.

    `count_fields` maps a hash field to a companion count field, so a changed genome
    set can report "12 -> 15" instead of two opaque hashes.
    """
    if not old:
        return ["no previous run state was found"]

    differences = []
    for key, label in field_labels.items():
        old_val = old.get(key)
        new_val = new.get(key)

        if key in deferred_fields and (old_val is None or new_val is None):
            continue

        if old_val == new_val:
            continue

        count_key = _COUNT_COMPANIONS.get(key)
        if count_key:
            differences.append(
                f"{label} changed ({old.get(count_key)} -> {new.get(count_key)})")
        elif key.endswith("_sha256"):
            # a raw hash diff tells the user nothing useful
            differences.append(f"{label} changed")
        else:
            differences.append(f"{label} changed ({old_val!r} -> {new_val!r})")

    return differences


# hash field -> the count field that makes its diff readable
_COUNT_COMPANIONS = {
    "accessions_sha256": "num_accessions",
    "local_genomes_sha256": "num_local_genomes",
    "targets_sha256": "num_targets",
}


def state_path(work_dir, filename=STATE_FILENAME):
    return os.path.join(work_dir, filename)


def load_state(work_dir, filename=STATE_FILENAME):
    """
    Read the stored run state, returning None if absent or unreadable.
    """
    path = state_path(work_dir, filename)
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


def save_state(work_dir, state, filename=STATE_FILENAME):
    """
    Write the run state atomically.
    """
    path = state_path(work_dir, filename)
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
    Record a stage as done, along with the artifacts it produced.

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
    Express `path` relative to `work_dir` when it sits inside it.
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
    """
    Turn a stored artifact key back into a usable path.
    """
    if os.path.isabs(key) or not work_dir:
        return key
    return os.path.join(work_dir, key)


def stage_is_reusable(state, stage, work_dir=None):
    """
    True only if the stage was marked complete AND every artifact it recorded is still
    present at the recorded size.

    The size check matters: the pipelines write atomically, so a file at its final path
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


def invalidate_from(state, stage, stage_order):
    """
    Drop the given stage and everything downstream of it.

    Used when a stage has to be re-run: anything computed from its output is no longer
    trustworthy, even if it was marked complete.
    """
    if stage not in stage_order:
        return state

    start = stage_order.index(stage)
    completed = state.setdefault("completed", {})
    for name in stage_order[start:]:
        completed.pop(name, None)
    return state


def save_sidecar(work_dir, filename, payload):
    """
    Write a stage sidecar (arbitrary JSON that a stage needs to hand downstream)
    atomically.
    """
    path = os.path.join(work_dir, filename)
    tmp_path = path + ".part"
    try:
        with open(tmp_path, "w") as f:
            json.dump(payload, f)
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise
    return path


def load_sidecar(work_dir, filename):
    path = os.path.join(work_dir, filename)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


################################################################################
# ResumeProfile
################################################################################

class ResumeProfile:
    """
    One program's resume configuration, with the generic machinery bound to it.

    Constructed once per program, near that program's argument definitions:

        RESUME = ResumeProfile(
            name="gen-scg-hmms",
            field_labels={"accessions_sha256": "the set of target genomes", ...},
            stages=["genomes", "pfams", "search"],
            deferred_fields=("pfam_version",),
        )

    Then `RESUME.compare(old, new)`, `RESUME.load(work_dir)`, and so on. Binding the
    config to the methods is what lets the per-program modules drop to a declaration
    plus a `build_fingerprint` -- there's no re-export layer to maintain.

    `stages` is optional. A profile without it supports the fingerprint half only, and
    the state-file methods raise rather than silently writing a state file the program
    never reads.
    """

    __slots__ = ("name", "field_labels", "stages", "deferred_fields", "state_filename")

    def __init__(self, name, field_labels, stages=None, deferred_fields=(),
                 state_filename=STATE_FILENAME):
        self.name = name
        self.field_labels = dict(field_labels)
        self.stages = list(stages) if stages else []
        self.deferred_fields = tuple(deferred_fields)
        self.state_filename = state_filename

        unknown = set(self.deferred_fields) - set(self.field_labels)
        if unknown:
            # a deferred field that isn't compared at all is almost certainly a typo,
            # and it would silently weaken the check rather than fail
            raise ValueError(
                f"{name}: deferred_fields not present in field_labels: {sorted(unknown)}")

    def __repr__(self):
        return f"<ResumeProfile {self.name!r} stages={self.stages}>"

    # --- fingerprints ---------------------------------------------------------

    def compare(self, old, new):
        """
        Return human-readable descriptions of what differs, empty if they match.
        """
        return compare_fingerprints(old, new, self.field_labels,
                                    deferred_fields=self.deferred_fields)

    def refusal_message(self, differences, flag="`--resume`"):
        """
        Format a fingerprint mismatch into the message the CLI layer shows.

        Kept here so all three programs refuse in the same voice, rather than each
        assembling its own bullet list.
        """
        return (f"{flag} was specified, but this run doesn't match the previous one in "
                "that directory:\n        - " + "\n        - ".join(differences)
                + "\n\n      Resuming would mix results from two different runs. Use a "
                "new output directory, or start over with `-F`.")

    # --- state file -----------------------------------------------------------

    def _require_stages(self):
        if not self.stages:
            raise RuntimeError(
                f"{self.name}: this profile has no stages, so it doesn't use a state "
                "file -- it tracks progress elsewhere.")

    def _check_stage(self, stage):
        self._require_stages()
        if stage not in self.stages:
            raise ValueError(
                f"{self.name}: unknown stage {stage!r}; expected one of {self.stages}")

    def state_path(self, work_dir):
        self._require_stages()
        return state_path(work_dir, self.state_filename)

    def load(self, work_dir):
        self._require_stages()
        return load_state(work_dir, self.state_filename)

    def save(self, work_dir, state):
        self._require_stages()
        return save_state(work_dir, state, self.state_filename)

    def new(self, fingerprint):
        return new_state(fingerprint)

    def mark_complete(self, state, stage, artifacts=None, work_dir=None):
        self._check_stage(stage)
        return mark_stage_complete(state, stage, artifacts=artifacts, work_dir=work_dir)

    def is_reusable(self, state, stage, work_dir=None):
        self._check_stage(stage)
        return stage_is_reusable(state, stage, work_dir=work_dir)

    def invalidate_from(self, state, stage):
        self._check_stage(stage)
        return invalidate_from(state, stage, self.stages)

    # --- sidecars -------------------------------------------------------------

    def save_sidecar(self, work_dir, filename, payload):
        return save_sidecar(work_dir, filename, payload)

    def load_sidecar(self, work_dir, filename):
        return load_sidecar(work_dir, filename)
