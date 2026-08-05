"""
Phase timing and peak-memory instrumentation

This is deliberately *not* exposed as a CLI flag, it's enabled by setting the
environment variable `GTT_DEBUG_TIMING` to a truthy value, e.g.:

    GTT_DEBUG_TIMING=1 gtt gen-scg-hmms -W nitrospirota ...

At initial writin git's just for gen-scg-hmms, but i'm putting it here and trying
to make it somewhat generic in case it's helpful elsewhere at some point

When unset (the normal case), every entry point here is a single boolean check
and nothing is recorded, printed, or written.

This is to help decide where tuning might help, and once a
memory-derived heuristic exists anywhere in the pipeline, we can just ask a
user to run things with this on to help debug.

Usage: `begin(label)` at each phase boundary closes out the previous phase and
opens a new one, so callers only ever make one call per phase. `finish()` closes
the last open phase and is safe to call more than once.

NOTE ON `ru_maxrss`: this is a high-water mark for the *whole process*, and it
only ever increases. "Peak RSS at end of phase N" is therefore not "memory used
by phase N" -- it's the ceiling the whole process had reached by that point. Both the
starting and ending values are recorded per phase so that a phase which set a
new high is distinguishable from one that merely inherited it, and there is a `delta`
column.
"""

import os
import sys
import time
import resource

_OFF_VALUES = {"", "0", "false", "no", "off"}

def debug_timing_enabled(environ=None):
    """True if the GTT_DEBUG_TIMING switch is set to anything meaning "on"."""
    raw = (os.environ if environ is None else environ).get("GTT_DEBUG_TIMING", "")
    return str(raw).strip().lower() not in _OFF_VALUES

# read once at import; this is a developer switch, not runtime-configurable
DEBUG_TIMING = debug_timing_enabled()

STATS_FILENAME = "phase-stats.tsv"

# (label, elapsed_seconds, rss_start_mb, rss_end_mb, cur_rss_end_mb)
_completed = []
_current = None


def _peak_rss_mb():
    """
    Process peak RSS in MB.

    `ru_maxrss` units differ by platform: Linux reports KiB, macOS/BSD report
    bytes. Getting this backwards understates by 1024x, which would silently
    validate whatever chunk size happened to be picked, so it's handled
    explicitly rather than assumed.
    """
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 ** 2)
    return rss / 1024


def _current_rss_mb():
    """
    Process RSS *right now*, in MB -- the value that can go down again.

    This is the necessary companion to `_peak_rss_mb`. Because `ru_maxrss` only ever
    rises, it cannot answer the question that actually matters when deciding how much
    memory a later phase may safely use: is an earlier phase's memory still resident, or
    was it handed back? A phase's `rss_start` of 1.4 GB is equally consistent with 1.4 GB
    still being held and with 100 MB being held after an earlier peak of 1.4 GB, and
    those two worlds imply very different headroom.

    There's no stdlib call for this, and psutil isn't a dependency, so it's read from
    /proc on Linux and from `ps` on macOS. Returns None if it can't be determined, which
    callers render as "n/a" rather than guessing.
    """
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/self/statm") as f:
                pages = int(f.read().split()[1])
            return pages * os.sysconf("SC_PAGE_SIZE") / (1024 ** 2)

        if sys.platform == "darwin":
            import subprocess
            out = subprocess.run(
                ["ps", "-o", "rss=", "-p", str(os.getpid())],
                capture_output=True, text=True, timeout=5)
            value = out.stdout.strip()
            if value:
                return int(value) / 1024      # ps reports KiB on macOS
    except Exception:
        # instrumentation must never take down a run
        return None

    return None


def checkpoint(label):
    """
    Print current RSS at a point *inside* a phase, when GTT_DEBUG_TIMING is on.

    The phase table shows where memory ended up per phase, which is enough to see that a
    phase retained something but not enough to see *which step* did. Dropping a couple of
    these into a phase localizes it: the step whose checkpoint jumps is the one holding
    the memory. No-op unless the debug switch is set, so these can be left in place.
    """
    if not DEBUG_TIMING:
        return

    cur = _current_rss_mb()
    shown = "n/a" if cur is None else f"{cur:.1f} MB"
    print(f"      [GTT_DEBUG_TIMING] {label}: current RSS {shown}", file=sys.stderr)


def begin(label):
    """Close out the current phase (if any) and start timing a new one."""
    if not DEBUG_TIMING:
        return

    finish()

    global _current
    _current = (label, time.monotonic(), _peak_rss_mb())


def finish():
    """
    Close out the currently open phase, if there is one. Safe to call repeatedly
    and safe to call when no phase is open, so it can go in an exception path
    without needing to know how far the run got.
    """
    if not DEBUG_TIMING:
        return

    global _current
    if _current is None:
        return

    label, started, rss_start = _current
    _completed.append((label, time.monotonic() - started, rss_start, _peak_rss_mb(),
                       _current_rss_mb()))
    _current = None


def reset():
    """Drop all recorded state (used by the tests)."""
    global _current
    _completed.clear()
    _current = None


def recorded():
    """Return the completed phase records, for tests and for the writer."""
    return list(_completed)


def _rows():
    """Completed phases plus a total row, as display-ready string tuples."""
    rows = []
    for label, elapsed, rss_start, rss_end, cur_end in _completed:
        rows.append((label,
                     f"{elapsed:.2f}",
                     f"{rss_start:.1f}",
                     f"{rss_end:.1f}",
                     f"{rss_end - rss_start:+.1f}",
                     "n/a" if cur_end is None else f"{cur_end:.1f}"))

    if _completed:
        total_elapsed = sum(r[1] for r in _completed)
        peak = max(r[3] for r in _completed)
        rows.append(("TOTAL", f"{total_elapsed:.2f}", "", f"{peak:.1f}", "", ""))

    return rows


def report(indent="      "):
    """Print the phase table to stderr, so it never pollutes parsable stdout."""
    if not DEBUG_TIMING:
        return

    finish()
    rows = _rows()
    if not rows:
        return

    header = ("phase", "seconds", "rss_start_mb", "rss_end_mb", "rss_delta_mb",
              "cur_rss_end_mb")
    widths = [max(len(header[i]), max(len(r[i]) for r in rows))
              for i in range(len(header))]

    def fmt(cells):
        return indent + "  ".join(c.ljust(widths[i]) for i, c in enumerate(cells)).rstrip()

    print(file=sys.stderr)
    print(f"{indent}[GTT_DEBUG_TIMING] phase stats", file=sys.stderr)
    print(fmt(header), file=sys.stderr)
    print(indent + "  ".join("-" * w for w in widths), file=sys.stderr)
    for row in rows:
        print(fmt(row), file=sys.stderr)
    print(file=sys.stderr)


def write_tsv(out_dir):
    """
    Write the phase table into `out_dir` so runs can be diffed against each
    other rather than compared by scrollback. Never raises: instrumentation
    failing must not take down a run that otherwise succeeded.
    """
    if not DEBUG_TIMING or not out_dir:
        return None

    finish()
    rows = _rows()
    if not rows:
        return None

    path = os.path.join(out_dir, STATS_FILENAME)
    tmp_path = path + ".part"
    try:
        with open(tmp_path, "w") as out:
            out.write("phase\tseconds\trss_start_mb\trss_end_mb\trss_delta_mb"
                      "\tcur_rss_end_mb\n")
            for row in rows:
                out.write("\t".join(row) + "\n")
        os.replace(tmp_path, path)
    except OSError:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass
        return None

    return path
