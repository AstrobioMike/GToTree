import importlib
import os
import pytest # type: ignore
from gtotree.utils import phase_stats


@pytest.fixture
def enabled(monkeypatch):
    """
    DEBUG_TIMING is read once at import (deliberately -- it's a developer switch,
    not runtime-configurable), so toggling it for a test means setting it on the
    module rather than only in the environment.
    """
    monkeypatch.setattr(phase_stats, "DEBUG_TIMING", True)
    phase_stats.reset()
    yield phase_stats
    phase_stats.reset()


@pytest.fixture
def disabled(monkeypatch):
    monkeypatch.setattr(phase_stats, "DEBUG_TIMING", False)
    phase_stats.reset()
    yield phase_stats
    phase_stats.reset()


def test_disabled_records_nothing(disabled, tmp_path, capsys):
    disabled.begin("Phase 1")
    disabled.finish()
    disabled.report()

    assert disabled.recorded() == []
    assert disabled.write_tsv(str(tmp_path)) is None
    assert not os.path.exists(tmp_path / phase_stats.STATS_FILENAME)
    assert capsys.readouterr().err == ""


def test_env_var_default_is_off():
    """Absent the env var, the module must come up disabled."""
    if "GTT_DEBUG_TIMING" in os.environ:
        pytest.skip("GTT_DEBUG_TIMING set in this environment")
    assert importlib.reload(phase_stats).DEBUG_TIMING is False


def test_begin_closes_previous_phase(enabled):
    enabled.begin("Phase 1")
    enabled.begin("Phase 2")

    # Phase 1 closed when Phase 2 opened; Phase 2 still open
    assert [r[0] for r in enabled.recorded()] == ["Phase 1"]

    enabled.finish()
    assert [r[0] for r in enabled.recorded()] == ["Phase 1", "Phase 2"]


def test_finish_is_idempotent(enabled):
    enabled.begin("Phase 1")
    enabled.finish()
    enabled.finish()
    enabled.finish()

    assert len(enabled.recorded()) == 1


def test_finish_with_no_open_phase_is_safe(enabled):
    enabled.finish()
    assert enabled.recorded() == []


def test_records_elapsed_and_rss(enabled):
    enabled.begin("Phase 1")
    enabled.finish()

    label, elapsed, rss_start, rss_end, cur_end = enabled.recorded()[0]
    assert label == "Phase 1"
    assert elapsed >= 0
    assert rss_start > 0
    # ru_maxrss is a high-water mark, so it can never go down
    assert rss_end >= rss_start
    # current RSS is measured separately and may be unavailable on some platforms
    assert cur_end is None or cur_end > 0


def test_write_tsv_contents(enabled, tmp_path):
    enabled.begin("Phase 1")
    enabled.begin("Phase 2")

    path = enabled.write_tsv(str(tmp_path))
    assert path is not None

    lines = open(path).read().strip().split("\n")
    header = lines[0].split("\t")
    assert header == ["phase", "seconds", "rss_start_mb", "rss_end_mb", "rss_delta_mb",
                      "cur_rss_end_mb"]

    # write_tsv closes the open phase, so both land, plus a TOTAL row
    labels = [line.split("\t")[0] for line in lines[1:]]
    assert labels == ["Phase 1", "Phase 2", "TOTAL"]


def test_write_tsv_is_atomic(enabled, tmp_path):
    enabled.begin("Phase 1")
    enabled.write_tsv(str(tmp_path))

    # no .part left behind at the destination
    assert not list(tmp_path.glob("*.part"))


def test_write_tsv_survives_unwritable_dir(enabled, tmp_path):
    """Instrumentation must never take down a run that otherwise succeeded."""
    enabled.begin("Phase 1")
    assert enabled.write_tsv(str(tmp_path / "does-not-exist")) is None


def test_write_tsv_with_no_phases(enabled, tmp_path):
    assert enabled.write_tsv(str(tmp_path)) is None


def test_report_goes_to_stderr(enabled, capsys):
    enabled.begin("Phase 1")
    enabled.report()

    captured = capsys.readouterr()
    assert captured.out == ""          # must not pollute parsable stdout
    assert "Phase 1" in captured.err
    assert "GTT_DEBUG_TIMING" in captured.err


def test_total_row_uses_peak_not_sum(enabled):
    enabled.begin("Phase 1")
    enabled.begin("Phase 2")
    enabled.finish()

    rows = enabled._rows()
    total = rows[-1]
    assert total[0] == "TOTAL"

    peak = max(r[3] for r in enabled.recorded())
    assert float(total[3]) == pytest.approx(peak, abs=0.1)

def test_current_rss_can_fall_while_peak_cannot(enabled):
    """
    The reason the current-RSS column exists: `ru_maxrss` only ever rises, so it cannot
    distinguish "an earlier phase's memory is still held" from "it was handed back after
    an earlier peak". Current RSS can fall, and that difference is what determines how
    much headroom a later phase actually has.
    """
    enabled.begin("Phase 1")
    blob = bytearray(200 * 1024 * 1024)
    blob[::4096] = b"\x01" * len(blob[::4096])   # touch pages so they're truly resident
    enabled.begin("Phase 2")
    del blob
    enabled.finish()

    (_, _, _, p1_peak, p1_cur), (_, _, p2_start, p2_peak, p2_cur) = enabled.recorded()

    # the peak never goes down between phases
    assert p2_peak >= p1_peak

    if p1_cur is not None and p2_cur is not None:
        # ...but current RSS is allowed to, and that's the whole point of the column
        assert p2_cur <= p1_cur + 50


def test_checkpoint_is_silent_when_disabled(disabled, capsys):
    disabled.checkpoint("somewhere")
    assert capsys.readouterr().err == ""


def test_checkpoint_reports_to_stderr(enabled, capsys):
    """Sub-phase checkpoints go to stderr, like the phase table, so stdout stays clean."""
    enabled.checkpoint("after the big read")
    err = capsys.readouterr().err
    assert "after the big read" in err
    assert "current RSS" in err
