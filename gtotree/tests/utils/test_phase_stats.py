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

    label, elapsed, rss_start, rss_end = enabled.recorded()[0]
    assert label == "Phase 1"
    assert elapsed >= 0
    assert rss_start > 0
    # ru_maxrss is a high-water mark, so it can never go down
    assert rss_end >= rss_start


def test_rss_delta_attributed_to_allocating_phase(enabled):
    """
    The point of keeping both start and end RSS: a phase that inherits a high
    water mark set by an earlier phase should show a ~zero delta, not the peak.
    """
    enabled.begin("Phase 1")
    blob = bytearray(120 * 1024 * 1024)
    enabled.begin("Phase 2")
    del blob
    enabled.finish()

    (_, _, _, p1_end), (_, _, p2_start, p2_end) = enabled.recorded()

    assert p1_end - enabled.recorded()[0][2] > 50    # Phase 1 grew
    assert abs(p2_end - p2_start) < 50               # Phase 2 inherited, didn't grow


def test_write_tsv_contents(enabled, tmp_path):
    enabled.begin("Phase 1")
    enabled.begin("Phase 2")

    path = enabled.write_tsv(str(tmp_path))
    assert path is not None

    lines = open(path).read().strip().split("\n")
    header = lines[0].split("\t")
    assert header == ["phase", "seconds", "rss_start_mb", "rss_end_mb", "rss_delta_mb"]

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
