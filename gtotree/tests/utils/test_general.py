import socket
import urllib.error
import urllib.request
import pytest # type: ignore

from gtotree.utils import general
from gtotree.utils.general import download_with_tqdm


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

class _FakeResp:
    """Context-manager stand-in for a urlopen response with a chunked body."""
    def __init__(self, headers=None, chunks=None):
        self.headers = headers or {}
        self._chunks = list(chunks if chunks is not None else [])

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self, n=-1):
        if not self._chunks:
            return b""
        return self._chunks.pop(0)


def _http_error(url="http://x/y", code=500):
    return urllib.error.HTTPError(url, code, f"err{code}", {}, None)


def _body(data=b"payload"):
    """A fresh _FakeResp streaming `data` then EOF (headers carry the size)."""
    return _FakeResp(headers={"Content-Length": str(len(data))}, chunks=[data])


@pytest.fixture
def no_sleep(monkeypatch):
    """Make retry waits instant."""
    monkeypatch.setattr(general.time, "sleep", lambda *_a, **_k: None)


# --------------------------------------------------------------------------- #
# success / retry / error behavior (gate off, the default)
# --------------------------------------------------------------------------- #

def test_success_single_attempt(tmp_path, monkeypatch):
    """Healthy download: streamed once, file written, no retries."""
    opens = {"n": 0}

    def fake_urlopen(target, *a, **k):
        opens["n"] += 1
        return _body(b"payload")

    monkeypatch.setattr(general.urllib.request, "urlopen", fake_urlopen)

    out = tmp_path / "out.bin"
    download_with_tqdm("http://x/y", "label", str(out))

    assert out.read_bytes() == b"payload"


def test_retries_then_succeeds(tmp_path, monkeypatch, no_sleep):
    """Transient failures are retried; eventual success writes the file."""
    calls = {"n": 0}

    def flaky_urlopen(target, *a, **k):
        # size-probe call passes a plain str url; stream call passes a Request.
        # Only fail on the streaming opens so we count transfer attempts.
        is_stream = isinstance(target, urllib.request.Request)
        if not is_stream:
            # size probe: no Content-Length, succeeds
            return _FakeResp(headers={})
        calls["n"] += 1
        if calls["n"] < 3:
            raise ConnectionResetError("transient")
        return _FakeResp(headers={}, chunks=[b"ok"])

    monkeypatch.setattr(general.urllib.request, "urlopen", flaky_urlopen)

    out = tmp_path / "out.bin"
    download_with_tqdm("http://x/y", "label", str(out),
                       attempts=5, retry_wait=0)

    assert calls["n"] == 3          # failed twice, succeeded on the third
    assert out.read_bytes() == b"ok"


def test_exhausts_attempts_then_raises(tmp_path, monkeypatch, no_sleep):
    """If every attempt fails transiently, the last error propagates."""
    calls = {"n": 0}

    def always_fail(target, *a, **k):
        if not isinstance(target, urllib.request.Request):
            return _FakeResp(headers={})       # size probe ok
        calls["n"] += 1
        raise socket.timeout("nope")

    monkeypatch.setattr(general.urllib.request, "urlopen", always_fail)

    with pytest.raises(socket.timeout):
        download_with_tqdm("http://x/y", "label", str(tmp_path / "out.bin"),
                           attempts=3, retry_wait=0)

    assert calls["n"] == 3          # tried exactly `attempts` times


def test_404_not_retried(tmp_path, monkeypatch, no_sleep):
    """A 404 is definitive: raised immediately, never retried."""
    calls = {"n": 0}

    def not_found(target, *a, **k):
        if not isinstance(target, urllib.request.Request):
            # let the size probe fail cleanly too -> swallowed, total=None
            raise _http_error(code=404)
        calls["n"] += 1
        raise _http_error(code=404)

    monkeypatch.setattr(general.urllib.request, "urlopen", not_found)

    with pytest.raises(urllib.error.HTTPError) as ei:
        download_with_tqdm("http://x/y", "label", str(tmp_path / "out.bin"),
                           attempts=5, retry_wait=0)

    assert ei.value.code == 404
    assert calls["n"] == 1          # exactly one streaming attempt, no retries


def test_non_404_http_error_is_retried(tmp_path, monkeypatch, no_sleep):
    """A non-404 HTTP error (e.g. 503) is transient and retried."""
    calls = {"n": 0}

    def flaky_http(target, *a, **k):
        if not isinstance(target, urllib.request.Request):
            return _FakeResp(headers={})
        calls["n"] += 1
        if calls["n"] < 2:
            raise _http_error(code=503)
        return _FakeResp(headers={}, chunks=[b"ok"])

    monkeypatch.setattr(general.urllib.request, "urlopen", flaky_http)

    out = tmp_path / "out.bin"
    download_with_tqdm("http://x/y", "label", str(out), attempts=5, retry_wait=0)

    assert calls["n"] == 2
    assert out.read_bytes() == b"ok"


def test_retries_disabled_single_shot(tmp_path, monkeypatch, no_sleep):
    """retries=False collapses to a single attempt that raises on failure."""
    calls = {"n": 0}

    def fail_once(target, *a, **k):
        if not isinstance(target, urllib.request.Request):
            return _FakeResp(headers={})
        calls["n"] += 1
        raise ConnectionResetError("transient")

    monkeypatch.setattr(general.urllib.request, "urlopen", fail_once)

    with pytest.raises(ConnectionResetError):
        download_with_tqdm("http://x/y", "label", str(tmp_path / "out.bin"),
                           retries=False)

    assert calls["n"] == 1


# --------------------------------------------------------------------------- #
# urlopen passthrough
# --------------------------------------------------------------------------- #

def test_urlopen_passthrough(monkeypatch):
    """urlopen=True returns the response object directly, no file written."""
    sentinel = _FakeResp(headers={"Content-Length": "5"})
    monkeypatch.setattr(general.urllib.request, "urlopen",
                        lambda *a, **k: sentinel)

    result = download_with_tqdm("http://x/y", "label", urlopen=True)
    assert result is sentinel


# --------------------------------------------------------------------------- #
# speed-gated route rerolling
# --------------------------------------------------------------------------- #

def test_speed_gate_slow_reroll_then_accept(tmp_path, monkeypatch, no_sleep):
    """
    A persistently slow stream is rerolled on non-final attempts, then accepted
    on the final attempt so the download still completes. 'Slow' is simulated by
    freezing the clock so the probe window is immediately 'elapsed' while only a
    tiny number of bytes have arrived -> measured rate below the floor.
    """
    ticks = {"t": 0.0}
    def fake_monotonic():
        ticks["t"] += 100.0     # each call jumps 100s -> always past probe window
        return ticks["t"]
    monkeypatch.setattr(general.time, "monotonic", fake_monotonic)

    opens = {"stream": 0}
    def fake_urlopen(target, *a, **k):
        if not isinstance(target, urllib.request.Request):
            return _FakeResp(headers={})           # size probe
        opens["stream"] += 1
        # tiny first chunk -> rate far below 2 MB/s when probe fires
        return _FakeResp(headers={}, chunks=[b"x" * 1024])

    monkeypatch.setattr(general.urllib.request, "urlopen", fake_urlopen)

    out = tmp_path / "out.bin"
    download_with_tqdm("http://x/y", "label", str(out),
                       speed_gate=True, min_mbps=2.0, probe_seconds=1.0,
                       attempts=3, retry_wait=0)

    # first two attempts reroll (too slow), third (final) accepts and finishes
    assert opens["stream"] == 3
    assert out.read_bytes() == b"x" * 1024


def test_speed_gate_fast_completes_first_try(tmp_path, monkeypatch, no_sleep):
    """A fast stream passes the probe and completes on the first attempt."""
    # clock stays put (0 elapsed) so a big first chunk => very high rate.
    monkeypatch.setattr(general.time, "monotonic", lambda: 0.0)

    opens = {"stream": 0}
    def fake_urlopen(target, *a, **k):
        if not isinstance(target, urllib.request.Request):
            return _FakeResp(headers={})
        opens["stream"] += 1
        return _FakeResp(headers={}, chunks=[b"x" * (4 * 1024 * 1024)])

    monkeypatch.setattr(general.urllib.request, "urlopen", fake_urlopen)

    out = tmp_path / "out.bin"
    download_with_tqdm("http://x/y", "label", str(out),
                       speed_gate=True, min_mbps=2.0, probe_seconds=1.0,
                       attempts=3, retry_wait=0)

    assert opens["stream"] == 1     # passed on the first attempt, no reroll
    assert out.read_bytes() == b"x" * (4 * 1024 * 1024)


# --------------------------------------------------------------------------- #
# atomic write: a failed/interrupted download must not leave a file at the
# destination (which os.path.isfile gates would later mistake for complete),
# nor leave a stray .part alongside it.
# --------------------------------------------------------------------------- #

def test_failed_download_leaves_no_file_at_destination(tmp_path, monkeypatch, no_sleep):
    """Every attempt fails: no partial file at the final path, no .part left."""
    def failing_urlopen(target, *a, **k):
        is_stream = isinstance(target, urllib.request.Request)
        if not is_stream:
            return _FakeResp(headers={})          # size probe succeeds
        raise ConnectionResetError("boom")        # every stream attempt fails

    monkeypatch.setattr(general.urllib.request, "urlopen", failing_urlopen)

    out = tmp_path / "out.bin"
    with pytest.raises(ConnectionResetError):
        download_with_tqdm("http://x/y", "label", str(out),
                           attempts=3, retry_wait=0)

    assert not out.exists()                        # no truncated file to be trusted
    assert not (tmp_path / "out.bin.part").exists()  # temp cleaned up


def test_404_leaves_no_partial_file(tmp_path, monkeypatch):
    """A definitive 404 (never retried) also leaves nothing behind."""
    def urlopen_404(target, *a, **k):
        is_stream = isinstance(target, urllib.request.Request)
        if not is_stream:
            return _FakeResp(headers={})
        raise _http_error(code=404)

    monkeypatch.setattr(general.urllib.request, "urlopen", urlopen_404)

    out = tmp_path / "out.bin"
    with pytest.raises(urllib.error.HTTPError):
        download_with_tqdm("http://x/y", "label", str(out))

    assert not out.exists()
    assert not (tmp_path / "out.bin.part").exists()
