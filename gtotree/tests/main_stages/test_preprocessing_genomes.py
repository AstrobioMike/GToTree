import gzip
import http.server
import io
import socket
import threading
import time
import urllib.error
from unittest.mock import patch

import pytest  # type: ignore

import gtotree.main_stages.preprocessing_genomes as P
from gtotree.main_stages.preprocessing_genomes import (
    _sleep_backoff,
    download_and_unzip_accession,
    NCBI_MAX_BACKOFF,
    NCBI_MAX_RETRY_AFTER,
)


def _http_error(code, headers=None):
    return urllib.error.HTTPError("http://example/x", code, "msg",
                                  headers or {}, io.BytesIO(b""))


# ---------------------------------------------------------------------------
# _sleep_backoff -- policy selection
# ---------------------------------------------------------------------------

def test_sawtooth_resets_each_cycle():
    # non-throttle sleeps cycle 1, 2, 4, 8, 16 then start over, so a straggler never
    # parks a worker thread on an ever-growing sleep
    seen = []
    with patch.object(P.time, "sleep", side_effect=seen.append), \
         patch.object(P.random, "uniform", return_value=0.0):
        for attempt in range(1, 12):
            _sleep_backoff(attempt)
    assert seen == [1.0, 2.0, 4.0, 8.0, 16.0, 1.0, 2.0, 4.0, 8.0, 16.0, 1.0]


def test_throttled_without_header_is_exponential_and_capped():
    seen = []
    with patch.object(P.time, "sleep", side_effect=seen.append), \
         patch.object(P.random, "uniform", return_value=0.0):
        for attempt in range(1, 9):
            _sleep_backoff(attempt, throttled=True)
    assert seen[:5] == [1.0, 2.0, 4.0, 8.0, 16.0]
    assert all(s == NCBI_MAX_BACKOFF for s in seen[5:])


def test_throttled_honors_retry_after():
    with patch.object(P.time, "sleep") as mock_sleep:
        _sleep_backoff(1, err=_http_error(429, {"Retry-After": "12"}), throttled=True)
    mock_sleep.assert_called_once_with(12.0)


def test_retry_after_honored_beyond_max_backoff():
    # Retry-After gets the larger ceiling: a server saying "come back in 120s" is
    # obeyed rather than clamped to NCBI_MAX_BACKOFF, which would just burn retries
    # on requests we already know will be refused
    with patch.object(P.time, "sleep") as mock_sleep:
        _sleep_backoff(1, err=_http_error(429, {"Retry-After": "120"}), throttled=True)
    mock_sleep.assert_called_once_with(120.0)


def test_absurd_retry_after_is_capped():
    with patch.object(P.time, "sleep") as mock_sleep:
        _sleep_backoff(1, err=_http_error(429, {"Retry-After": "99999"}), throttled=True)
    mock_sleep.assert_called_once_with(NCBI_MAX_RETRY_AFTER)


def test_invalid_retry_after_falls_back_to_exponential():
    with patch.object(P.time, "sleep") as mock_sleep, \
         patch.object(P.random, "uniform", return_value=0.0):
        _sleep_backoff(3, err=_http_error(429, {"Retry-After": "not-a-number"}),
                       throttled=True)
    mock_sleep.assert_called_once_with(4.0)      # 2 ** (3 - 1)


def test_retry_after_ignored_when_not_throttled():
    # a plain 5xx/connection blip takes the sawtooth path, which does not consult
    # Retry-After at all (only an explicit throttle does)
    with patch.object(P.time, "sleep") as mock_sleep, \
         patch.object(P.random, "uniform", return_value=0.0):
        _sleep_backoff(1, err=_http_error(503, {"Retry-After": "600"}))
    mock_sleep.assert_called_once_with(1.0)


# ---------------------------------------------------------------------------
# download_and_unzip_accession -- classification of failures
# ---------------------------------------------------------------------------

def _classify(exc, max_retries=3):
    """run the retry loop against a always-failing fetch, capturing throttled flags

    Patches `_fetch_to_file`, which replaced the bare `urllib.request.urlretrieve` call
    (urlretrieve accepts no timeout, so a stalled server hung a pool thread forever).
    Same seam, one layer in.
    """
    seen = []
    def fake_backoff(attempt, err=None, throttled=False):
        seen.append(throttled)
    with patch.object(P, "_fetch_to_file", side_effect=exc), \
         patch.object(P, "_sleep_backoff", fake_backoff):
        try:
            download_and_unzip_accession("http://example/x", "/tmp/gtt_backoff_test",
                                         max_retries=max_retries)
        except Exception:
            pass
    return seen


def test_bare_429_is_treated_as_throttle():
    assert _classify(_http_error(429)) == [True, True]


def test_bare_5xx_is_not_treated_as_throttle():
    assert _classify(_http_error(503)) == [False, False]


def test_5xx_with_retry_after_is_treated_as_throttle():
    # a server that bothered to say when to come back is throttling us, even if it
    # didn't use a 429 to say so
    assert _classify(_http_error(503, {"Retry-After": "7"})) == [True, True]


def test_connection_error_is_not_treated_as_throttle():
    assert _classify(urllib.error.URLError("connection reset")) == [False, False]


def test_404_fails_fast_without_retrying():
    seen = _classify(_http_error(404), max_retries=5)
    assert seen == []        # permanent -> raised immediately, no backoff at all


# ---------------------------------------------------------------------------
# download timeouts
# ---------------------------------------------------------------------------

def _loopback_is_usable(timeout=2.0):
    """
    Can this host connect to a socket it just bound on 127.0.0.1?
    """
    srv = socket.socket()
    try:
        srv.bind(("127.0.0.1", 0))
        srv.listen(1)
        with socket.create_connection(srv.getsockname(), timeout=timeout):
            return True
    except OSError:
        return False
    finally:
        srv.close()


needs_loopback = pytest.mark.skipif(
    not _loopback_is_usable(),
    reason="cannot connect to a local listening socket -- loopback TCP is blocked on "
           "this host (endpoint security agent, packet filter, or lo0 misconfigured)",
)


class _TimeoutHandler(http.server.BaseHTTPRequestHandler):
    """Serves /ok as gzip; /stall sends headers then goes quiet."""

    payload = b">seq1\nMSEQVENCE\n"
    stall_seconds = 30

    def log_message(self, *a):
        pass

    def do_GET(self):
        if self.path == "/stall":
            self.send_response(200)
            self.send_header("Content-Length", "1000000")
            self.end_headers()
            threading.Event().wait(self.stall_seconds)
            return

        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            gz.write(self.payload)
        body = buf.getvalue()
        self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


@pytest.fixture
def stalling_server():
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _TimeoutHandler)
    srv.daemon_threads = True
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()
    srv.server_close()


@needs_loopback
def test_download_fetches_and_gunzips(stalling_server, tmp_path, monkeypatch):
    monkeypatch.setattr(P, "NCBI_DOWNLOAD_TIMEOUT", 5)

    dest = tmp_path / "genome.faa"
    P.download_and_unzip_accession(f"{stalling_server}/ok", str(dest), max_retries=1)

    assert dest.read_bytes() == _TimeoutHandler.payload
    assert sorted(p.name for p in tmp_path.iterdir()) == ["genome.faa"]


@needs_loopback
def test_an_unresponsive_server_times_out_and_is_retried(stalling_server, tmp_path,
                                                         monkeypatch):
    """
    A server that accepts the connection then stops sending must raise, so the retry
    loop can act on it, rather than holding the thread indefinitely.
    """
    monkeypatch.setattr(P, "NCBI_DOWNLOAD_TIMEOUT", 1)
    backoffs = []
    monkeypatch.setattr(P, "_sleep_backoff", lambda attempt, **k: backoffs.append(attempt))

    dest = tmp_path / "genome.faa"
    started = time.monotonic()
    with pytest.raises(OSError):
        P.download_and_unzip_accession(f"{stalling_server}/stall", str(dest),
                                       max_retries=3)
    elapsed = time.monotonic() - started

    assert elapsed < 15, f"took {elapsed:.1f}s -- the timeout is not being applied"
    assert backoffs == [1, 2], "a timeout should be retried like any transient failure"
    assert sorted(p.name for p in tmp_path.iterdir()) == [], "partial files left behind"
