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
        # Every open is now a streaming Request open (the old up-front size
        # probe was removed; total is read off the stream response instead).
        # The non-Request branch is retained only as a defensive no-op.
        is_stream = isinstance(target, urllib.request.Request)
        if not is_stream:
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
        raise TimeoutError("nope")

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


# --------------------------------------------------------------------------- #
# _stream_once: the bar's total is resolved from the download response's own
# Content-Length header (no separate probe request). Whether or not the server
# reports a length, the full payload must land on disk.
# --------------------------------------------------------------------------- #

def test_stream_once_writes_full_payload_known_length(tmp_path, monkeypatch):
    """A Content-Length response: total is read off the response, the bar is
    promoted to bounded, and every byte is streamed to disk."""
    payload = b"z" * 300_000  # larger than a single 256 KB read
    dest = str(tmp_path / "out.bin")

    def fake_urlopen(target, *a, **k):
        return _FakeResp(headers={"Content-Length": str(len(payload))},
                         chunks=[payload[:256 * 1024], payload[256 * 1024:]])

    monkeypatch.setattr(general.urllib.request, "urlopen", fake_urlopen)

    general._stream_once("http://x/y", dest, "label", leave=False,
                         floor_bytes_per_s=0.0, probe_seconds=5.0)

    from pathlib import Path
    assert Path(dest).read_bytes() == payload
    assert not (tmp_path / "out.bin.part").exists()


def test_stream_once_writes_full_payload_unknown_length(tmp_path, monkeypatch):
    """No Content-Length (total stays None): the whole payload is still written
    via the len(buf) bar-update branch, and the bar stays indeterminate."""
    payload = b"q" * 300_000
    dest = str(tmp_path / "out.bin")

    def fake_urlopen(target, *a, **k):
        # headers carry NO Content-Length
        return _FakeResp(headers={},
                         chunks=[payload[:256 * 1024], payload[256 * 1024:]])

    monkeypatch.setattr(general.urllib.request, "urlopen", fake_urlopen)

    general._stream_once("http://x/y", dest, "label", leave=False,
                         floor_bytes_per_s=0.0, probe_seconds=5.0)

    from pathlib import Path
    assert Path(dest).read_bytes() == payload
    assert not (tmp_path / "out.bin.part").exists()


# ---------------------------------------------------------------------------
# atomic_write_text
# ---------------------------------------------------------------------------

class TestAtomicWriteText:
    """Writes via a .part file that is only moved into place once the write succeeds."""

    def test_writes_the_content(self, tmp_path):
        target = tmp_path / "out.txt"
        general.atomic_write_text(str(target), lambda f: f.write("hello"))
        assert target.read_text() == "hello"

    def test_a_failed_write_creates_no_target_and_leaves_no_part(self, tmp_path):
        target = tmp_path / "out.txt"

        def boom(f):
            f.write("half a file")
            raise RuntimeError("interrupted mid-write")

        with pytest.raises(RuntimeError):
            general.atomic_write_text(str(target), boom)

        assert not target.exists()
        assert not (tmp_path / "out.txt.part").exists()

    def test_a_failed_write_preserves_the_previous_version(self, tmp_path):
        target = tmp_path / "out.txt"
        general.atomic_write_text(str(target), lambda f: f.write("good version"))

        def boom(f):
            f.write("corrupt")
            raise RuntimeError("interrupted")

        with pytest.raises(RuntimeError):
            general.atomic_write_text(str(target), boom)

        assert target.read_text() == "good version"

    def test_cleans_up_after_a_keyboardinterrupt(self, tmp_path):
        target = tmp_path / "out.txt"

        def boom(f):
            raise KeyboardInterrupt

        with pytest.raises(KeyboardInterrupt):
            general.atomic_write_text(str(target), boom)
        assert not (tmp_path / "out.txt.part").exists()


# ---------------------------------------------------------------------------
# run data persistence
# ---------------------------------------------------------------------------

class TestRunDataPersistence:

    def _run_data(self, tmp_path):
        rd = general.RunData()
        rd.run_data_path = str(tmp_path / "run-data.json")
        rd.ncbi_accs = [general.GenomeData.from_acc("GCF_000000001.1")]
        rd.update_all_input_genomes()
        rd.fingerprint = {"hmm": "deadbeef"}
        return rd

    def test_round_trips(self, tmp_path):
        rd = self._run_data(tmp_path)
        general.write_run_data(rd)
        back = general.read_run_data(rd.run_data_path)
        assert back.fingerprint == {"hmm": "deadbeef"}
        assert [g.id for g in back.ncbi_accs] == ["GCF_000000001.1"]

    def test_reading_a_missing_file_returns_none(self, tmp_path):
        assert general.read_run_data(str(tmp_path / "nope.json")) is None

    def test_reading_a_damaged_file_raises_corrupt_run_data(self, tmp_path):
        path = tmp_path / "run-data.json"
        path.write_text('{"ncbi_accs": [{"id": "GCF_1", "source": "acc"')
        with pytest.raises(general.CorruptRunData):
            general.read_run_data(str(path))

    def test_an_interrupted_write_preserves_the_loadable_previous_state(
            self, tmp_path, monkeypatch):
        """The resume checkpoint is written atomically."""
        rd = self._run_data(tmp_path)
        general.write_run_data(rd)

        def exploding_dump(obj, fh, **kw):
            fh.write('{"partial": ')
            raise RuntimeError("interrupted mid-dump")

        monkeypatch.setattr(general.json, "dump", exploding_dump)
        rd.fingerprint = {"hmm": "bad"}
        with pytest.raises(RuntimeError):
            general.write_run_data(rd)
        monkeypatch.undo()

        assert not (tmp_path / "run-data.json.part").exists()
        assert general.read_run_data(rd.run_data_path).fingerprint == {"hmm": "deadbeef"}


@pytest.mark.parametrize("value,expected", [
    (b"PF00001", "PF00001"),
    ("PF00001", "PF00001"),
    (bytearray(b"MOCK"), "MOCK"),
    (None, None),
])
def test_decode_pyhmmer_text_normalizes_to_str(value, expected):
    """pyhmmer returns name/accession as bytes or str depending on version."""
    assert general.decode_pyhmmer_text(value) == expected
