import dataclasses
import json
import socket
import urllib.error
import urllib.request
from pathlib import Path
import pytest # type: ignore

from gtotree.utils.misc import general
from gtotree.utils.misc.general import download_with_tqdm


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

    def test_derived_fields_are_not_written(self, tmp_path):
        """
        all_input_genomes aliases the source lists, so writing it would put every
        genome in the file twice.
        """
        rd = self._run_data(tmp_path)
        general.write_run_data(rd)
        on_disk = json.loads(Path(rd.run_data_path).read_text())
        for name in general.DERIVED_RUN_DATA_FIELDS:
            assert name not in on_disk
        assert [gd["id"] for gd in on_disk["ncbi_accs"]] == ["GCF_000000001.1"]

    def test_every_non_derived_field_is_written(self, tmp_path):
        """Skipping derived fields must not silently drop anything else."""
        rd = self._run_data(tmp_path)
        general.write_run_data(rd)
        on_disk = json.loads(Path(rd.run_data_path).read_text())
        expected = {f.name for f in dataclasses.fields(rd)} - general.DERIVED_RUN_DATA_FIELDS
        assert set(on_disk) == expected

    def test_all_input_genomes_is_rebuilt_on_read(self, tmp_path):
        rd = self._run_data(tmp_path)
        rd.amino_acid_files = [general.GenomeData.from_path("/in/mock.faa", "amino-acid-fasta")]
        rd.update_all_input_genomes()
        general.write_run_data(rd)

        back = general.read_run_data(rd.run_data_path)
        assert [gd.id for gd in back.all_input_genomes] == ["GCF_000000001.1", "mock"]

    def test_rebuilt_all_input_genomes_aliases_the_source_lists(self, tmp_path):
        """
        Mutating a genome through either view has to be visible through the other,
        or per-genome progress recorded by a stage wouldn't reach the next write.
        """
        rd = self._run_data(tmp_path)
        general.write_run_data(rd)

        back = general.read_run_data(rd.run_data_path)
        assert back.all_input_genomes[0] is back.ncbi_accs[0]

        back.all_input_genomes[0].mark_hmm_search_done()
        assert back.ncbi_accs[0].hmm_search_done is True

    def test_a_legacy_file_carrying_derived_fields_still_loads(self, tmp_path):
        """run-data.json written before derived fields were dropped must still resume."""
        rd = self._run_data(tmp_path)
        general.write_run_data(rd)

        on_disk = json.loads(Path(rd.run_data_path).read_text())
        # a stale duplicate, deliberately disagreeing with the source list
        on_disk["all_input_genomes"] = [
            {"id": "STALE", "source": "ncbi-accession", "full_path": None,
             "provided_path": None, "basename": "STALE"}
        ]
        Path(rd.run_data_path).write_text(json.dumps(on_disk))

        back = general.read_run_data(rd.run_data_path)
        assert [gd.id for gd in back.all_input_genomes] == ["GCF_000000001.1"]

    def test_writing_does_not_mutate_the_run_data(self, tmp_path):
        """run_data_as_dict references rather than deep-copies; it must stay read-only."""
        rd = self._run_data(tmp_path)
        rd.tax_info_dict = {"GCF_000000001.1": {"domain": "Bacteria"}}
        general.write_run_data(rd)
        assert rd.tax_info_dict == {"GCF_000000001.1": {"domain": "Bacteria"}}
        assert rd.all_input_genomes[0] is rd.ncbi_accs[0]

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


# ---------------------------------------------------------------------------
# get_path_rel_to_outdir
# ---------------------------------------------------------------------------

def test_get_path_rel_to_outdir_keeps_outdir_when_basename_repeats_it(tmp_path, monkeypatch):
    """
    The final tree is named after the output directory (`<outdir>/<outdir>.tre`), so the
    directory name appears twice in its path. Re-expressing it must still keep the
    directory, not collapse to the bare filename.
    """
    from types import SimpleNamespace
    from gtotree.utils.misc.messaging import get_path_rel_to_outdir

    # args.output_dir is the user's (usually relative) spelling; run_data.output_dir is
    # its abspath, so the two only line up from the directory the run started in
    monkeypatch.chdir(tmp_path)
    args = SimpleNamespace(output_dir="GToTree-output")
    outdir = tmp_path / "GToTree-output"

    tree = str(outdir / "GToTree-output.tre")
    assert get_path_rel_to_outdir(tree, args) == "GToTree-output/GToTree-output.tre"

    # a path with no repeated component behaves the same way
    aln = str(outdir / "aligned-SCGs.faa")
    assert get_path_rel_to_outdir(aln, args) == "GToTree-output/aligned-SCGs.faa"
