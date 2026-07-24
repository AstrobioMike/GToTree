import io
import urllib.error
from unittest.mock import patch

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
    """run the retry loop against a always-failing fetch, capturing throttled flags"""
    seen = []
    def fake_backoff(attempt, err=None, throttled=False):
        seen.append(throttled)
    with patch.object(P.urllib.request, "urlretrieve", side_effect=exc), \
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
