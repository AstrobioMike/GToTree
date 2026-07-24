"""
Tests for the shared `run_pooled_stage` helper.

All four preprocessing stages (genbank, fasta, amino-acid, NCBI accessions) and
gen-scg-hmms route through this

  * `worker` runs concurrently and must not raise
  * `apply_result` runs single-threaded on the calling thread
  * `max_workers_cap` clamps below --num-jobs
  * a KeyboardInterrupt drops queued work rather than draining the backlog
"""

import argparse
import threading
import time

import pytest # type: ignore

from gtotree.utils.general import run_pooled_stage


def _args(num_jobs):
    return argparse.Namespace(num_jobs=num_jobs)


def test_applies_every_item():
    applied = []
    run_pooled_stage(list(range(50)),
                     lambda item, rd: item * 2,
                     lambda item, result, rd: applied.append(result),
                     _args(8), None)
    assert sorted(applied) == [i * 2 for i in range(50)]


def test_returns_run_data():
    sentinel = {"state": "mine"}
    out = run_pooled_stage([1, 2], lambda i, rd: i,
                           lambda i, r, rd: None, _args(2), sentinel)
    assert out is sentinel


def test_run_data_is_passed_to_both_callbacks():
    seen = {"worker": [], "apply": []}
    rd = {"key": "value"}
    run_pooled_stage([1],
                     lambda i, d: seen["worker"].append(d) or i,
                     lambda i, r, d: seen["apply"].append(d),
                     _args(1), rd)
    assert seen["worker"] == [rd]
    assert seen["apply"] == [rd]


def test_runs_concurrently():
    live = [0]
    peak = [0]
    lock = threading.Lock()

    def worker(item, rd):
        with lock:
            live[0] += 1
            peak[0] = max(peak[0], live[0])
        time.sleep(0.02)
        with lock:
            live[0] -= 1
        return item

    run_pooled_stage(list(range(16)), worker, lambda i, r, d: None, _args(8), None)
    assert peak[0] > 1


def test_max_workers_cap_clamps_below_num_jobs():
    live = [0]
    peak = [0]
    lock = threading.Lock()

    def worker(item, rd):
        with lock:
            live[0] += 1
            peak[0] = max(peak[0], live[0])
        time.sleep(0.02)
        with lock:
            live[0] -= 1
        return item

    run_pooled_stage(list(range(40)), worker, lambda i, r, d: None,
                     _args(100), None, max_workers_cap=4)
    assert peak[0] <= 4


def test_apply_result_runs_on_a_single_thread():
    """
    Shared appends (the combined fasta, GenomeData bookkeeping) happen in
    apply_result, so it must never run concurrently.
    """
    threads = set()

    def worker(item, rd):
        time.sleep(0.005)
        return item

    def apply_result(item, result, rd):
        threads.add(threading.current_thread().name)

    run_pooled_stage(list(range(24)), worker, apply_result, _args(8), None)
    assert len(threads) == 1


def test_worker_exception_propagates():
    """
    Documented contract: workers must not raise. If one does, the stage aborts partway
    through rather than silently dropping the item -- which is why every real worker
    wraps its body and returns a status object instead.
    """
    def worker(item, rd):
        if item == 3:
            raise RuntimeError("worker blew up")
        return item

    with pytest.raises(RuntimeError, match="worker blew up"):
        run_pooled_stage(list(range(6)), worker, lambda i, r, d: None, _args(2), None)


def test_num_jobs_floor_of_one():
    applied = []
    run_pooled_stage([1, 2, 3], lambda i, rd: i,
                     lambda i, r, d: applied.append(i), _args(0), None)
    assert sorted(applied) == [1, 2, 3]


def test_empty_item_list():
    applied = []
    run_pooled_stage([], lambda i, rd: i, lambda i, r, d: applied.append(i),
                     _args(4), None)
    assert applied == []


def test_keyboard_interrupt_cancels_queued_work():
    """
    Regression test for the interrupt fix. Without cancel_futures, a ctrl-c waits for
    EVERY queued item to finish -- at GToTree's 10k-30k genome scale that's a very long
    wait to quit. Only the in-flight workers should have to finish.
    """
    workers = 2
    items = 40
    per_item = 0.05
    drain_time = items * per_item / workers      # ~1.0s if it drains everything

    def worker(item, rd):
        time.sleep(per_item)
        return item

    def interrupt():
        time.sleep(per_item * 2)
        raise KeyboardInterrupt

    applied = []

    def apply_result(item, result, rd):
        applied.append(item)
        if len(applied) == 2:
            raise KeyboardInterrupt

    start = time.time()
    with pytest.raises(KeyboardInterrupt):
        run_pooled_stage(list(range(items)), worker, apply_result,
                         _args(workers), None)
    elapsed = time.time() - start

    assert elapsed < drain_time * 0.6, (
        f"took {elapsed:.2f}s; looks like the queued backlog was drained "
        f"(full drain would be ~{drain_time:.2f}s)")
