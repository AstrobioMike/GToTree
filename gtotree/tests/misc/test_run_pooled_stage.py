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

from gtotree.utils.misc.general import run_pooled_stage


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


# ---------------------------------------------------------------------------
# bounded submission
# ---------------------------------------------------------------------------

def test_submission_is_windowed_not_all_up_front():
    """
    Work is submitted in a rolling window, so queued futures and their results stay
    O(window) rather than growing with the number of items.
    """
    import concurrent.futures as cf
    from gtotree.utils.misc.general import WORKER_QUEUE_DEPTH

    workers = 2
    items = 500
    window = workers * WORKER_QUEUE_DEPTH

    submitted = [0]
    at_first_apply = []
    real_submit = cf.ThreadPoolExecutor.submit

    def counting_submit(self, fn, *a, **k):
        submitted[0] += 1
        return real_submit(self, fn, *a, **k)

    def apply_result(item, result, rd):
        if not at_first_apply:
            at_first_apply.append(submitted[0])

    cf.ThreadPoolExecutor.submit = counting_submit
    try:
        run_pooled_stage(list(range(items)), lambda i, rd: i, apply_result,
                         _args(workers), None)
    finally:
        cf.ThreadPoolExecutor.submit = real_submit

    assert at_first_apply, "apply_result never ran"
    assert at_first_apply[0] <= window, (
        f"{at_first_apply[0]} items were submitted before the first result was applied; "
        f"the window is {window} -- submission is not bounded")
    # and everything must still get submitted eventually
    assert submitted[0] == items


def test_every_item_is_processed_exactly_once():
    """Topping up the window neither drops nor duplicates work."""
    items = list(range(997))          # deliberately not a multiple of the window
    applied = []

    def apply_result(item, result, rd):
        applied.append(result)

    run_pooled_stage(items, lambda i, rd: i, apply_result, _args(4), None)

    assert sorted(applied) == items


def test_a_slow_item_does_not_stall_the_window():
    """One long-running item does not hold back the rest."""
    def worker(item, rd):
        time.sleep(0.25 if item == 0 else 0.001)
        return item

    applied = []
    start = time.time()
    run_pooled_stage(list(range(60)), worker, lambda i, r, d: applied.append(i),
                     _args(4), None)
    elapsed = time.time() - start

    assert len(applied) == 60
    # if the slow item blocked top-ups, the other 59 would queue behind it
    assert elapsed < 1.0, f"took {elapsed:.2f}s -- the window appears to stall"


# ---------------------------------------------------------------------------
# context propagation into workers
# ---------------------------------------------------------------------------

def test_workers_see_the_callers_contextvars():
    """
    A bare thread starts with an empty context, so a ContextVar set on the calling
    thread would read back as its default inside a worker. `log_file_var` is the one
    that matters: preflight sets it to an absolute path inside the output directory,
    and a worker seeing the relative default would send any reporter output to the
    user's cwd instead.
    """
    from gtotree.utils.misc.context import log_file_var

    wanted = "/tmp/some-output-dir/gtotree-runlog.txt"
    token = log_file_var.set(wanted)
    try:
        seen = []
        run_pooled_stage(list(range(20)),
                         lambda item, rd: log_file_var.get(),
                         lambda item, result, rd: seen.append(result),
                         _args(4), None)
    finally:
        log_file_var.reset(token)

    assert set(seen) == {wanted}, \
        "workers fell back to the ContextVar default instead of the caller's value"


def test_worker_contextvar_writes_do_not_leak_back_to_the_caller():
    """Each worker gets its own copy, so a stray set() inside one stays contained."""
    from gtotree.utils.misc.context import log_file_var

    token = log_file_var.set("caller-value")
    try:
        run_pooled_stage(list(range(10)),
                         lambda item, rd: log_file_var.set(f"worker-{item}"),
                         lambda item, result, rd: None,
                         _args(4), None)
        assert log_file_var.get() == "caller-value"
    finally:
        log_file_var.reset(token)


def test_concurrent_workers_do_not_share_one_context():
    """
    A Context cannot be entered re-entrantly, so reusing a single copy across
    concurrent submissions would raise RuntimeError rather than run.
    """
    barrier = threading.Barrier(4, timeout=5)

    def worker(item, rd):
        barrier.wait()          # forces genuine overlap
        return item

    applied = []
    run_pooled_stage(list(range(16)), worker,
                     lambda i, r, rd: applied.append(r), _args(4), None)

    assert sorted(applied) == list(range(16))
