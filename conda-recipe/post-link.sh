#!/usr/bin/env bash
# Intentionally a no-op. This used to silence GNU parallel's citation notice
# (`parallel --citation`), but GToTree no longer depends on parallel -- all
# parallelism is handled in-process with Python's ThreadPoolExecutor -- so there is
# nothing to do at post-link time. Kept as an empty script to avoid changing the
# recipe's file layout; safe to `git rm` if you'd rather drop it entirely.
exit 0
