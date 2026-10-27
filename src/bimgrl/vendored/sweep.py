# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/sweep.py):
# only bisect_k; the rest of the sweep is rewritten in bimgrl/sweep.py.
from __future__ import annotations

from collections.abc import Callable


def bisect_k(passes_k: Callable[[int], bool], k_pass: int, k_fail: int, allowed: Callable[[int], bool]) -> int:
    """The smallest allowed k in (k_fail, k_pass] that passes, for dt = base / k, assuming that
    every k above the first that passes passes too. `k_pass` passes and `k_fail` does not"""
    candidates = [k for k in range(k_fail + 1, k_pass) if allowed(k)] + [k_pass]
    lo, hi = 0, len(candidates) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if passes_k(candidates[mid]):
            hi = mid
        else:
            lo = mid + 1
    return candidates[lo]
