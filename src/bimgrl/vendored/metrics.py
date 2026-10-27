# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/metrics.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""Errors of a run against the reference, and whether the run is physical"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

#: How far a gate may go outside [0, 1], and a concentration below 0, in its own unit
GATE_TOL = 1e-6
CONC_TOL = 1e-9


def rrms_sd(V: np.ndarray, Vref: np.ndarray) -> float:
    """The RRMS of Spiteri & Dean (2010) as printed (p. 3597):
    sqrt((1/N) sum (V - Vref)^2 / sum Vref^2). The 1/N inside the root makes it
    the relative RMS divided by sqrt(N)"""
    V, Vref = np.asarray(V), np.asarray(Vref)
    return float(np.sqrt(np.sum((V - Vref) ** 2) / np.sum(Vref**2) / len(Vref)))


def relrms(V: np.ndarray, Vref: np.ndarray) -> float:
    """sqrt(sum (V - Vref)^2 / sum Vref^2)"""
    V, Vref = np.asarray(V), np.asarray(Vref)
    return float(np.sqrt(np.sum((V - Vref) ** 2) / np.sum(Vref**2)))


def max_abs_dv(V: np.ndarray, Vref: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(V) - np.asarray(Vref))))


def sd_points(t0: float, T: float, n: int = 100) -> np.ndarray:
    """The `n` equally spaced times of S&D's RRMS: t0 + i T / n, i = 1..n"""
    return t0 + np.arange(1, n + 1) * T / n


def state_kinds(
    ode, voltage: str, exempt: Sequence[str] = (), order: Sequence[str] | None = None
) -> list[str]:
    """``voltage``, ``gate`` (dimensionless), ``conc`` (a unit with M or mol) or ``other`` for
    each state, in the given `order` of names (the columns of the generated code, which
    are not the order of ``ode.states``); the `exempt` states are ``other``"""
    by_name = {state.name: state for state in ode.states}
    kinds = []
    for name in order if order is not None else by_name:
        state = by_name[name]
        unit = (state.unit_str or "").strip()
        if state.name == voltage:
            kinds.append("voltage")
        elif state.name in exempt:
            kinds.append("other")
        elif unit in ("", "1", "dimensionless"):
            kinds.append("gate")
        elif "M" in unit or "mol" in unit:
            kinds.append("conc")
        else:
            kinds.append("other")
    return kinds


def bounds(
    kinds: Sequence[str], ref_min: np.ndarray, ref_max: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The physical range of every state: [0, 1] for a gate and [0, inf) for a concentration,
    widened to what the reference itself visits (a published initial value can be outside,
    e.g. X_to1 = 2.63 in Winslow 1999), with the tolerances; the rest is not checked"""
    kinds = np.asarray(kinds)
    lo = np.full(len(kinds), -np.inf)
    hi = np.full(len(kinds), np.inf)
    gate, conc = kinds == "gate", kinds == "conc"
    lo[gate] = np.minimum(0.0, ref_min[gate]) - GATE_TOL
    hi[gate] = np.maximum(1.0, ref_max[gate]) + GATE_TOL
    lo[conc] = np.minimum(0.0, ref_min[conc]) - CONC_TOL
    return lo, hi


def classify(
    states: np.ndarray,
    kinds: Sequence[str],
    lower: np.ndarray | None = None,
    upper: np.ndarray | None = None,
) -> str:
    """``nan`` if a state is not finite, ``unphysical`` if a state is outside its range (by
    default a gate outside [0, 1] or a negative concentration, beyond the tolerances, see
    :func:`bounds`), and ``ok`` otherwise

    `states` has one row per time and one column per state
    """
    states = np.atleast_2d(states)
    if not np.all(np.isfinite(states)):
        return "nan"
    if lower is None or upper is None:
        n = states.shape[1]
        lower, upper = bounds(kinds, np.zeros(n), np.zeros(n))
    if np.any(states < lower) or np.any(states > upper):
        return "unphysical"
    return "ok"
