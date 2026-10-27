# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/features.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""Biomarkers of an action potential, the same for the reference and the schemes

``vpeak``, ``cat`` and ``tact`` are those of
``benchmarks/block_rush_larsen/features.py``, imported unchanged, so that the
survey and the benchmark of one model measure them alike, and so is
``apd90_grid``, its APD90 from grid point to grid point. The grid is 0.1 ms, or
the time step when that is longer, and at a time step of 1 ms the grid APD90
moves in steps of 1 ms, so the survey's ``apd90`` and ``apd50`` interpolate
both ends linearly instead, as the draft's Table 2 does:

* ``apd50``, ``apd90``: from ``tact`` to the first time after the peak of V where
  it has repolarized by 50 % or 90 % of ``vpeak - vrest``, linearly interpolated

* ``vrest``: V at the start of the trace (the beat starts in diastole)
* ``dvdt_max``: the largest forward difference of V
* ``ca_dia``: Ca at the start of the trace
* ``cat90``: from ``tact`` to the first time after the peak of Ca where it has
  decayed by 90 % of ``cat``, linearly interpolated
* ``nai_dia``: Na_i at the start of the trace

For autorhythmic models, :func:`autorhythmic_features` takes the last complete
cycle instead, with the maximum diastolic potential ``mdp`` in place of ``vrest``,
and the mean ``cycle_length``.
"""

from __future__ import annotations

import numpy as np
from .features_base import features as _brl_features

NAN = float("nan")


def _crossing_down(t: np.ndarray, y: np.ndarray, start: int, level: float) -> float:
    """The first time after index `start` where `y` falls below `level`, linearly interpolated"""
    after = np.where(y[start:] < level)[0]
    if not len(after):
        return NAN
    i = start + after[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def _calcium(t: np.ndarray, Ca: np.ndarray | None, tact: float, base: float | None = None) -> dict:
    if Ca is None:
        return dict(cat=NAN, ca_dia=NAN, cat90=NAN)
    base = float(Ca[0]) if base is None else base
    i_peak = int(np.nanargmax(Ca))
    cat = float(Ca[i_peak] - base)
    level = float(Ca[i_peak]) - 0.9 * cat
    return dict(cat=cat, ca_dia=base, cat90=_crossing_down(t, Ca, i_peak, level) - tact)


def paced_features(
    t: np.ndarray, V: np.ndarray, Ca: np.ndarray | None = None, Na: np.ndarray | None = None
) -> dict[str, float]:
    """The biomarkers of one paced beat that starts in diastole"""
    base = _brl_features(t, V, Ca if Ca is not None else np.zeros_like(V))
    i_peak = int(np.nanargmax(V))
    vpeak, vrest = float(V[i_peak]), float(V[0])
    out = dict(
        vpeak=base["vpeak"],
        vrest=vrest,
        dvdt_max=float(np.nanmax(np.diff(V) / np.diff(t))),
        tact=base["tact"],
        apd50=_crossing_down(t, V, i_peak, vpeak - 0.5 * (vpeak - vrest)) - base["tact"],
        apd90=_crossing_down(t, V, i_peak, vpeak - 0.9 * (vpeak - vrest)) - base["tact"],
        apd90_grid=base["apd90"],
    )
    out.update(_calcium(t, Ca, base["tact"]))
    if Ca is not None:
        out["cat"] = base["cat"]
    out["nai_dia"] = float(Na[0]) if Na is not None else NAN
    return out


def _upstrokes(t: np.ndarray, V: np.ndarray, level: float) -> np.ndarray:
    """The times of the upward crossings of `level`, linearly interpolated"""
    i = np.where((V[:-1] < level) & (V[1:] >= level))[0] + 1
    return t[i - 1] + (level - V[i - 1]) * (t[i] - t[i - 1]) / (V[i] - V[i - 1])


def autorhythmic_features(
    t: np.ndarray, V: np.ndarray, Ca: np.ndarray | None = None, Na: np.ndarray | None = None
) -> dict[str, float]:
    """The biomarkers of the last complete cycle of a spontaneously beating trace

    The cycle is the one between the second to last and the last upstroke (the
    upward crossing of the midpoint between the minimum and the maximum of V),
    and its MDP is the minimum of V since the upstroke before it
    """
    vmid = 0.5 * (float(np.nanmin(V)) + float(np.nanmax(V)))
    ups = _upstrokes(t, V, vmid)
    out = dict(n_up=float(len(ups)), cycle_length=NAN)
    keys = (
        "vpeak",
        "mdp",
        "dvdt_max",
        "tact",
        "apd50",
        "apd90",
        "cat",
        "ca_dia",
        "cat90",
        "nai_dia",
    )
    out.update({key: NAN for key in keys})
    if len(ups) < 2:
        return out
    out["cycle_length"] = float(np.mean(np.diff(ups)))
    t_up, t_next = ups[-2], ups[-1]
    t_prev = ups[-3] if len(ups) >= 3 else t[0]
    before = (t >= t_prev) & (t < t_up)
    cycle = (t >= t_up) & (t < t_next)
    i0 = int(np.argmax(cycle))
    tc, Vc = t[cycle], V[cycle]
    mdp = float(np.min(V[before])) if before.any() else float(V[0])
    i_peak = int(np.argmax(Vc))
    vpeak = float(Vc[i_peak])
    slope = np.diff(V) / np.diff(t)
    out.update(
        vpeak=vpeak,
        mdp=mdp,
        # the steepest step of the cycle, the one into its first point included
        dvdt_max=float(np.max(slope[max(i0 - 1, 0) : i0 + len(tc) - 1])),
        tact=float(t_up),
        apd50=_crossing_down(tc, Vc, i_peak, vpeak - 0.5 * (vpeak - mdp)) - t_up,
        apd90=_crossing_down(tc, Vc, i_peak, vpeak - 0.9 * (vpeak - mdp)) - t_up,
    )
    if Ca is not None:
        out.update(_calcium(tc, Ca[cycle], float(t_up), base=float(np.min(Ca[cycle]))))
    if Na is not None:
        out["nai_dia"] = float(np.interp(t_up, t, Na))
    return out


def feature_errors(actual: dict[str, float], expected: dict[str, float]) -> dict[str, float]:
    """The errors of `actual` against `expected`: absolute, but relative for ``cat``"""
    errors = {}
    for key in expected:
        if key not in actual or key == "n_up":
            continue
        if key == "cat":
            errors[key] = (actual[key] - expected[key]) / expected[key]
        else:
            errors[key] = actual[key] - expected[key]
    return errors
