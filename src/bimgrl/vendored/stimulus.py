# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/stimulus.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""The stimulus of the survey: a half-open interval in every period

Every paced model gets the same stimulus in place of its own, on the
half-open interval ``[start + k period, start + k period + duration)``::

    phase = t - stim_start
    phase_mod = phase - stim_period * floor((phase + stim_eps) / stim_period)
    i_stim = piecewise(phase >= -stim_eps and phase_mod < stim_duration - stim_eps,
                       stim_amplitude, 0)

The time of step ``n`` is ``n * dt`` from an integer counter, and the time
steps divide the start, the duration and the period, so that a step is
stimulated if and only if it starts inside the interval, and exactly
``duration / dt`` steps are stimulated in every period. ``stim_eps`` absorbs
the rounding of ``n * dt``. A closed interval with an accumulated time stimulates
one step too many in some periods, which is the error that TWorld.jl found.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np

from ..meta import Stim

STIM_EPS = 1e-9

#: The time steps of the sweep, in ms: 1e-4 to 1 as in the brief, extended down to 1e-5 for the
#: models where a scheme still fails at 1e-4 (GRL in Winslow 1999 loses the occupancy of its
#: Markov chains); the budget of a run limits how far down a scheme gets
DT_GRID = (
    0.00001,
    0.00002,
    0.00005,
    0.0001,
    0.0002,
    0.0005,
    0.001,
    0.002,
    0.005,
    0.01,
    0.02,
    0.05,
    0.1,
    0.2,
    0.5,
    1.0,
)

#: The names of the parameters of the stimulus
PARAMETERS = ("stim_start", "stim_duration", "stim_period", "stim_amplitude", "stim_eps")


def rewrite_stimulus(model, variable: str, stim: Stim) -> list[str]:
    """Replace the stimulus current `variable` of a myokit model with the half-open one

    Parameters
    ----------
    model : myokit.Model
        The model, in ms. It is changed in place
    variable : str
        The qualified name of the stimulus current, such as ``"membrane.i_stim"``
    stim : Stim
        The stimulus. The amplitude is the value of the current when it is on

    Returns
    -------
    list[str]
        The edits, for ``meta.toml`` and ``NOTES.md``
    """
    import myokit

    var = model.get(variable)
    component = var.parent(myokit.Component)
    time = model.time()
    time_unit = time.unit()
    old = var.rhs().code()

    values = dict(
        stim_start=(stim.start, time_unit),
        stim_duration=(stim.duration, time_unit),
        stim_period=(stim.period, time_unit),
        stim_amplitude=(stim.amplitude, var.unit()),
        stim_eps=(STIM_EPS, time_unit),
    )
    for name, (value, unit) in values.items():
        if component.can_add_variable(name):
            new = component.add_variable(name)
        else:
            new = component.get(name)
        # a Number, not a string: "-52.0" would be parsed as an expression, which
        # gotranx would write as an intermediate instead of a parameter
        new.set_rhs(myokit.Number(float(value), unit))
        new.set_unit(unit)
    for name in ("stim_phase", "stim_phase_mod"):
        if component.can_add_variable(name):
            component.add_variable(name).set_unit(time_unit)
    t = time.qname(component)
    component.get("stim_phase").set_rhs(f"{t} - stim_start")
    component.get("stim_phase_mod").set_rhs(
        "stim_phase - stim_period * floor((stim_phase + stim_eps) / stim_period)"
    )
    var.set_rhs(
        "piecewise(stim_phase >= -stim_eps and stim_phase_mod < stim_duration - stim_eps, "
        "stim_amplitude, 0)"
    )
    return [
        f"stimulus {variable} rewritten to the half-open interval "
        f"[start + k period, start + k period + duration) with start {stim.start:g} ms, "
        f"duration {stim.duration:g} ms, period {stim.period:g} ms, amplitude {stim.amplitude:g} "
        f"(was: {old})"
    ]


def aligned(dt: float, *values: float) -> bool:
    """Whether `dt` divides all of `values` (up to rounding, relative to the number of steps)"""
    return all(abs(x / dt - round(x / dt)) < 1e-9 * max(1.0, abs(x / dt)) for x in values)


def dt_grid(stim: Stim | None, window: float, grid: Sequence[float] = DT_GRID) -> list[float]:
    """The time steps of `grid` that divide the stimulus (start, duration, period) and the window"""
    values = [window] if stim is None else [window, stim.start, stim.duration, stim.period]
    return [dt for dt in grid if aligned(dt, *values)]


def breakpoints(stim: Stim | None, t0: float, t1: float) -> list[float]:
    """The times in (t0, t1) where the stimulus switches on or off, sorted"""
    if stim is None or stim.amplitude == 0:
        return []
    first = math.floor((t0 - stim.start) / stim.period) - 1
    times = []
    k = first
    while stim.start + k * stim.period < t1:
        for t in (stim.start + k * stim.period, stim.start + k * stim.period + stim.duration):
            if t0 < t < t1:
                times.append(t)
        k += 1
    return sorted(times)


def forced_parameters(p: np.ndarray, index: Callable[[str], int], on: bool) -> np.ndarray:
    """A copy of the parameters with the stimulus always on (`on`) or always off

    Used for the reference, which is integrated segment by segment between the
    breakpoints, so that the conditional never switches inside a segment
    """
    p = np.array(p, dtype=float, copy=True)
    if on:
        p[index("stim_start")] = -1e9
        p[index("stim_duration")] = 1e12
        p[index("stim_period")] = 1e13
    else:
        p[index("stim_amplitude")] = 0.0
    return p
