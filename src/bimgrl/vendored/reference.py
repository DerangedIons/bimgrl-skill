# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/reference.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""The "true" solution: scipy LSODA at tight tolerances, segment by segment

The interval is cut at every time where the stimulus switches on or off, and in
each segment the stimulus is forced on or off through its parameters (see
:func:`survey.stimulus.forced_parameters`), so that the solver never sees the
discontinuity. Every segment keeps its dense output, and the reference is
evaluated at exactly the times where a scheme is sampled. A time on a
breakpoint belongs to the segment that starts there.

The dense reference is not saved: it is recomputed once per model and protocol
in the process that needs it. :func:`save_traces` writes the signals that the
features and the figures need, on a 0.01 ms grid.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..meta import Stim
from .stimulus import STIM_EPS, breakpoints, forced_parameters


def numpy_model(ode) -> dict:
    """The generated numpy code of `ode`: ``rhs``, the initial values and the indices"""
    from gotranx.codegen import PythonCodeGenerator
    from gotranx.codegen.python import Format

    codegen = PythonCodeGenerator(ode, format=Format.none)
    source = "\n".join(
        [
            codegen.imports(),
            codegen.parameter_index(),
            codegen.state_index(),
            codegen.initial_parameter_values(),
            codegen.initial_state_values(),
            codegen.rhs(),
        ]
    )
    namespace: dict = {}
    exec(compile(source, f"<{ode.name}>", "exec"), namespace)
    return namespace


def stim_active(stim: Stim | None, t: float) -> bool:
    """Whether the half-open stimulus is on at time `t`"""
    if stim is None or stim.amplitude == 0:
        return False
    phase = t - stim.start
    phase_mod = phase - stim.period * np.floor((phase + STIM_EPS) / stim.period)
    return bool(phase >= -STIM_EPS and phase_mod < stim.duration - STIM_EPS)


@dataclass
class Segment:
    t0: float
    t1: float
    on: bool
    sol: Callable[[np.ndarray], np.ndarray]


@dataclass
class Reference:
    t0: float
    t1: float
    segments: list[Segment]
    rtol: float
    atol: float
    nfev: int
    njev: int
    wall: float
    method: str
    y0: np.ndarray = field(repr=False, default=None)

    def __call__(self, t: np.ndarray) -> np.ndarray:
        """The states at the times `t`, one column per time"""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        starts = np.array([s.t0 for s in self.segments])
        index = np.clip(np.searchsorted(starts, t, side="right") - 1, 0, len(starts) - 1)
        out = np.empty((len(self.y0), len(t)))
        for i in np.unique(index):
            mask = index == i
            out[:, mask] = self.segments[i].sol(t[mask])
        return out

    @property
    def end_state(self) -> np.ndarray:
        return self(np.array([self.t1]))[:, 0]


def compute_reference(
    ns: dict,
    p: np.ndarray,
    y0: np.ndarray,
    t0: float,
    t1: float,
    stim: Stim | None,
    method: str = "LSODA",
    rtol: float = 1e-10,
    atol: float = 1e-12,
    jac: Callable | None = None,
) -> Reference:
    """Integrate from `t0` to `t1`, cut at the stimulus breakpoints

    Parameters
    ----------
    ns : dict
        The generated numpy code, see :func:`numpy_model`
    p, y0 : np.ndarray
        The parameters (with the stimulus of the protocol) and the initial state
    t0, t1 : float
        The interval in ms
    stim : Stim | None
        The stimulus that `p` holds, for the breakpoints. None for no stimulus
    method, rtol, atol : optional
        Passed to ``scipy.integrate.solve_ivp``
    jac : Callable | None, optional
        ``jac(t, y, p)``, the exact Jacobian, for the implicit solvers
    """
    from scipy.integrate import solve_ivp

    rhs = ns["rhs"]
    cuts = [t0, *breakpoints(stim, t0, t1), t1]
    segments = []
    y = np.array(y0, dtype=float)
    nfev = njev = 0
    start = time.perf_counter()
    for a, b in zip(cuts[:-1], cuts[1:]):
        on = stim_active(stim, a)
        if stim is None or stim.amplitude == 0:
            pa = p
        else:
            pa = forced_parameters(p, ns["parameter_index"], on=on)
        kwargs = {}
        if jac is not None:
            kwargs["jac"] = lambda t, y, pa=pa: jac(t, y, pa)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            sol = solve_ivp(
                lambda t, y, pa=pa: rhs(t, y, pa),
                (a, b),
                y,
                method=method,
                rtol=rtol,
                atol=atol,
                dense_output=True,
                **kwargs,
            )
        if not sol.success:
            raise RuntimeError(f"{method} failed on [{a:g}, {b:g}] ms: {sol.message}")
        nfev += sol.nfev
        njev += sol.njev
        segments.append(Segment(a, b, on, sol.sol))
        y = sol.y[:, -1]
    return Reference(
        t0=t0,
        t1=t1,
        segments=segments,
        rtol=rtol,
        atol=atol,
        nfev=nfev,
        njev=njev,
        wall=time.perf_counter() - start,
        method=method,
        y0=np.array(y0, dtype=float),
    )


def self_convergence(
    ref: Reference,
    ns: dict,
    p: np.ndarray,
    y0: np.ndarray,
    stim: Stim | None,
    v_index: int,
    t: np.ndarray,
    method: str | None = None,
    factor: float = 10.0,
    jac: Callable | None = None,
) -> dict:
    """Recompute with tolerances `factor` times tighter (or with another `method`), and
    return the largest difference in the membrane potential at the times `t`"""
    rtol, atol = (ref.rtol, ref.atol) if method else (ref.rtol / factor, ref.atol / factor)
    other = compute_reference(
        ns, p, y0, ref.t0, ref.t1, stim, method=method or ref.method, rtol=rtol, atol=atol, jac=jac
    )
    dv = np.abs(other(t)[v_index] - ref(t)[v_index])
    return dict(
        max_abs_dv=float(np.max(dv)),
        rtol=rtol,
        atol=atol,
        method=other.method,
        nfev=other.nfev,
        wall=other.wall,
        other=other,
    )


def save_traces(
    ref: Reference,
    path: Path,
    signals: dict[str, int],
    all_states_every_ms: float = 1.0,
    grid_ms: float = 0.01,
) -> None:
    """Write the `signals` (name: state index) on a `grid_ms` grid, and all the states every
    `all_states_every_ms`, to an ``.npz`` file"""
    t = ref.t0 + np.arange(round((ref.t1 - ref.t0) / grid_ms) + 1) * grid_ms
    y = ref(t)
    t_states = (
        ref.t0 + np.arange(round((ref.t1 - ref.t0) / all_states_every_ms) + 1) * all_states_every_ms
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        t=t,
        t_states=t_states,
        states=ref(t_states).T,
        **{name: y[index] for name, index in signals.items()},
    )
