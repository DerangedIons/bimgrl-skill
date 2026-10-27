"""Numeric access to a model that is not (yet) in the canonical .ode form

Every adapter gives batched ``rhs(t, Y)`` and ``jac(t, Y)`` for ``Y`` of shape (m, n), the state
names, the initial state and the parameters. Three kinds:

* :class:`OdeAdapter`: the gotranx numpy model (in-process; exact symbolic Jacobian).
* :class:`SubprocessAdapter`: a long-running process in any language speaking JSON lines on
  stdin/stdout (``references/contracts.md``): ``{"op": "info"}``, ``{"op": "rhs", "t": [...],
  "y": [[...], ...]}`` -> ``{"dy": [[...], ...]}``, ``{"op": "jac", ...}`` -> ``{"J": [[[...]]]}``
  (dense per point) or ``{"jac": false}`` in `info`, in which case the harness takes central
  differences through ``rhs``; ``{"op": "set_params", "parameters": {...}}``; ``{"op": "quit"}``.
  One process start, batched requests: a 13-checkpoint finite-difference Jacobian of a
  93-state model is 13 requests, not 2500.
* :class:`CtypesAdapter`: a shared library with ``int rhs(double t, const double* y, double* dy,
  const double* p)``.

The Julia helper ``scripts/jl/adapter.jl`` serves the protocol for any ``f(du, u, p, t)`` with
ForwardDiff for the Jacobian.
"""

from __future__ import annotations

import ctypes
import json
import subprocess
from pathlib import Path

import numpy as np


class Adapter:
    names: list[str]
    y0: np.ndarray
    parameters: dict
    units: dict | None = None
    has_jac: bool = False

    def rhs(self, t, Y: np.ndarray) -> np.ndarray:  # (m, n) -> (m, n)
        raise NotImplementedError

    def jac(self, t, Y: np.ndarray, scales: np.ndarray | None = None) -> np.ndarray:  # (m, n) -> (m, n, n)
        return fd_jacobian(self, t, Y, scales)

    def close(self) -> None:
        pass


def fd_jacobian(adapter: Adapter, t, Y: np.ndarray, scales: np.ndarray | None = None) -> np.ndarray:
    """Central differences through the batched rhs: for every point, 2n evaluations in one call"""
    Y = np.atleast_2d(np.asarray(Y, dtype=float))
    m, n = Y.shape
    scales = np.ones(n) if scales is None else np.asarray(scales, dtype=float)
    h = np.sqrt(np.finfo(float).eps) * np.maximum(np.abs(Y), scales[None, :])  # (m, n)
    batch = np.empty((m * 2 * n, n))
    tt = np.empty(m * 2 * n)
    t = np.broadcast_to(np.asarray(t, dtype=float), (m,))
    k = 0
    for i in range(m):
        for j in range(n):
            up, down = Y[i].copy(), Y[i].copy()
            up[j] += h[i, j]
            down[j] -= h[i, j]
            batch[k], batch[k + 1] = up, down
            tt[k] = tt[k + 1] = t[i]
            k += 2
    F = adapter.rhs(tt, batch)
    J = np.empty((m, n, n))
    k = 0
    for i in range(m):
        for j in range(n):
            J[i, :, j] = (F[k] - F[k + 1]) / (2 * h[i, j])
            k += 2
    return J


class OdeAdapter(Adapter):
    def __init__(self, model, p: np.ndarray | None = None):
        self.model = model
        self.names = list(model.names)
        self.y0 = model.y0
        self.p = model.p0 if p is None else np.asarray(p, dtype=float)
        self.parameters = dict(zip(model.parameter_names, self.p))
        by_name = {s.name: s for s in model.ode.states}
        self.units = {n: (by_name[n].unit_str or "") for n in self.names}
        self.has_jac = True
        self._jac = model.jacobian()

    def rhs(self, t, Y):
        Y = np.atleast_2d(Y)
        t = np.broadcast_to(np.asarray(t, dtype=float), (len(Y),))
        return np.vstack([self.model.rhs(ti, y, self.p) for ti, y in zip(t, Y)])

    def jac(self, t, Y, scales=None):
        Y = np.atleast_2d(Y)
        t = np.broadcast_to(np.asarray(t, dtype=float), (len(Y),))
        return np.stack([np.asarray(self._jac(ti, y, self.p), dtype=float) for ti, y in zip(t, Y)])


class SubprocessAdapter(Adapter):
    def __init__(self, command: list[str], cwd: Path | None = None):
        self.proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1, cwd=cwd)
        info = self._call({"op": "info"})
        self.names = list(info["states"])
        self.y0 = np.asarray(info["y0"], dtype=float)
        self.parameters = dict(info.get("parameters", {}))
        self.units = info.get("units")
        self.has_jac = bool(info.get("jac", False))
        self.info = info

    def _call(self, request: dict) -> dict:
        assert self.proc.stdin and self.proc.stdout
        self.proc.stdin.write(json.dumps(request) + "\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline()
        if not line:
            err = self.proc.stderr.read() if self.proc.stderr else ""
            raise SystemExit(f"the adapter exited without answering {request.get('op')}: {err[-2000:]}")
        reply = json.loads(line)
        if "error" in reply:
            raise SystemExit(f"adapter error on {request.get('op')}: {reply['error']}")
        return reply

    def rhs(self, t, Y):
        Y = np.atleast_2d(np.asarray(Y, dtype=float))
        t = np.broadcast_to(np.asarray(t, dtype=float), (len(Y),))
        out = np.asarray(self._call({"op": "rhs", "t": t.tolist(), "y": Y.tolist()})["dy"], dtype=float)
        if not np.all(np.isfinite(out)):
            raise SystemExit("the adapter returned a non-finite derivative")
        return out

    def jac(self, t, Y, scales=None):
        if not self.has_jac:
            return fd_jacobian(self, t, Y, scales)
        Y = np.atleast_2d(np.asarray(Y, dtype=float))
        t = np.broadcast_to(np.asarray(t, dtype=float), (len(Y),))
        return np.asarray(self._call({"op": "jac", "t": t.tolist(), "y": Y.tolist()})["J"], dtype=float)

    def set_params(self, parameters: dict) -> None:
        self._call({"op": "set_params", "parameters": parameters})
        self.parameters.update(parameters)

    def close(self):
        try:
            self._call({"op": "quit"})
        except Exception:
            pass
        self.proc.terminate()


class CtypesAdapter(Adapter):
    def __init__(self, library: Path, names: list[str], y0: np.ndarray, p: np.ndarray, units: dict | None = None):
        self.lib = ctypes.CDLL(str(library))
        dbl_p = ctypes.POINTER(ctypes.c_double)
        self.lib.rhs.argtypes = [ctypes.c_double, dbl_p, dbl_p, dbl_p]
        self.lib.rhs.restype = ctypes.c_int
        self.names = list(names)
        self.y0 = np.asarray(y0, dtype=float)
        self.p = np.ascontiguousarray(p, dtype=float)
        self.parameters = {f"p[{i}]": float(v) for i, v in enumerate(self.p)}
        self.units = units

    def rhs(self, t, Y):
        Y = np.atleast_2d(np.asarray(Y, dtype=float))
        t = np.broadcast_to(np.asarray(t, dtype=float), (len(Y),))
        out = np.empty_like(Y)
        dbl_p = ctypes.POINTER(ctypes.c_double)
        for i, (ti, y) in enumerate(zip(t, Y)):
            y = np.ascontiguousarray(y)
            dy = np.empty(len(y))
            rc = self.lib.rhs(ti, y.ctypes.data_as(dbl_p), dy.ctypes.data_as(dbl_p), self.p.ctypes.data_as(dbl_p))
            if rc != 0:
                raise SystemExit(f"rhs returned {rc}")
            out[i] = dy
        return out


def open_adapter(spec: str, ws=None) -> Adapter:
    """``ode`` (the workspace's model), ``cmd:<command line>`` (subprocess), ``so:<lib>:<meta.json>``"""
    if spec == "ode":
        return OdeAdapter(ws.model)
    if spec.startswith("cmd:"):
        import shlex

        return SubprocessAdapter(shlex.split(spec[4:]), cwd=ws.model_dir if ws else None)
    if spec.startswith("so:"):
        _, lib, meta = spec.split(":", 2)
        info = json.loads(Path(meta).read_text())
        return CtypesAdapter(Path(lib), info["states"], info["y0"], info["parameters"], info.get("units"))
    raise SystemExit(f"unknown adapter spec {spec!r}: use ode, cmd:<command>, or so:<library>:<meta.json>")
