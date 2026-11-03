"""The reference solution: LSODA at rtol 1e-10 / atol 1e-12, cut at the stimulus edges, written
in the trace contract on the 0.1 ms grid (protocol A: the window; protocol B: the last beat,
plus the state at the start of every beat), with a self-convergence check"""

from __future__ import annotations

import time

import numpy as np

from .model import Model
from .protocol import SAMPLE_MS, Protocol, beat_starts, protocols
from .trace import RunMeta, write_trace
from .vendored.reference import Reference, self_convergence
from .workspace import Workspace, md_table


def grid(protocol: Protocol) -> np.ndarray:
    """The 0.1 ms grid of the compared window: the whole of A, the last beat of B"""
    start = 0.0 if protocol.period is None else protocol.t1 - protocol.period
    n = round((protocol.t1 - start) / SAMPLE_MS)
    return start + np.arange(n + 1) * SAMPLE_MS


def compute(model: Model, protocol: Protocol, rtol: float, atol: float) -> tuple[Reference, dict]:
    start = time.perf_counter()
    ref = model.reference(protocol, rtol=rtol, atol=atol)
    wall = time.perf_counter() - start
    t = grid(protocol)
    conv = self_convergence(ref, model.ns, model.parameters(protocol), model.y0, protocol.stim, model.v_index, t)
    info = dict(
        method=ref.method, rtol=rtol, atol=atol, nfev=ref.nfev, wall_s=wall,
        self_convergence_max_abs_dv=conv["max_abs_dv"], self_convergence_rtol=conv["rtol"],
    )
    return ref, info


def write(ws: Workspace, protocol: Protocol, ref: Reference, info: dict) -> None:
    model = ws.model
    t = grid(protocol)
    y = ref(t).T
    tb = np.array(beat_starts(protocol))
    yb = ref(tb).T
    meta = RunMeta(
        model=model.meta.id, scheme="reference", dt=None, reference=True,
        solver=dict(name=info["method"], rtol=info["rtol"], atol=info["atol"]),
        protocol=protocol.name, beats=protocol.beats, bcl=protocol.period,
        stimulus=None if protocol.stim is None else vars(protocol.stim),
        status="ok", wall_s_total=info["wall_s"],
    )
    write_trace(ws.out / f"reference_{protocol.name}", t, y, model.names, meta, beats=(tb, yb))


def run(ws: Workspace, protocols_: list[str] = ("A", "B"), beats: int = 10, rtol: float = 1e-10, atol: float = 1e-12) -> dict:
    model = ws.model
    protos = protocols(model.meta, beats)
    rows, out = [], {}
    for name in protocols_:
        protocol = protos[name]
        ref, info = compute(model, protocol, rtol, atol)
        write(ws, protocol, ref, info)
        out[name] = info
        rows.append([name, protocol.beats, protocol.t1, info["method"], info["nfev"], info["wall_s"], info["self_convergence_max_abs_dv"]])
        print(f"reference {name}: {protocol.t1:g} ms, {info['nfev']} rhs calls, {info['wall_s']:.1f} s, "
              f"self-convergence max|dV| {info['self_convergence_max_abs_dv']:.2e} mV")
    ws.write_json("reference.json", out)
    ws.write_md("reference.md", "# Reference solutions\n\n" + md_table(
        ["protocol", "beats", "t_end (ms)", "method", "rhs calls", "wall (s)", "self-convergence max |dV| (mV)"], rows))
    return out
