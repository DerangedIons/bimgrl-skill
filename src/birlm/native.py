"""The native path: stiffness and invariants through an adapter, for a model that is not in
.ode form. Checkpoints come from a trace the user's own code produced (any scheme; the
upstroke, plateau and diastole must be covered)."""

from __future__ import annotations

import numpy as np

from .adapters import Adapter
from .trace import read_trace
from .vendored.stiffness import OWNER, fe_limit, p2_checkpoint, participation
from .workspace import Workspace, md_table


def checkpoints_from_trace(trace, n_dense: int = 60) -> tuple[np.ndarray, np.ndarray]:
    """The paper's offsets after the upstroke (max dV/dt) plus a spread of dense samples"""
    t = trace.t
    V = trace.column(trace.meta.voltage) if getattr(trace.meta, "voltage", None) else None
    picks = np.linspace(0, len(t) - 1, n_dense).astype(int)
    return t[picks], trace.y[picks]


def stiffness(ws: Workspace, adapter: Adapter, trace_path, threshold: float = -20.0, voltage: str | None = None) -> dict:
    tr = read_trace(trace_path)
    names = adapter.names
    if tr.names != names:
        raise SystemExit("the trace's state columns do not match the adapter's state names")
    v = voltage or ws.meta.voltage
    V = tr.column(v)
    up = int(np.argmax(np.gradient(V, tr.t)))
    offsets = (-1.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 300.0, 500.0)
    want = [tr.t[up] + o for o in offsets]
    idx = sorted({int(np.argmin(np.abs(tr.t - w))) for w in want if tr.t[0] <= w <= tr.t[-1]} | set(np.linspace(0, len(tr.t) - 1, 40).astype(int)))
    times, Y = tr.t[idx], tr.y[idx]
    scales = np.array([1.0 if k == "voltage" else 1e-3 if k == "gate" else 1e-4 for k in _kinds(ws, names)])
    J = adapter.jac(times, Y, scales)
    rows, per_t = [], []
    for ti, Ji in zip(times, J):
        if not np.all(np.isfinite(Ji)):
            continue
        lam, p = participation(Ji)
        fe = {names[k] for k in range(len(names)) if abs(Ji[k, k]) < 1e-8}
        modes = []
        for j in np.argsort(lam.real):
            if lam[j].real > threshold:
                break
            order = np.argsort(-p[:, j])
            owners = [(names[k], float(p[k, j])) for k in order if p[k, j] > OWNER]
            modes.append(dict(lam_re=float(lam[j].real), lam_im=float(lam[j].imag), fe_limit_ms=fe_limit(lam[j]), owners=owners,
                              fe_share=float(sum(p[k, j] for k in range(len(names)) if names[k] in fe))))
        per_t.append(dict(t=float(ti), modes=modes))
    table = {}
    for e in per_t:
        for m in e["modes"]:
            key = ";".join(sorted(n for n, _ in m["owners"]))
            r = table.setdefault(key, dict(owners=m["owners"], min_lam_re=0.0, n=0, fe_share=0.0))
            r["n"] += 1
            r["fe_share"] = max(r["fe_share"], m["fe_share"])
            if m["lam_re"] < r["min_lam_re"]:
                r.update(min_lam_re=m["lam_re"], t=e["t"], owners=m["owners"], fe_limit_ms=m["fe_limit_ms"])
    rows = sorted(table.values(), key=lambda r: r["min_lam_re"])
    p2 = min(p2_checkpoint(Ji, "FE") for Ji in J if np.all(np.isfinite(Ji)))
    result = dict(model=ws.meta.id, source=str(trace_path), n_checkpoints=len(per_t), threshold=threshold, modes=rows, p2_fe_ms=float(p2),
                  note="update type per state is not known for a native integrator: read it off the user's code (exponential vs forward Euler) "
                       "and compare with the fe_share (share of the mode on states whose Jacobian diagonal is zero)")
    ws.write_json("stiffness_native.json", result)
    md = [f"# Stiffness (native adapter) of `{ws.meta.id}`", "", f"{len(per_t)} checkpoints from `{trace_path}`; forward-Euler P2 limit {p2:.4g} ms.", "",
          md_table(["owners (participation at min)", "min Re λ (/ms)", "FE limit (µs)", "t (ms)", "checkpoints", "share on zero-diagonal states"],
                   [[", ".join(f"{n} {p:.2f}" for n, p in r["owners"]), r["min_lam_re"], 1e3 * r["fe_limit_ms"], r["t"], r["n"], f"{r['fe_share']:.2f}"] for r in rows]),
          "", result["note"]]
    ws.write_md("stiffness_native.md", "\n".join(md))
    print("\n".join(md))
    return result


def _kinds(ws: Workspace, names):
    meta = ws.meta
    out = []
    for n in names:
        if n == meta.voltage:
            out.append("voltage")
        elif n in meta.kinds:
            out.append(meta.kinds[n])
        elif any(x in n.lower() for x in ("ca", "na", "k_", "ki", "cl", "mg")) and not n.lower().startswith(("x", "y", "f", "d", "m", "h", "j")):
            out.append("conc")
        else:
            out.append("gate")
    return out


def invariants_numeric(ws: Workspace, adapter: Adapter, trace_path, groups: list[list[str]] | None = None, tol: float = 1e-10) -> dict:
    """Left null vectors shared by the Jacobians at the checkpoints: candidate linear invariants
    cᵀx. Restricted to the given groups of states (else all), rationalised to small integers
    or to ratios of the adapter's parameters when they match"""
    tr = read_trace(trace_path)
    names = adapter.names
    idx = np.linspace(0, len(tr.t) - 1, 25).astype(int)
    J = adapter.jac(tr.t[idx], tr.y[idx])
    J = J[np.all(np.isfinite(J), axis=(1, 2))]
    n = len(names)
    groups = groups or [names]
    found = []
    for g in groups:
        rows = [names.index(s) for s in g]
        M = np.concatenate([Ji[rows, :] for Ji in J], axis=1)  # (|g|, K n)
        U, S, _ = np.linalg.svd(M, full_matrices=True)
        small = S / S[0] < tol if len(S) else np.array([])
        nulls = [U[:, k] for k in range(len(S)) if small[k]] + [U[:, k] for k in range(len(S), U.shape[1])]
        for c in nulls:
            c = c / c[np.argmax(np.abs(c))]
            found.append(dict(members=[g[i] for i in range(len(g)) if abs(c[i]) > 1e-8], weights={g[i]: round(float(c[i]), 10) for i in range(len(g)) if abs(c[i]) > 1e-8},
                              residual=float(np.max(np.abs(c @ M))), rationalised=_rationalise(c, g, adapter.parameters)))
    result = dict(model=ws.meta.id, source=str(trace_path), n_checkpoints=len(J), tol=tol, invariants=found,
                  note="numeric left null space: an invariant conserved by the ODE up to round-off; conservation up to boundary fluxes cannot be told from an open system here")
    ws.write_json("invariants_native.json", result)
    md = [f"# Linear invariants (native adapter) of `{ws.meta.id}`", "", md_table(["members", "weights", "max |cᵀJ|", "as parameter ratios"],
          [[", ".join(f["members"]), ", ".join(f"{k}: {v:g}" for k, v in f["weights"].items()), f["residual"], f["rationalised"]] for f in found]) if found else "none found", "", result["note"]]
    ws.write_md("invariants_native.md", "\n".join(md))
    print("\n".join(md))
    return result


def _rationalise(c: np.ndarray, names: list[str], params: dict) -> str:
    cand = {"1": 1.0}
    vols = {k: float(v) for k, v in params.items() if isinstance(v, (int, float)) and v > 0 and any(x in k.lower() for x in ("v", "vol"))}
    for a, va in vols.items():
        for b, vb in vols.items():
            if a != b:
                cand[f"{a}/{b}"] = va / vb
    out = []
    for i, x in enumerate(c):
        if abs(x) < 1e-8:
            continue
        best = min(cand.items(), key=lambda kv: abs(abs(x) - kv[1]) / kv[1])
        if abs(abs(x) - best[1]) / best[1] < 1e-6:
            out.append(f"{names[i]}: {'-' if x < 0 else ''}{best[0]}")
        else:
            out.append(f"{names[i]}: {x:.6g}")
    return ", ".join(out)
