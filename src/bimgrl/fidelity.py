"""Prove that a translated model.ode reproduces the original code

1. parameters: every .ode default equals params.json after the declared scaling (1e-12 rel);
2. pointwise RHS: at every row of rhs_samples.csv (t, y_*, dy_* from the original code), the
   .ode right-hand side matches to 1e-9 relative to max(|dy|, scale); this is the proof, it is
   integrator-free and localises a mismatch to a state;
3. trajectory: the .ode integrated by LSODA from the first row of the user's trace against
   that trace, with tight thresholds when the trace came from an adaptive solver (rtol <= 1e-8),
   moderate thresholds plus a mandatory pass of check 2 otherwise.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .trace import read_trace
from .workspace import Workspace, md_table


def run(ws: Workspace, samples: Path | None, params: Path | None, trace: Path | None, map_file: Path | None, rhs_tol: float = 1e-9) -> dict:
    model = ws.model
    names = model.names
    tmap = json.loads(Path(map_file).read_text()) if map_file else {}
    smap = tmap.get("states", {})      # original name -> ode name
    pmap = tmap.get("parameters", {})  # original name -> ode name
    scale = tmap.get("scale", {})      # ode name -> {factor}
    report: dict = {}
    md = ["# Fidelity of the translation", ""]

    # 1. parameters
    if params:
        orig = json.loads(Path(params).read_text())
        index = model.ns["parameter_index"]
        p0 = model.p0
        bad = []
        for oname, value in orig.items():
            ode_name = pmap.get(oname, oname)
            try:
                i = index(ode_name)
            except Exception:
                bad.append((oname, "not a parameter of the .ode"))
                continue
            factor = scale.get(ode_name, {}).get("factor", 1.0)
            if abs(p0[i] - value * factor) > 1e-12 * max(abs(value * factor), 1e-300):
                bad.append((oname, f".ode {p0[i]!r} vs original {value!r} x {factor}"))
        report["parameters"] = dict(checked=len(orig), mismatches=bad)
        md += [f"## 1. Parameters: {len(orig)} checked, {len(bad)} mismatches", ""] + [f"- {n}: {why}" for n, why in bad[:30]]

    # 2. pointwise RHS
    if samples:
        with open(samples) as f:
            header = f.readline().strip().split(",")
        data = np.loadtxt(samples, delimiter=",", skiprows=1, ndmin=2)
        ycols = [h[2:] for h in header if h.startswith("y_")]
        dcols = [h[3:] for h in header if h.startswith("dy_")]
        if set(ycols) != set(dcols):
            raise SystemExit("rhs_samples.csv needs matching y_<state> and dy_<state> columns")
        col = {h: i for i, h in enumerate(header)}
        p = model.p0
        worst = {n: 0.0 for n in names}
        n_rows = len(data)
        for row in data:
            t = row[col["t"]] * tmap.get("time_factor", 1.0)
            y = np.zeros(len(names))
            dy_orig = np.zeros(len(names))
            for oname in ycols:
                ode_name = smap.get(oname, oname)
                i = names.index(ode_name)
                f = scale.get(ode_name, {}).get("factor", 1.0)
                y[i] = row[col[f"y_{oname}"]] * f
                dy_orig[i] = row[col[f"dy_{oname}"]] * f * tmap.get("time_factor", 1.0) ** -1
            dy = model.rhs(t, y, p)
            for i, n in enumerate(names):
                sc = scale.get(n, {}).get("rhs_scale", 1.0)
                err = abs(dy[i] - dy_orig[i]) / max(abs(dy_orig[i]), sc)
                worst[n] = max(worst[n], err)
        failing = {n: e for n, e in worst.items() if e > rhs_tol}
        report["rhs"] = dict(rows=n_rows, tol=rhs_tol, worst=worst, failing=failing)
        md += ["", f"## 2. Pointwise RHS: {n_rows} samples, tolerance {rhs_tol:g}: {'PASS' if not failing else 'FAIL'}", ""]
        md.append(md_table(["state", "worst relative error"], [[n, f"{e:.2e}"] for n, e in sorted(worst.items(), key=lambda kv: -kv[1])[:15]]))

    # 3. trajectory
    if trace:
        from .meta import Stim
        from .protocol import Protocol
        from .validate import compare

        tr = read_trace(trace)
        if tr.names != names:
            md.append("\n## 3. Trajectory: skipped, the trace columns are not the .ode state names (apply the map first)")
        else:
            adaptive = tr.meta.reference or (tr.meta.solver and float(tr.meta.solver.get("rtol", 1)) <= 1e-8)
            stim = Stim(**tr.meta.stimulus) if tr.meta.stimulus and tr.meta.stimulus.get("amplitude") is not None else None
            protocol = Protocol(tr.meta.protocol or "A", 1, None, float(tr.t[-1]), stim)
            ref = model.reference(protocol, y0=tr.y[0])
            yref = ref(tr.t).T
            res = compare(model.meta, names, protocol, tr.t, yref, tr.y)  # errors of the .ode against the user's trace
            level = "tight" if adaptive else "moderate"
            ok = res.get(f"pass_{level}", False) and (not samples or not report["rhs"]["failing"])
            report["trajectory"] = dict(level=level, adaptive=bool(adaptive), pass_=bool(ok), errors={k: v for k, v in res.items() if k.startswith("err_")})
            md += ["", f"## 3. Trajectory against `{Path(trace).name}` ({'adaptive solver: tight' if adaptive else 'fixed-step trace: moderate + RHS check'}): {'PASS' if ok else 'FAIL'}", ""]
            md.append(md_table(["biomarker", "error"], [[k[4:], f"{v:+.4g}"] for k, v in res.items() if k.startswith("err_")]))
            if samples and not report["rhs"]["failing"] and not res.get(f"pass_{level}"):
                md.append("\nThe equations match pointwise but the trajectories differ: look at the protocol (stimulus timing, amplitude), the initial state, or the solver tolerance of the trace, not at the equations.")
    ws.write_json("fidelity.json", report)
    ws.write_md("fidelity.md", "\n".join(md))
    print("\n".join(md))
    return report
