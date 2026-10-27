"""`bimgrl validate`: a trace file from any language against a reference trace, with the same
biomarkers, thresholds, physical-range and invariant-drift checks as the sweep"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .protocol import Protocol
from .trace import read_trace
from .validate import compare, violations
from .vendored.metrics import bounds, classify
from .workspace import Workspace, md_table


def _interp(t_to: np.ndarray, t_from: np.ndarray, y_from: np.ndarray) -> np.ndarray:
    return np.column_stack([np.interp(t_to, t_from, y_from[:, j]) for j in range(y_from.shape[1])])


def run(ws: Workspace, test: Path, reference: Path | None, level: str = "moderate", name: str | None = None) -> dict:
    meta = ws.meta
    tr = read_trace(test)
    if reference is None:
        proto = tr.meta.protocol or "B"
        reference = ws.out / f"reference_{proto}.csv"
    ref = read_trace(reference)
    if not ref.meta.reference or (ref.meta.solver and float(ref.meta.solver.get("rtol", 1)) > 1e-8):
        raise SystemExit(f"{reference} is not a reference trace (run_meta.reference must be true and solver.rtol <= 1e-8)")
    if tr.names != ref.names:
        missing = set(ref.names) - set(tr.names)
        raise SystemExit(f"state columns differ from the reference; missing in the test: {sorted(missing)[:10]}")
    names = tr.names
    # the compared window: where the test and the reference overlap
    t0, t1 = max(tr.t[0], ref.t[0]), min(tr.t[-1], ref.t[-1])
    mask = (tr.t >= t0 - 1e-9) & (tr.t <= t1 + 1e-9)
    t, y = tr.t[mask], tr.y[mask]
    yref = _interp(t, ref.t, ref.y)
    kinds = ws.model.kinds if ws.ode_path.exists() else _kinds_from_meta(meta, names)
    lo, hi = bounds(kinds, ref.y.min(axis=0), ref.y.max(axis=0))
    status = classify(tr.y, kinds, lo, hi)
    stim = tr.meta.stimulus
    from .meta import Stim

    protocol = Protocol(tr.meta.protocol or "B", tr.meta.beats or 1, tr.meta.bcl, t1, Stim(**stim) if stim else None)
    out = dict(test=str(test), reference=str(reference), status=status, scheme=tr.meta.scheme, dt=tr.meta.dt)
    if status != "ok":
        out["violations"] = violations(names, tr.t, tr.y, lo, hi)
    else:
        invariants = ws.read_json("invariants.json")["invariants"] if (ws.out / "invariants.json").exists() else []
        tb = yb = yrb = None
        if tr.beats is not None and ref.beats is not None:
            tb, yb = tr.beats
            trb, yrb_all = ref.beats
            yrb = _interp(tb, trb, yrb_all)
        out.update(compare(meta, names, protocol, t, y, yref, tb, yb, yrb, invariants))
    stem = name or Path(test).stem
    ws.write_json(f"validation_{stem}.json", out)
    rows = []
    for key, unit in (("vpeak", "mV"), ("apd90", "ms"), ("cat", "rel"), ("tact", "ms"), ("cycle_length", "ms")):
        if f"err_{key}" in out:
            rows.append([key, out.get(f"ref_{key}"), out.get(key), out[f"err_{key}"], unit])
    md = [f"# Validation of `{Path(test).name}` against `{Path(reference).name}`", "", f"Status: **{status}**; scheme {tr.meta.scheme}, dt {tr.meta.dt}", ""]
    if rows:
        md += [md_table(["biomarker", "reference", "test", "error", "unit"], rows), "",
               f"Moderate thresholds: {'pass' if out.get('pass_moderate') else 'FAIL'}; tight: {'pass' if out.get('pass_tight') else 'FAIL'}"]
        if out.get("invariant_drift_flags"):
            md += ["", f"Invariant drift flags: {out['invariant_drift_flags']}"]
    if out.get("violations"):
        md += ["", "First violations of the physical range:"] + [f"- t = {v['t']:.3g}: {v['state']} = {v['value']:.4g} outside [{v['lo']:.3g}, {v['hi']:.3g}]" for v in out["violations"]]
    ws.write_md(f"validation_{stem}.md", "\n".join(md))
    print("\n".join(md))
    return out


def _kinds_from_meta(meta, names):
    out = []
    for n in names:
        if n == meta.voltage:
            out.append("voltage")
        elif n in meta.kinds:
            out.append(meta.kinds[n])
        else:
            out.append("other")
    return out
