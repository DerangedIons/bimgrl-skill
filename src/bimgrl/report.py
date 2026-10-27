"""The report: for every scheme and threshold set the cheapest passing run, the speed-ups, the
Table-S1 row, the stability limits against the predictions, and an error-vs-cost figure"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from .validate import passes
from .workspace import Workspace, md_table


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def load_runs(ws: Workspace) -> list[dict]:
    path = ws.out / "runs.csv"
    if not path.exists():
        raise SystemExit(f"{path} not found: run `bimgrl sweep` first")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def cost(row: dict) -> float:
    """Seconds per beat: measured if timed, else the calibration estimate"""
    if row.get("time_per_beat_s"):
        return _f(row["time_per_beat_s"])
    n = _f(row.get("n_steps")) / max(_f(row.get("beats")), 1)
    return _f(row.get("est_step_us")) * 1e-6 * n


def best(rows: list[dict], meta, protocol: str, scheme: str, set_name: str) -> dict | None:
    kind, has_ca = meta.protocol, bool(meta.calcium)
    ok = [r for r in rows if r["scheme"] == scheme and r["protocol"] == protocol and passes(r, set_name, kind, has_ca)]
    ok = [r for r in ok if not math.isnan(cost(r))]
    return min(ok, key=cost) if ok else None


def run(ws: Workspace, protocol: str = "B", baseline: str = "GRL") -> dict:
    meta = ws.meta
    rows = load_runs(ws)
    schemes = sorted({r["scheme"] for r in rows if r["protocol"] == protocol}, key=lambda s: ["FE", "RL", "GRL", "BIGRL", "BIMGRL", "BIMGRL-auto"].index(s) if s in ("FE", "RL", "GRL", "BIGRL", "BIMGRL", "BIMGRL-auto") else 99)
    table = []
    chosen = {}
    for scheme in schemes:
        for set_name in ("moderate", "tight"):
            b = best(rows, meta, protocol, scheme, set_name)
            chosen[(scheme, set_name)] = b
            table.append([scheme, set_name, f"{_f(b['dt']):g}" if b else "none", f"{1e3 * cost(b):.3g}" if b else "-",
                          "measured" if b and b.get("time_per_beat_s") else ("estimated" if b else "-")])
    speed = []
    for set_name in ("moderate", "tight"):
        base = chosen.get((baseline, set_name))
        for scheme in schemes:
            if scheme == baseline:
                continue
            b = chosen.get((scheme, set_name))
            if base and b:
                speed.append([set_name, scheme, f"{cost(base) / cost(b):.3g}", f"{_f(base['dt']):g} -> {_f(b['dt']):g}"])
            elif b and not base:
                speed.append([set_name, scheme, f"{baseline} fails", f"- -> {_f(b['dt']):g}"])
            else:
                speed.append([set_name, scheme, "-", "-"])
    status_table = []
    for r in sorted((r for r in rows if r["protocol"] == protocol), key=lambda r: (r["scheme"], -_f(r["dt"]))):
        def g(key, fmt="%+.3f"):
            return fmt % _f(r[key]) if r.get(key) else ""

        status_table.append([r["scheme"], f"{_f(r['dt']):g}", r["status"], g("err_vpeak", "%+.2f"), g("err_apd90", "%+.2f"),
                             g("err_cat", "%+.2f").replace(" ", "") and f"{100 * _f(r['err_cat']):+.2f}" if r.get("err_cat") else "", g("err_tact"), r.get("pass_moderate", ""), r.get("pass_tight", ""),
                             f"{1e3 * cost(r):.3g}" if not math.isnan(cost(r)) else "", r.get("invariant_drift_flags", "")])
    stab = {}
    if (ws.out / "stability.json").exists():
        stab = json.loads((ws.out / "stability.json").read_text())
    pred = {}
    if (ws.out / "stiffness.json").exists():
        pred = json.loads((ws.out / "stiffness.json").read_text())["prediction"]
    out = [f"# Report: `{meta.id}`, protocol {protocol}", ""]
    out += ["## Cheapest run within each threshold set", "", md_table(["scheme", "thresholds", "dt (ms)", "cost (ms per beat)", "cost source"], table)]
    out += ["", f"## Speed-up against {baseline} at matched accuracy", "", md_table(["thresholds", "scheme", "speed-up", "dt (base -> scheme)"], speed)]
    if stab or pred:
        out += ["", "## Stability limits: measured against predicted", ""]
        for scheme, d in stab.items():
            out.append(f"- {scheme}: largest step that runs {d['limit_ms']:g} ms" + (f" (fails at {d['fails_at_ms']:g} ms)" if d.get("fails_at_ms") else ""))
        if pred:
            p1g, p1r, p2g = pred["p1_grl"], pred["p1_rl_fe"], pred["p2_grl"]
            out.append(f"- predicted: GRL P1 {p1g['limit_ms'] if p1g['limit_ms'] is None else f'{p1g['limit_ms']:.4g} ms'}; GRL P2 {p2g['limit_ms']:.4g} ms; RL+FE P1 {p1r['limit_ms'] if p1r['limit_ms'] is None else f'{p1r['limit_ms']:.4g} ms'}")
    out += ["", "## All runs", "", md_table(["scheme", "dt", "status", "dVpeak (mV)", "dAPD90 (ms)", "dCaT (%)", "dt_act (ms)", "moderate", "tight", "cost (ms/beat)", "drift flags"], status_table)]
    # Table-S1 row
    bm, bt = chosen.get(("BIMGRL", "moderate")), chosen.get(("BIMGRL", "tight"))
    gm, gt = chosen.get((baseline, "moderate")), chosen.get((baseline, "tight"))
    n_imp = n_sub = 0
    blocks = ws.model_dir / "blocks.toml"
    if blocks.exists():
        text = blocks.read_text()
        n_imp, n_sub = text.count("[[implicit_blocks]]"), text.count("[[substep_blocks]]")
    row = [meta.name or meta.id, f"{meta.species}, {meta.cell_type}".strip(", "), len(ws.model.names) if ws.ode_path.exists() else "", f"{n_imp} / {n_sub}",
           f"{_f(gm['dt']):g} / {1e3 * cost(gm):.3g}" if gm else "none", f"{_f(bm['dt']):g} / {1e3 * cost(bm):.3g}" if bm else "none",
           f"{cost(gm) / cost(bm):.3g}" if gm and bm else ("GRL fails" if bm else "-"), f"{cost(gt) / cost(bt):.3g}" if gt and bt else ("GRL fails" if bt else "-")]
    out += ["", "## Supplementary Table S1 row", "", md_table(["Model", "Species, cell", "States", "Blocks (implicit / sub-step)", "GRL dt (ms) / cost (ms per beat)", "BIMGRL dt (ms) / cost (ms per beat)", "Speed-up, moderate", "Speed-up, tight"], [row])]
    fig = figure(ws, rows, protocol)
    if fig:
        out += ["", f"Figure: `{fig.relative_to(ws.model_dir)}` (error in units of the moderate threshold against cost per beat)."]
    text = "\n".join(out)
    ws.write_md("report.md", text)
    ws.write_json("report.json", dict(chosen={f"{k[0]}:{k[1]}": v for k, v in chosen.items()}, table_s1=row))
    print(text)
    return dict(chosen=chosen, row=row)


def figure(ws: Workspace, rows: list[dict], protocol: str) -> Path | None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    from .vendored.select import BIOMARKER_SETS

    tol = BIOMARKER_SETS["moderate"]
    fig, ax = plt.subplots(figsize=(5, 4))
    for scheme, marker in (("GRL", "o"), ("BIGRL", "s"), ("BIMGRL", "^"), ("BIMGRL-auto", "D"), ("RL", "v"), ("FE", "x")):
        pts = []
        for r in rows:
            if r["scheme"] != scheme or r["protocol"] != protocol or r["status"] != "ok":
                continue
            errs = [abs(_f(r.get("err_vpeak"))) / tol["vpeak"], abs(_f(r.get("err_apd90"))) / tol["apd90"], abs(_f(r.get("err_tact"))) / tol["tact"]]
            if r.get("err_cat"):
                errs.append(abs(_f(r["err_cat"])) / tol["cat"])
            e = max(x for x in errs if not math.isnan(x)) if any(not math.isnan(x) for x in errs) else math.nan
            c = cost(r)
            if not math.isnan(e) and not math.isnan(c):
                pts.append((c * 1e3, max(e, 1e-3), _f(r["dt"])))
        if pts:
            pts.sort(key=lambda p: -p[2])
            ax.loglog([p[0] for p in pts], [p[1] for p in pts], marker=marker, label=scheme)
    ax.axhline(1.0, ls="--", color="k", lw=0.8)
    ax.set_xlabel("cost (ms per beat)")
    ax.set_ylabel("largest biomarker error / moderate threshold")
    ax.set_title(f"{ws.meta.id}, protocol {protocol}")
    ax.legend()
    path = ws.out / "error_vs_cost.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
