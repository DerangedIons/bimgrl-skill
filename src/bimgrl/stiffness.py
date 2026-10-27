"""Stiffness along the reference: eigenmodes, participation factors, owners, how each scheme
advances the owners, and the predicted stability limits

Checkpoints: the paper's 13 (1 ms before the stimulus and 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200,
300, 500 ms after it) in the last beat of protocol B, plus every ``every_ms`` through that beat
(the survey's sampling). A mode with Re λ ≤ threshold is stiff; its owners are the states with
participation above 0.1. Two scheme views: GRL as gotranx generates it (exponential wherever the
Jacobian diagonal is non-zero) and RL+FE (exponential for the HH gates only, as TWorld.jl).

Predictions: P1 = 2|Re λ|/|λ|² of the stiffest mode owned (≥ 50 % participation) by states
advanced by forward Euler; P2 = the largest h at which the linearised scheme step does not grow
faster than the model (catches coupled modes that GRL damps poorly).
"""

from __future__ import annotations

import numpy as np

from .model import Model
from .protocol import protocols
from .vendored.stiffness import OWNER, fe_limit, p2_checkpoint, participation
from .workspace import Workspace, md_table

PAPER_OFFSETS = (-1.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 300.0, 500.0)


def checkpoints(protocol, every_ms: float) -> tuple[np.ndarray, np.ndarray]:
    start = 0.0 if protocol.period is None else protocol.t1 - protocol.period
    length = protocol.t1 - start
    stim_t = start + (protocol.stim.start if protocol.stim else 0.0)
    paper = np.array([stim_t + o for o in PAPER_OFFSETS])
    paper = paper[(paper >= start) & (paper <= protocol.t1)]
    dense = start + np.arange(0.0, length + 1e-9, every_ms)
    return np.unique(np.round(np.concatenate([paper, dense]), 6)), paper


def modes_at(J: np.ndarray, names, fe_grl: set, fe_rl: set, threshold: float) -> list[dict]:
    lam, p = participation(J)
    out = []
    for j in np.argsort(lam.real):
        if lam[j].real > threshold:
            break
        order = np.argsort(-p[:, j])
        owners = [(names[k], float(p[k, j])) for k in order if p[k, j] > OWNER]
        out.append(
            dict(
                lam_re=float(lam[j].real), lam_im=float(lam[j].imag), tau_us=1e3 / abs(lam[j].real),
                fe_limit_ms=fe_limit(lam[j]), owners=owners,
                owner_key=";".join(sorted(n for n, _ in owners)),
                fe_share_grl=float(sum(p[k, j] for k in range(len(names)) if names[k] in fe_grl)),
                fe_share_rl=float(sum(p[k, j] for k in range(len(names)) if names[k] in fe_rl)),
            )
        )
    return out


def analyse(model: Model, structure: dict, protocol, every_ms: float, threshold: float) -> dict:
    ref = model.reference(protocol)
    times, paper = checkpoints(protocol, every_ms)
    jac = model.jacobian()
    p = model.parameters(protocol)
    names = model.names
    upd_grl = {s["name"]: s["update_grl"] for s in structure["states"]}
    upd_rl = {s["name"]: s["update_rl"] for s in structure["states"]}
    fe_grl = {n for n, u in upd_grl.items() if u == "FE"}
    fe_rl = {n for n, u in upd_rl.items() if u == "FE"}
    Y = ref(times).T
    per_t, Js, skipped = [], [], []
    with np.errstate(all="ignore"):
        for t, y in zip(times, Y):
            J = np.asarray(jac(t, y, p), dtype=float)
            if not np.all(np.isfinite(J)):
                skipped.append(float(t))
                continue
            # a numerically zero diagonal at this point also means forward Euler in gotranx's GRL
            fe_here = fe_grl | {names[k] for k in range(len(names)) if abs(J[k, k]) < 1e-8}
            per_t.append(dict(t=float(t), paper=bool(np.any(np.isclose(paper, t))), modes=modes_at(J, names, fe_here, fe_rl, threshold)))
            Js.append(J)
    # Table-1 view: one row per owner set
    table: dict[str, dict] = {}
    for entry in per_t:
        for m in entry["modes"]:
            row = table.setdefault(m["owner_key"], dict(owners=m["owners"], min_lam_re=0.0, t=None, n=0, part=[], fe_share_grl=0.0, fe_share_rl=0.0))
            row["n"] += 1
            row["part"].append(dict(m["owners"]))
            row["fe_share_grl"] = max(row["fe_share_grl"], m["fe_share_grl"])
            row["fe_share_rl"] = max(row["fe_share_rl"], m["fe_share_rl"])
            if m["lam_re"] < row["min_lam_re"]:
                row.update(min_lam_re=m["lam_re"], lam_im=m["lam_im"], t=entry["t"], owners=m["owners"])
    rows = []
    for key, row in sorted(table.items(), key=lambda kv: kv[1]["min_lam_re"]):
        states = sorted({n for d in row["part"] for n in d})
        mean_part = {n: float(np.mean([d.get(n, 0.0) for d in row["part"]])) for n in states}
        rows.append(
            dict(
                owners=states, min_lam_re=row["min_lam_re"], lam_im=row.get("lam_im", 0.0), t_at_min=row["t"],
                tau_us=1e3 / abs(row["min_lam_re"]), fe_limit_ms=fe_limit(complex(row["min_lam_re"], row.get("lam_im", 0.0))),
                n_checkpoints=row["n"], mean_participation=mean_part,
                update_grl={n: upd_grl[n] for n in states}, update_rl={n: upd_rl[n] for n in states},
                fe_share_grl=row["fe_share_grl"], fe_share_rl=row["fe_share_rl"],
            )
        )
    # predictions
    def p1(share_key):
        lim = [r["fe_limit_ms"] for r in rows if r[share_key] >= 0.5]
        which = [r for r in rows if r[share_key] >= 0.5]
        if not lim:
            return dict(limit_ms=None, mode=None)
        i = int(np.argmin(lim))
        return dict(limit_ms=lim[i], mode=which[i]["owners"], lam_re=which[i]["min_lam_re"])

    # P2 costs ~60 eigendecompositions per checkpoint and scheme: the paper's checkpoints plus
    # every tenth dense one are enough to find the minimum
    sel = [i for i, e in enumerate(per_t) if e["paper"] or i % 10 == 0]
    p2 = {}
    for scheme in ("GRL", "FE"):
        limits = [p2_checkpoint(Js[i], scheme) for i in sel]
        k = int(np.argmin(limits)) if limits else None
        p2[scheme] = dict(limit_ms=float(limits[k]) if limits else None, t=per_t[sel[k]]["t"] if limits else None, n_checkpoints=len(sel))
    all_eigs = np.concatenate([np.linalg.eigvals(J) for J in Js]) if Js else np.array([])
    return dict(
        model=model.meta.id, protocol=protocol.name, threshold=threshold, every_ms=every_ms,
        n_checkpoints=len(per_t), skipped=skipped, paper_checkpoints=[float(x) for x in paper],
        min_lam_re=float(all_eigs.real.min()) if len(all_eigs) else None,
        modes=rows, checkpoints=per_t,
        prediction=dict(
            p1_grl=p1("fe_share_grl"), p1_rl_fe=p1("fe_share_rl"),
            p2_grl=p2["GRL"], p2_fe=p2["FE"],
        ),
    )


def to_markdown(r: dict) -> str:
    out = [f"# Stiffness of `{r['model']}` (protocol {r['protocol']})", ""]
    out.append(f"{r['n_checkpoints']} checkpoints (the paper's 13 plus every {r['every_ms']:g} ms); modes with Re λ ≤ {r['threshold']:g} /ms. "
               f"Most negative eigenvalue over all checkpoints: {r['min_lam_re']:.4g} /ms.")
    out += ["", "## Stiff modes by owner set (Table-1 view)", ""]
    out.append(md_table(
        ["owners (mean participation)", "min Re λ (/ms)", "τ (µs)", "FE limit (µs)", "t at min (ms)", "checkpoints", "GRL update of owners", "RL+FE update of owners", "FE share GRL / RL+FE"],
        [[", ".join(f"{n} {p:.2f}" for n, p in m["mean_participation"].items()), m["min_lam_re"], m["tau_us"], 1e3 * m["fe_limit_ms"], m["t_at_min"], m["n_checkpoints"],
          ", ".join(f"{n}: {u}" for n, u in m["update_grl"].items()), ", ".join(f"{n}: {u}" for n, u in m["update_rl"].items()),
          f"{m['fe_share_grl']:.2f} / {m['fe_share_rl']:.2f}"] for m in r["modes"]],
    ))
    p = r["prediction"]
    out += ["", "## Predicted stability limits", ""]

    def fmt(x):
        return "none (no FE-owned stiff mode)" if x.get("limit_ms") is None else f"{x['limit_ms']:.4g} ms (mode {', '.join(x['mode'])}, Re λ {x['lam_re']:.4g} /ms)"

    out.append(f"- GRL as generated (exponential wherever the diagonal is non-zero): P1 {fmt(p['p1_grl'])}; P2 (linearised step) {p['p2_grl']['limit_ms']:.4g} ms at t = {p['p2_grl']['t']:g} ms")
    out.append(f"- RL+FE (exponential for HH gates only): P1 {fmt(p['p1_rl_fe'])}")
    out.append(f"- forward Euler everywhere: P2 {p['p2_fe']['limit_ms']:.4g} ms")
    out += ["", "P1 is the forward Euler limit 2|Re λ|/|λ|² of the stiffest mode whose owners the scheme advances explicitly. "
            "P2 is the largest step at which the linearised scheme step grows no faster than the model, which also sees coupled modes "
            "that GRL updates exponentially one state at a time and damps poorly. Verify with `bimgrl sweep --stability`."]
    return "\n".join(out)


def run(ws: Workspace, protocol: str = "B", beats: int = 2, every_ms: float = 1.0, threshold: float = -20.0) -> dict:
    model = ws.model
    structure = ws.read_json("structure.json")
    proto = protocols(model.meta, beats)[protocol]
    result = analyse(model, structure, proto, every_ms, threshold)
    ws.write_json("stiffness.json", result)
    ws.write_md("stiffness.md", to_markdown(result))
    return result
