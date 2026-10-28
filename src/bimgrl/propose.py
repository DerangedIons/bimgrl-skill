"""Propose the blocks from the structure, the invariants, the stiffness and the compartments,
by the rules of the paper ("Choice of blocks"), each block with the finding that produced it;
then compare with gotranx's own ``suggest-blocks`` and list the judgment calls

Rules:
1. Every Markov chain of an ion channel or of the release channel (a closed or exactly
   conserving chain that feeds V, or that the compartments name as the release channel) is an
   implicit block. A chain that does not feed V (contraction) is listed, not blocked.
2. Every binding cluster (a Ca pool with its buffer states) that owns a stiff mode is an
   implicit block; the non-stiff ones (cytosolic TnC, CaM) stay explicit.
3. The SR and the junctional pool, coupled by release, share one block, so that the release
   conserves Ca; if the junction has no buffer states it joins the SR block as a bare state.
4. V with the INa (and INaL) gates: a proportional sub-step block, rate trigger on V with
   atol 0.5 mV, h_min 0.001 ms.
5. The dyad: junctional Ca with its buffers, the release channel states, and the SR with its
   buffers: a binary sub-step block, h_min 0.005 ms, rate triggers on the junctional Ca
   (rtol 0.02) and on each open state of the release channel (atol 1e-3).
6. An implicit block that shares states with a sub-step block is placed entirely inside it.
"""

from __future__ import annotations

from pathlib import Path

from . import stiffness
from .workspace import Workspace, md_table

V_ATOL = 0.5
V_H_MIN = 0.001
DYAD_H_MIN = 0.005
DYAD_RTOL = 0.02
OPEN_ATOL = 1e-3


def _toml_list(xs):
    return "[" + ", ".join(f'"{x}"' for x in xs) + "]"


def propose(ws: Workspace, dt: float = 0.1, stiff_threshold: float = -20.0, nest_v_markov: bool = False) -> dict:
    meta = ws.meta
    comps = meta.compartments
    structure = ws.read_json("structure.json")
    inv = ws.read_json("invariants.json")
    # the default file holds every mode at or below the default threshold, so a stricter threshold
    # filters it; a looser one needs the probe written by `bimgrl stiffness --threshold` at that value
    stiff_file = f"{stiffness.artifact(max(stiff_threshold, stiffness.THRESHOLD))}.json"
    stiff = ws.read_json(stiff_file) if (ws.out / stiff_file).exists() else None
    feeds_v = set(structure["feeds_v"]) | {structure["voltage"]}
    kind_of = {s["name"]: s["kind"] for s in structure["states"]}

    stiff_owners: dict[str, float] = {}
    if stiff:
        for m in stiff["modes"]:
            if m["min_lam_re"] <= stiff_threshold:
                for n in m["owners"]:
                    stiff_owners[n] = min(stiff_owners.get(n, 0.0), m["min_lam_re"])

    implicit: list[dict] = []
    notes: list[str] = []
    judgment: list[str] = []
    release = set(comps.release_channel)

    # 1. Markov chains
    for g in inv["invariants"]:
        if g["level"] != "subgroup" or g["structure"] != "markov":
            continue
        members = g["members"]
        if g["class"] == "improper":
            judgment.append(f"chain {', '.join(members)} is flagged as not conserving (see invariants.md); fix the model before blocking it")
            continue
        is_channel = bool(set(members) & feeds_v) or bool(set(members) & release)
        if not is_channel:
            notes.append(f"chain {', '.join(members)} does not feed V (contraction?): left explicit")
            continue
        implicit.append(dict(states=members, why=f"Markov chain ({g['class']}); rule 1"))

    # 2. binding clusters that own a stiff mode
    clusters = [g for g in inv["invariants"] if g["level"] == "subgroup" and g["structure"] == "binding"]
    cluster_of = {}
    for g in clusters:
        members = g["members"]
        owners = [n for n in members if n in stiff_owners]
        if g["class"] == "improper":
            judgment.append(f"binding cluster {', '.join(members)} is flagged as not conserving; fix the model before blocking it")
            continue
        if owners:
            lam = min(stiff_owners[n] for n in owners)
            block = dict(states=list(members), why=f"pool with its buffers, stiff mode Re λ {lam:.4g} /ms owned by {', '.join(owners)}; rule 2")
            implicit.append(block)
            for n in members:
                cluster_of[n] = block
        else:
            notes.append(f"binding cluster {', '.join(members)} owns no stiff mode (threshold {stiff_threshold:g} /ms): left explicit")
    if stiff is None:
        probe = "" if stiff_threshold <= stiffness.THRESHOLD else f" --threshold {stiff_threshold:g}"
        judgment.append(f"no {stiff_file}: every binding cluster was treated as non-stiff; run `bimgrl stiffness{probe}` first")

    # 3. SR + junction share a block
    sr, junc = comps.sr, comps.junction
    if sr and junc:
        b_sr, b_j = cluster_of.get(sr), cluster_of.get(junc)
        if b_sr is not None and b_j is not None and b_sr is not b_j:
            merged = dict(states=b_sr["states"] + b_j["states"], why=f"{b_sr['why']} + {b_j['why']}; joined so that the release conserves Ca (rule 3)")
            implicit = [b for b in implicit if b is not b_sr and b is not b_j] + [merged]
            for n in merged["states"]:
                cluster_of[n] = merged
        elif b_sr is not None and b_j is None:
            b_sr["states"] = b_sr["states"] + [junc]
            b_sr["why"] += f"; {junc} joined as a bare pool so that the release conserves Ca (rule 3)"
            cluster_of[junc] = b_sr
        elif b_j is not None and b_sr is None:
            b_j["states"] = b_j["states"] + [sr]
            b_j["why"] += f"; {sr} joined as a bare pool so that the release conserves Ca (rule 3)"
            cluster_of[sr] = b_j
    elif junc and not sr:
        judgment.append("a junctional pool but no SR was named in meta.toml [compartments]: check")

    substep: list[dict] = []
    # 4. V block
    v = meta.voltage
    na = [g for g in comps.na_gates if kind_of.get(g) == "gate"]
    if na:
        substep.append(dict(states=[v, *na], h_min=V_H_MIN, mode="proportional",
                            rate_triggers=[dict(state=v, atol=V_ATOL)], why="V with the INa/INaL gates, refined in proportion to |dV/dt| (rule 4)"))
    else:
        judgment.append("no INa gates named in meta.toml [compartments].na_gates: no V block proposed")

    # 5. dyad
    if junc:
        dyad = [junc]
        if junc in cluster_of:
            dyad = list(cluster_of[junc]["states"])
        chain_blocks = [b for b in implicit if set(b["states"]) & release]
        for b in chain_blocks:
            dyad += [s for s in b["states"] if s not in dyad]
        dyad += [s for s in release if s not in dyad]
        if sr and sr not in dyad:
            dyad.append(sr)
            if sr in cluster_of:
                dyad += [s for s in cluster_of[sr]["states"] if s not in dyad]
            judgment.append(f"the SR ({sr}) is in the dyad block by default; drop it if the Ca transient does not improve (the paper kept it where it helped)")
        triggers = [dict(state=junc, rtol=DYAD_RTOL)] + [dict(state=o, atol=OPEN_ATOL) for o in comps.release_open if o in dyad]
        if not comps.release_open:
            judgment.append("no open states of the release channel named: the dyad block triggers only on the junctional Ca")
        substep.append(dict(states=dyad, h_min=DYAD_H_MIN, mode="binary", rate_triggers=triggers,
                            why="junctional Ca with its buffers, the release channel and the SR, refined while release is active (rule 5)"))
    else:
        notes.append("no junctional/dyadic pool: no dyad block (as for models whose Ca enters the cytosol directly)")

    # 6. nesting: an implicit block that shares states with a sub-step block goes inside it
    for b in implicit:
        for s in substep:
            shared = set(b["states"]) & set(s["states"])
            if shared and not set(b["states"]) <= set(s["states"]):
                s["states"] += [x for x in b["states"] if x not in s["states"]]
                s["why"] += f"; absorbed the implicit block {', '.join(b['states'])} (rule 6)"
    if nest_v_markov:
        vblock = next((s for s in substep if s["mode"] == "proportional"), None)
        if vblock:
            for b in implicit:
                if set(b["states"]) & feeds_v and not set(b["states"]) <= set(vblock["states"]):
                    vblock["states"] += b["states"]
                    vblock["why"] += f"; nested the voltage-dependent chain {', '.join(b['states'])}"
    else:
        for b in implicit:
            if "Markov" in b["why"] and set(b["states"]) & feeds_v and not any(set(b["states"]) <= set(s["states"]) for s in substep):
                judgment.append(f"the voltage-dependent chain {', '.join(b['states'])} is an implicit block outside the V block; nest it (--nest-v-markov) if its states go negative at large steps (Winslow 1999 L-type)")

    # gotranx's own suggestion for comparison
    suggested, suggest_msg = None, ""
    try:
        from gotranx.block_detection import DetectionSettings, blocks_to_toml, detect_blocks

        det = detect_blocks(ws.model.ode, DetectionSettings(dt=dt, membrane_potential=v))
        suggested = det
        (ws.out / f"blocks_suggested_dt{dt:g}.toml").write_text(blocks_to_toml(det))
    except Exception as e:  # noqa: BLE001
        suggest_msg = f"gotranx suggest-blocks failed: {type(e).__name__}: {e}"

    toml = render_toml(implicit, substep, meta.id, round_=1)
    # validate
    try:
        from gotranx.blocks import parse_blocks, validate_blocks
        import tomllib

        validate_blocks(ws.model.ode, parse_blocks(tomllib.loads(toml), strict=True))
        valid = "valid (gotranx validate_blocks)"
    except Exception as e:  # noqa: BLE001
        valid = f"INVALID: {e}"
        judgment.append(valid)

    result = dict(implicit=implicit, substep=substep, notes=notes, judgment=judgment, valid=valid, toml=toml, dt=dt,
                  suggested=_describe_detection(suggested) if suggested else None, suggest_msg=suggest_msg)
    return result


def _describe_detection(det) -> dict:
    return dict(
        implicit=[list(b.states) for b in det.blocks.implicit],
        substep=[dict(states=list(b.states), mode=b.mode, h_min=b.h_min,
                      rate_triggers=[dict(state=t.state, atol=t.atol, rtol=t.rtol) for t in b.rate_triggers]) for b in det.blocks.substep],
        reasons={",".join(k): v for k, v in det.reasons.items()},
    )


def render_toml(implicit, substep, model_id: str, round_: int) -> str:
    lines = [f"# Blocks for {model_id} (round {round_}; see NOTES.md), proposed by bimgrl propose from the", "# paper's rules; every block carries the finding that produced it.", ""]
    for b in implicit:
        lines += [f"# {b['why']}", "[[implicit_blocks]]", f"states = {_toml_list(b['states'])}", ""]
    for s in substep:
        lines += [f"# {s['why']}", "[[substep_blocks]]", f"states = {_toml_list(s['states'])}", f"h_min = {s['h_min']}", f'mode = "{s["mode"]}"']
        trig = ", ".join("{ " + ", ".join(f"{k} = {v!r}" if k != "state" else f'state = "{v}"' for k, v in t.items()) + " }" for t in s["rate_triggers"])
        lines += [f"rate_triggers = [{trig}]", ""]
    return "\n".join(lines)


def to_markdown(r: dict) -> str:
    out = ["# Proposed blocks", "", f"Validation: {r['valid']}", ""]
    out += ["## Implicit blocks", ""] + ([f"- {', '.join(b['states'])}: {b['why']}" for b in r["implicit"]] or ["- none"])
    out += ["", "## Sub-step blocks", ""] + ([f"- {', '.join(s['states'])} ({s['mode']}, h_min {s['h_min']}, triggers {s['rate_triggers']}): {s['why']}" for s in r["substep"]] or ["- none"])
    out += ["", "## Left explicit (by rule)", ""] + ([f"- {n}" for n in r["notes"]] or ["- nothing"])
    out += ["", "## Judgment calls for the user", ""] + ([f"- {j}" for j in r["judgment"]] or ["- none"])
    out += ["", f"## gotranx suggest-blocks at dt {r['dt']:g}, for comparison", ""]
    if r["suggested"]:
        s = r["suggested"]
        out += [f"- implicit: {b}" for b in s["implicit"]] or ["- implicit: none"]
        out += [f"- sub-step: {b['states']} ({b['mode']}, h_min {b['h_min']}, {b['rate_triggers']})" for b in s["substep"]] or ["- sub-step: none"]
        out += ["", "Reasons: " + "; ".join(f"{k}: {v}" for k, v in s["reasons"].items())]
        mine_i = {tuple(sorted(b["states"])) for b in r["implicit"]}
        theirs_i = {tuple(sorted(b)) for b in s["implicit"]}
        out += ["", "Differences: " + (
            "; ".join([f"only here: {list(b)}" for b in mine_i - theirs_i] + [f"only gotranx: {list(b)}" for b in theirs_i - mine_i]) or "none in the implicit blocks")]
        out.append("gotranx does not detect the dyad; its V block comes from the numeric fast-gate rule, which is the cross-check of the INa gate list.")
    else:
        out.append(r["suggest_msg"] or "not run")
    return "\n".join(out)


def run(ws: Workspace, dt: float = 0.1, stiff_threshold: float = -20.0, nest_v_markov: bool = False, write: bool = True) -> dict:
    r = propose(ws, dt, stiff_threshold, nest_v_markov)
    ws.write_json("propose.json", {k: v for k, v in r.items() if k != "toml"})
    ws.write_md("propose.md", to_markdown(r))
    target = ws.model_dir / "blocks.toml"
    path = ws.out / "blocks_proposed.toml"
    path.write_text(r["toml"])
    if write and not target.exists():
        target.write_text(r["toml"])
        print(f"wrote {target}")
    else:
        print(f"wrote {path} ({target} exists and was left alone)" if target.exists() else f"wrote {path}")
    return r
