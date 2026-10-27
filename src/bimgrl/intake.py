"""Draft ``meta.toml`` from the equations alone, and say what the user must confirm

Heuristics (each is reported with its evidence, nothing is silently assumed):

* membrane potential: the state in mV, else the state that most derivatives depend on;
* cytosolic Ca / Na / K: a concentration state whose name, lower-cased without underscores,
  is ``cai`` / ``nai`` / ``ki`` (or ``cac``, ``nac``, ``kc``);
* SR Ca: the Ca concentration with the largest initial value (SR Ca is 100-10000x cytosolic);
* junctional / dyadic Ca: a Ca concentration other than the cytosol that the SR exchanges
  with directly (J[junction, sr] and J[sr, junction] both non-zero);
* release channel: dimensionless states coupled both ways to the junctional Ca and depending on
  SR Ca (RyR gates and Markov states; the Ca-dependent ICaL gates fail the SR test);
* INa gates: HH gates named m*, h*, j* or containing "na";
* stimulus: the half-open ``stim_*`` parameters if present, else any symbol with "stim" in its
  name is reported and the model is marked as needing the stimulus rewrite.
"""

from __future__ import annotations

import re
from pathlib import Path

from .meta import Compartments, ModelMeta, Stim, dumps_meta, save_meta
from .model import Model
from .structure import analyse as analyse_structure
from .workspace import Workspace, md_table

_NORMAL = {
    "calcium": ("cai", "cac", "cacyt", "cacytosol", "cai_free", "caifree"),
    "sodium": ("nai", "nac", "nacyt"),
    "potassium": ("ki", "kc", "kcyt"),
}


def _norm(name: str) -> str:
    return name.lower().replace("_", "")


def _is_conc(unit: str) -> bool:
    unit = unit.strip()
    return "M" in unit or "mol" in unit


_BOUND = re.compile(r"buf|trap|bound|csqn|calm|cam\b|camk|casig|trpn|tnc|myo|sig", re.I)


def _is_ca(name: str) -> bool:
    """A free Ca pool: named with "ca" and not like a buffer or a signalling trap"""
    return "ca" in name.lower() and not _BOUND.search(name)


def draft_meta(model_dir: Path, voltage: str | None = None) -> tuple[ModelMeta, dict]:
    """The draft meta and the evidence for each guess"""
    import gotranx
    from gotranx.block_detection import _membrane_potential
    from gotranx.linearization import jacobian

    model_dir = Path(model_dir)
    ode = gotranx.load_ode(model_dir / "model.ode")
    names = [s.name for s in ode.states]
    units = {s.name: (s.unit_str or "").strip() for s in ode.states}
    init = {s.name: float(s.value) for s in ode.states}
    jac = jacobian(ode)
    evidence: dict = {}

    # membrane potential
    if voltage is None:
        mv = [n for n in names if units[n] in ("mV", "millivolt")]
        if len(mv) == 1:
            voltage = mv[0]
            evidence["voltage"] = f"the only state in mV"
        else:
            voltage = _membrane_potential(None, names, jac)
            evidence["voltage"] = "the state that most other derivatives depend on (no unique mV state)"
    else:
        evidence["voltage"] = "given on the command line"

    conc = [n for n in names if _is_conc(units[n])]
    meta = ModelMeta(id=model_dir.name, voltage=voltage, language="ode")
    for key, patterns in _NORMAL.items():
        hits = [n for n in conc if _norm(n) in patterns]
        if hits:
            setattr(meta, key, hits[0])
            evidence[key] = f"name {hits[0]!r} matches {patterns[:2]}"
        else:
            evidence[key] = "no concentration state named like " + "/".join(patterns[:2])

    # compartments
    comps = Compartments()
    comps.cytosol = meta.calcium or ""
    ca_conc = [n for n in conc if _is_ca(n)]
    deps = {}
    for (row, col) in jac:
        deps.setdefault(row, set()).add(col)

    def coupled(a, b):
        return b in deps.get(a, ()) and a in deps.get(b, ())

    def exchanges(j, s):
        """`j` receives from `s` directly, or through a release state r that both depend on
        (O'Hara-Rudy writes its release as a flux state between the junction and the JSR)"""
        if coupled(j, s):
            return True
        return any(s in deps.get(r, ()) and r in deps.get(s, ()) for r in deps.get(j, ()) if r not in (j, s))

    if ca_conc:
        low_ref = init[comps.cytosol] if comps.cytosol else min(init[n] for n in ca_conc)
        high = [n for n in ca_conc if init[n] > 5 * low_ref]
        low = [n for n in ca_conc if n not in high and n != comps.cytosol]
        if high:
            # the junction is a low pool that exchanges directly with a high (SR) pool; the SR
            # is the high pool it exchanges with (a second high pool, NSR, is only noted)
            pairs = [(j, s) for j in low for s in high if exchanges(j, s)]
            if pairs:
                comps.junction, comps.sr = pairs[0]
                evidence["junction"] = f"low-Ca pool exchanging directly with a high-Ca (SR) pool: {pairs}"
                evidence["sr"] = f"high-Ca pool ({init[comps.sr]:.3g} mM) exchanging with {comps.junction}"
                if len(pairs) > 1:
                    evidence["junction"] += " (several: first taken, confirm)"
            else:
                comps.sr = max(high, key=lambda n: init[n])
                evidence["sr"] = f"largest initial Ca ({init[comps.sr]:.3g} mM)"
                evidence["junction"] = "no low-Ca pool exchanges directly with the SR: no dyadic pool (release goes to the cytosol)"
            other_high = [n for n in high if n != comps.sr]
            if other_high:
                evidence["sr"] += f"; other high-Ca pools (network SR?): {other_high}"
        else:
            evidence["sr"] = "no Ca state is 5x above the cytosolic value (no SR?)"
        remaining = [n for n in low if n != comps.junction]
        if remaining:
            comps.subsarcolemma = remaining[0]
            evidence["subsarcolemma"] = f"remaining low-Ca pool(s): {remaining}"
    dimless = [n for n in names if units[n] in ("", "1", "dimensionless")]
    if comps.junction:
        j, sr = comps.junction, comps.sr
        # release-channel states: coupled both ways to the junctional Ca and depending on the SR
        # Ca; any unit (ORd writes its release as flux states in mM/ms), but no Ca pool
        rc = [
            n for n in names
            if n not in ca_conc and not _BOUND.search(n) and n in deps.get(j, ()) and sr in deps.get(n, ())
            and (j in deps.get(n, ()) or n in deps.get(sr, ()))
        ]
        comps.release_channel = rc
        evidence["release_channel"] = (
            f"dimensionless states coupled both ways to {j} and depending on {sr}" if rc
            else f"no dimensionless state is coupled both ways to {j} and depends on {sr}"
        )
        # open states: the ones whose name says so, else all
        opens = [n for n in rc if re.search(r"(^|_)o(_|$)|open", n, re.I)]
        comps.release_open = opens or rc
    # INa gates
    gates = _hh_gates(ode, jac)
    na = [g for g in gates if re.match(r"^[mhj]([A-Z_]|$|p)", g) or "na" in g.lower()]
    comps.na_gates = na
    evidence["na_gates"] = f"HH gates named m*/h*/j* or *na*: {na}" if na else "no HH gate named like an INa gate"
    comps.status = "auto"
    meta.compartments = comps

    # stimulus
    pnames = {p.name for p in ode.parameters}
    if {"stim_start", "stim_duration", "stim_period", "stim_amplitude"} <= pnames:
        pv = {p.name: float(p.value) for p in ode.parameters}
        meta.stim = Stim(pv["stim_start"], pv["stim_duration"], pv["stim_period"], pv["stim_amplitude"])
        meta.bcl = pv["stim_period"]
        meta.protocol = "paced"
        meta.window_a = pv["stim_period"]
        evidence["stimulus"] = "half-open stim_* parameters found"
    else:
        stimmy = sorted({a.name for a in (*ode.parameters, *ode.intermediates) if "stim" in a.name.lower()})
        if stimmy:
            meta.protocol = "paced"
            evidence["stimulus"] = (
                f"symbols named like a stimulus: {stimmy}; the stimulus must be rewritten to the "
                "half-open stim_* parameters (bimgrl convert --stimulus or references/contracts.md)"
            )
        else:
            meta.protocol = "autorhythmic"
            meta.window_b = 10 * meta.bcl
            evidence["stimulus"] = "no stimulus symbol: assumed autorhythmic (confirm)"
    meta.kinds = {}
    meta.status = "draft"
    return meta, dict(evidence=evidence, units=units, init=init, gates=gates)


def _hh_gates(ode, jac) -> list[str]:
    from gotranx.block_detection import state_dependencies
    from gotranx.linearization import diagonal_jacobian

    from .structure import hh_gates

    return hh_gates(ode, state_dependencies(ode), diagonal_jacobian(ode))


def to_markdown(meta: ModelMeta, info: dict) -> str:
    ev = info["evidence"]
    out = [f"# Intake of `{meta.id}`", ""]
    out.append("Draft `meta.toml` written. Confirm every row marked *confirm* before `bimgrl propose`.")
    out += ["", "## Guesses and their evidence", ""]
    rows = [
        ["voltage", meta.voltage, ev.get("voltage", "")],
        ["calcium (cytosol)", meta.calcium or "-", ev.get("calcium", "")],
        ["sodium", meta.sodium or "-", ev.get("sodium", "")],
        ["potassium", meta.potassium or "-", ev.get("potassium", "")],
        ["sr", meta.compartments.sr or "-", ev.get("sr", "")],
        ["junction", meta.compartments.junction or "-", ev.get("junction", "") + " *confirm*"],
        ["subsarcolemma", meta.compartments.subsarcolemma or "-", ev.get("subsarcolemma", "")],
        ["release_channel", ", ".join(meta.compartments.release_channel) or "-", ev.get("release_channel", "") + " *confirm*"],
        ["release_open", ", ".join(meta.compartments.release_open) or "-", "name contains O/open, else all *confirm*"],
        ["na_gates", ", ".join(meta.compartments.na_gates) or "-", ev.get("na_gates", "") + " *confirm*"],
        ["stimulus", repr(meta.stim) if meta.stim else "-", ev.get("stimulus", "")],
        ["protocol", meta.protocol, "paced if a stimulus exists"],
    ]
    out.append(md_table(["field", "value", "evidence"], rows))
    out += ["", "## States", ""]
    out.append(
        md_table(
            ["state", "unit", "init", "HH gate"],
            [[n, info["units"][n], info["init"][n], "yes" if n in info["gates"] else ""] for n in info["units"]],
        )
    )
    return "\n".join(out)


def run(ws: Workspace, voltage: str | None = None, force: bool = False) -> ModelMeta:
    path = ws.model_dir / "meta.toml"
    if path.exists() and not force:
        raise SystemExit(f"{path} exists; use --force to overwrite the draft")
    meta, info = draft_meta(ws.model_dir, voltage)
    save_meta(meta, path)
    ws.write_md("intake.md", to_markdown(meta, info))
    ws.write_json("intake.json", dict(meta=dumps_meta(meta), **info))
    return meta
