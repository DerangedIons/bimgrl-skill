"""Accuracy of a run against the reference: biomarkers, thresholds, the physical range, and
the drift of every invariant. Works on arrays, so the in-process sweep and the trace-file
path (`bimgrl validate`) measure exactly alike."""

from __future__ import annotations

import math

import numpy as np

from .invariants import evaluate as evaluate_invariants
from .meta import ModelMeta
from .protocol import SAMPLE_MS, Protocol
from .vendored.features import autorhythmic_features, feature_errors, paced_features
from .vendored.metrics import bounds, classify, max_abs_dv, relrms, rrms_sd, sd_points
from .vendored.select import BIOMARKER_SETS, CYCLE_LENGTH

#: relative drift per beat above which an invariant is reported as not conserved. An exact
#: invariant (a closed Markov chain) is measured against its own initial value, so round-off
#: is the only legitimate change; a boundary one (total Ca up to membrane currents) against
#: the reference, where the integration error of the currents is legitimate, so the floor is
#: the level at which it starts to matter for the biomarkers
DRIFT_FLOOR = {"exact": 1e-7, "boundary": 1e-3}


def invariant_drift(invariants: list[dict], names: list[str], y_beats: np.ndarray, yref_beats: np.ndarray | None) -> dict:
    """Per invariant: the relative change over the run (``drift[name]``), and a flag string
    naming the invariants whose per-beat drift exceeds their floor. Closed chains (class
    closure) get a bound check on their sum instead; sub-groups that exchange with the rest
    of the model are not checked (they are not conserved alone)"""
    out: dict = {}
    q = evaluate_invariants(invariants, names, y_beats)
    qref = evaluate_invariants(invariants, names, yref_beats) if yref_beats is not None else {}
    flags = []
    for inv in invariants:
        name, cls = inv["name"], inv.get("class")
        if name not in q:
            continue
        if cls == "closure":
            total = float(inv.get("closure_total", 1.0))
            excess = float(np.max(q[name])) - total
            out[f"drift[{name}]"] = excess
            if excess > 1e-6:
                flags.append(f"{name}: sum exceeds {total:g} by {excess:.2e}")
            continue
        if cls == "exact":
            rel = (q[name] - q[name][0]) / max(abs(float(q[name][0])), 1e-300)
        elif cls == "boundary" and name in qref:
            rel = (q[name] - qref[name]) / max(abs(float(qref[name][0])), 1e-300)
        else:
            continue
        per_beat = np.diff(rel) if len(rel) > 1 else rel
        out[f"drift[{name}]"] = float(rel[-1])
        if np.max(np.abs(per_beat)) > DRIFT_FLOOR[cls]:
            flags.append(f"{name}: {np.max(np.abs(per_beat)):.2e}/beat")
    out["invariant_drift_flags"] = "; ".join(flags)
    return out


def features(meta: ModelMeta, names: list[str], t: np.ndarray, y: np.ndarray, autorhythmic: bool) -> dict:
    V = y[:, names.index(meta.voltage)]
    Ca = y[:, names.index(meta.calcium)] if meta.calcium else None
    Na = y[:, names.index(meta.sodium)] if meta.sodium else None
    if autorhythmic:
        return autorhythmic_features(t, V, Ca, Na)
    return paced_features(t, V, Ca, Na)


def passes(row: dict, set_name: str, kind: str, has_ca: bool) -> bool:
    if row.get("status") != "ok":
        return False
    tol = BIOMARKER_SETS[set_name]

    def within(key, limit):
        v = row.get(key)
        return v is not None and not (isinstance(v, float) and math.isnan(v)) and abs(float(v)) <= limit + 1e-9

    checks = [within("err_vpeak", tol["vpeak"]), within("err_apd90", tol["apd90"])]
    if has_ca:
        checks.append(within("err_cat", tol["cat"]))
    if kind == "autorhythmic":
        checks.append(within("err_cycle_length", CYCLE_LENGTH[set_name]))
    else:
        checks.append(within("err_tact", tol["tact"]))
    return all(checks)


def compare(meta: ModelMeta, names: list[str], protocol: Protocol, t: np.ndarray, y: np.ndarray,
            yref: np.ndarray, t_beats=None, y_beats=None, yref_beats=None, invariants: list[dict] | None = None) -> dict:
    """Metrics of a finished run (`t`, `y` on the compared window; `yref` at the same times)"""
    iv = names.index(meta.voltage)
    start = t[0]
    length = t[-1] - t[0]
    grid = t[1] - t[0]
    out: dict = {}
    points = sd_points(start, length)
    idx = (points - t[0]) / grid
    if np.all(np.abs(idx - np.round(idx)) < 1e-6):
        idx = np.round(idx).astype(int)
        idx = idx[idx < len(t)]
        out["rrms_sd"] = rrms_sd(y[idx, iv], yref[idx, iv])
    out["relrms"] = relrms(y[:, iv], yref[:, iv])
    out["max_abs_dv"] = max_abs_dv(y[:, iv], yref[:, iv])
    auto = meta.protocol == "autorhythmic" or protocol.stim is None
    actual, expected = features(meta, names, t - start, y, auto), features(meta, names, t - start, yref, auto)
    out.update(actual)
    out.update({f"ref_{k}": v for k, v in expected.items()})
    out.update({f"err_{k}": v for k, v in feature_errors(actual, expected).items()})
    kind = meta.protocol
    for set_name in ("tight", "moderate"):
        out[f"pass_{set_name}"] = passes(dict(out, status="ok"), set_name, kind, bool(meta.calcium))
    # drift of the slow states and of the invariants, at the start of every beat
    if t_beats is not None and len(t_beats) > 1 and yref_beats is not None:
        for key, state in (("nai", meta.sodium), ("ki", meta.potassium)):
            if state:
                d = y_beats[:, names.index(state)] - yref_beats[:, names.index(state)]
                out[f"drift_{key}"], out[f"drift_{key}_max"] = float(d[-1]), float(np.max(np.abs(d)))
        if invariants:
            out.update(invariant_drift(invariants, names, y_beats, yref_beats))
    return out


def status_of(y_all: np.ndarray, kinds: list[str], lo, hi, fail: int) -> str:
    status = classify(y_all, kinds, lo, hi)
    return "nan" if fail >= 0 else status


def violations(names, t, y, lo, hi, n=5):
    found = []
    for i, name in enumerate(names):
        col = y[:, i]
        bad = np.where(~np.isfinite(col) | (col < lo[i]) | (col > hi[i]))[0]
        if len(bad):
            k = bad[0]
            found.append(dict(t=float(t[k]), state=name, value=float(col[k]), lo=float(lo[i]), hi=float(hi[i])))
    return sorted(found, key=lambda d: d["t"])[:n]
