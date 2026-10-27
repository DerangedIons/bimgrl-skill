"""The sweep: a scheme at every allowed time step, coarsest first, against the reference

Rows go to ``bimgrl/runs.csv`` one at a time; a row that is already there (same scheme,
configuration, protocol, dt) is not run again, so an interrupted sweep resumes. A failure is a
row, never an exception. Protocol B integrates its first N-1 beats saving the state at the
start of every beat (drift of the slow states and of the invariants), then the last beat on
the 0.1 ms grid, which is compared with the reference.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import time
import traceback
from pathlib import Path

import numpy as np

from .model import Model
from .protocol import Protocol, allowed_dt, dt_base, protocols, save_every
from .stepper import SCHEMES, CStepper, SchemeConfig, compiler, scheme_config
from .trace import RunMeta, write_trace
from .validate import compare, passes, status_of, violations
from .vendored.stimulus import DT_GRID
from .vendored.sweep import bisect_k  # noqa: E402
from .workspace import Workspace, md_table

RUN_KEY = ("scheme", "config_key", "protocol", "beats", "dt")
COLUMNS = [
    *RUN_KEY, "blocks_file", "refined", "status", "message", "n_steps", "wall_s", "est_step_us",
    "time_per_beat_s", "pass_tight", "pass_moderate",
    "vpeak", "apd90", "cat", "tact", "cycle_length", "err_vpeak", "err_apd90", "err_cat", "err_tact", "err_cycle_length",
    "rrms_sd", "max_abs_dv", "drift_nai", "drift_ki", "invariant_drift_flags", "violations", "trace",
]


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, (float, np.floating)):
        return "" if math.isnan(v) else f"{float(v):.12g}"
    return str(v)


class RunStore:
    def __init__(self, path: Path):
        self.path = Path(path)
        self._rows: dict[tuple, dict] = {}
        if self.path.exists() and self.path.stat().st_size:
            with open(self.path, newline="") as f:
                for row in csv.DictReader(f):
                    self._rows[self.key(row)] = row

    @staticmethod
    def key(row: dict) -> tuple:
        return tuple(_fmt(row.get(k)) for k in RUN_KEY)

    def get(self, key):
        return self._rows.get(key)

    def rows(self) -> list[dict]:
        return list(self._rows.values())

    def append(self, row: dict) -> None:
        line = {k: _fmt(row.get(k)) for k in COLUMNS}
        extra = {k: _fmt(v) for k, v in row.items() if k.startswith("drift[")}
        # invariant drift columns are open-ended: keep them in a side file per key
        new_file = not self.path.exists() or self.path.stat().st_size == 0
        with open(self.path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS)
            if new_file:
                w.writeheader()
            w.writerow(line)
        if extra:
            side = self.path.with_name("runs_drift.csv")
            new_side = not side.exists() or side.stat().st_size == 0
            with open(side, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=[*RUN_KEY, "invariant", "drift"])
                if new_side:
                    w.writeheader()
                for k, v in extra.items():
                    w.writerow({**{kk: line[kk] for kk in RUN_KEY}, "invariant": k[6:-1], "drift": v})
        self._rows[self.key(line)] = line


def simulate(stepper: CStepper, protocol: Protocol, dt: float, p: np.ndarray, u0=None) -> dict:
    every = save_every(dt)
    if protocol.period is None:
        n = round(protocol.t1 / dt)
        t, y, fail = stepper.run(dt, n, every, p, u0=u0)
        return dict(t=t, y=y, fail=fail, t_beats=t[:1], y_beats=y[:1], n_steps=n)
    per_beat = round(protocol.period / dt)
    t_b, y_b, fail = stepper.run(dt, per_beat * (protocol.beats - 1), per_beat, p, u0=u0)
    if fail >= 0 or not np.all(np.isfinite(y_b[-1])):
        return dict(t=t_b, y=y_b, fail=max(fail, 0), t_beats=t_b, y_beats=y_b, n_steps=per_beat * protocol.beats)
    t, y, fail = stepper.run(dt, per_beat, every, p, u0=y_b[-1])
    t_beats = np.concatenate([t_b, [t_b[-1] + protocol.period]])
    y_beats = np.vstack([y_b, y[-1:]])
    return dict(t=t + t_b[-1], y=y, fail=fail, t_beats=t_beats, y_beats=y_beats, n_steps=per_beat * protocol.beats)


class Runner:
    """Everything one sweep shares: the model, the reference, the compiled steppers"""

    def __init__(self, ws: Workspace, blocks_file: Path | None = None, suggested_file: Path | None = None, save_traces: bool = False):
        self.ws = ws
        self.model: Model = ws.model
        self.blocks_file = Path(blocks_file) if blocks_file else ws.model_dir / "blocks.toml"
        self.suggested_file = suggested_file
        self.build_dir = ws.out / "build"
        self.cc = compiler()
        self.store = RunStore(ws.out / "runs.csv")
        self.save_traces = save_traces
        self._steppers: dict[str, CStepper] = {}
        inv_path = ws.out / "invariants.json"
        self.invariants = ws.read_json("invariants.json")["invariants"] if inv_path.exists() else []

    def config(self, label: str) -> SchemeConfig | None:
        return scheme_config(label, self.model, self.blocks_file, self.suggested_file)

    def stepper(self, config: SchemeConfig) -> CStepper:
        if config.key not in self._steppers:
            self._steppers[config.key] = CStepper(self.model.ode, config, self.build_dir, self.cc)
        return self._steppers[config.key]

    def run_one(self, protocol: Protocol, label: str, dt: float, budget: float = 600.0, refined: str = "", force: bool = False) -> dict:
        model = self.model
        row = dict(scheme=label, protocol=protocol.name, beats=protocol.beats, dt=dt, refined=refined, status="error", message="")
        start = time.perf_counter()
        try:
            config = self.config(label)
            if config is None:
                row.update(status="n/a", message="no implicit blocks", config_key="")
                return row
            row["config_key"] = config.key
            row["blocks_file"] = config.blocks_file or ""
            existing = self.store.get(self.store.key(row))
            if existing is not None and existing["status"] not in ("error",) and not force:
                return existing
            stepper = self.stepper(config)
            p = model.parameters(protocol, stepper.parameters)
            # cost estimate from a short calibration, to respect the budget
            n_cal = min(2000, round(protocol.t1 / dt))
            est_step = stepper.time(dt, n_cal, repeats=1, parameters=p) / n_cal
            row["est_step_us"] = 1e6 * est_step
            n_total = round(protocol.t1 / dt)
            if est_step * n_total > budget:
                row.update(status="budget", message=f"estimated {est_step * n_total:.0f} s > budget {budget:.0f} s", n_steps=n_total)
                return row
            sim = simulate(stepper, protocol, dt, p)
            row["n_steps"] = sim["n_steps"]
            lo, hi = model.bounds(protocol)
            y_all = np.vstack([sim["y_beats"], sim["y"]])
            status = status_of(y_all, model.kinds, lo, hi, sim["fail"])
            row["status"] = status
            if status != "ok":
                v = violations(model.names, np.concatenate([sim["t_beats"], sim["t"]]), y_all, lo, hi)
                row["violations"] = "; ".join(f"t={d['t']:.3g} {d['state']}={d['value']:.4g} outside [{d['lo']:.3g},{d['hi']:.3g}]" for d in v)
            else:
                ref = model.reference(protocol)
                yref = ref(sim["t"]).T
                yref_b = ref(sim["t_beats"]).T
                row.update(compare(model.meta, model.names, protocol, sim["t"], sim["y"], yref, sim["t_beats"], sim["y_beats"], yref_b, self.invariants))
                if self.save_traces:
                    tdir = self.ws.out / "traces"
                    meta = RunMeta(model=model.meta.id, scheme=label, dt=dt, protocol=protocol.name, beats=protocol.beats,
                                   bcl=protocol.period, stimulus=None if protocol.stim is None else vars(protocol.stim),
                                   status=status, n_steps=sim["n_steps"], blocks=config.blocks_file, language="c")
                    path = write_trace(tdir / f"{label}_{protocol.name}{protocol.beats}_dt{dt:g}", sim["t"], sim["y"], model.names, meta,
                                       beats=(sim["t_beats"], sim["y_beats"]))
                    row["trace"] = str(path.relative_to(self.ws.model_dir))
        except SystemExit as e:
            row.update(status="error", message=str(e)[-500:])
        except Exception as e:  # noqa: BLE001 - a failure is a row
            row.update(status="error", message=f"{type(e).__name__}: {e}"[-500:])
            traceback.print_exc()
        row["wall_s"] = time.perf_counter() - start
        row.setdefault("config_key", "")
        self.store.append(row)
        return row

    def time_row(self, protocol: Protocol, label: str, dt: float, beats: int = 1, repeats: int = 5) -> float:
        config = self.config(label)
        stepper = self.stepper(config)
        p = self.model.parameters(protocol, stepper.parameters)
        length = protocol.period if protocol.period is not None else protocol.t1
        n = round(length / dt) * beats
        return stepper.time(dt, n, repeats=repeats, parameters=p) / beats


def describe(row: dict) -> str:
    s = f"{row['scheme']:12s} dt {float(row['dt']):<8g} {row['status']:10s}"
    if row.get("est_step_us") not in (None, ""):
        s += f" {float(row['est_step_us']):7.2f} us/step"
    if row["status"] == "ok":
        def f(k, fmt):
            v = row.get(k)
            return fmt % float(v) if v not in (None, "") else "   -  "
        s += f"  dVpeak {f('err_vpeak', '%+6.2f')} dAPD90 {f('err_apd90', '%+6.2f')} dCaT {f('err_cat', '%+7.4f')} dtact {f('err_tact', '%+6.3f')}"
        s += f"  tight {row.get('pass_tight')} moderate {row.get('pass_moderate')}"
        if row.get("invariant_drift_flags"):
            s += f"  DRIFT {row['invariant_drift_flags']}"
    elif row.get("violations"):
        s += f"  {row['violations']}"
    elif row.get("message"):
        s += f"  {row['message']}"
    return s


def sweep(ws: Workspace, labels=("GRL", "BIMGRL"), protocol_name: str = "B", beats: int = 10, dts=None, budget: float = 600.0,
          blocks_file=None, suggested_file=None, stability: bool = False, save_traces: bool = False, force: bool = False) -> list[dict]:
    runner = Runner(ws, blocks_file, suggested_file, save_traces)
    protocol = protocols(ws.meta, beats)[protocol_name]
    grid = sorted((dt for dt in (dts or DT_GRID) if allowed_dt(dt, protocol)), reverse=True)
    if dts and len(grid) < len(dts):
        print(f"dropped time steps that do not line up with the stimulus/grid: {sorted(set(dts) - set(grid))}")
    rows = []
    kind, has_ca = ws.meta.protocol, bool(ws.meta.calcium)
    for label in labels:
        streak = 0
        for dt in grid:
            row = runner.run_one(protocol, label, dt, budget, force=force)
            rows.append(row)
            print(describe(row), flush=True)
            if row["status"] in ("budget", "n/a"):
                break
            ok = row["status"] == "ok" and passes(row, "tight", kind, has_ca)
            streak = streak + 1 if ok else 0
            if streak >= 2 and not dts:
                break
        if stability and label in ("GRL", "RL", "FE"):
            stability_limit(runner, protocol, label, [r for r in rows if r["scheme"] == label], budget)
    return rows


def stability_limit(runner: Runner, protocol: Protocol, label: str, rows: list[dict], budget: float) -> float | None:
    """Bisect the largest stable time step of an explicit scheme between the largest ok and the
    smallest failing time step already run"""
    tested = sorted((r for r in rows if r["status"] in ("ok", "nan", "unphysical")), key=lambda r: float(r["dt"]))
    ok = [float(r["dt"]) for r in tested if r["status"] == "ok"]
    bad = [float(r["dt"]) for r in tested if r["status"] != "ok" and ok and float(r["dt"]) > max(ok)]
    if not ok or not bad:
        print(f"{label}: no bracket for the stability limit (ok {ok}, failing {bad})")
        return None
    base = dt_base(protocol)
    k_pass, k_fail = round(base / max(ok)), round(base / min(bad))

    def passes_k(k: int) -> bool:
        row = runner.run_one(protocol, label, base / k, budget, refined="stability")
        print(describe(row), flush=True)
        return row["status"] == "ok"

    k = bisect_k(passes_k, k_pass, k_fail, lambda k: allowed_dt(base / k, protocol))
    limit = base / k
    print(f"{label} stability limit: {limit:g} ms (largest allowed step that runs; {base / (k - 1) if k > 1 else float('nan'):g} ms fails)")
    out = runner.ws.out / "stability.json"
    data = {}
    if out.exists():
        import json
        data = json.loads(out.read_text())
    data[label] = dict(limit_ms=limit, protocol=protocol.name, beats=protocol.beats, fails_at_ms=base / (k - 1) if k > 1 else None)
    runner.ws.write_json("stability.json", data)
    return limit


def time_rows(ws: Workspace, labels, protocol_name: str, beats_timed: int = 1, repeats: int = 5, blocks_file=None, max_one_beat_s: float = 120.0):
    """Measure the time per beat of every ok row of the sweep, serially; writes runs.csv in place"""
    runner = Runner(ws, blocks_file)
    protocol = protocols(ws.meta, 10)[protocol_name]
    rows = runner.store.rows()
    for row in rows:
        if row["scheme"] not in labels or row["status"] != "ok" or row["protocol"] != protocol_name:
            continue
        if row.get("time_per_beat_s"):
            continue
        dt = float(row["dt"])
        est = float(row["est_step_us"]) * 1e-6 * round((protocol.period or protocol.t1) / dt)
        if est > max_one_beat_s:
            row["time_per_beat_s"] = est
            print(f"{row['scheme']} dt {dt:g}: untimed, estimated {est:.1f} s per beat")
        else:
            row["time_per_beat_s"] = runner.time_row(protocol, row["scheme"], dt, beats_timed, repeats)
            print(f"{row['scheme']} dt {dt:g}: {row['time_per_beat_s'] * 1e3:.3f} ms per beat (best of {repeats})")
    _rewrite(runner.store, rows)


def _rewrite(store: RunStore, rows: list[dict]) -> None:
    with open(store.path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in COLUMNS})
