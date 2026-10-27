"""The trace contract: what every scheme in every language writes, and what validation reads

* ``<name>.csv``: header ``t,<state names in model order>``; one row per sample; full precision.
  Optional extra columns ``aux_<intermediate>``.
* ``<name>.json`` (run meta): model, language, scheme, dt, reference (bool), solver {name, rtol,
  atol}, t_unit, beats, bcl, stimulus, status (ok|nan|unphysical|timeout|error), wall_s_total,
  wall_s_per_beat, n_steps, parameters (overrides), float_type, created.
* ``<name>_beats.csv`` (optional, protocol B): the state at the start of every beat, for drift.

See references/contracts.md.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class RunMeta:
    model: str
    scheme: str
    dt: float | None = None
    reference: bool = False
    solver: dict = field(default_factory=dict)
    language: str = "ode"
    t_unit: str = "ms"
    protocol: str = "A"
    beats: int = 1
    bcl: float | None = None
    stimulus: dict | None = None
    status: str = "ok"
    message: str = ""
    n_steps: int | None = None
    wall_s_total: float | None = None
    wall_s_per_beat: float | None = None
    parameters: dict = field(default_factory=dict)
    float_type: str = "float64"
    blocks: str | None = None
    created: str = field(default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%S"))


def write_trace(path: Path, t: np.ndarray, y: np.ndarray, names: list[str], meta: RunMeta,
                aux: dict[str, np.ndarray] | None = None, beats: tuple[np.ndarray, np.ndarray] | None = None) -> Path:
    path = Path(path)
    stem = str(path)[:-4] if str(path).endswith(".csv") else str(path)  # a stem may contain dots (dt0.25)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["t", *names] + [f"aux_{k}" for k in (aux or {})]
    data = np.column_stack([t, np.asarray(y)] + [np.asarray(v) for v in (aux or {}).values()])
    np.savetxt(stem + ".csv", data, delimiter=",", header=",".join(cols), comments="", fmt="%.17g")
    Path(stem + ".json").write_text(json.dumps(asdict(meta), indent=2) + "\n")
    if beats is not None:
        tb, yb = beats
        np.savetxt(stem + "_beats.csv", np.column_stack([tb, yb]), delimiter=",",
                   header=",".join(["t", *names]), comments="", fmt="%.17g")
    return Path(stem + ".csv")


@dataclass
class Trace:
    t: np.ndarray
    y: np.ndarray
    names: list[str]
    meta: RunMeta
    aux: dict[str, np.ndarray]
    beats: tuple[np.ndarray, np.ndarray] | None = None

    def column(self, name: str) -> np.ndarray:
        return self.y[:, self.names.index(name)]


def read_trace(path: Path) -> Trace:
    path = Path(path)
    stem = str(path)[:-4] if str(path).endswith(".csv") else str(path)
    csv = Path(stem + ".csv")
    with open(csv) as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(csv, delimiter=",", skiprows=1, ndmin=2)
    if header[0] != "t":
        raise SystemExit(f"{csv}: the first column must be t, got {header[0]!r}")
    names = [h for h in header[1:] if not h.startswith("aux_")]
    aux_names = [h for h in header[1:] if h.startswith("aux_")]
    y = data[:, 1 : 1 + len(names)]
    aux = {h[4:]: data[:, 1 + len(names) + i] for i, h in enumerate(aux_names)}
    meta_path = Path(stem + ".json")
    meta = RunMeta(**json.loads(meta_path.read_text())) if meta_path.exists() else RunMeta(model="?", scheme="?")
    beats = None
    bpath = Path(stem + "_beats.csv")
    if bpath.exists():
        b = np.loadtxt(bpath, delimiter=",", skiprows=1, ndmin=2)
        beats = (b[:, 0], b[:, 1:])
    return Trace(t=data[:, 0], y=y, names=names, meta=meta, aux=aux, beats=beats)
