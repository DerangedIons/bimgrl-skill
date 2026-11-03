"""``meta.toml``: what the skill knows about a model that is not in its equations

The schema follows the ``ModelMeta`` of the cell-model survey (blocked-rush-larsen @ 58af805),
extended with the compartment naming that the block rules need and that no heuristic gets right
every time (``[compartments]``), which the intake gate asks the user to confirm.
"""

from __future__ import annotations

import json
import math
import tomllib
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Literal

Protocol = Literal["paced", "autorhythmic", "dimensionless"]


@dataclass(frozen=True)
class Stim:
    """A periodic stimulus on ``[start + k period, start + k period + duration)``, in ms; the
    amplitude is in the unit of the model's stimulus current"""

    start: float
    duration: float
    period: float
    amplitude: float


@dataclass
class Compartments:
    """The Ca compartments and the release channel, by state name. Empty strings mean "none"."""

    #: cytosolic free Ca (the biomarker Ca)
    cytosol: str = ""
    #: junctional / dyadic / subspace free Ca, the pool the release channel feeds
    junction: str = ""
    #: sub-sarcolemmal free Ca (Shannon/Bers/Grandi lineage), if separate from the cytosol
    subsarcolemma: str = ""
    #: SR free Ca
    sr: str = ""
    #: states of the release channel (RyR) whose open fraction gates release; Markov chains
    #: list all their states, a single gate lists itself
    release_channel: list[str] = field(default_factory=list)
    #: the open states of the release channel (for the dyad trigger), subset of release_channel
    release_open: list[str] = field(default_factory=list)
    #: gates of the fast Na current (and INaL): the V sub-step block
    na_gates: list[str] = field(default_factory=list)
    #: how the compartments were filled: "auto", "confirmed", or free text
    status: str = "auto"


@dataclass
class ModelMeta:
    id: str
    name: str = ""
    citation: str = ""
    doi: str = ""
    species: str = ""
    cell_type: str = ""
    source: str = ""
    language: str = "ode"
    time_unit: str = "ms"
    voltage_unit: str = "mV"
    protocol: Protocol = "paced"
    voltage: str = "V"
    calcium: str | None = None
    sodium: str | None = None
    potassium: str | None = None
    #: the qualified myokit name of the stimulus current (CellML path), e.g. "membrane.i_stim"
    stim_variable: str | None = None
    stim: Stim | None = None
    bcl: float = 1000.0
    window_a: float = 0.0
    window_b: float = 0.0
    unphysical_exempt: list[str] = field(default_factory=list)
    remove_singularities: bool = True
    jacobian: str = "exact"
    rl_gates: list[str] | None = None
    #: state kinds that override the unit-based rule: name -> voltage|gate|markov|conc|other
    kinds: dict[str, str] = field(default_factory=dict)
    compartments: Compartments = field(default_factory=Compartments)
    edits: list[str] = field(default_factory=list)
    status: str = "draft"
    notes: str = ""

    @property
    def has_ca(self) -> bool:
        return bool(self.calcium)


def _from_dict(data: dict[str, Any]) -> ModelMeta:
    data = dict(data)
    known = {f.name for f in fields(ModelMeta)}
    unknown = set(data) - known
    if unknown:
        raise ValueError(f"Unknown keys in meta.toml: {sorted(unknown)}")
    if data.get("stim") is not None:
        data["stim"] = Stim(**data["stim"])
    if data.get("compartments") is not None:
        data["compartments"] = Compartments(**data["compartments"])
    return ModelMeta(**data)


def load_meta(path: Path) -> ModelMeta:
    path = Path(path)
    if path.is_dir():
        path = path / "meta.toml"
    with open(path, "rb") as f:
        return _from_dict(tomllib.load(f))


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isinf(value) or math.isnan(value):
            raise ValueError(f"Cannot write {value} to meta.toml")
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, dict):
        return "{ " + ", ".join(f"{k} = {_toml_value(v)}" for k, v in value.items()) + " }"
    if isinstance(value, (list, tuple)):
        if any(isinstance(v, dict) for v in value) or sum(len(str(v)) for v in value) > 60:
            return "[\n" + "".join(f"    {_toml_value(v)},\n" for v in value) + "]"
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    raise TypeError(f"Cannot write {value!r} to meta.toml")


def dumps_meta(meta: ModelMeta) -> str:
    data = asdict(meta)
    stim = data.pop("stim")
    comps = data.pop("compartments")
    kinds = data.pop("kinds")
    lines = [f"{key} = {_toml_value(value)}" for key, value in data.items() if value is not None]
    if kinds:
        lines += ["", "[kinds]"] + [f"{k} = {_toml_value(v)}" for k, v in kinds.items()]
    if stim is not None:
        lines += ["", "[stim]"] + [f"{key} = {_toml_value(value)}" for key, value in stim.items()]
    lines += ["", "[compartments]"] + [f"{k} = {_toml_value(v)}" for k, v in comps.items()]
    return "\n".join(lines) + "\n"


def save_meta(meta: ModelMeta, path: Path) -> None:
    path = Path(path)
    if path.is_dir():
        path = path / "meta.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps_meta(meta))
