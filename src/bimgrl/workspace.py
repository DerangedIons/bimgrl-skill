"""Where a model's files live: ``<model_dir>/model.ode``, ``meta.toml``, ``blocks.toml``,
``NOTES.md`` and the analysis outputs under ``<model_dir>/bimgrl/``"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .meta import ModelMeta, load_meta
from .model import Model


class _Encoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            v = float(o)
            return None if math.isnan(v) else v
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (set, frozenset, tuple)):
            return sorted(o) if isinstance(o, (set, frozenset)) else list(o)
        if isinstance(o, complex):
            return {"re": o.real, "im": o.imag}
        if isinstance(o, Path):
            return str(o)
        return str(o)


def _clean(x):
    """NaN and inf become null (JSON has no spelling for them); containers are walked"""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        v = float(x)
        return None if (math.isnan(v) or math.isinf(v)) else v
    return x


class Workspace:
    def __init__(self, model_dir: Path):
        self.model_dir = Path(model_dir).resolve()
        self.out = self.model_dir / "bimgrl"
        self.out.mkdir(parents=True, exist_ok=True)
        self._model: Model | None = None
        self._meta: ModelMeta | None = None

    @property
    def ode_path(self) -> Path:
        return self.model_dir / "model.ode"

    @property
    def meta(self) -> ModelMeta:
        if self._meta is None:
            path = self.model_dir / "meta.toml"
            if not path.exists():
                raise SystemExit(f"{path} not found: run `bimgrl intake {self.model_dir}` first")
            self._meta = load_meta(path)
        return self._meta

    @property
    def model(self) -> Model:
        if self._model is None:
            if not self.ode_path.exists():
                raise SystemExit(f"{self.ode_path} not found (canonical form is a gotranx .ode)")
            self._model = Model.load(self.ode_path, self.meta)
        return self._model

    def write_json(self, name: str, data: Any) -> Path:
        path = self.out / name
        path.write_text(json.dumps(_clean(data), indent=2, cls=_Encoder, allow_nan=False) + "\n")
        return path

    def read_json(self, name: str) -> Any:
        path = self.out / name
        if not path.exists():
            raise SystemExit(f"{path} not found: run the step that produces it first")
        return json.loads(path.read_text())

    def write_md(self, name: str, text: str) -> Path:
        path = self.out / name
        path.write_text(text.rstrip() + "\n")
        return path


def md_table(headers: list[str], rows: list[list[Any]]) -> str:
    def cell(v):
        if isinstance(v, float):
            if math.isnan(v):
                return "nan"
            return f"{v:.4g}"
        return str(v).replace("|", "\\|")

    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(cell(v) for v in row) + " |" for row in rows]
    return "\n".join(out)
