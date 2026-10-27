"""Generate the scheme code with gotranx for the user's target (C, Python/numpy, Julia,
CytoZoo), with the blocks of blocks.toml"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from .stepper import GOTRANX_SCHEME
from .workspace import Workspace

BACKENDS = {"c": "ode2c", "python": "ode2py", "julia": "ode2julia", "cytozoo": "ode2cytozoo"}


def run(ws: Workspace, backend: str, label: str, out: Path | None, blocks: Path | None = None, extra: list[str] = ()) -> Path:
    scheme = GOTRANX_SCHEME[label]
    blocks = Path(blocks) if blocks else ws.model_dir / "blocks.toml"
    out = Path(out) if out else ws.out / "generated" / f"{ws.meta.id}_{label.lower().replace('-', '_')}"
    out.parent.mkdir(parents=True, exist_ok=True)
    gotranx = Path(sys.executable).with_name("gotranx")
    cmd = [str(gotranx), BACKENDS[backend], str(ws.ode_path), "--scheme", scheme]
    if scheme == "block_rush_larsen":
        cmd += ["--blocks", str(blocks)]
    if backend == "cytozoo":
        cmd += ["--v-name", ws.meta.voltage]
    cmd += ["-o", str(out), *extra]
    print(" ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(result.stdout[-2000:] + result.stderr[-3000:])
    produced = sorted(out.parent.glob(out.name + "*"))
    print("wrote " + ", ".join(str(p) for p in produced))
    return out
