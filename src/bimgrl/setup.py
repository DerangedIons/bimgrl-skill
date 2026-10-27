"""Install gotranx into the skill's environment from config.toml (or GOTRANX_SRC)"""

from __future__ import annotations

import os
import subprocess
import sys
import tomllib
from pathlib import Path


def skill_dir() -> Path:
    return Path(__file__).resolve().parents[2]


def run(root: Path | None = None) -> None:
    root = Path(root) if root else skill_dir()
    cfg = tomllib.loads((root / "config.toml").read_text())
    spec = os.environ.get("GOTRANX_SRC") or cfg["gotranx_spec"]
    python = sys.executable
    print(f"installing gotranx from {spec} into {python}")
    cmd = ["uv", "pip", "install", "--python", python, spec]
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        subprocess.run([python, "-m", "pip", "install", spec], check=True)
    except subprocess.CalledProcessError as e:
        raise SystemExit(
            f"gotranx could not be installed from {spec} (exit {e.returncode}). The block_rush_larsen "
            "scheme lives in DerangedIons/gotranx-private until it is upstreamed: you need access to "
            "that repository, or set GOTRANX_SRC to a local clone at the pinned commit."
        )
    out = subprocess.run([str(Path(python).with_name("gotranx")), "--version"], capture_output=True, text=True).stdout.strip()
    want = cfg.get("gotranx_version")
    print(out + ("" if not want or want in out else f" (config.toml expects {want})"))
