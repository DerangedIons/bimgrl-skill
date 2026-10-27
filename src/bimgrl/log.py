"""NOTES.md: the model's log, in the survey's section convention, and the round header of
blocks.toml. Every tuning round is appended as a structured entry so that the record of what
was changed and why survives the session."""

from __future__ import annotations

import re
import time
from pathlib import Path

from .workspace import Workspace

TEMPLATE = """# {name}

## Source

- (where the model came from: repository, file, commit, conversion)

## Edits

- (changes made to the equations before analysis, e.g. the stimulus rewrite)

## Protocol and sanity checks

- (protocol, reference biomarkers, self-convergence)

## Structure and invariants

- (see bimgrl/structure.md and bimgrl/invariants.md: chains, binding clusters, flags)

## Stiffness

- (see bimgrl/stiffness.md: the stiff modes, their owners, the predicted and measured limits)

## Blocks

- (the blocks and the rule or finding behind each; see bimgrl/propose.md)

## Tuning log

## Promising dt

- (GRL: ...; BIMGRL: ...)
"""


def ensure_notes(ws: Workspace) -> Path:
    path = ws.model_dir / "NOTES.md"
    if not path.exists():
        path.write_text(TEMPLATE.format(name=ws.meta.name or ws.meta.id))
    return path


def append_round(ws: Workspace, round_: int, change: str, result: str, decision: str, blocks_file: str = "blocks.toml") -> Path:
    path = ensure_notes(ws)
    text = path.read_text()
    stamp = time.strftime("%Y-%m-%d")
    entry = f"- Round {round_} ({stamp}, {blocks_file}): {change}. Result: {result}. Decision: {decision}\n"
    if "## Tuning log" in text:
        head, tail = text.split("## Tuning log", 1)
        # append at the end of the section (before the next heading)
        m = re.search(r"\n## ", tail[1:])
        if m:
            cut = m.start() + 1
            tail = tail[:cut].rstrip("\n") + "\n" + entry + "\n" + tail[cut:]
        else:
            tail = tail.rstrip("\n") + "\n" + entry
        text = head + "## Tuning log" + tail
    else:
        text = text.rstrip("\n") + "\n\n## Tuning log\n\n" + entry
    path.write_text(text)
    return path


def set_round_header(ws: Workspace, round_: int, blocks_file: str = "blocks.toml") -> Path:
    path = ws.model_dir / blocks_file
    text = path.read_text()
    header = f"# Blocks for {ws.meta.id} (round {round_}; see NOTES.md)"
    lines = text.splitlines()
    if lines and lines[0].startswith("# Blocks for"):
        lines[0] = header
    else:
        lines.insert(0, header)
    path.write_text("\n".join(lines) + "\n")
    return path


def add_section_line(ws: Workspace, section: str, line: str) -> Path:
    path = ensure_notes(ws)
    text = path.read_text()
    marker = f"## {section}"
    if marker not in text:
        text = text.rstrip("\n") + f"\n\n{marker}\n\n- {line}\n"
    else:
        head, tail = text.split(marker, 1)
        m = re.search(r"\n## ", tail[1:])
        if m:
            cut = m.start() + 1
            tail = tail[:cut].rstrip("\n") + f"\n- {line}\n\n" + tail[cut:]
        else:
            tail = tail.rstrip("\n") + f"\n- {line}\n"
        text = head + marker + tail
    path.write_text(text)
    return path
