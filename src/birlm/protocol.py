"""Protocols A (one window) and B (N beats, the last compared), as in the survey"""

from __future__ import annotations

from dataclasses import dataclass

from .meta import ModelMeta, Stim
from .vendored.stimulus import aligned

#: the grid that features and errors are computed on, or the time step if it is longer
SAMPLE_MS = 0.1


@dataclass(frozen=True)
class Protocol:
    name: str
    beats: int
    #: the length of a beat (paced B), None if the whole run is one piece
    period: float | None
    t1: float
    stim: Stim | None


def protocols(meta: ModelMeta, beats_b: int = 10) -> dict[str, Protocol]:
    out = {}
    stim_a = meta.stim if meta.protocol != "autorhythmic" else None
    window_a = meta.window_a or (meta.bcl if meta.stim else 1000.0)
    out["A"] = Protocol("A", 1, None, window_a, stim_a)
    if meta.protocol == "autorhythmic" or meta.stim is None:
        window_b = meta.window_b or 10 * meta.bcl
        out["B"] = Protocol("B", beats_b, None, window_b * beats_b / 10, None)
    else:
        stim = Stim(meta.stim.start, meta.stim.duration, meta.bcl, meta.stim.amplitude)
        out["B"] = Protocol("B", beats_b, meta.bcl, beats_b * meta.bcl, stim)
    return out


def allowed_dt(dt: float, protocol: Protocol) -> bool:
    """Whether `dt` lines up with the stimulus, the window and the sampling grid"""
    values = [protocol.t1] + ([] if protocol.period is None else [protocol.period])
    if protocol.stim is not None:
        values += [protocol.stim.start, protocol.stim.duration, protocol.stim.period]
    grid = SAMPLE_MS if dt <= SAMPLE_MS + 1e-12 else 1.0
    return aligned(dt, *values) and (dt >= 1.0 or aligned(dt, grid))


def dt_base(protocol: Protocol) -> float:
    return protocol.stim.duration if protocol.stim is not None else 1.0


def save_every(dt: float) -> int:
    return max(1, round(SAMPLE_MS / dt)) if dt <= SAMPLE_MS + 1e-12 else 1


def beat_starts(protocol: Protocol) -> list[float]:
    """The start of every beat of protocol B (0, bcl, 2 bcl, ...), or [0] for one window"""
    if protocol.period is None:
        return [0.0]
    return [k * protocol.period for k in range(protocol.beats + 1)]
