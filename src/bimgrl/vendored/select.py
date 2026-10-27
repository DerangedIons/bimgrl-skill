# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/select.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""The fastest configuration of each scheme that meets each set of thresholds

The sets are those of S&D (``sd5``, ``sd1``: their RRMS of V below 5 % and 1 %, in the
standard form sqrt(sum (V - Vref)^2 / sum Vref^2) on their 100 points, which is what their Table 3
needs: RL reaches 5 % at dt 0.1-0.2 ms for TT06 and 0.5-1 ms for BR77, as S&D's 0.152 and
0.849 ms; the formula printed in the paper has 1/N inside the root and would pass every stable
run of these models at 1 ms. It is ``10 * rrms_sd``, since ``rrms_sd`` is the printed formula on
the same 100 points)
and the biomarker sets of ``benchmarks/block_rush_larsen`` (``tight``,
``moderate``). The biomarker sets leave out the calcium transient for a model
without a calcium state, and use the cycle length in place of the activation
time for an autorhythmic model, whose upstroke time accumulates the error of
every cycle before it. A dimensionless model has the S&D sets only.

The choice is by measured time per beat, not by the time step, since a lot of
sub-stepping can make a larger time step slower.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Mapping

#: The thresholds on the errors: mV, ms, relative, ms
BIOMARKER_SETS = {
    "tight": dict(vpeak=1.0, apd90=1.0, cat=0.005, tact=0.02),
    "moderate": dict(vpeak=2.0, apd90=2.0, cat=0.01, tact=0.05),
}
#: The cycle length thresholds of autorhythmic models, in ms
CYCLE_LENGTH = {"tight": 1.0, "moderate": 2.0}
SD_SETS = {"sd5": 0.05, "sd1": 0.01}
SETS = ("sd5", "sd1", "tight", "moderate")
BASELINES = ("FE", "RL", "GRL")


def _f(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _within(value, tol: float) -> bool:
    return abs(_f(value)) <= tol + 1e-9  # a NaN is not within


def sd_rrms(row: Mapping) -> float:
    """S&D's RRMS of V in the standard form on their 100 points: ``sqrt(N) * rrms_sd``, N = 100"""
    return 10.0 * _f(row.get("rrms_sd"))


def passes(row: Mapping, set_name: str, kind: str, has_ca: bool) -> bool:
    """Whether a run meets the set of thresholds `set_name`"""
    if row.get("status") != "ok":
        return False
    if set_name in SD_SETS:
        return sd_rrms(row) < SD_SETS[set_name]
    if kind == "dimensionless":
        return False
    tol = BIOMARKER_SETS[set_name]
    checks = [
        _within(row.get("err_vpeak"), tol["vpeak"]),
        _within(row.get("err_apd90"), tol["apd90"]),
    ]
    if has_ca:
        checks.append(_within(row.get("err_cat"), tol["cat"]))
    if kind == "autorhythmic":
        checks.append(_within(row.get("err_cycle_length"), CYCLE_LENGTH[set_name]))
    else:
        checks.append(_within(row.get("err_tact"), tol["tact"]))
    return all(checks)


def best(rows: Iterable[Mapping], models: Mapping[str, tuple[str, bool]]) -> list[dict]:
    """For every (model, protocol, beats, scheme, set) the passing run with the least
    ``time_per_beat_s``. `models` gives (protocol kind, has calcium) by model id"""
    groups: dict[tuple, list[Mapping]] = defaultdict(list)
    for row in rows:
        if math.isnan(_f(row.get("time_per_beat_s"))):
            continue
        groups[(row["model"], row["protocol"], str(row["beats"]), row["scheme"])].append(row)
    chosen = []
    for (model, protocol, beats, scheme), group in sorted(groups.items()):
        kind, has_ca = models[model]
        for set_name in SETS:
            ok = [r for r in group if passes(r, set_name, kind, has_ca)]
            if ok:
                row = min(ok, key=lambda r: _f(r["time_per_beat_s"]))
                chosen.append(dict(row, set=set_name))
    return chosen


def speedups(chosen: Iterable[Mapping]) -> list[dict]:
    """The chosen rows with ``speedup_vs_<baseline>`` = time(best baseline) / time(this row)"""
    chosen = list(chosen)
    times = {
        (r["model"], r["protocol"], str(r["beats"]), r["set"], r["scheme"]): _f(
            r["time_per_beat_s"]
        )
        for r in chosen
    }
    out = []
    for r in chosen:
        row = dict(r)
        for baseline in BASELINES:
            base = times.get(
                (r["model"], r["protocol"], str(r["beats"]), r["set"], baseline), math.nan
            )
            row[f"speedup_vs_{baseline}"] = base / _f(r["time_per_beat_s"])
        out.append(row)
    return out
