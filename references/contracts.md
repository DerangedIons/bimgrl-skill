# Contracts: the .ode form, meta.toml, the RHS adapter, the trace

## The canonical `.ode`

gotranx's ODE format, with three conventions the pipeline relies on.

**Units.** Time in ms, membrane potential in mV, concentrations in mM (any unit containing `M` or `mol` is a concentration for the physical-range check; `1` or an empty unit is a gate). `bimgrl convert` converts CellML time and voltage with myokit; a hand translation must do the same.

**The stimulus.** Every paced model carries its stimulus as five parameters and a half-open interval, so that a step is stimulated if and only if it starts inside the interval and exactly `duration / dt` steps are stimulated per period:

```
parameters("membrane",
stim_start=ScalarParam(0.0, unit="ms"),
stim_duration=ScalarParam(1.0, unit="ms"),
stim_period=ScalarParam(1000.0, unit="ms"),
stim_amplitude=ScalarParam(-53.0, unit="A/F"),
stim_eps=ScalarParam(1e-9, unit="ms"))

expressions("membrane")
stim_phase = t - stim_start
stim_phase_mod = stim_phase - stim_period * floor((stim_phase + stim_eps) / stim_period)
Istim = Conditional(And(Ge(stim_eps, -stim_phase), Gt(stim_duration, stim_eps + stim_phase_mod)), stim_amplitude, 0.0)
```

A closed interval with an accumulated time stimulates one step too many in some periods. The reference integrator cuts the interval at every switch and forces the stimulus on or off through these parameters, so it never sees the discontinuity; the sweep sets them per protocol.

**Names.** Keep the model's own state and parameter names. The generated code orders states by name (uppercase first), so never assume V is first.

## `meta.toml`

Written as a draft by `bimgrl intake`, confirmed by the user. Keys (from `src/bimgrl/meta.py`):

| key | meaning |
|---|---|
| `id`, `name`, `citation`, `doi`, `species`, `cell_type`, `source`, `language` | provenance; `language` is `ode`, `cellml`, `julia`, `python`, `c`, `matlab` |
| `time_unit`, `voltage_unit` | of the original model (the `.ode` is in ms and mV) |
| `protocol` | `paced`, `autorhythmic` or `dimensionless` |
| `voltage`, `calcium`, `sodium`, `potassium` | state names (cytosolic concentrations) |
| `stim_variable` | the qualified myokit name of the stimulus current, for CellML conversion |
| `[stim]` | `start`, `duration`, `period`, `amplitude` |
| `bcl`, `window_a`, `window_b` | protocol B cycle length; the window of protocol A; the window of ten beats for an autorhythmic model |
| `unphysical_exempt` | states left out of the physical-range check (contraction states that may go negative) |
| `remove_singularities`, `jacobian` | gotranx load option; `exact` (symbolic) or `numeric` (central differences) |
| `rl_gates` | override of the HH-gate rule for the RL scheme |
| `[kinds]` | overrides of the unit-based state kind: `name = "voltage" | "gate" | "markov" | "conc" | "other"` |
| `[compartments]` | see below |
| `edits`, `status`, `notes` | the edits made to the equations; `draft` until confirmed |

```toml
[compartments]
cytosol = "cai"                # the biomarker Ca
junction = "caj"               # the pool the release channel feeds (dyad, subspace)
subsarcolemma = "casl"         # Shannon/Bers/Grandi lineage; empty otherwise
sr = "casr"
release_channel = ["ryr_R", "ryr_O", "ryr_I", "ryr_CaRI"]   # all states of the RyR chain, or the single gate
release_open = ["ryr_O"]       # the open states, for the dyad trigger
na_gates = ["m", "h", "j", "mL", "hL"]
status = "confirmed"           # "auto" until a human has checked it
```

## The RHS adapter (native path)

A model in any language is served by one long-running process that reads one JSON object per line on stdin and answers one per line on stdout. Batching matters: a finite-difference Jacobian of a 93-state model at 13 checkpoints is 13 requests of 186 points each, not 2500 round trips, and one process start amortises Julia or MATLAB latency.

| request | reply |
|---|---|
| `{"op": "info"}` | `{"n": 93, "states": [...names in order...], "y0": [...], "parameters": {"name": value, ...}, "units": {"state": "mM", ...} or null, "jac": true or false}` |
| `{"op": "rhs", "t": [t_1, ..., t_m], "y": [[...n...], ...m rows...]}` | `{"dy": [[...n...], ...m rows...]}` |
| `{"op": "jac", "t": [...], "y": [[...], ...]}` | `{"J": [[[...n×n...]], ...m matrices...]}` (only if `info` said `"jac": true`; else the harness takes central differences through `rhs`) |
| `{"op": "set_params", "parameters": {"name": value}}` | `{"ok": true}` |
| `{"op": "quit"}` | the process exits |
| any error | `{"error": "message"}`, which aborts the step with that message |

Finite-difference steps are $h_j = \sqrt{\epsilon}\, \max(|y_j|, s_j)$ with $s_j$ 1 mV for the voltage, 1e-3 for a gate and 1e-4 mM for a concentration. A non-finite derivative in any reply fails the call and names the checkpoint. Pass the adapter as `--adapter cmd:"julia --project=... adapter.jl mymodel.jl"`; the command runs with the model directory as its working directory.

The second mechanism is a shared library, `--adapter so:<library>:<meta.json>`, with the C signature

```c
int rhs(double t, const double *y, double *dy, const double *p);   /* return 0 on success */
```

and `meta.json` holding `states`, `y0`, `parameters` (an array in the order the library expects) and optionally `units`.

The Julia helper `scripts/jl/adapter.jl` serves the protocol for any `f(du, u, p, t)` with ForwardDiff for the Jacobian. A minimal Python adapter:

```python
#!/usr/bin/env python3
import json, sys
import numpy as np
from mymodel import NAMES, Y0, P, rhs   # rhs(t, y, p) -> dy

for line in sys.stdin:
    req = json.loads(line)
    try:
        if req["op"] == "info":
            out = {"n": len(NAMES), "states": NAMES, "y0": list(Y0), "parameters": P, "units": None, "jac": False}
        elif req["op"] == "rhs":
            out = {"dy": [rhs(t, np.array(y), P).tolist() for t, y in zip(req["t"], req["y"])]}
        elif req["op"] == "set_params":
            P.update(req["parameters"]); out = {"ok": True}
        elif req["op"] == "quit":
            break
        else:
            out = {"error": f"unknown op {req['op']}"}
    except Exception as e:
        out = {"error": repr(e)}
    print(json.dumps(out), flush=True)
```

A MATLAB or C adapter has the same shape: a loop reading a line (`fgetl(stdin)` / `fgets`), decoding JSON (`jsondecode` / a small JSON library), evaluating the model's RHS for every row of `y`, and printing one encoded line with a flush. Only `info`, `rhs` and `quit` are required.

## The trace contract

What every scheme in every language writes for validation, and what `bimgrl reference` writes for the reference.

**`<name>.csv`**: header `t,<state names in model order>`, one row per sample, full precision (`%.17g`). Optional extra columns `aux_<intermediate>` for currents or fluxes.

**`<name>.json`** (the run meta; fields of `RunMeta` in `src/bimgrl/trace.py`):

```json
{
  "model": "tworld", "scheme": "BIMGRL", "dt": 0.1, "reference": false,
  "solver": {}, "language": "julia", "t_unit": "ms",
  "protocol": "B", "beats": 10, "bcl": 1000.0,
  "stimulus": {"start": 0.0, "duration": 1.0, "period": 1000.0, "amplitude": -53.0},
  "status": "ok", "message": "", "n_steps": 100000,
  "wall_s_total": 3.3, "wall_s_per_beat": 0.33,
  "parameters": {}, "float_type": "float64", "blocks": "blocks.toml",
  "created": "2026-10-27T14:00:00"
}
```

A reference trace has `"reference": true` and `"solver": {"name": "LSODA", "rtol": 1e-10, "atol": 1e-12}`; `bimgrl validate` refuses a reference whose `rtol` is above 1e-8. `status` is one of `ok`, `nan`, `unphysical`, `timeout`, `error`.

**`<name>_beats.csv`** (optional, protocol B): header `t,<states>`, the state at the start of every beat, for the drift of the slow states and of the invariants.

**`bimgrl validate <model_dir> --test trace.csv [--reference ref.csv] [--level moderate|tight] [--name NAME]`** reads both, takes the overlap of the two time ranges, interpolates the reference at the test's times, classifies the test's physical range against the reference's own range, computes the biomarkers of `thresholds.md`, the pass/fail per set, the Na⁺ and K⁺ drift and every invariant's drift from the `_beats.csv` files, and writes `bimgrl/validation_<name>.md` and `.json`. Protocol B traces should carry the last beat on a 0.1 ms grid (or the scheme's own step if coarser), as `bimgrl reference` does; a trace written with `%.6g` cannot pass the invariant floors, so write full precision.
