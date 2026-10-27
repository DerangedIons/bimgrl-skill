# Translating a model to the canonical `.ode`

The deterministic pipeline keys off gotranx's `.ode` form: it gives the symbolic Jacobian, the stoichiometric atoms, the generated schemes and the compiled sweep. A model written in Julia, Python, C or MATLAB is translated into it by the agent and the translation is then proven, not trusted, by `bimgrl fidelity`. A model in CellML needs no translation: `bimgrl convert` does it with myokit.

## When to translate and when to go native

Translate when the goal is a generated BIMGRL scheme for the model (C, Python, Julia or CytoZoo) and the model is a plain right-hand side. Stay native (the adapter of `contracts.md` plus the recipe of `implementing-blocks.md`) when the integrator itself is the product, hand-optimised and must keep its own code, as for TWorld.jl, or when the model uses constructs the `.ode` grammar has no spelling for (data-dependent loops, lookup tables, external inputs). The analyses run on both; codegen and the automated sweep run only on the canonical path.

## `.ode` essentials

```
states("membrane",
V=ScalarParam(-85.23, unit="mV", description="membrane potential"))

states("calcium",
Ca_i=ScalarParam(1.2e-4, unit="mM"),
Ca_SR=ScalarParam(3.64, unit="mM"))

parameters("membrane",
Cm=ScalarParam(0.185, unit="uF"),
F=ScalarParam(96485.3415, unit="C/mmol"),
stim_start=ScalarParam(10.0, unit="ms"),
stim_duration=ScalarParam(1.0, unit="ms"),
stim_period=ScalarParam(1000.0, unit="ms"),
stim_amplitude=ScalarParam(-52.0, unit="A/F"),
stim_eps=ScalarParam(1e-9, unit="ms"))

expressions("membrane")
stim_phase = t - stim_start
stim_phase_mod = stim_phase - stim_period * floor((stim_phase + stim_eps) / stim_period)
i_Stim = Conditional(And(Ge(stim_eps, -stim_phase), Gt(stim_duration, stim_eps + stim_phase_mod)), stim_amplitude, 0.0)
E_K = (R*T/F) * log(K_o/K_i)
i_K1 = g_K1 * (V - E_K) # A/F
dV_dt = -(i_K1 + i_Stim + ...) # mV
```

States and parameters are grouped by component; every state carries a unit (concentrations in a unit containing `M`, gates `1`) and an initial value. Expressions follow, one per line, intermediates before the `dX_dt` derivative of every state. Conditionals are `Conditional(cond, a, b)` with `Ge`, `Gt`, `Le`, `Lt`, `And`, `Or`; `floor`, `exp`, `log`, `sqrt`, `pow` are available. Keep the original names for states, parameters and currents so that the stoichiometric report and the user's code speak the same language, and keep every current as its own intermediate (not expanded into the derivative) so the bookkeeping can see it. Time is `t`, in ms; voltage in mV; rewrite the stimulus to the half-open form above.

## `translation_map.json`

The map from the original code to the `.ode`, written by the agent alongside `model.ode`:

```json
{
  "states": {"u[1]": "V", "u[2]": "Ca_i"},
  "parameters": {"p.gK1": "g_K1"},
  "scale": {"V": {"factor": 1000.0, "from_unit": "V", "to_unit": "mV"}, "t": {"factor": 1000.0, "from_unit": "s", "to_unit": "ms"}},
  "time_unit": "s",
  "stimulus": {"amplitude": "p.Istim", "start": "p.tstim", "duration": "p.dur", "period": "p.bcl"}
}
```

## What the original code must produce

Two files from the user's own implementation, untouched by the translation:

`rhs_samples.csv`: header `t, y_<state>..., dy_<state>...` in the original names and units, at least 50 rows spread over one beat (upstroke, plateau, repolarisation, diastole), each row a state and the derivative the original code returns at it. The rows of the user's own trace, with the RHS evaluated at each, are the natural source.

`params.json`: every parameter the original code uses, by original name, with its value.

Optionally a trace of the original code in the trace contract (`contracts.md`), with its run meta saying which solver produced it.

## `bimgrl fidelity`

```
bimgrl fidelity <model_dir> --samples rhs_samples.csv --params params.json [--trace trace.csv] [--map translation_map.json] [--rhs-tol 1e-9]
```

Three checks, in order, written to `bimgrl/fidelity.md`:

1. **Parameters.** Every `.ode` default equals the original value after the map's scaling, to 1e-12 relative. Differences are listed by name.
2. **Pointwise right-hand side.** The `.ode` is compiled to numpy and evaluated at every row of the samples after mapping and scaling; the requirement is $|dy_{ode} - dy_{orig}| \le \text{rhs-tol} \cdot \max(|dy_{orig}|, s_i)$ per state, with $s_i$ the kind scale (1 mV, 1e-3 for a gate, 1e-4 mM). This is the proof. It is independent of any integrator, and a failure localises to a state, and through the atom table to an intermediate.
3. **Trajectory.** The `.ode` is integrated with LSODA from the first row of the user's trace with the declared stimulus and compared with the trace: tight thresholds if the trace came from an adaptive solver with `rtol ≤ 1e-8`, moderate thresholds plus a mandatory pass of check 2 otherwise, because a trace from a fixed-step scheme carries its own error and cannot prove fidelity by itself.

A pass of check 2 with a fail of check 3 is reported as a protocol or initial-state mismatch, not as an equation error. A Float32 model cannot pass the default `--rhs-tol`; relax it to 1e-5 and say so in NOTES.md.
