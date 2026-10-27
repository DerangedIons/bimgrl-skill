---
name: bimgrl
description: Accelerate the fixed-step integration of a cardiac cell model with the block-implicit multirate generalized Rush-Larsen scheme (BIMGRL) from the paper "Targeted implicit integration accelerates tissue simulations of advanced cardiac cell models". Deterministic scripts do the analysis (Jacobian eigenmodes with participation factors, conservation-law and charge-bookkeeping checks, dependency structure, block proposal by the paper's rules, gotranx code generation, dt sweeps against an LSODA reference with biomarker thresholds); the agent confirms compartments, reviews blocks and tunes tolerances. Use when asked to speed up, stabilise, or find the stiffness or time-step limit of a cardiac (or any stiff biochemical) cell model in CellML, gotranx .ode, Julia, Python, C or MATLAB, to apply BIMGRL or GRL blocks, to check whether a model's Markov chains and buffers conserve what they should, or to reproduce the paper's survey for a new model.
license: MIT
allowed-tools: Bash, Read, Write, Edit, Glob, Grep, AskUserQuestion
metadata:
  version: "0.1.0"
  paper: "Targeted implicit integration accelerates tissue simulations of advanced cardiac cell models"
---

# bimgrl: apply BIMGRL to a cardiac cell model

## What this is

Generalized Rush-Larsen (GRL) advances every state alone. That is fine while each fast time scale belongs to one state, as it does for a gating variable, and it breaks when a fast mode belongs to a *group*: a Ca pool with its buffers, or the states of a Markov chain. The group's mode has no diagonal entry to be exponential in, so GRL either steps it with forward Euler (and the buffer mode of T-World limits the step to 0.0077 ms) or damps it poorly and leaks the group's invariant (a Markov chain whose probabilities sum to 1.25 by beat 10, a 0.19 % Ca leak per release). BIMGRL gives each such group one linearly implicit step in increment form, which keeps the invariant to round-off, and sub-steps the two briefly fast subsystems (V with the INa gates, the dyad during release) only while they are active. T-World then runs at 0.1 ms instead of 0.005 ms, 19 times faster, with the same biomarkers, alternans, EADs and DADs.

This skill turns the paper's workflow into scripts. The scripts do every measurement; you do the judgment at the gates marked below, and you never type a number the scripts can produce.

## Setup (once per machine)

```
cd <skill-dir>
uv venv .venv && uv pip install --python .venv/bin/python -e .
uv run --directory <skill-dir> bimgrl setup      # gotranx from config.toml (or GOTRANX_SRC=/path/to/clone)
```

`<skill-dir>` is the directory of this file. Every command below is `uv run --directory <skill-dir> bimgrl <command> <model_dir> [flags]`; `<model_dir>` holds `model.ode`, `meta.toml`, `blocks.toml`, `NOTES.md`, and the outputs under `<model_dir>/bimgrl/`. A C compiler is needed for the sweep (`CC`, else gcc); Julia 1.10+ only for the Julia on-ramp (`scripts/jl`).

## The two paths

- **Canonical path**: the model is, or becomes, a gotranx `.ode` file. Everything below runs. CellML converts deterministically (`convert`); a Julia, Python, C or MATLAB model is translated by you and *proven* with `fidelity` (pointwise right-hand-side agreement to 1e-9) before any analysis is trusted. Read `references/translation.md` before translating.
- **Native path**: the integrator must stay in its own code (a hand-optimised Julia model). Run `stiffness --adapter` and `invariants --adapter` through an RHS adapter (JSON lines; `scripts/jl/adapter.jl` serves any `f(du,u,p,t)` with ForwardDiff), implement the blocks yourself with `references/implementing-blocks.md`, and validate every scheme through the trace contract with `validate`. Read `references/contracts.md` first.

Prefer the canonical path whenever translation is feasible: it is the only one where the blocks are generated rather than hand-written.

## Phases and gates

Run the phases in order; each prints a Markdown report and leaves a JSON for the next one. Stop at every **Gate** and put the question to the user with their options; do not continue on an assumption.

### 0. Intake

```
bimgrl intake <model_dir> [--voltage V]
```

Drafts `meta.toml` from the equations: membrane potential, cytosolic Ca/Na/K, the Ca compartments (SR, junction/dyad, sub-sarcolemma), the release-channel states and their open states, the INa gates, the stimulus. Every guess carries its evidence.

**Gate 0** (always): show the table of guesses and have the user confirm the rows marked *confirm*. The compartment naming is the one input no heuristic gets right every time (two SR pools, a release flux written as a state, Ca-dependent ICaL gates that look like release states). Set `compartments.status = "confirmed"` in `meta.toml` when done. If the stimulus is not in the half-open `stim_*` form, stop: the model needs `convert` or the rewrite in `references/contracts.md` first.

### 1. Structure

```
bimgrl structure <model_dir>
```

Dependency graph from the symbolic Jacobian sparsity: strongly connected components, what feeds V, dead and observer states, duplicate dependency rows (phosphorylated copies), structurally linear subsystems (chains, buffer pairs), and how GRL (exponential wherever the diagonal is non-zero) and RL+FE (gates only) advance each state. No gate; read `structure.md` once so the later reports make sense.

### 2. Invariants and conservation checks

```
bimgrl invariants <model_dir>
```

From the stoichiometry of the equations (not the Jacobian): every Markov chain and binding cluster with its weights and class (`exact`, `conserved-internally`, `closure`, `improper`), the total-ion invariants with their volume-ratio weights (`exact`, `boundary`), the algebraic closures, and the charge bookkeeping of every current in dV/dt. Read `references/invariants.md` to interpret it.

**Gate 2** (when the Flags section is not empty): an `improper` chain or cluster, a weight-not-1 warning inside a compartment, or an `unaccounted` / `miscounted` current means the model is written inconsistently (a rate missing from one equation, a volume ratio in one but not the other, a stimulus that drifts K⁺). Show the flagged atoms and columns to the user and ask whether to fix the model, proceed as is, or stop. Never block a chain that does not conserve: the implicit step would conserve the wrong thing exactly.

### 3. Reference

```
bimgrl reference <model_dir> [--beats 10] [--protocol A,B]
```

LSODA at rtol 1e-10 / atol 1e-12, cut at the stimulus edges, with a self-convergence check; written in the trace contract. If `self-convergence max |dV|` is above 1e-5 mV, say so: every error below is measured against this.

### 4. Stiffness

```
bimgrl stiffness <model_dir> [--every 1] [--threshold -20]
```

Eigenmodes and participation factors at the paper's 13 checkpoints plus every 1 ms through a beat; the Table-1 view by owner set with the update type of each owner under GRL and under RL+FE; the predicted limits P1 (forward-Euler limit of the stiffest explicitly stepped mode) and P2 (linearised step). Report the three or four stiffest owner sets in words: which modes bind, and whether their owners are stepped explicitly. The prediction is checked against the measured limit in phase 6.

### 5. Baseline

```
bimgrl sweep <model_dir> --scheme GRL --stability
```

GRL from the coarsest allowed step down, with the stability limit bisected, biomarkers against the reference, the physical-range check, and the drift of every invariant per beat. Three things to read off: the measured limit against the phase-4 prediction (agreement within a factor of ~1.5 confirms the mode diagnosis); a run that is `ok` but carries a `DRIFT` flag (a chain or a Ca total leaking), which is the "runs but gives the wrong answer" case that needs an implicit block even if nothing crashes; and the step at which GRL meets the moderate thresholds, the number BIMGRL has to beat.

### 6. Propose blocks

```
bimgrl propose <model_dir> [--dt 0.1] [--nest-v-markov]
```

Applies the paper's rules (`references/block-rules.md`) to the findings: chains of channels and of the release channel, binding clusters that own a stiff mode, SR joined with the junction so the release conserves, V with the INa gates (proportional, atol 0.5 mV), the dyad (binary, triggers on the junctional Ca and the open states), nesting. Every block carries the finding that produced it; gotranx's own `suggest-blocks` is run for comparison; `validate_blocks` checks the file.

**Gate 6** (always): walk the user through the proposed `blocks.toml` block by block, with the reason for each, then the "Judgment calls" list (keep the SR in the dyad? add a Ca-dependent ICaL gate such as TT06's fCass? nest a voltage-dependent chain in the V block?). Record the decisions. This is the paper's "refinements restricted to membership and tolerances": nothing else about the scheme is ever changed.

### 7. Sweep and validate

```
bimgrl sweep <model_dir> --scheme GRL,BIGRL,BIMGRL [--dt 0.5,0.2,0.1,0.05,0.02,0.01,0.005] [--save-traces]
```

Every scheme at every allowed step (1-2-5 per decade, dividing the stimulus edges and the 0.1 ms grid), coarsest first, stopping after two consecutive tight passes; status `ok | nan | unphysical | budget`, errors in Vpeak, APD90, CaT amplitude, activation time, pass/fail for the moderate and tight sets, invariant drift, and the first state to leave its range when a run fails. BIGRL (implicit blocks only) isolates what the implicit step bought; BIMGRL adds the sub-steps. For a scheme implemented elsewhere (native path), run its trace through `bimgrl validate <model_dir> --test trace.csv` instead; the metrics are the same code.

### 8. Tune (at most four rounds)

Read the failures: a Vpeak overshoot at large steps is the under-resolved upstroke (V-block trigger atol, 0.5 → 0.1); a CaT error that grows with the step and a release leak that cancels it at one step (the T-World round-2 story) means the SR belongs in the dyad; a chain state going negative at the largest step means the chain nests in the V block; a `DRIFT` flag on a scheme with blocks means a block boundary cuts an exchange. Propose one change, have the user confirm, edit `blocks.toml`, and record it:

```
bimgrl log <model_dir> --round N --change "..." --result "..." --decision "adopted|rejected: why" [--freeze]
bimgrl sweep <model_dir> --scheme BIMGRL --force
```

**Gate 8** (every round): the change and its reason before the edit; the numbers before the decision. Stop after four rounds or when the moderate thresholds hold at a step at least ten times GRL's; `--freeze` stamps the round into `blocks.toml`. Blocks are frozen before any timing.

### 9. Time and report

```
bimgrl sweep <model_dir> --scheme GRL,BIMGRL --beats 100 --dt <chosen steps>
bimgrl time <model_dir> --scheme GRL,BIGRL,BIMGRL
bimgrl report <model_dir>
bimgrl codegen <model_dir> --backend c|python|julia|cytozoo --scheme BIMGRL
```

The 100-beat run confirms the configuration chosen at beat 10 still holds (the survey lost 13 of 39 models here, mostly by 1 to 2 % in CaT). `time` measures the cost per beat in C, best of five, on one core; `report` picks the cheapest passing run of each scheme (not the largest step), the speed-ups at both threshold sets, the measured against the predicted stability limits, the Table-S1 row and an error-vs-cost figure. `codegen` emits the chosen scheme for the user's target. Hand back `report.md`, the frozen `blocks.toml`, `NOTES.md`, and a paragraph in the paper's style (which modes limited the step, which blocks removed them, the speed-up at matched accuracy).

## What you decide and what you do not

You decide: compartment naming, whether a flagged conservation law is a model bug or a convention, block membership and trigger tolerances within the rules, when a round is "good enough", and how to describe the result. You do not decide: any biomarker, error, limit, drift, cost or pass/fail (the scripts print them; quote them), the thresholds (`references/thresholds.md`), the form of the implicit step, or the reference solver and its tolerances. If a script fails on the model (a sympy timeout, a compile error), report the message and the file it wrote, and fix the cause rather than working around the script.

## Reference files

| File | Load when |
|---|---|
| `references/workflow.md` | the full procedure with every flag and output, and the native-path variant |
| `references/block-rules.md` | at Gate 6 and in tuning: the rules, the blocks.toml schema, what each knob trades off |
| `references/invariants.md` | at Gate 2: how to read the classes, weights, flags and the drift monitor |
| `references/thresholds.md` | when a run passes or fails and the user asks why |
| `references/contracts.md` | before a native-path run or a trace from another language: adapter protocol, trace format, meta.toml |
| `references/translation.md` | before translating a non-.ode model, and to run `fidelity` |
| `references/implementing-blocks.md` | native path: the increment-form step, the joint SR-junction solve, the two wrong fixes |
| `references/notes-template.md` | writing NOTES.md and the tuning log |
