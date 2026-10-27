# The workflow, phase by phase

Everything here is driven from one command, run from the skill directory against a model directory that holds `model.ode` and `meta.toml`:

```
uv run --directory <skill-dir> bimgrl <command> <model_dir> [flags]
```

Every phase writes a JSON (for the next script) and a Markdown report (for the agent and the user) under `<model_dir>/bimgrl/`. Scripts never raise on a model failure: a run that diverges, goes unphysical or exceeds its budget is a row with a status, so the agent can reason about it. The agent's job is the judgment at the gates below, never the arithmetic.

## The mental model (from the T-World implicit tutorial)

1. Note which variables misbehave and how. Sign-flipping growth means an explicit instability; a smooth drift away from the reference means a conservation leak.
2. Compute the Jacobian at snapshots through a beat, diagonalise, list the fastest eigenvalues.
3. Use participation factors, not raw eigenvector components, to say which states own a mode. They are unit-free, so a concentration and a gate can be ranked against each other.
4. Cross the mode list with how each state is integrated. The binding limit is the fastest mode whose owners are stepped explicitly (forward Euler).
5. If the fast modes live in small blocks, treat only those blocks implicitly with an analytic sub-Jacobian in increment form, and check the conservation of the block's invariant.
6. Validate against an adaptive reference, first at the old step, then at the new one, then on the threshold phenomena (alternans, EADs, DADs).

The pipeline below is this list made deterministic, with the invariant check promoted to its own phase because it catches the failures that "runs and looks like an AP" misses.

## Two paths

**Canonical path.** The model is a gotranx `.ode` file (converted from CellML by `bimgrl convert`, or translated by the agent and proven by `bimgrl fidelity`, see `translation.md`). Every phase runs, gotranx generates the schemes, and the sweep measures cost in compiled C.

**Native path.** The model must stay in its own code (a hand-optimised Julia integrator, a legacy C model). The analysis phases run through an RHS adapter (`contracts.md`), the agent implements the blocks in the user's integrator by the recipe in `implementing-blocks.md`, and validation goes through the trace contract with `bimgrl validate`. There is no `propose`, `codegen` or `sweep` on this path; the block rules in `block-rules.md` are applied by hand.

## Phase 0: setup

```
uv run --directory <skill-dir> bimgrl setup [--skill-dir DIR]
```

Installs gotranx from the spec in `config.toml` (or `GOTRANX_SRC=/path/to/clone`) into the skill's environment. The `block_rush_larsen` scheme, `suggest-blocks` and the symbolic Jacobian live in the DerangedIons fork until they are upstreamed, so this needs access to that repository or a local clone at the pinned commit. Nothing else in the skill works without it on the canonical path; the native path needs only numpy and scipy.

## Phase 1: intake

```
bimgrl intake <model_dir> [--voltage NAME] [--force]
```

Reads `model.ode`, writes a draft `meta.toml` and `bimgrl/intake.md`. The draft names the membrane potential (the only mV state, else the state most derivatives depend on), the cytosolic Ca, Na and K (by normalised name), the stimulus (the half-open `stim_*` parameters if present, else a warning that the stimulus must be rewritten), and the `[compartments]` table: the SR (the free Ca pool with the largest initial value), the junction (the pool exchanging directly with the SR), the sub-sarcolemmal pool, the release channel (dimensionless states coupled both ways to the junction and depending on the SR), its open states, and the INa gates (HH gates named m*, h*, j* or containing "na").

**Gate.** The agent shows the user every row marked *confirm* in `intake.md` and edits `meta.toml` accordingly. Compartment naming is the one thing no heuristic gets right every time, and every block rule downstream depends on it. The agent decides alone: nothing here.

## Phase 2: structure

```
bimgrl structure <model_dir>
```

From the symbolic Jacobian sparsity: the dependency graph, strongly connected components, the cone of influence of V, dead states (f ≡ 0), observer states (nothing depends on them, typically contraction), identical dependency rows (symmetry candidates such as phosphorylated copies), the structurally linear subsystems (Markov chains, linked buffers), the HH gates, and for every state how GRL (exponential wherever the Jacobian diagonal is non-zero) and RL+FE (exponential for HH gates only) advance it. Writes `structure.json`, `structure.md`, `structure.dot`.

**Gate.** None. The agent reads `structure.md` to know the model's shape before interpreting the stiffness table.

## Phase 3: invariants

```
bimgrl invariants <model_dir> [--no-numeric]
bimgrl invariants <model_dir> --adapter cmd:"..." --trace trace.csv [--group a,b,c]   # native path
```

The stoichiometric analysis of `invariants.md`: Markov chains and binding clusters with their conservation class, the total-ion invariants with their volume-ratio weights, the closures, and the charge bookkeeping of every current in dV/dt. Writes `invariants.json`, `invariants.md`.

**Gate.** Hard stop on any flag: a chain or cluster of class `improper`, a binding cluster with weights other than 1, or a current that is `unaccounted` or `miscounted`. The agent shows the user the atom and the columns the report names and asks whether the model should be corrected before going on. A model that does not conserve cannot be rescued by a scheme that does, and an unaccounted stimulus drifts K⁺ over hundreds of beats, which the 100-beat check would otherwise blame on the scheme. The agent decides alone: which flags are known properties of the published model (documented in NOTES.md) rather than errors.

## Phase 4: reference

```
bimgrl reference <model_dir> [--protocol A,B] [--beats 10] [--rtol 1e-10] [--atol 1e-12]
```

LSODA at tight tolerances, cut at every stimulus edge so the solver never sees the discontinuity, with a self-convergence check at ten times tighter tolerances. Protocol A is one window from the model's initial state; protocol B is N beats at the basic cycle length with the last beat compared and the state at every beat start saved for drift. Writes `reference_A.csv/json`, `reference_B.csv/json`, `reference_B_beats.csv`, `reference.md`.

**Gate.** The agent checks the self-convergence line (max |ΔV| should be well below 1e-5 mV) and the reference biomarkers against the model's published values when it knows them. A reference that disagrees with the paper the model came from usually means a protocol or unit problem, not a numerics problem.

## Phase 5: stiffness

```
bimgrl stiffness <model_dir> [--protocol B] [--beats 2] [--every 1.0] [--threshold -20]
bimgrl stiffness <model_dir> --adapter cmd:"..." --trace trace.csv   # native path
```

The Jacobian at the paper's 13 checkpoints plus every millisecond through the last beat; eigenvalues; participation factors $p_{kj} = |V_{kj} W_{jk}|$ normalised per mode; the owners of every mode with Re λ ≤ threshold; the owners' update type under GRL and under RL+FE; the forward Euler limit $2|\mathrm{Re}\,\lambda| / |\lambda|^2$; the P1 prediction (the stiffest mode owned by explicitly advanced states) and the P2 prediction (the largest step at which the linearised scheme step grows no faster than the model, which also sees coupled modes that GRL updates one state at a time). Writes `stiffness.json`, `stiffness.md`.

**Gate.** The agent reads the Table-1 view and names, in words, what limits the step and why: which modes are buffer binding, which are gates, which are chains. It decides alone whether the threshold should be relaxed to show more modes. It does not yet decide blocks.

## Phase 6: baseline sweep with the drift monitor

```
bimgrl sweep <model_dir> --scheme GRL --stability [--dt 0.1,0.05,0.02,0.01,0.005]
```

GRL at every allowed step from the coarsest down, each run classified ok / nan / unphysical / budget, compared with the reference on the last beat, and every invariant of phase 3 evaluated at the start of every beat. The `--stability` bisection finds the largest step that runs. Writes rows to `runs.csv` (and `runs_drift.csv` for the per-invariant drift) and `stability.json`.

**Gate.** The agent compares the measured limit with the P1 and P2 predictions of phase 5 and reports whether they agree (the paper's 0.0077 ms predicted against 0.0079 ms failing). A `DRIFT` flag on a run that otherwise passes the thresholds is the Gauthier case, a scheme that runs but gives the wrong answer; the agent reports it as a finding of its own.

## Phase 7: propose blocks

```
bimgrl propose <model_dir> [--dt 0.1] [--threshold -20] [--nest-v-markov] [--no-write]
```

The rule engine of `block-rules.md` over the outputs of phases 2, 3 and 5 and the compartments of `meta.toml`: implicit blocks for every channel or release Markov chain and every binding cluster that owns a stiff mode, the SR and the junction joined, the V/INa block, the dyad block, straddling blocks nested. Every block carries the finding that produced it. gotranx's own `suggest-blocks` runs at the target step for comparison and the differences are listed. Writes `<model_dir>/blocks.toml` (unless it exists), `bimgrl/blocks_proposed.toml`, `propose.md`, `propose.json`.

**Gate.** The agent walks the user through `propose.md`: each block and its reason, the "left explicit" list, and the "judgment calls" list (the SR in the dyad, Ca-dependent gates such as fCass, a voltage-dependent chain outside the V block). Changes are made in `blocks.toml` and noted in NOTES.md with `bimgrl log --section Blocks --line "..."`.

## Phase 8: sweep and validate

```
bimgrl sweep <model_dir> --scheme GRL,BIGRL,BIMGRL [--protocol B] [--beats 10] [--dt ...] [--budget 600] [--blocks FILE] [--suggested FILE] [--stability] [--save-traces] [--time] [--force]
bimgrl codegen <model_dir> [--backend c|python|julia|cytozoo] [--scheme GRL|BIGRL|BIMGRL] [--blocks FILE] [-o OUT]
bimgrl validate <model_dir> --test trace.csv [--reference ref.csv] [--level moderate|tight] [--name NAME]
```

The sweep generates and compiles each scheme once per configuration, runs the dt grid coarsest first and stops after two consecutive steps pass the tight thresholds (or at the budget), and records biomarker errors, pass/fail per threshold set, the physical-range check with the first violating states, the invariant drift, and a cost estimate. `codegen` writes the same schemes for the user's own target. `validate` applies the same checks to a trace produced by any code, against `reference_<protocol>.csv`.

**Gate.** The agent reads the rows as the tuning log of `thresholds.md` reads them: which step passes which set, where it fails and on which biomarker, whether a failure is a stability loss (nan, unphysical) or an accuracy loss. It decides alone which steps to add with `--dt`.

## Phase 9: tune, at most four rounds

```
bimgrl log <model_dir> --round N --change "..." --result "..." --decision "..." [--blocks FILE] [--freeze]
```

Only two things may change between rounds: block membership and trigger tolerances. Each round is one entry in NOTES.md with what changed, the numbers, and the decision. The blocks are frozen (`--freeze` stamps the round into the `blocks.toml` header) before anything is timed.

**Gate.** Every change is proposed by the agent with its reason and the number it expects to move, and confirmed by the user before the rerun. The agent never changes a block without logging it.

## Phase 10: time and report

```
bimgrl time <model_dir> [--scheme GRL,BIMGRL] [--protocol B] [--beats 1] [--repeats 5] [--blocks FILE]
bimgrl sweep <model_dir> --beats 100 --scheme GRL,BIMGRL --dt <chosen steps>
bimgrl report <model_dir> [--protocol B] [--baseline GRL]
```

Timing is serial, in C, best of N, one beat; it must run on an idle machine. The 100-beat sweep at the chosen steps confirms that the configuration picked at beat 10 still holds. The report picks, per scheme and threshold set, the cheapest passing run (not the largest step: heavy sub-stepping can make a larger step slower), the speed-ups at matched accuracy, the measured stability limits against the predictions, the Table-S1 row, and an error-versus-cost figure. Writes `report.md`, `report.json`, `error_vs_cost.png`.

**Gate.** The agent hands the user `report.md`, the frozen `blocks.toml`, the completed `NOTES.md` and the generated code, and states plainly which thresholds were met at which step and what the speed-up is. On the native path the deliverable is instead the recipe checklist of `implementing-blocks.md` with the validation results.

## Outputs at a glance

| phase | files under `bimgrl/` |
|---|---|
| intake | `intake.md`, `intake.json`, and `../meta.toml` |
| structure | `structure.json`, `structure.md`, `structure.dot` |
| invariants | `invariants.json`, `invariants.md` (native: `invariants_native.*`) |
| reference | `reference_A.*`, `reference_B.*`, `reference_B_beats.csv`, `reference.md` |
| stiffness | `stiffness.json`, `stiffness.md` (native: `stiffness_native.*`) |
| sweep | `runs.csv`, `runs_drift.csv`, `stability.json`, `build/`, `traces/` |
| propose | `blocks_proposed.toml`, `blocks_suggested_dt*.toml`, `propose.md`, `propose.json`, and `../blocks.toml` |
| codegen | `generated/` |
| validate | `validation_<name>.md`, `validation_<name>.json` |
| log | `../NOTES.md` |
| report | `report.md`, `report.json`, `error_vs_cost.png` |
