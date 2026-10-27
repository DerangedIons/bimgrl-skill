# bimgrl

A Claude Code and Codex skill that applies the block-implicit multirate generalised Rush–Larsen scheme (BIMGRL) to a cardiac cell model, with deterministic scripts for everything that is arithmetic and the agent reserved for the judgment calls. It is the workflow of the paper *Targeted implicit integration accelerates tissue simulations of advanced cardiac cell models*, distributed as a Claude and Codex skill.

The problem it solves: generalised Rush–Larsen advances each state alone, exponentially where the state has a linear self-term and by forward Euler otherwise, which is fine as long as every fast time scale belongs to one state. In detailed models the fastest modes belong to groups: a Ca pool with its buffers (binding at hundreds per millisecond), the occupancies of a Markov chain. Forward Euler on those states caps the step at $2/|\lambda|$, a few microseconds, and even the exponential update, applied one state at a time, breaks the group's conservation law, so chain sums and Ca mass drift. An *implicit block* solves such a group together with one cheap linearly implicit step in increment form, which removes the stability limit and keeps the invariant exact. What remains at large steps is accuracy loss in two briefly fast processes, the upstroke and the Ca release; a *sub-step block* refines only those states, only while they are active, decided per cell from its own state. The skill finds the groups, proves what they conserve, proposes the blocks by the paper's rules, generates the schemes with gotranx, and sweeps the time step against a tight reference.

## Install

```sh
git clone <this repository> ~/dev/bimgrl-skill
cd ~/dev/bimgrl-skill
uv venv && uv pip install -e .          # numpy, scipy, sympy, myokit, lxml
uv run bimgrl setup                      # gotranx from config.toml (or GOTRANX_SRC=/path/to/clone)
ln -s ~/dev/bimgrl-skill ~/dotfiles/llm/skills/bimgrl   # the folder both Claude Code and Codex read
```

`bimgrl setup` installs gotranx from the spec in `config.toml`. The `block_rush_larsen` scheme, `suggest-blocks` and the symbolic Jacobian live in the DerangedIons fork until they are upstreamed, so this needs access to that repository or a local clone at the pinned commit. A C compiler (`gcc`, or `CC`) is needed for the sweep. The Julia helpers under `scripts/jl/` need Julia 1.10 or later and instantiate their own project on first use.

## Commands, in pipeline order

All of the form `uv run --directory <skill-dir> bimgrl <command> <model_dir> [flags]`; outputs land in `<model_dir>/bimgrl/`.

| command | does |
|---|---|
| `setup` | install gotranx into the skill environment |
| `convert` | CellML to `model.ode` (ms, mV, half-open stimulus) with a draft `meta.toml` |
| `intake` | draft `meta.toml` from `model.ode`; report what the user must confirm (compartments, gates, stimulus) |
| `fidelity` | prove a hand-translated `model.ode` against the original code (parameters, pointwise RHS, trajectory) |
| `structure` | dependency graph, strongly connected components, cone of influence of V, dead and observer states, linear subsystems, how each scheme advances each state |
| `invariants` | stoichiometric conservation analysis: chains, binding clusters, total-ion invariants with volume-ratio weights, "written improperly" flags, charge bookkeeping of every current |
| `reference` | LSODA reference traces at rtol 1e-10, cut at the stimulus edges, with self-convergence |
| `stiffness` | eigenmodes, participation factors, owners and their update type, forward-Euler and GRL limits (predicted) |
| `sweep` | run schemes over the dt grid against the reference, coarsest first, with status, biomarker errors, invariant drift, and the stability-limit bisection |
| `propose` | `blocks.toml` by the paper's rules, each block with the finding behind it, compared with gotranx `suggest-blocks` |
| `codegen` | the schemes for your target: C, Python, Julia, CytoZoo |
| `validate` | the same checks on a trace from any language, against the reference trace |
| `log` | append a tuning round to `NOTES.md`; freeze the round into `blocks.toml` |
| `time` | time the passing runs in C, best of five |
| `report` | cheapest passing run per threshold set, speed-ups, measured against predicted limits, the Table-S1 row, an error-versus-cost figure |

`references/workflow.md` walks the phases with their agent gates; `references/block-rules.md` has the rules; `references/contracts.md` the file formats and the adapter protocol for models in other languages.

## Layout

```
SKILL.md                 the skill: phases, gates, commands
agents/openai.yaml       Codex metadata
config.toml              where gotranx is installed from
pyproject.toml           the Python package (src/bimgrl) and its dependencies
src/bimgrl/              the commands; vendored/ holds the survey's biomarker, metric and stiffness code with provenance
scripts/jl/              Julia helpers: adapter (ForwardDiff Jacobian over JSON lines), Rodas5P reference, sparsity
references/              workflow, block rules, implementing blocks by hand, contracts, invariants, translation, thresholds, notes template
assets/                  templates for meta.toml and blocks.toml
evals/                   prompts and fixtures for testing the skill
```

## What the scripts decide and what the agent decides

The scripts decide every number: which states own which mode, what the model conserves and with which weights, whether a chain is written consistently, which step is stable, which biomarker fails at which step, what a run costs. The agent decides what no heuristic gets right every time, and only after showing the user the evidence: which Ca pool is the junction and which states form the release channel, whether a conservation flag is an error of the model or a documented property, whether the SR belongs in the dyad block, which tolerance to change in a tuning round, and how to read a failure (stability or accuracy). Every change to the blocks is logged with its numbers before the next run, and the blocks are frozen before anything is timed.
