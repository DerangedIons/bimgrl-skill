# Choice of blocks

Two kinds of block exist because detailed cell models break the assumption behind GRL in two different ways. GRL advances each state alone, exponentially where it has a linear self-term and by forward Euler otherwise, so it is fine as long as every fast time scale belongs to a single state. A Ca pool with its buffers, or a Markov chain, has fast modes that belong to the group and to no single state; an implicit block solves the group together. The upstroke and the Ca release are fast only briefly; a sub-step block refines only those states, only while they are active.

## The rules

1. **Each Ca compartment with fast buffers is an implicit block with its buffer states.** Mass-action binding $k_{on} c (B_{max} - B) - k_{off} B$ is a stiff exchange between the free and the bound form; an increment-form implicit step keeps $c + B$ exact at any step. Buffers whose binding is not stiff (cytosolic troponin, calmodulin in T-World, every buffer that owns no mode below the threshold) stay explicit, because the block would cost and gain nothing.
2. **Compartments coupled by the RyR release, the SR and the junctional or subspace pool, share one implicit block.** Solved separately, the SR loses $h(J_{rel} + g\,\Delta c_{SR})$ while the junction gains $h(J_{rel} - g\,\Delta c_j)$; the mismatch is small per step but accumulates over release-heavy protocols (0.19 % of the release in T-World, which overloaded the SR). Together they see the same linearised release. If the junction has no buffer states it joins the SR block as a bare pool.
3. **Every Markov chain of a voltage- or ligand-gated channel and of the RyR is an implicit block.** The exponential update treats each occupancy alone, so the probability of leaving one state does not equal the probability of arriving in the next, and the sum drifts (T-World's IKr at 1.0454 after prepacing with GRL; the Gauthier CaRU at 1.25 by beat 10 and 41 by beat 600). One implicit step on the chain keeps the sum to round-off.
4. **The membrane potential and the gates of the fast Na⁺ current, with the late Na⁺ current where present, form a proportional sub-step block.** The number of sub-steps is set so that V changes by about the tolerance per sub-step: a rate trigger on V with atol 0.5 mV, tightened to as low as 0.1 mV where the activation time error demands it, and `h_min` 0.001 ms. In practice the block engages for less than 2 ms per beat.
5. **The dyad, with the junctional Ca, its buffers, the release channel and the SR where its inclusion improves the transient, forms a binary sub-step block** with `h_min` 0.005 ms (0.0025 ms for Hund–Rudy 2004), triggered by the rate of the junctional Ca (rtol 0.02) and by the open occupancy of the release channel (atol 1e-3). It engages for roughly the first 65 ms of a beat at 1 Hz.
6. **An implicit block that shares states with a sub-step block lies entirely inside it.** gotranx refuses a straddling block; the sub-step block absorbs the rest of the implicit one.
7. **A voltage-dependent Markov chain in an implicit block is nested in the V block when its states go negative at large steps.** Evaluating the chain's rates across sub-steps while holding a Jacobian from the start of the step drove the L-type states of Winslow 1999 negative at 0.1 ms; nesting removes the inconsistency. T-World's IKr did not need it.

The blocks are chosen on the same protocol they are tested on, refined in at most four rounds that change only membership and tolerances, and frozen before timing. The same file is used at every time step.

## blocks.toml

```toml
# Blocks for <model> (round 2; see NOTES.md)
[[implicit_blocks]]
states = ["C0", "C1", "C2", "I_", "O_"]

[[implicit_blocks]]
states = ["caj", "buffers_SLLj", "buffers_SLHj", "casr", "buffers_Csqn"]

[[substep_blocks]]
states = ["v", "m", "h", "j", "mL", "hL"]
h_min = 0.001
mode = "proportional"
rate_triggers = [{ state = "v", atol = 0.5 }]

[[substep_blocks]]
states = ["caj", "buffers_SLLj", "buffers_SLHj", "casr", "buffers_Csqn", "ryr_R", "ryr_O", "ryr_I", "ryr_CaRI"]
h_min = 0.005
mode = "binary"
rate_triggers = [{ state = "caj", rtol = 0.02 }, { state = "ryr_O", atol = 1e-3 }]
```

| table | key | meaning |
|---|---|---|
| `[[implicit_blocks]]` | `states` | state names advanced together by one linearly implicit step |
| `[[substep_blocks]]` | `states` | state names advanced in sub-steps while the rest is held |
| | `h_min` | the shortest sub-step; at most `ceil(dt / h_min)` sub-steps are taken |
| | `mode` | `binary` (all sub-steps or one) or `proportional` (as many as the triggers ask, `ceil(r)` clamped to `[1, n_max]`) |
| | `rate_triggers` | list of `{ state, atol, rtol }`; activity $r = |f_x|\,dt / (atol + rtol\,|x|)$, fires when $r > 1$ |
| | `value_triggers` | list of `{ name, threshold }` on a state or intermediate; fires when `value / threshold > 1` (a negative threshold means "below") |
| | `freeze` | intermediates held at their start-of-step values during the sub-steps |
| | `average` | pass the sub-step-averaged derivatives to the slow states (default true; this is what keeps Na⁺ and Ca⁺⁺ conserved at the block boundary) |
| | `advance_time` | sub-steps see `t + k h` (default true) |

The same tables may live under `[tool.gotranx]` in a `pyproject.toml`; a file with that section must not also have tables at the top level. `gotranx.blocks.validate_blocks` checks that every name is a state, that no state is in two implicit or two sub-step blocks, that no implicit block straddles a sub-step block, that every trigger state is in its block, that `atol` and `rtol` are non-negative with a positive sum, and that a frozen or trigger intermediate survives `remove_unused`. An implicit block whose states do not depend on each other is only a warning.

## Defaults and what they trade

| knob | default | smaller means | larger means |
|---|---|---|---|
| V block `atol` | 0.5 mV | more sub-steps in the upstroke, better activation time, more cost in tissue at the wavefront | fewer sub-steps, peak overshoot, later activation |
| V block `h_min` | 0.001 ms | finer upstroke at the largest steps | caps the refinement |
| dyad `rtol` on junctional Ca | 0.02 | engages earlier and longer in the release, better CaT | misses the onset, CaT amplitude error |
| dyad open-state `atol` | 1e-3 | engages at the first channel openings | engages late |
| dyad `h_min` | 0.005 ms | finer release, more cost per beat | release under-resolved at coarse dt |
| stiffness threshold for blocking a buffer cluster | −20 /ms | more clusters implicit, more cost, no accuracy gain if the mode was not limiting | a stiff cluster left explicit, which brings the stability limit back |

## What `bimgrl propose` automates and what it leaves to you

Automated, each with its justification in `propose.md`: the channel and release chains (rule 3) from the invariants' sub-groups, the stiff binding clusters (rule 1) from the stiffness owners, the SR–junction join (rule 2) from the compartments, the V block (rule 4) from the INa gates named in `meta.toml`, the dyad (rule 5) from the junction, its cluster, the release chain and the SR, and the nesting (rule 6). gotranx `suggest-blocks` runs alongside; its V block is the numeric cross-check of the gate list and its κ-based stiff pairs the cross-check of the clusters.

Left to judgment, and listed under "Judgment calls": whether the SR stays in the dyad (the paper kept it where the transient improved; T-World needed it from round 2), whether a Ca-dependent gate of ICaL belongs in the dyad (TT06's `fCass` is in the paper's block, not in the rule), whether a voltage-dependent chain is nested in the V block (`--nest-v-markov`), and every tolerance change in the tuning rounds. Rule 7 is the one rule the engine does not apply by default, because the paper applied it only where the chain went negative.
