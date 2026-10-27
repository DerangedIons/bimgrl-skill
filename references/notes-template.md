# NOTES.md and the tuning log

Every model directory keeps a `NOTES.md` in the survey's section convention. `bimgrl log` creates it from the skeleton below if it is missing, appends tuning rounds, and adds lines to sections; nothing in it is overwritten.

## Skeleton

```markdown
# <model name>

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
```

## Tuning-log entries

```
bimgrl log <model_dir> --round N --change "<what changed>" --result "<the numbers>" --decision "<adopted or rejected, and why>" [--freeze]
```

writes one bullet under `## Tuning log`:

```
- Round N (YYYY-MM-DD, blocks.toml): <change>. Result: <numbers>. Decision: <decision>
```

Two entries in the style of the survey's T-World notes, which tell the story the log is for:

```
- Round 1 (2026-10-01, blocks.toml): implicit blocks for the IKr chain, the SL pool, the junction pool and the SR, each with its own buffers; V/INa proportional (atol 0.5, h_min 0.001); dyad binary (h_min 0.005) with caj rtol 0.02, ryr_O and ryr_O_p atol 1e-3. Result: protocol B, dt 0.5 / 0.2 / 0.1 / 0.05 all ok; dt 0.2: Vpeak +0.53, APD90 +0.99, CaT -0.19 %, t_act +0.008 (within moderate); dt 0.1: +0.38, +0.55, CaT +1.35 %. GRL ok at 0.005 (+0.33, +0.03, -0.28 %), unphysical at 0.01. Decision: adopted as the baseline; the positive CaT error at 0.1 is suspicious (a leak would be positive).
- Round 2 (2026-10-14, blocks.toml): casr and buffers_Csqn moved into the junction block and into the dyad sub-step block, so that the release is solved jointly in both pools. Result: protocol B, dt 0.5 / 0.2 / 0.1 / 0.05: CaT -6.29 / -1.21 / +0.43 / -0.06 % (round 1: -2.92 / -0.19 / +1.35 / +0.37 %); Vpeak, APD90, t_act unchanged within 0.05 mV / 0.05 ms / 0.001 ms; about 5 % more time per step at 0.1. Decision: adopted. The error is now monotone in dt; round 1 met the moderate CaT threshold at 0.2 only because the release leak (+1.35 % at 0.1) cancelled a negative error that grows with dt. Moderate thresholds now met at dt 0.1, not 0.2.
```

The second entry is the reason the log exists: a configuration that passed the thresholds by cancellation would have been frozen without it.

## The blocks.toml header

The first line of `blocks.toml` names the round it belongs to, so that a file and a log entry can always be matched:

```
# Blocks for tworld (round 2; see NOTES.md)
```

`bimgrl log --freeze` stamps the current round into that header; `bimgrl propose` writes round 1. Old rows in `runs.csv` keep the configuration key of the blocks they were run with, so a re-run after a change never mixes rounds.
