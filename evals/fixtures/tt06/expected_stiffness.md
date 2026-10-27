# Stiffness of `tt06` (protocol B)

1002 checkpoints (the paper's 13 plus every 1 ms); modes with Re λ ≤ -20 /ms. Most negative eigenvalue over all checkpoints: -989.4 /ms.

## Stiff modes by owner set (Table-1 view)

| owners (mean participation) | min Re λ (/ms) | τ (µs) | FE limit (µs) | t at min (ms) | checkpoints | GRL update of owners | RL+FE update of owners | FE share GRL / RL+FE |
|---|---|---|---|---|---|---|---|---|
| m 1.00 | -989.4 | 1.011 | 2.022 | 2000 | 697 | m: exp | m: exp | 0.00 / 0.00 |

## Predicted stability limits

- GRL as generated (exponential wherever the diagonal is non-zero): P1 none (no FE-owned stiff mode); P2 (linearised step) 0.4187 ms at t = 1011 ms
- RL+FE (exponential for HH gates only): P1 none (no FE-owned stiff mode)
- forward Euler everywhere: P2 0.002021 ms

P1 is the forward Euler limit 2|Re λ|/|λ|² of the stiffest mode whose owners the scheme advances explicitly. P2 is the largest step at which the linearised scheme step grows no faster than the model, which also sees coupled modes that GRL updates exponentially one state at a time and damps poorly. Verify with `bimgrl sweep --stability`.
