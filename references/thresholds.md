# Biomarkers, thresholds, the dt grid and the cost

## Biomarkers (`src/bimgrl/vendored/features.py`)

Computed on the compared window: all of protocol A, the last beat of protocol B, on a 0.1 ms grid (or the scheme's step if coarser). The reference is evaluated at the same times, so the reference and the scheme are measured alike.

| biomarker | definition |
|---|---|
| `vpeak` | maximum of V |
| `tact` | the upward crossing of the midpoint between rest (V at the start of the window) and peak, linearly interpolated |
| `apd90`, `apd50` | from `tact` to the first time after the peak where V has repolarised by 90 % (50 %) of `vpeak − vrest`, linearly interpolated at both ends |
| `cat` | the Ca transient amplitude, max Ca minus Ca at the start of the window; its error is relative |
| `cat90` | from `tact` to 90 % decay of the transient |
| `vrest`, `dvdt_max`, `ca_dia`, `nai_dia` | as named |
| `cycle_length` | autorhythmic models only: the mean interval between upstrokes (crossings of the midpoint of V); the last complete cycle gives the other biomarkers, with the maximum diastolic potential in place of rest |

## Threshold sets

| set | Vpeak | APD90 | CaT amplitude | activation time | cycle length (autorhythmic, replaces activation time) |
|---|---|---|---|---|---|
| moderate | 2 mV | 2 ms | 1 % | 0.05 ms | 2 ms |
| tight | 1 mV | 1 ms | 0.5 % | 0.02 ms | 1 ms |

The Ca metric is omitted for a model without a Ca state. A run passes a set when its status is `ok` and every error is within the set. The moderate set is what the survey calls "within the thresholds"; the tight set is where the sweep's descent stops (two consecutive passing steps).

## The physical range

A run is `nan` when any saved state is not finite, `unphysical` when a gate leaves [0, 1] by more than 1e-6 or a concentration goes below 0 by more than 1e-9, after widening the range to what the reference itself visits (a published initial value can be outside: `X_to1 = 2.63` in Winslow 1999). States in `unphysical_exempt` are not checked. Kinds come from the units (`1` or empty is a gate, a unit containing `M` or `mol` is a concentration) with `[kinds]` in `meta.toml` overriding. The sweep reports the first states to leave the range and when.

## The dt grid

Steps are 1, 2 and 5 per decade from 1e-5 to 1 ms, kept only where the step divides the stimulus start, duration and period, the window, and the 0.1 ms sampling grid (above 0.1 ms, the 1 ms grid). Then step n is at `n dt` from an integer counter and exactly `duration / dt` steps are stimulated per period. The sweep runs from the largest step down and stops after two consecutive steps pass the tight set, or at the first run whose cost estimate exceeds the budget (smaller steps cost more). `--dt` overrides the grid with an explicit list and runs all of it.

The stability limit of an explicit scheme (`--stability`) is bisected between the largest step that runs and the smallest that fails, as `dt = duration / k` over allowed `k`, assuming monotonicity.

## Cost

Cost is the time per beat spent in the stepping loop in compiled C (`gcc -O3 -march=native`), timed with `clock_gettime`, one thread, best of 5 (`bimgrl time`), on an idle machine. The sweep's rows carry a calibration estimate (seconds per step from 2000 steps) until they are timed; the report says which it used.

The report chooses, per scheme and threshold set, the **cheapest passing run, not the largest passing step**: heavy sub-stepping can make a larger step slower than a smaller one. The speed-up is the ratio of the baseline's cost to the scheme's at the same threshold set; "GRL fails" means the baseline met the set at no step.
