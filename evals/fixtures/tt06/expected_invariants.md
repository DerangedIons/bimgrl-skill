# Invariants and conservation checks of `tt06`

43 atoms after inlining 4 affine intermediates.

Equations with a common factor divided out (rapid buffering): `Ca_SR`: Ca_sr_bufsr; `Ca_i`: Ca_i_bufc; `Ca_ss`: Ca_ss_bufss; `Xr1`: 1/tau_xr1; `Xr2`: 1/tau_xr2; `Xs`: 1/tau_xs; `d`: 1/tau_d; `f2`: 1/tau_f2; `fCass`: 1/tau_fCass; `f`: 1/tau_f; `h`: 1/tau_h; `j`: 1/tau_j; `m`: 1/tau_m; `r`: 1/tau_r; `s`: 1/tau_s

Their conserved quantity enters through the buffered total: `Ca_SR` -> `1.0*(-Buf_sr*K_buf_sr + Ca_SR*(Ca_SR + K_buf_sr))/(Ca_SR + K_buf_sr)`; `Ca_i` -> `1.0*(-Buf_c*K_buf_c + Ca_i*(Ca_i + K_buf_c))/(Ca_i + K_buf_c)`; `Ca_ss` -> `1.0*(-Buf_ss*K_buf_ss + Ca_ss*(Ca_ss + K_buf_ss))/(Ca_ss + K_buf_ss)`

## Flags

- none: every chain and binding cluster conserves, every current is accounted for

## Sub-groups (candidate implicit blocks)

| name | structure | class | members | weights | external atoms | numeric check |
|---|---|---|---|---|---|---|

## Networks (total-ion invariants)

| name | class | ions | members | weights | residual (d/dt of the total) | numeric check |
|---|---|---|---|---|---|---|
| network:Ca:Ca_SR | boundary | Ca | Ca_SR, Ca_i, Ca_ss | Ca_SR: 1, Ca_i: V_c/V_sr, Ca_ss: V_ss/V_sr | 1.0*Cm/(F*V_sr) * i_NaCa; -0.5*Cm/(F*V_sr) * i_b_Ca; -0.5*Cm/(F*V_sr) * i_p_Ca; -0.5*Cm/(F*V_sr) * i_CaL | - |

## Algebraic closures (a state written as const minus the others)

- none

## Charge bookkeeping of the currents in dV/dt

| current | ion counts (charge units) | total | status |
|---|---|---|---|
| i_K1 | K_i: 1 | 1 | ok |
| i_Ks | K_i: 1 | 1 | ok |
| i_b_Ca | Ca_i: 1 | 1 | ok |
| i_p_Ca | Ca_i: 1 | 1 | ok |
| i_NaCa | Ca_i: -2, Na_i: 3 | 1 | ok |
| i_NaK | Na_i: 3, K_i: -2 | 1 | ok |
| i_p_K | K_i: 1 | 1 | ok |
| i_Na | Na_i: 1 | 1 | ok |
| i_Kr | K_i: 1 | 1 | ok |
| i_to | K_i: 1 | 1 | ok |
| i_b_Na | Na_i: 1 | 1 | ok |
| i_CaL | Ca_ss: 1 | 1 | ok |
