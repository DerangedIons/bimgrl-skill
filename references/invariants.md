# Invariants and conservation checks

`bimgrl invariants` answers two questions a Jacobian cannot: what the model conserves by construction, and whether what should be conserved actually is. Both come from the stoichiometry of the equations rather than from their linearisation.

## Why not the Jacobian's null space

A linear invariant $c^\top x$ needs $c^\top f(x) \equiv 0$ for every $x$. A left null vector of $J(x)$ only says $c^\top f$ is stationary at that one point; a vector found at diastole is not an invariant in the upstroke. And the right null space, which TWorld's `conservation_law_analysis.jl` computed, is the set of directions $f$ is insensitive to, which has nothing to do with conserved combinations of equations.

## The stoichiometric view

Write every derivative as $f_i = \sum_k S_{ik}(p)\, v_k(x, p)$. The **atoms** $v_k$ are the monomials of $f_i$ in the states and the non-affine intermediates: a current, a flux, a product rate × state, a binding term. $S$ holds only parameters. A weight vector with $c^\top S = 0$ then gives $c^\top f \equiv 0$ for every state and every rate function, which is why a Markov chain with voltage-dependent rates conserves its occupancy by construction: the rate is inside the atom, not in $S$.

To make the same current the same atom in every equation, intermediates that are affine sums of other intermediates (`ICa_tot_junc = ICaL_junc + ICab_junc + ...`, `I_tot = ...`) are inlined first; parameter-only intermediates (`vjunc = 0.00539 * vcell`) stay as symbols so that the weights read `vjunc/vmyo` rather than a float. An equation whose every term shares a common factor, the rapid-buffering factor $g(c) = 1/(1 + dB/dc)$ of a pool without buffer states (TT06's `Ca_i_bufc`), has that factor divided out and recorded: $f = g(c)\,(\ldots)$ means $\frac{d}{dt}[c + B(c)] = (\ldots)$, and the report prints the buffered total $c + B(c)$ with $B = \int (1/g - 1)\,dc$ when sympy can integrate it. Atoms are **shared** (in two or more equations), **own** (in one), or **boundary** (a current in dV/dt, or anything containing the stimulus).

## Sub-groups: what an implicit block should contain

A **Markov chain** is a structurally linear subsystem of dimensionless states (gotranx's `linear_subsystems`). A **binding cluster** is a free concentration with every bound state it shares a binding atom with. For each sub-group the *internal* atoms (shared by its members only) must have zero column sums under the weights $c$, found as the null space of the internal sub-matrix over the parameter field (with a timeout, then exact rationals at the default parameter values).

| class | meaning |
|---|---|
| `exact` | no external atom: the weighted sum is a constant of the ODE (T-World's IKr `C0 + C1 + C2 + I_ + O_`) |
| `conserved-internally` | the binding or the chain conserves, and the group exchanges with the rest through its external atoms (every pool with its buffers: release, diffusion, uptake, currents) |
| `closure` | no null vector among the ODE states because the chain is closed by an intermediate `const - sum of states` (T-World's `RIcleft = 1 - R - O - I - CaRI`); conserved by construction, and the implicit block takes the ODE states only; the monitor checks the bound `sum <= const` |
| `improper` | structurally a chain or a cluster, yet no full-support weights exist: the report prints, atom by atom, the columns and their weighted sum, and for a cluster the single atom whose removal restores conservation |

A binding cluster whose weights are not all 1 gets a **warning**: a pool and its own buffers share a compartment and a volume, so a weight of 0.5 means one binding equation is scaled differently from the other (a stray factor, or a volume ratio on one side only).

## Networks: total-ion invariants

Joining every state linked by a shared atom (release, diffusion, uptake, binding) gives the exchange networks. Their null vectors are the total-ion invariants, with the volume ratios as weights:

| class | meaning |
|---|---|
| `exact` | nothing enters or leaves (a closed chain, again) |
| `boundary` | only membrane currents change the total: the report lists them with their coefficients, which is the statement $\frac{d}{dt}(\text{total Ca}) = -\frac{C_m}{2 F v_{myo}}(I_{CaL} + I_{Cab} + I_{pCa} - 2 I_{NaCa})$ |
| `open` | an atom of the network also feeds something that carries no mass, or an own atom remains: not conserved, hidden from the table, kept in the JSON |

For T-World the Ca network comes out as `cai + buffers(cytosol) + (vjunc/vmyo)(caj + SLLj + SLHj) + (vsl/vmyo)(casl + SLLsl + SLHsl) + (vsr/vmyo)(casr + Csqn)`, which is exactly the `total_ca` the survey wrote by hand with 0.000539/0.65 and the other ratios; the Na network likewise with the Na buffers. For TT06 the buffered totals of the three rapid-buffered pools combine into `Ca_SR + (V_c/V_sr) Ca_i + (V_ss/V_sr) Ca_ss`, boundary class.

## Charge bookkeeping

Every current $I_k$ in $dV/dt$ moves ions. For every concentration equation the coefficient of $I_k$ divided by its coefficient in $dV/dt$ is a multiple of one compartment constant $C_m/(z F\, \mathrm{vol})$; the smallest such multiple in each equation defines one charge unit, so every current gets integer counts per ion in charge units. A properly written current sums to 1 (INaK: Na⁺ 3, K⁺ −2; INaCa: Na⁺ 3, Ca⁺⁺ −2; Cl⁻ currents count with the anion sign). Statuses:

| status | meaning |
|---|---|
| `ok` | the counts sum to 1 |
| `unaccounted` | the current is in dV/dt and in no concentration equation: a slow drift source. For the stimulus the report estimates the K⁺ drift per beat (amplitude × duration × the K⁺ compartment constant) and names the remedy of O'Hara–Rudy 2011, which assigns the stimulus to K⁺ |
| `miscounted` | the counts sum to something else, or are not integers: a wrong valence, compartment volume or factor |
| `charge-neutral` | sums to 0: an exchanger written without its net charge |

A clamped ion (a concentration whose derivative is identically zero, or a parameter) is skipped and said so.

## The drift monitor

`validate.py` evaluates every invariant with numeric weights at the start of every beat of a run and flags it by class: an `exact` invariant is measured against its own initial value, where round-off is the only legitimate change, floor 1e-7 relative per beat; a `boundary` invariant is measured against the reference's value at the same beat, where the integration error of the currents is legitimate, floor 1e-3 relative per beat, the level at which it starts to move the biomarkers; a `closure` chain gets the bound check on its sum. Sub-groups that exchange with the rest of the model are not checked alone. GRL's IKr drift in T-World (the prepaced state sums to 1.0454) and the 0.19 % release leak of two separate SR and junction blocks are far above these floors; LSODA at rtol 1e-10 is far below.

## Worked readings

**T-World, as published.** IKr `exact`; both RyR chains `closure` through `RIcleft` and `RIcleft_p`; the contraction chain `conserved-internally` with the Land closure `XU` (not an ion channel, left explicit by the block rules); binding clusters for casl, caj, casr, nai/naj/nasl with weights 1; the Ca and Na networks `boundary`; every current `ok` including `Istim` into `ki`; no flags.

**T-World with `dC0_dt = -C0*IKr_alpha + 1.25*C1*IKr_beta`.** The IKr chain is `improper`: "atom `C1*IKr_beta`: columns {C0: 1.25, C1: −1} sum to 0.25 under weights all 1". The atom and the equation it enters wrong are named; the fix is one coefficient.

**T-World with `dbuffers_SLLj_dt = 2.0*d_buffers_SLLj`.** The junction cluster still has a null vector, `caj + 0.5 buffers_SLLj + buffers_SLHj`, so it is not `improper`, but the warning fires: "weights other than 1 inside one compartment: {buffers_SLLj: 0.5}". A weight is a volume ratio, and the free and the bound form of one buffer in one compartment cannot have different volumes.

**Gauthier-type failure.** A chain reported `exact` by the equations whose sum nevertheless drifts under GRL in the baseline sweep is the scheme's fault, not the model's; the drift flag on the GRL rows is the finding, and the chain's implicit block is the fix.
