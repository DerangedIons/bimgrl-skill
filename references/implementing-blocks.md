# Implementing the blocks by hand (native path)

When the model must stay in its own integrator, the agent implements the two block kinds there. This is the recipe the T-World work converged on after two wrong fixes; the anti-patterns at the end are the ones that passed by eye and failed against the reference.

## The implicit block for a pool with its buffers

Write the step for the increments, not for the new values. For a pool $c$ with buffers $B_i$ evaluate at the start of the step the binding rates $d_i = k_{on,i}\, c\,(B_{max,i} - B_i) - k_{off,i} B_i$ and their derivatives $a_i = \partial d_i / \partial c = k_{on,i}(B_{max,i} - B_i)$ and $b_i = -\partial d_i / \partial B_i = k_{on,i} c + k_{off,i}$. The linearly implicit (Rosenbrock–Euler) step is

$$\Delta B_i = h\,(d_i + a_i \Delta c - b_i \Delta B_i), \qquad \Delta c = h\,\Big(f_c - \sum_i (a_i \Delta c - b_i \Delta B_i) - k\,\Delta c\Big),$$

where $f_c$ is the full derivative of $c$ including binding and $k$ collects the pool's own linear sinks that are also taken implicitly (diffusion out of the junction and the sub-sarcolemmal space, the release conductance for the SR). Each buffer exchanges only with its pool, so solve the first equation for $\Delta B_i = h\,(d_i + a_i \Delta c)/(1 + h b_i)$ and substitute it into the second. The binding contribution to the pool becomes $-\sum_i (a_i \Delta c - b_i \Delta B_i) = -\sum_i \frac{a_i \Delta c - h b_i d_i}{1 + h b_i}$, which leaves one scalar equation per pool,

$$\Delta c = \frac{h\Big(f_c + \sum_i \frac{h\, b_i d_i}{1 + h b_i}\Big)}{1 + h k + \sum_i \frac{h a_i}{1 + h b_i}},$$

with no matrix to assemble or factorise. Because the binding term enters the pool and buffer equations with opposite sign and the same increments, $c + \sum_i B_i$ changes only by the non-binding fluxes, at any $h$. The binding mode is multiplied by $1/(1 + h(a + b))$, between 0 and 1 for every $h > 0$, so it is always damped.

## The joint SR–junction solve

The release $J_{rel} = g\,(c_{SR} - c_j)$ with $g = k_s\, O$ (the open fraction at the start of the step, times the volume ratio on the junction side) is the sink of the SR and the source of the junction. Take it implicitly in both: with the buffers eliminated as above, the two pools form a 2×2 system

$$\begin{pmatrix} D_{SR} + h g & -h g \\ -h g\, r & D_j + h g\, r \end{pmatrix} \begin{pmatrix} \Delta c_{SR} \\ \Delta c_j \end{pmatrix} = \begin{pmatrix} N_{SR} \\ N_j \end{pmatrix},$$

with $N$ and $D$ the numerator and denominator of each pool's scalar equation and $r = v_{SR}/v_j$. Solve it in closed form; the determinant is at least $1 + h g (1 + r)$, so it is well posed for any step. Both pools then see the same linearised release $h[J_{rel} + g(\Delta c_{SR} - \Delta c_j)]$ and the release neither creates nor removes Ca.

## The implicit block for a Markov chain

For a chain with occupancies $x$ and rate matrix $A(V, c)$ evaluated at the start of the step, solve $(I - h A)\,\Delta x = h A x$ by Gaussian elimination without pivoting. Pivoting is unnecessary because the off-diagonal entries of $A$ are non-negative rates and every column sums to zero, so $I - hA$ is column diagonally dominant; the elimination then stays stable, and the column-sum property means $\sum_i \Delta x_i = 0$ exactly, which keeps the occupancy sum. In exact arithmetic the step also preserves non-negativity. Recompute $A$ from the current voltage and Ca at every (sub-)step; never carry a Jacobian from the start of a macro step across sub-steps (that inconsistency is what sent Winslow's L-type states negative).

## Sub-step blocks

Within a step of length $h$, advance the fast states in $n$ sub-steps while every other state is held at its start-of-step value and the currents from outside the block are frozen. Decide $n$ per cell from its own state at every step: for the V block $n = \mathrm{clamp}(\lceil |dV/dt|\, h / \Delta V_{tol} \rceil, 1, \lceil h / h_{min} \rceil)$ with $\Delta V_{tol}$ about 1 mV (the `atol` of the rate trigger); for the dyad all $\lceil h / h_{min} \rceil$ sub-steps when the junctional Ca changes by more than 2 % per step or an open fraction exceeds 1e-4, else one. Inside the sub-steps recompute only what the block owns: $I_{Na}$ and $I_{NaL}$ from the sub-step voltage and gates; the release, the RyR rates and the junction–SR implicit step from the sub-step Ca.

Pass the exchange fluxes back averaged over the sub-steps: the junctional and sub-sarcolemmal Na⁺ pools receive $I_{Na} + I_{NaL}$ averaged, so the Na⁺ entering the cell matches the charge that moved the membrane potential; the sub-sarcolemmal pool receives the junction-to-SL diffusion averaged, so no Ca is lost at the block boundary. Without the averaging the slow states see the start-of-step flux and the block leaks mass in proportion to how much the fast states moved. This is what gotranx's `average = true` does.

## Anti-patterns that passed by eye

1. **Per-variable exponential update on an exchange pair.** Applying GRL's exponential update separately to $c$ and to $B$ multiplies the two halves of one binding exchange by different factors $\varphi(a h) \neq \varphi(b h)$, so the exchange no longer cancels and mass leaks. At dt 0.005 ms the SR nadir moved from 0.459 to 0.504 mM. The action potential looked fine.
2. **Backward Euler in absolute form.** Writing the implicit step for the new values rather than the increments and linearising the product $k_{on}\, c\,(B_{max} - B)$ in both equations counts the product term twice. CaT came out at 582 instead of 327 nM.

Neither would have been noticed by "does it run and look like an AP". Both were caught by comparing SR load and Ca-transient metrics with Rodas5P at reltol 1e-8 on the same right-hand side. The lesson is the validation checklist.

## Validation checklist

1. **Bit-for-bit no-op.** With every block switched off, the new integrator must reproduce the old one exactly at the old step. With the blocks on at a step small enough that no sub-step is ever taken, the implicit blocks alone must change the result only at round-off level of what the old step already resolved.
2. **At the old step.** Run the blocks at the old dt and compare SR load, CaT amplitude and diastolic Ca with the adaptive reference (Rodas5P or LSODA, rtol ≤ 1e-8); differences here are implementation errors, not step-size errors.
3. **At the new step.** Sweep dt upward and apply the moderate and tight thresholds of `thresholds.md` on beat 10; run 100 to 600 beats and watch every invariant (chain sums, total Ca, Na⁺ and K⁺ drift) with `bimgrl validate` on your trace.
4. **Threshold phenomena.** Restitution, alternans onset, EAD threshold and the beat of the first EAD, DAD threshold. These test the dynamics, not the end state, and are where a scheme that is "almost right" shows.
