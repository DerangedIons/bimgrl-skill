# Tissue, GPUs and other kinds of models

A block configuration chosen on a single cell is a starting point for anything bigger. Two things change in tissue: the cell step is no longer the whole cost, and cells at different phases share one time step. Both change which blocks pay. The numbers below are from upstream `DerangedIons/blocked-rush-larsen` (`AGENTS.md`, `examples/tissue_monodomain.py`).

## Tissue simulations (2D and 3D)

The pattern is operator splitting (Godunov): a PDE step for diffusion, then the cell step at every node. In fenicsx-beat:

```python
import beat, beat.c_backend
from gotranx.blocks import load_blocks

blocks = load_blocks("blocks.toml")
kwargs = dict(implicit_blocks=blocks.implicit, substep_blocks=blocks.substep,
              base_scheme="rush_larsen", stiff_states=blocks.stiff_states, freeze_boundary_fluxes=True)
fun = beat.c_backend.from_ode(ode, "block_rush_larsen", codegen_kwargs=kwargs, comm=comm)
odesolver = beat.odesolver.DolfinODESolver(v_ode=..., v_pde=pde.state, fun=fun, ...)
solver = beat.MonodomainSplittingSolver(pde=pde, ode=odesolver)
```

- **Use the C backend** (`beat.c_backend.from_ode`, optionally `openmp=True`). It loops over the nodes, so only the nodes that need sub-steps take them. The numpy code steps all the nodes as one array, at the largest number of sub-steps any node needs.
- **The splitting limits the time step, not the cell model.** In the 2D example (TT06 epi, h = 0.25 mm) the activation-time error against RL at dt 0.01 ms grows with dt for both schemes: RL 0.85 ms at 0.05 and 1.8 ms at 0.1; BIRL-M 0.95 ms at 0.1 and 3.9 ms at 0.2, although BIRL-M is accurate to 0.5 ms in a single cell. Choose dt by the error in tissue (activation and repolarisation times, an ECG), not by single-cell biomarkers.
- **Drop sub-step blocks whose cost does not shrink with dt.** In tissue some node is always in release, so a dyad block costs at every step: the cell step of BIRL-M took 1.85 s at dt 0.1 ms with it and 0.44 s without (RL: 0.28 s). Within 1 ms of activation error, BIRL-M with the V and calcium blocks was 1.6× faster than RL; with the dyad block, slower than RL.
- **The PDE step does not get cheaper.** The gain is bounded by the share of the cell step in the total. With the GRL-base scheme (tag `grl-edition`), the PDE step was about two thirds of the time with TT06 in the Niederer slab (speed-up over GRL about 4× at 1 ms activation error); with T-World in a 903,010-node biventricular mesh, where the cell step dominates, BIMGRL at dt 0.1 ms was 17× faster than GRL at 0.005 ms.
- **Solve the PDE step tightly** (`ksp_rtol` 1e-10). With PETSc's default 1e-5 the error accumulates over the steps; in the first study a patch of tissue fired by itself 6 ms early.
- **Prepare the initial states** by pacing a single cell with the same scheme and time step.

Tissue with the RL base has only been tried in the 2D example; the larger tissue studies used the GRL base.

## GPUs and other backends

Every gotranx code generator takes `--scheme block_rush_larsen --blocks <file>`; `bimgrl codegen --backend` wraps the ones marked.

| backend | command | sub-step blocks |
|---|---|---|
| C | `ode2c` (`codegen --backend c`); fenicsx-beat `beat.c_backend` | yes, one loop per cell |
| C++ | `ode2cpp` (`--backend cpp`; `--float-type float` for single precision) | yes |
| CUDA | `ode2cuda` (`--backend cuda`): a header of kernels, one thread per cell | yes, see below |
| numpy | `ode2py` (`--backend python`); arrays of cells | yes, all cells take the largest number |
| JAX | `ode2py --backend jax` | yes, see below |
| Julia | `ode2julia`, `ode2cytozoo` (`--backend julia`, `cytozoo`) | yes |
| UFL | `ode2ufl` | no, implicit blocks only |

The threads of a GPU warp run in lockstep, which splits the two block kinds:

- **Implicit blocks** are branch-free, the same work for every cell, and suit a GPU.
- **Sub-step blocks** make each thread loop over its own number of sub-steps, so a warp takes as long as its slowest cell; near a wavefront the cells differ. Keep sub-step blocks few and small, and compare with a version without them (more implicit blocks and a smaller step).
- **JAX** runs both branches of a conditional under `vmap`, so a batch pays for its largest number of sub-steps; `jax.grad` does not go through the sub-step loop, `jax.jacfwd` does.
- **Single precision** has not been tested with implicit blocks, which solve without pivoting; check against the double-precision results.

gotranx's tests cover the generated CUDA and C++ code, but it has not been run on a GPU or timed, and neither has JAX.

## Other kinds of models

Nothing in the scheme is specific to the heart; what carries over depends on the model.

- **Models with Hodgkin–Huxley gates** (neurons, other excitable cells) fit the RL base directly: the gates go in `stiff_states`, and the membrane potential with its fast gates makes a V block (point `suggest-blocks` at it with `--membrane-potential`).
- **Models without gates** (chemical kinetics, signalling, pharmacology): RL is plain forward Euler, so the blocks do all the work. Put each stiff, linearly coupled group in an implicit block, and sub-step a fast transient with a `value_triggers` entry on the flux that drives it. `suggest-blocks` finds no V block in a model without a membrane potential, which is fine.
- **Check the assumption of the implicit solve**: column diagonal dominance of $I - hJ$ (pools, buffers, chains; `block-rules.md`). For another kind of block, check against the reference.
- **The biomarkers are cardiac** (stimulus protocols, action-potential biomarkers). For another model, write the reference and the error measure for what matters in it before choosing blocks.

## Limitations

- The scheme is first order, like RL and GRL.
- With the RL base and frozen boundary fluxes, what crosses the boundary of an implicit block is conserved to round-off (total Ca in a closed T-World cell changes by at most 1e-13). A flux between a sub-step block and a state outside it is not: the outside state sees the flux averaged over the starts of the sub-steps.
- The survey checks single cells over 10 and 100 beats at 1 Hz (rodent models at their paper's rate); alternans, early afterdepolarisations and other pacing rates have not been tested with the RL base.
