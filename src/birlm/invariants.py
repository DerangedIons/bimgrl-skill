"""Linear invariants and conservation-law checks, from the stoichiometry of the equations

Write every state derivative as ``f_i = sum_k S_ik(p) * v_k(x, p)``: the *atoms* ``v_k`` are the
monomials of ``f_i`` in the states and the non-affine intermediates (a current, a flux, a
``rate * state`` product, a binding term), and ``S`` holds only parameters. A weight vector
``c`` with ``cᵀS = 0`` then gives ``cᵀf ≡ 0`` for every state and every rate function, so a
Markov chain with voltage-dependent rates conserves its occupancy by construction, and a pool
with its buffers conserves mass with volume-ratio weights. The left null space of the
Jacobian would only hold at one point (and TWorld's attempt took the right one).

Two levels:

* **sub-groups** (what an implicit block should contain): each structurally linear subsystem
  (Markov chain) and each binding cluster (a pool with the buffers it shares a binding atom
  with). Their *internal* atoms (shared by members only) must have zero column sums under the
  weights; a sub-group that is structurally a chain or a binding cluster but has no such
  weights is flagged ``improper`` with the offending atoms, which is where a typo, a missing
  volume ratio or an inconsistent rate shows up.
* **networks** (what is conserved overall): the groups joined by any shared atom (release,
  diffusion, uptake). Their null vectors are the total-ion invariants: ``exact`` when nothing
  else enters, ``boundary`` when only membrane currents do (``d(total)/dt = currents``).

Plus the **charge bookkeeping**: every current in ``dV/dt`` should move exactly one charge
unit's worth of ions in the concentration equations; a current that moves none is a drift
source (the stimulus, background currents), one that moves another amount is miscounted.
"""

from __future__ import annotations

import re
import signal
from collections import defaultdict

import numpy as np
import sympy

from .model import Model
from .workspace import Workspace, md_table

_BOUND = re.compile(r"buf|trap|bound|csqn|calm|cam\b|camk|casig|trpn|tnc|myo|sig", re.I)


class _Timeout(Exception):
    pass


def _ev(e: sympy.Basic) -> sympy.Basic:
    """Rebuild an expression with evaluation: gotranx builds nested, unevaluated Adds and Muls,
    on which expand() and coeff() do nothing"""
    if not getattr(e, "args", None):
        return e
    return e.func(*[_ev(a) for a in e.args])


def _with_timeout(seconds: int, fn, *args):
    def handler(signum, frame):
        raise _Timeout()

    old = signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)
    try:
        return fn(*args)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


# ----------------------------------------------------------------------------------------
# atoms
# ----------------------------------------------------------------------------------------


def _affine_in(expr: sympy.Expr, gens: set[sympy.Symbol]) -> bool:
    """Whether `expr` is ``sum coef_j * g_j + const`` with parameter-only coefficients"""
    for term in sympy.Add.make_args(expr):
        coef, rest = term.as_independent(*gens, as_Add=False)
        if rest != 1 and not (isinstance(rest, sympy.Symbol) and rest in gens):
            return False
    return True


def extract(model: Model) -> dict:
    """The stoichiometry: atoms, coefficients, and which atoms are shared, own or boundary"""
    ode = model.ode
    state_syms = {s.name: s.symbol for s in ode.states}
    sym_state = {v: k for k, v in state_syms.items()}
    gens_all = set(state_syms.values()) | {i.symbol for i in ode.intermediates}
    # anything that is neither a state, a parameter nor an intermediate (the time, a missing
    # variable) is a generator too, so that a time-dependent expression is never a constant
    known = gens_all | {q.symbol for q in ode.parameters}
    for a in ode.sorted_assignments():
        gens_all |= a.expr.free_symbols - known

    # inline the affine intermediates (sums of currents and fluxes), in dependency order
    inlined: dict[sympy.Symbol, sympy.Expr] = {}
    affine_names = []
    constant_names = []
    for a in ode.sorted_assignments():
        if a.__class__.__name__ != "Intermediate":
            continue
        expr = _ev(a.expr)
        if not (expr.free_symbols & gens_all):
            constant_names.append(a.name)  # a parameter-only intermediate: stays a symbol
            gens_all.discard(a.symbol)
            continue
        expr = expr.xreplace(inlined) if inlined else expr
        if _affine_in(expr, gens_all):
            inlined[a.symbol] = sympy.expand(expr)
            affine_names.append(a.name)
    constant_syms = {i.symbol for i in ode.intermediates if i.name in constant_names}
    leaf_syms = {i.symbol for i in ode.intermediates if i.symbol not in inlined and i.symbol not in constant_syms}
    gens = set(state_syms.values()) | leaf_syms
    sym_leaf = {i.symbol: i.name for i in ode.intermediates}
    leaf_expr = {i.symbol: _ev(i.expr) for i in ode.intermediates if i.symbol not in inlined}

    # the dependence of every leaf intermediate on the states, through other intermediates
    from gotranx.block_detection import state_dependencies

    deps = state_dependencies(ode)

    stim_syms = {s for s in gens_all if "stim" in str(s).lower()}

    # atoms per equation. A factor common to every term of an equation (the rapid-buffering
    # factor 1/(1 + dB/dc) of a Ca pool without buffer states) is divided out and recorded, so
    # that the same current or flux is the same atom in every equation it enters
    atoms: dict[sympy.Expr, dict[str, sympy.Expr]] = defaultdict(dict)  # atom -> {state: coef}
    row_factor: dict[str, sympy.Expr] = {}
    for d in ode.state_derivatives:
        expr = _ev(d.expr)
        expr = sympy.expand(expr.xreplace(inlined) if inlined else expr)
        terms = []
        for term in sympy.Add.make_args(expr):
            coef, atom = term.as_independent(*gens, as_Add=False)
            if atom == 0 or coef == 0:
                continue
            terms.append((coef, atom))
        if not terms:
            continue
        common = terms[0][1] if len(terms) >= 2 else sympy.Integer(1)
        try:
            for _, atom in terms[1:]:
                common = sympy.gcd(common, atom)
        except Exception:  # a Piecewise (guard, stimulus) in a term: no common factor
            common = sympy.Integer(1)
        if common != 1 and common.free_symbols and not (common.free_symbols & set(state_syms.values())):
            # a common leaf-intermediate factor (no state symbol directly)
            row_factor[d.state.name] = common
            terms = [(coef, sympy.cancel(atom / common)) for coef, atom in terms]
        for coef, atom in terms:
            row = atoms[atom]
            row[d.state.name] = row.get(d.state.name, 0) + coef
    # drop cancelled entries
    for atom in list(atoms):
        atoms[atom] = {s: sympy.cancel(c) for s, c in atoms[atom].items() if sympy.cancel(c) != 0}
        if not atoms[atom]:
            del atoms[atom]

    v = model.meta.voltage
    v_atoms = {atom for atom, row in atoms.items() if v in row and atom != 1}

    def atom_states(atom) -> set[str]:
        names = set()
        for s in atom.free_symbols:
            if s in sym_state:
                names.add(sym_state[s])
            elif s in leaf_expr:
                names |= set(deps.get(s, frozenset()))
        return names

    def classify(atom) -> str:
        if atom in v_atoms or (atom.free_symbols & stim_syms):
            return "boundary"
        if atom == 1:
            return "constant"
        if len(atoms[atom]) >= 2:
            return "shared"
        return "own"

    kinds = {atom: classify(atom) for atom in atoms}

    # a row factor g(x_i) that depends only on its own state is a rapid-buffering factor:
    # f_i = g_i * (...) means d/dt [x_i + B_i(x_i)] = (...) with B_i = int (1/g_i - 1) dx_i
    buffered_total: dict[str, str] = {}
    for name, g in row_factor.items():
        x = state_syms[name]
        g_expr = g.xreplace({s: leaf_expr[s] for s in g.free_symbols if s in leaf_expr})
        # resolve nested intermediates a few levels
        for _ in range(5):
            nxt = g_expr.xreplace({s: leaf_expr[s] for s in g_expr.free_symbols if s in leaf_expr})
            if nxt == g_expr:
                break
            g_expr = nxt
        if g_expr.free_symbols & set(state_syms.values()) <= {x}:
            try:
                B = _with_timeout(20, lambda e: sympy.integrate(sympy.cancel(1 / e - 1), x), g_expr)
                if not B.has(sympy.Integral):
                    B0 = B.subs(x, 0)
                    if B0.is_finite:
                        B = B - B0  # bound buffer is zero at zero free concentration
                    buffered_total[name] = str(sympy.simplify(x + B))
            except Exception:
                pass
    return dict(
        atoms=atoms,
        kinds=kinds,
        atom_states={atom: atom_states(atom) for atom in atoms},
        affine_inlined=affine_names,
        constant_intermediates=constant_names,
        leaf_names=sym_leaf,
        state_syms=state_syms,
        v_atoms=v_atoms,
        row_factor={k: _label(v, sym_leaf) for k, v in row_factor.items()},
        buffered_total=buffered_total,
    )


def _label(atom, leaf_names) -> str:
    return "1" if atom == 1 else str(atom.xreplace({s: sympy.Symbol(n) for s, n in leaf_names.items()}))


# ----------------------------------------------------------------------------------------
# weights
# ----------------------------------------------------------------------------------------


def _numeric_params(model: Model, constant_names=()) -> dict[sympy.Symbol, sympy.Rational]:
    p = model.p0
    index = model.ns["parameter_index"]
    out = {q.symbol: sympy.Rational(repr(float(p[index(q.name)]))) for q in model.ode.parameters}
    names = set(constant_names)
    for a in model.ode.sorted_assignments():
        if a.name in names:
            try:
                out[a.symbol] = sympy.Rational(repr(float(a.expr.xreplace(out).evalf())))
            except (TypeError, ValueError):
                pass
    return out


def nullspace(matrix: sympy.Matrix, params: dict, timeout: int = 30) -> tuple[list, str]:
    """Null vectors of ``matrixᵀ`` (rows of the matrix are states): symbolic over the parameter
    field with a timeout, else exact-rational with the default parameter values"""
    if matrix.rows == 0 or matrix.cols == 0:
        return [sympy.Matrix([1])] if matrix.rows == 1 else [], "trivial"
    try:
        vecs = _with_timeout(timeout, lambda m: m.T.nullspace(), matrix)
        return [sympy.Matrix([sympy.cancel(x) for x in v]) for v in vecs], "symbolic"
    except (_Timeout, Exception):
        numeric = matrix.xreplace(params)
        return numeric.T.nullspace(), "numeric-rational"


def _normalise(vec: sympy.Matrix, members: list[str]) -> dict[str, sympy.Expr]:
    pivot = next((x for x in vec if x != 0), None)
    if pivot is None:
        return {}
    return {m: sympy.cancel(x / pivot) for m, x in zip(members, vec) if sympy.cancel(x / pivot) != 0}


def _full_support(vecs, n: int):
    """A combination of null vectors with no zero entry, if the span has one"""
    if not vecs:
        return None
    for v in vecs:
        if all(x != 0 for x in v):
            return v
    # try the sum (generic combination)
    total = sum(vecs[1:], vecs[0])
    if all(x != 0 for x in total):
        return total
    return None


# ----------------------------------------------------------------------------------------
# groups
# ----------------------------------------------------------------------------------------


def _union_find(names, edges):
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    groups = defaultdict(list)
    for n in names:
        groups[find(n)].append(n)
    return [sorted(g) for g in groups.values()]


def _ion(name: str) -> str | None:
    n = name.lower().replace("_", "")
    if _BOUND.search(name):
        return None
    for prefix, ion in (("ca", "Ca"), ("na", "Na"), ("cl", "Cl"), ("k", "K"), ("mg", "Mg")):
        if n.startswith(prefix):
            return ion
    return None


def _closures(model: Model, state_syms: dict) -> list[dict]:
    """Intermediates of the form ``const - sum of states`` (the closing state of a chain)"""
    out = []
    syms = set(state_syms.values())
    for i in model.ode.intermediates:
        e = sympy.expand(_ev(i.expr))
        if not (e.free_symbols & syms) or e.free_symbols - syms:
            continue
        coeffs = {str(x): e.coeff(x) for x in e.free_symbols & syms}
        rest = sympy.expand(e - sum(e.coeff(x) * x for x in e.free_symbols & syms))
        if len(coeffs) >= 2 and all(c == -1 for c in coeffs.values()) and rest.is_number and rest > 0:
            out.append(dict(name=i.name, expr=str(i.expr), states=sorted(coeffs), total=float(rest)))
    return out


def analyse(model: Model, numeric_check: bool = True) -> dict:
    ex = extract(model)
    closures = _closures(model, ex["state_syms"])
    closed_by = {}
    for c in closures:
        for st in c["states"]:
            closed_by[st] = c["name"]
    atoms, kinds, atom_states, leaf_names = ex["atoms"], ex["kinds"], ex["atom_states"], ex["leaf_names"]
    params = _numeric_params(model, ex["constant_intermediates"])
    names = model.names
    kind_of = dict(zip(names, model.kinds))
    bound = {n for n in names if _BOUND.search(n) and kind_of[n] == "conc"}
    free_conc = {n for n in names if kind_of[n] == "conc" and n not in bound}

    def matrix(members, atom_list):
        return sympy.Matrix([[atoms[a].get(m, 0) for a in atom_list] for m in members])

    def label_cols(atom_list):
        return [_label(a, leaf_names) for a in atom_list]

    results = []

    # ---- sub-groups: Markov chains and binding clusters ---------------------------------
    from gotranx.block_detection import linear_subsystems

    dimless = {n for n in names if kind_of[n] in ("gate", "markov")}
    chains = [list(g) for g in linear_subsystems(model.ode) if set(g) <= dimless and len(g) >= 2]
    # binding pairs: a free pool and a bound state that share an atom depending on both
    edges = []
    for atom, row in atoms.items():
        if kinds[atom] != "shared":
            continue
        touched = set(row)
        pools = touched & free_conc
        bounds_ = touched & bound
        for p in pools:
            for b in bounds_:
                if {p, b} <= atom_states[atom]:
                    edges.append((p, b))
    clusters = [g for g in _union_find(sorted(free_conc | bound), edges) if len(g) >= 2]

    subgroups = [("markov", g) for g in chains] + [("binding", g) for g in clusters]
    for kind, members in subgroups:
        mset = set(members)
        internal = [a for a, row in atoms.items() if kinds[a] == "shared" and set(row) <= mset]
        external = [a for a, row in atoms.items() if (set(row) & mset) and a not in internal]
        M = matrix(members, internal)
        vecs, how = nullspace(M, params)
        vec = _full_support(vecs, len(members))
        entry = dict(
            name=f"{kind}:{members[0]}" if kind == "markov" else f"binding:{members[0]}",
            level="subgroup",
            structure=kind,
            members=members,
            internal_atoms=label_cols(internal),
            external_atoms=[
                dict(atom=_label(a, leaf_names), kind=kinds[a], states=sorted(set(atoms[a]) & mset))
                for a in external
            ],
            solver=how,
            suggested_block=dict(implicit=members),
        )
        closing = {closed_by[m] for m in members if m in closed_by}
        if vec is not None:
            w = _normalise(vec, members)
            entry["class"] = "exact" if not external else "conserved-internally"
            entry["weights"] = {k: str(v) for k, v in w.items()}
            entry["weights_numeric"] = {k: float(v.xreplace(params)) for k, v in w.items()}
            if kind == "binding":
                # a pool and its own buffers share a volume, so their weights should all be 1;
                # another weight means the bound and the free form are scaled differently
                odd = {k: v for k, v in entry["weights_numeric"].items() if abs(v - 1) > 1e-9}
                if odd:
                    entry["warning"] = (
                        f"weights other than 1 inside one compartment: {odd}; the binding terms "
                        "of these states are scaled differently from their pool (a volume ratio "
                        "or a factor in one equation but not the other?)"
                    )
        elif len(closing) == 1 and set(members) <= set(next(c for c in closures if c["name"] in closing)["states"]):
            c = next(c for c in closures if c["name"] in closing)
            entry["class"] = "closure"
            entry["closure"] = c["name"]
            entry["closure_total"] = c["total"]
            entry["note"] = (
                f"closed by the intermediate {c['name']} = {c['total']:g} - sum of the states: the "
                f"group conserves by construction; the implicit block takes the {len(members)} ODE states"
            )
            entry["weights"] = {m: "1" for m in members}
            entry["weights_numeric"] = {m: 1.0 for m in members}
            entry["bound_check"] = f"sum of {', '.join(members)} must stay <= {c['total']:g}"
        else:
            entry["class"] = "improper"
            entry.update(_diagnose(members, internal, M, params, atoms, leaf_names, kind))
        results.append(entry)

    # ---- networks: everything joined by shared atoms --------------------------------------
    net_edges = []
    for atom, row in atoms.items():
        if kinds[atom] == "shared":
            r = sorted(row)
            net_edges += [(r[0], x) for x in r[1:]]
    chain_members = {m for _, g in subgroups for m in g}
    eligible = {n for n in names if kind_of[n] == "conc" or n in chain_members}
    net_edges = []
    for atom, row in atoms.items():
        if kinds[atom] != "boundary":
            r = sorted(set(row) & eligible)
            net_edges += [(r[0], x) for x in r[1:]]
    nets = [g for g in _union_find(sorted(eligible), net_edges) if len(g) >= 2]
    for members in nets:
        mset = set(members)
        # rows outside the mass-carrying states (a gate that mirrors a buffer) do not count
        shared = [a for a, row in atoms.items() if kinds[a] != "boundary" and len(set(row) & eligible) >= 2 and (set(row) & eligible) <= mset]
        other = [a for a, row in atoms.items() if (set(row) & mset) and a not in shared]
        M = matrix(members, shared)
        vecs, how = nullspace(M, params, timeout=60)
        ions = sorted({i for i in (_ion(n) for n in members) if i})
        for vec in vecs:
            w = _normalise(vec, members)
            support = [m for m in members if m in w]
            if len(support) < 2:
                continue
            residual = []
            for a in other:
                total = sympy.cancel(sum(w.get(m, 0) * atoms[a].get(m, 0) for m in support))
                if total != 0:
                    residual.append(dict(atom=_label(a, leaf_names), kind=kinds[a], coefficient=str(total)))
            if not residual:
                cls = "exact"
            elif all(r["kind"] == "boundary" for r in residual):
                cls = "boundary"
            else:
                cls = "open"
            results.append(
                dict(
                    name=f"network:{'+'.join(ions) or members[0]}:{support[0]}",
                    level="network",
                    structure="exchange",
                    ions=ions,
                    members=support,
                    weights={k: str(v) for k, v in w.items()},
                    weights_numeric={k: float(v.xreplace(params)) for k, v in w.items()},
                    residual=residual,
                    solver=how,
                    **{"class": cls},
                )
            )

    closures = _closures(model, ex["state_syms"])

    # ---- numeric cross-check of the exact invariants ---------------------------------------
    if numeric_check:
        rng = np.random.default_rng(0)
        y0, p = model.y0, model.p0
        samples = [y0] + [y0 * rng.uniform(0.8, 1.2, len(y0)) for _ in range(3)]
        for r in results:
            if r.get("class") not in ("exact", "conserved-internally") or "weights_numeric" not in r:
                continue
            if r["level"] == "subgroup" and r["class"] in ("conserved-internally", "closure"):
                continue
            c = np.zeros(len(names))
            for k, v in r["weights_numeric"].items():
                c[names.index(k)] = v
            worst = 0.0
            for y in samples:
                f = model.rhs(0.0, y, p)
                worst = max(worst, abs(float(c @ f)) / max(float(np.max(np.abs(f))), 1e-300))
            r["numeric_check"] = dict(max_rel_cT_f=worst, ok=bool(worst < 1e-9))

    bookkeeping = charge_bookkeeping(model, ex, params)
    return dict(
        model=model.meta.id,
        affine_inlined=ex["affine_inlined"],
        n_atoms=len(atoms),
        row_factors=ex["row_factor"],
        buffered_totals=ex["buffered_total"],
        invariants=results,
        closures=closures,
        bookkeeping=bookkeeping,
    )


def _diagnose(members, internal, M, params, atoms, leaf_names, kind) -> dict:
    """For a sub-group that should conserve but does not: the intended weights and the atoms
    whose column sums break them, and whether removing one atom restores a null vector"""
    out = {}
    if kind == "markov":
        intended = {m: sympy.Integer(1) for m in members}
    else:
        intended = None
        # try dropping one atom at a time
        for j, a in enumerate(internal):
            Mj = M.copy()
            Mj.col_del(j)
            vecs, _ = nullspace(Mj, params, timeout=10)
            vec = _full_support(vecs, len(members))
            if vec is not None:
                intended = _normalise(vec, members)
                out["culprit_atom"] = _label(a, leaf_names)
                break
        if intended is None:
            intended = {m: sympy.Integer(1) for m in members}
    out["intended_weights"] = {k: str(v) for k, v in intended.items()}
    bad = []
    for j, a in enumerate(internal):
        col = {m: atoms[a].get(m, 0) for m in members}
        total = sympy.cancel(sum(intended.get(m, 0) * col[m] for m in members))
        if total != 0:
            bad.append(
                dict(
                    atom=_label(a, leaf_names),
                    column={m: str(c) for m, c in col.items() if c != 0},
                    weighted_sum=str(total),
                    numeric=float(total.xreplace(params)) if total.free_symbols else float(total),
                )
            )
    out["residual_atoms"] = bad
    return out


# ----------------------------------------------------------------------------------------
# charge / ion bookkeeping
# ----------------------------------------------------------------------------------------


def charge_bookkeeping(model: Model, ex: dict, params: dict) -> dict:
    atoms, kinds, leaf_names = ex["atoms"], ex["kinds"], ex["leaf_names"]
    v = model.meta.voltage
    names = model.names
    kind_of = dict(zip(names, model.kinds))
    conc = [n for n in names if kind_of[n] == "conc" and not _BOUND.search(n) and _ion(n)]
    currents = [a for a in ex["v_atoms"] if a != 1]
    if not currents:
        return dict(currents=[], compartments={}, note="no current atoms in dV/dt")

    a_v = {a: atoms[a][v] for a in currents}
    # kappa[n][current] = (coefficient of the current in d[n]/dt) / (its coefficient in dV/dt),
    # the current matched by its symbol so that a leftover factor (a buffering factor that was
    # not common to the whole equation) is kept and cancels in the ratios below
    kappa: dict[str, dict] = {n: {} for n in conc}
    for atom, row in atoms.items():
        if atom == 1:
            continue
        for cur in currents:
            if cur == atom or (cur.is_Symbol and cur in atom.free_symbols):
                left = sympy.cancel(atom / cur)
                for n in conc:
                    if n in row:
                        kappa[n][cur] = kappa[n].get(cur, 0) + sympy.cancel(row[n] * left / a_v[cur])
    # per compartment: the scale of one charge unit = the smallest |kappa| over its currents;
    # a kappa that is still state-dependent after cancelling is evaluated at the initial state
    y0 = model.y0
    state_vals = {model.ode.states[i].symbol: float(y0[model.index(model.ode.states[i].name)]) for i in range(len(model.ode.states))}

    def num(e):
        e = sympy.sympify(e).xreplace(params)
        if e.free_symbols:
            e = e.xreplace(state_vals)
        try:
            return float(e)
        except TypeError:
            return float("nan")

    for n in conc:
        kappa[n] = {cur: num(val) for cur, val in kappa[n].items()}
    scale = {}
    for n, ks in kappa.items():
        vals = [abs(x) for x in ks.values() if x == x and x != 0]
        scale[n] = min(vals) if vals else float("nan")
    rows = []
    for a in currents:
        label = _label(a, leaf_names)
        counts = {}
        total = 0.0
        integral = True
        for n in conc:
            if a in kappa[n] and scale[n] == scale[n]:
                q = kappa[n][a] / scale[n] * (-1 if _ion(n) == "Cl" else 1)
                counts[n] = round(q, 3)
                total += q
                if abs(q - round(q)) > 1e-6:
                    integral = False
        if not counts:
            status = "unaccounted"
        elif not integral:
            status = "miscounted (non-integer charge count)"
        elif abs(total - 1) < 1e-6:
            status = "ok"
        elif abs(total) < 1e-6:
            status = "charge-neutral (exchanger written without net charge?)"
        else:
            status = f"miscounted (charge count {total:.3g}, expected 1)"
        is_stim = "stim" in label.lower()
        rows.append(dict(current=label, counts=counts, charge_count=round(total, 6), status=status, stimulus=is_stim))
    # estimated drift of the stimulus into cytosolic K, if the stimulus is unaccounted
    note = ""
    k_state = model.meta.potassium
    stim = model.meta.stim
    if stim and k_state in scale and scale[k_state] == scale[k_state]:
        for r in rows:
            if r["stimulus"] and r["status"] == "unaccounted":
                r["est_drift_per_beat"] = dict(
                    state=k_state,
                    mM=abs(stim.amplitude) * stim.duration * scale[k_state],
                    note="the stimulus charge is not assigned to any ion; assigning it to K+ "
                    "(O'Hara-Rudy 2011 remedy) removes this drift",
                )
    return dict(currents=rows, compartments={n: scale[n] for n in conc}, note=note)


# ----------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------


def to_markdown(r: dict) -> str:
    out = [f"# Invariants and conservation checks of `{r['model']}`", ""]
    out.append(f"{r['n_atoms']} atoms after inlining {len(r['affine_inlined'])} affine intermediates.")
    if r.get("row_factors"):
        out += ["", "Equations with a common factor divided out (rapid buffering): "
                + "; ".join(f"`{k}`: {v}" for k, v in r["row_factors"].items())]
        if r.get("buffered_totals"):
            out += ["", "Their conserved quantity enters through the buffered total: "
                    + "; ".join(f"`{k}` -> `{v}`" for k, v in r["buffered_totals"].items())]
    subs = [x for x in r["invariants"] if x["level"] == "subgroup"]
    nets = [x for x in r["invariants"] if x["level"] == "network" and x["class"] != "open"]
    flags = [x for x in subs if x["class"] == "improper"]
    unacc = [c for c in r["bookkeeping"].get("currents", []) if c["status"] != "ok"]
    out += ["", "## Flags", ""]
    warned = [x for x in subs if x.get("warning")]
    if not flags and not unacc and not warned:
        out.append("- none: every chain and binding cluster conserves, every current is accounted for")
    for x in [x for x in subs if x.get("warning")]:
        out.append(f"- **{x['name']}** ({', '.join(x['members'])}): {x['warning']}")
    for x in flags:
        out.append(f"- **{x['name']}** ({x['structure']}, {', '.join(x['members'])}) does not conserve with any weights.")
        if x.get("culprit_atom"):
            out.append(f"  - removing the atom `{x['culprit_atom']}` restores conservation with weights {x['intended_weights']}")
        for b in x.get("residual_atoms", []):
            out.append(f"  - atom `{b['atom']}`: columns {b['column']} sum to `{b['weighted_sum']}` ({b['numeric']:.3g}) under weights {x['intended_weights']}")
    for c in unacc:
        line = f"- current `{c['current']}`: {c['status']}; charge counts {c['counts']}"
        if "est_drift_per_beat" in c:
            d = c["est_drift_per_beat"]
            line += f"; estimated drift {d['mM']:.3g} mM of {d['state']} per beat ({d['note']})"
        out.append(line)
    out += ["", "## Sub-groups (candidate implicit blocks)", ""]
    out.append(
        md_table(
            ["name", "structure", "class", "members", "weights", "external atoms", "numeric check"],
            [
                [
                    x["name"], x["structure"], x["class"], ", ".join(x["members"]),
                    ", ".join(f"{k}: {v}" for k, v in x.get("weights", {}).items()) or "-",
                    ", ".join(sorted({e["atom"] for e in x["external_atoms"]})) or "none",
                    ("ok" if x["numeric_check"]["ok"] else f"FAIL {x['numeric_check']['max_rel_cT_f']:.1e}") if "numeric_check" in x else "-",
                ]
                for x in subs
            ],
        )
    )
    out += ["", "## Networks (total-ion invariants)", ""]
    out.append(
        md_table(
            ["name", "class", "ions", "members", "weights", "residual (d/dt of the total)", "numeric check"],
            [
                [
                    x["name"], x["class"], "+".join(x["ions"]), ", ".join(x["members"]),
                    ", ".join(f"{k}: {v}" for k, v in x["weights"].items()),
                    "; ".join(f"{e['coefficient']} * {e['atom']}" for e in x["residual"]) or "0",
                    ("ok" if x["numeric_check"]["ok"] else f"FAIL {x['numeric_check']['max_rel_cT_f']:.1e}") if "numeric_check" in x else "-",
                ]
                for x in nets
            ],
        )
    )
    out += ["", "## Algebraic closures (a state written as const minus the others)", ""]
    out += [f"- `{c['name']} = {c['expr']}`" for c in r["closures"]] or ["- none"]
    out += ["", "## Charge bookkeeping of the currents in dV/dt", ""]
    out.append(
        md_table(
            ["current", "ion counts (charge units)", "total", "status"],
            [[c["current"], ", ".join(f"{k}: {v:g}" for k, v in c["counts"].items()) or "-", c["charge_count"], c["status"]] for c in r["bookkeeping"].get("currents", [])],
        )
    )
    return "\n".join(out)


def evaluate(invariants: list[dict], names: list[str], y: np.ndarray) -> dict[str, np.ndarray]:
    """``Q(t) = cᵀ y(t)`` for every invariant with numeric weights; `y` has one row per time"""
    out = {}
    for inv in invariants:
        w = inv.get("weights_numeric")
        if not w:
            continue
        c = np.zeros(len(names))
        for k, v in w.items():
            c[names.index(k)] = v
        out[inv["name"]] = np.asarray(y) @ c
    return out


def run(ws: Workspace, numeric_check: bool = True) -> dict:
    result = analyse(ws.model, numeric_check=numeric_check)
    ws.write_json("invariants.json", result)
    ws.write_md("invariants.md", to_markdown(result))
    return result
