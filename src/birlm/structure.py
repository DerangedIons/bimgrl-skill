"""The structure of the model from the Jacobian sparsity: the dependency graph, strongly
connected components, the cone of influence of V, dead and observer states, duplicate
dependency rows, the structurally linear subsystems, and how each scheme advances each state

This is the ODEGraph v1 slice on the symbolic Jacobian: an edge ``x_j -> x_i`` for every
structurally non-zero ``df_i/dx_j``.
"""

from __future__ import annotations

from collections import defaultdict

import sympy

from .model import Model
from .workspace import Workspace, md_table


def sparsity(jac: dict) -> dict[str, set[str]]:
    """``deps[row] = {columns}``: the states that f_row depends on"""
    deps: dict[str, set[str]] = defaultdict(set)
    for (row, col), _ in jac.items():
        deps[row].add(col)
    return deps


def tarjan_scc(names: list[str], deps: dict[str, set[str]]) -> list[list[str]]:
    """Strongly connected components (edges col -> row), in reverse topological order of the
    condensation, i.e. a component comes after everything it depends on... iterative Tarjan"""
    index = {}
    low = {}
    on_stack = set()
    stack: list[str] = []
    out: list[list[str]] = []
    counter = 0
    # successors of x: the rows whose derivative depends on x
    succ: dict[str, list[str]] = defaultdict(list)
    for row, cols in deps.items():
        for col in cols:
            succ[col].append(row)
    for root in names:
        if root in index:
            continue
        work = [(root, iter(succ[root]))]
        index[root] = low[root] = counter
        counter += 1
        stack.append(root)
        on_stack.add(root)
        while work:
            node, it = work[-1]
            advanced = False
            for nxt in it:
                if nxt not in index:
                    index[nxt] = low[nxt] = counter
                    counter += 1
                    stack.append(nxt)
                    on_stack.add(nxt)
                    work.append((nxt, iter(succ[nxt])))
                    advanced = True
                    break
                if nxt in on_stack:
                    low[node] = min(low[node], index[nxt])
            if advanced:
                continue
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[node])
            if low[node] == index[node]:
                comp = []
                while True:
                    x = stack.pop()
                    on_stack.discard(x)
                    comp.append(x)
                    if x == node:
                        break
                out.append(sorted(comp))
    return out


def reachable(start: set[str], graph: dict[str, set[str]]) -> set[str]:
    seen = set(start)
    todo = list(start)
    while todo:
        x = todo.pop()
        for y in graph.get(x, ()):
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def hh_gates(ode, deps_by_symbol, diagonal) -> list[str]:
    """The Hodgkin-Huxley gates: dimensionless states whose derivative is affine in themselves
    and depends on no other dimensionless state (which leaves out Markov chains)"""
    dimensionless = {
        s.name for s in ode.states if (s.unit_str or "").strip() in ("", "1", "dimensionless")
    }
    derivative = {d.state.name: d for d in ode.state_derivatives}
    gates = []
    for state in ode.states:
        name = state.name
        if name not in dimensionless:
            continue
        rhs_states = deps_by_symbol[derivative[name].symbol]
        diag_states = frozenset().union(
            *(deps_by_symbol.get(s, frozenset()) for s in diagonal[name].free_symbols)
        )
        if name in diag_states or diagonal[name].is_zero:
            continue
        if (rhs_states - {name}) & dimensionless:
            continue
        gates.append(name)
    return gates


def analyse(model: Model) -> dict:
    from gotranx.block_detection import linear_subsystems, state_dependencies
    from gotranx.linearization import diagonal_jacobian, jacobian

    ode = model.ode
    names = model.names
    jac = jacobian(ode)
    diag = diagonal_jacobian(ode)
    deps = sparsity(jac)
    deps_by_symbol = state_dependencies(ode)
    v = model.meta.voltage

    # the derivative expressions, for dead states
    dexpr = {d.state.name: d.expr for d in ode.state_derivatives}
    dead = [n for n in names if dexpr[n] == 0 or dexpr[n].is_zero]

    succ: dict[str, set[str]] = defaultdict(set)
    for row, cols in deps.items():
        for col in cols:
            if col != row:
                succ[col].add(row)
    observers = [n for n in names if not succ[n] and n not in dead]

    sccs = tarjan_scc(names, deps)
    scc_of = {n: i for i, comp in enumerate(sccs) for n in comp}
    feeds_v = sorted(reachable({v}, {r: {c for c in cs if c != r} for r, cs in deps.items()}) - {v})
    driven_by_v = sorted(reachable({v}, succ) - {v})
    not_feeding_v = [n for n in names if n != v and n not in feeds_v]

    # identical dependency rows (excluding self): symmetry candidates
    rows: dict[tuple, list[str]] = defaultdict(list)
    for n in names:
        rows[tuple(sorted(deps.get(n, set()) - {n}))].append(n)
    duplicates = [sorted(g) for g in rows.values() if len(g) > 1]

    linear = [list(g) for g in linear_subsystems(ode, jac)]
    gates = hh_gates(ode, deps_by_symbol, diag)

    # how GRL (gotranx: exponential wherever the diagonal is structurally non-zero) and
    # RL+FE (exponential only for the HH gates) advance each state
    update_grl = {n: ("FE" if diag[n].is_zero else "exp") for n in names}
    update_rl = {n: ("exp" if n in gates else "FE") for n in names}

    states = []
    by_name = {s.name: s for s in ode.states}
    for n in names:
        states.append(
            dict(
                name=n,
                unit=(by_name[n].unit_str or "").strip(),
                kind=model.kinds[model.index(n)],
                init=float(by_name[n].value),
                n_in=len(deps.get(n, set()) - {n}),
                n_out=len(succ[n]),
                scc=scc_of[n],
                diag_zero=bool(diag[n].is_zero),
                update_grl=update_grl[n],
                update_rl=update_rl[n],
                feeds_v=n in feeds_v or n == v,
            )
        )
    return dict(
        model=model.meta.id,
        n_states=len(names),
        n_edges=len(jac),
        voltage=v,
        states=states,
        deps={n: sorted(deps.get(n, set())) for n in names},
        sccs=[dict(index=i, states=c, size=len(c)) for i, c in enumerate(sccs)],
        n_nontrivial_sccs=sum(1 for c in sccs if len(c) > 1),
        feeds_v=feeds_v,
        not_feeding_v=not_feeding_v,
        driven_by_v=driven_by_v,
        dead_states=dead,
        observer_states=observers,
        duplicate_dependency_rows=duplicates,
        linear_subsystems=linear,
        hh_gates=gates,
    )


def to_dot(result: dict, blocks: dict | None = None) -> str:
    """Graphviz DOT of the dependency graph, one node per state, optionally boxed by SCC"""
    lines = ["digraph model {", "  rankdir=LR; node [shape=box, fontsize=10];"]
    color = {"voltage": "gold", "gate": "lightblue", "markov": "plum", "conc": "palegreen"}
    for s in result["states"]:
        lines.append(f'  "{s["name"]}" [style=filled, fillcolor={color.get(s["kind"], "white")}];')
    for comp in result["sccs"]:
        if comp["size"] > 1:
            inner = " ".join(f'"{n}";' for n in comp["states"])
            lines.append(f'  subgraph cluster_scc{comp["index"]} {{ label="SCC {comp["index"]}"; {inner} }}')
    for row, cols in result["deps"].items():
        for col in cols:
            if col != row:
                lines.append(f'  "{col}" -> "{row}";')
    lines.append("}")
    return "\n".join(lines)


def to_markdown(r: dict) -> str:
    out = [f"# Structure of `{r['model']}`", ""]
    out.append(
        f"{r['n_states']} states, {r['n_edges']} structural Jacobian entries, "
        f"{r['n_nontrivial_sccs']} non-trivial strongly connected components. "
        f"Membrane potential: `{r['voltage']}`."
    )
    out += ["", "## States", ""]
    out.append(
        md_table(
            ["state", "unit", "kind", "init", "in", "out", "SCC", "GRL update", "RL+FE update", "feeds V"],
            [
                [s["name"], s["unit"], s["kind"], s["init"], s["n_in"], s["n_out"], s["scc"],
                 s["update_grl"], s["update_rl"], "yes" if s["feeds_v"] else "no"]
                for s in r["states"]
            ],
        )
    )
    out += ["", "## Strongly connected components (size > 1)", ""]
    big = [c for c in r["sccs"] if c["size"] > 1]
    out += [f"- SCC {c['index']} ({c['size']}): {', '.join(c['states'])}" for c in big] or ["- none"]
    out += ["", "## Structurally linear subsystems (Markov chains, linked buffers)", ""]
    out += [f"- {', '.join(g)}" for g in r["linear_subsystems"]] or ["- none"]
    out += ["", "## Hodgkin-Huxley gates", "", "- " + (", ".join(r["hh_gates"]) or "none")]
    out += ["", "## States that do not feed V (observers of the electrophysiology)", ""]
    out += ["- " + (", ".join(r["not_feeding_v"]) or "none")]
    out += ["", "## Dead states (f ≡ 0)", "", "- " + (", ".join(r["dead_states"]) or "none")]
    out += ["", "## Observer states (nothing depends on them)", "", "- " + (", ".join(r["observer_states"]) or "none")]
    out += ["", "## Identical dependency rows (symmetry candidates)", ""]
    out += [f"- {', '.join(g)}" for g in r["duplicate_dependency_rows"]] or ["- none"]
    return "\n".join(out)


def run(ws: Workspace) -> dict:
    result = analyse(ws.model)
    ws.write_json("structure.json", result)
    ws.write_md("structure.md", to_markdown(result))
    (ws.out / "structure.dot").write_text(to_dot(result))
    return result
