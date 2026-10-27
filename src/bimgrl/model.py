"""A model in the canonical ``.ode`` form, loaded through gotranx

Everything downstream that is symbolic (structure, invariants, codegen) needs the
:class:`gotranx.ode.ODE`; everything numeric (reference, stiffness, sweep) needs the generated
numpy right-hand side. Both live here, built once per process.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .meta import ModelMeta
from .protocol import Protocol
from .vendored.metrics import bounds as _bounds
from .vendored.reference import Reference, compute_reference, numpy_model
from .vendored.stimulus import PARAMETERS as STIM_PARAMETERS

KINDS = ("voltage", "gate", "markov", "conc", "other")


def require_gotranx():
    try:
        import gotranx  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            "gotranx is not installed in the skill environment; run `bimgrl setup` "
            "(see config.toml for the source it installs from)"
        ) from e


def state_kinds(ode, meta: ModelMeta, order) -> list[str]:
    """voltage | gate | conc | other for every state, from the units, with meta.kinds overriding"""
    by_name = {state.name: state for state in ode.states}
    kinds = []
    for name in order:
        state = by_name[name]
        unit = (state.unit_str or "").strip()
        if name in meta.kinds:
            kinds.append(meta.kinds[name])
        elif name == meta.voltage:
            kinds.append("voltage")
        elif name in meta.unphysical_exempt:
            kinds.append("other")
        elif unit in ("", "1", "dimensionless"):
            kinds.append("gate")
        elif "M" in unit or "mol" in unit:
            kinds.append("conc")
        else:
            kinds.append("other")
    return kinds


@dataclass
class Model:
    """One ``.ode`` model with its metadata, ready for analysis"""

    meta: ModelMeta
    ode_path: Path
    ode: object = None
    ns: dict = None
    #: state names in the order of the generated code (the columns of every trace)
    names: list[str] = None
    kinds: list[str] = None
    _references: dict = field(default_factory=dict)
    _jacobian: object = None

    @classmethod
    def load(cls, ode_path: Path, meta: ModelMeta) -> Model:
        require_gotranx()
        import gotranx

        ode = gotranx.load_ode(Path(ode_path), remove_singularities=meta.remove_singularities)
        ns = numpy_model(ode)
        index = ns["state_index"]
        names = sorted((s.name for s in ode.states), key=index)
        model = cls(meta=meta, ode_path=Path(ode_path), ode=ode, ns=ns, names=names)
        model.kinds = state_kinds(ode, meta, names)
        return model

    # --- names and indices -------------------------------------------------------------

    def index(self, name: str | None) -> int | None:
        return None if name is None else self.ns["state_index"](name)

    @property
    def v_index(self) -> int:
        return self.index(self.meta.voltage)

    @property
    def parameter_names(self) -> list[str]:
        index = self.ns["parameter_index"]
        return sorted((p.name for p in self.ode.parameters), key=index)

    @property
    def has_stim_parameters(self) -> bool:
        names = {q.name for q in self.ode.parameters}
        return all(p in names for p in STIM_PARAMETERS)

    @property
    def y0(self) -> np.ndarray:
        return np.array(self.ns["init_state_values"](), dtype=float)

    @property
    def p0(self) -> np.ndarray:
        return np.array(self.ns["init_parameter_values"](), dtype=float)

    def parameters(self, protocol: Protocol, base: np.ndarray | None = None) -> np.ndarray:
        """The parameters with the stimulus of `protocol`"""
        p = np.array(self.p0 if base is None else base, dtype=float)
        if not self.has_stim_parameters:
            if protocol.stim is not None:
                raise ValueError(
                    f"{self.meta.id} has no stim_* parameters; the stimulus must be rewritten "
                    "to the half-open form first (bimgrl convert / references/contracts.md)"
                )
            return p
        index = self.ns["parameter_index"]
        if protocol.stim is None:
            p[index("stim_amplitude")] = 0.0
        else:
            for name in ("start", "duration", "period", "amplitude"):
                p[index(f"stim_{name}")] = getattr(protocol.stim, name)
        return p

    # --- numerics ---------------------------------------------------------------------

    def rhs(self, t: float, y: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.asarray(self.ns["rhs"](t, y, p), dtype=float)

    def jacobian(self):
        """``jacobian(t, y, p)``: gotranx's exact symbolic Jacobian compiled to numpy, or central
        differences when meta.jacobian == "numeric" """
        if self._jacobian is None:
            if self.meta.jacobian == "numeric":
                self._jacobian = numeric_jacobian(self.rhs)
            else:
                from gotranx.block_detection import jacobian_function

                self._jacobian = jacobian_function(self.ode)
        return self._jacobian

    def reference(
        self, protocol: Protocol, rtol: float = 1e-10, atol: float = 1e-12, y0=None
    ) -> Reference:
        key = (protocol, rtol, atol, None if y0 is None else tuple(np.round(y0, 15)))
        if key not in self._references:
            self._references[key] = compute_reference(
                self.ns,
                self.parameters(protocol),
                self.y0 if y0 is None else y0,
                0.0,
                protocol.t1,
                protocol.stim,
                rtol=rtol,
                atol=atol,
            )
        return self._references[key]

    def bounds(self, protocol: Protocol, sample_ms: float = 0.1):
        ref = self.reference(protocol)
        y = ref(np.linspace(0.0, protocol.t1, round(protocol.t1 / sample_ms) + 1))
        return _bounds(self.kinds, y.min(axis=1), y.max(axis=1))


def numeric_jacobian(rhs, rel: float = 1e-6, floor: float = 1e-9):
    """``jacobian(t, y, p)`` by central differences of ``rhs(t, y, p)``"""

    def jacobian(t, y, p):
        y = np.asarray(y, dtype=float)
        J = np.empty((len(y), len(y)))
        for j in range(len(y)):
            h = rel * max(abs(y[j]), floor)
            up, down = y.copy(), y.copy()
            up[j] += h
            down[j] -= h
            J[:, j] = (np.asarray(rhs(t, up, p)) - np.asarray(rhs(t, down, p))) / (2 * h)
        return J

    return jacobian
