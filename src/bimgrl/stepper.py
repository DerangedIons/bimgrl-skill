"""The generated C of one scheme configuration, compiled and loaded with ctypes

The driver (stepping loop with saves, timing with ``clock_gettime``) and the compile command are
those of blocked-rush-larsen ``block_rush_larsen/run.py`` @ 58af805. Step n is at t = n dt from
an integer counter, from t = 0.
"""

from __future__ import annotations

import ctypes
import hashlib
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

SCHEMES = ("FE", "RL", "GRL", "BIGRL", "BIMGRL", "BIMGRL-auto")
GOTRANX_SCHEME = {
    "FE": "explicit_euler",
    "RL": "hybrid_rush_larsen",
    "GRL": "generalized_rush_larsen",
    "BIGRL": "block_rush_larsen",
    "BIMGRL": "block_rush_larsen",
    "BIMGRL-auto": "block_rush_larsen",
}

DRIVER = r"""
#include <stdlib.h>
#include <time.h>

/* Take nsteps steps from u0 at t = 0 and write the state to out every
   save_every steps (out starts with u0). Returns -1, or the step after which
   a state was not finite, which is where it stops. */
long run_beat(int n, const double *u0, const double *p, double dt, long nsteps,
              long save_every, double *out) {
    double *cur = malloc(sizeof(double) * n);
    double *next = malloc(sizeof(double) * n);
    long fail = -1, isave = 0;
    for (int j = 0; j < n; j++) cur[j] = out[j] = u0[j];
    isave++;
    for (long it = 0; it < nsteps && fail < 0; it++) {
        SCHEME(cur, it * dt, dt, p, next);
        double *swap = cur; cur = next; next = swap;
        if ((it + 1) % save_every == 0) {
            for (int j = 0; j < n; j++) {
                out[isave * n + j] = cur[j];
                if (!(fabs(cur[j]) < 1e300)) fail = it;
            }
            isave++;
        }
    }
    free(cur);
    free(next);
    return fail;
}

/* The seconds it takes to take nsteps steps from u0, `beats` times over, and
   the state at the end of the last (which keeps the compiler from discarding the work) */
double time_beats(int n, const double *u0, const double *p, double dt, long nsteps,
                  long beats, double *last) {
    double *cur = malloc(sizeof(double) * n);
    double *next = malloc(sizeof(double) * n);
    struct timespec start, stop;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (long beat = 0; beat < beats; beat++) {
        for (int j = 0; j < n; j++) cur[j] = u0[j];
        for (long it = 0; it < nsteps; it++) {
            SCHEME(cur, it * dt, dt, p, next);
            double *swap = cur; cur = next; next = swap;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &stop);
    for (int j = 0; j < n; j++) last[j] = cur[j];
    free(cur);
    free(next);
    return (stop.tv_sec - start.tv_sec) + 1e-9 * (stop.tv_nsec - start.tv_nsec);
}
"""


def _pointer(array: np.ndarray) -> ctypes.c_void_p:
    return ctypes.c_void_p(array.ctypes.data)


def compiler() -> str:
    """``$CC``, else gcc, else cc"""
    for cc in (os.environ.get("CC"), "gcc", "cc", "clang"):
        if not cc:
            continue
        try:
            subprocess.run([cc, "--version"], capture_output=True, check=True)
            return cc
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue
    raise SystemExit("no C compiler found; set CC")


def _compile(cc: str, c_file: Path, library: Path) -> None:
    command = [cc, "-O3", "-march=native", "-shared", "-fPIC", "-o", str(library), str(c_file), "-lm"]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"{' '.join(command)} failed:\n{e.stderr[-3000:]}")


@dataclass(frozen=True)
class SchemeConfig:
    label: str
    gotranx_scheme: str
    kwargs: dict = field(hash=False)
    identity: str = ""
    blocks_file: str | None = None

    @property
    def key(self) -> str:
        digest = hashlib.sha1(f"{self.gotranx_scheme}\n{self.identity}".encode()).hexdigest()
        return f"{self.label}-{digest[:10]}"


def scheme_config(label: str, model, blocks_file: Path | None = None, suggested_file: Path | None = None) -> SchemeConfig | None:
    """The configuration `label`; None if it does not apply (BIGRL without implicit blocks)"""
    from gotranx.blocks import load_blocks

    scheme = GOTRANX_SCHEME[label]
    if label in ("FE", "GRL"):
        return SchemeConfig(label, scheme, {}, "")
    if label == "RL":
        gates = model.meta.rl_gates or ws_gates(model)
        return SchemeConfig(label, scheme, {"stiff_states": list(gates)}, repr(list(gates)))
    path = suggested_file if label == "BIMGRL-auto" else blocks_file
    if path is None or not Path(path).exists():
        raise SystemExit(f"{label} needs a blocks file ({path})")
    blocks = load_blocks(path)
    if label == "BIGRL":
        if not blocks.implicit:
            return None
        kwargs = dict(implicit_blocks=blocks.implicit, substep_blocks=())
    else:
        kwargs = dict(implicit_blocks=blocks.implicit, substep_blocks=blocks.substep)
    identity = repr((tuple(kwargs["implicit_blocks"]), tuple(kwargs["substep_blocks"])))
    return SchemeConfig(label, scheme, kwargs, identity, str(path))


def ws_gates(model) -> list[str]:
    from gotranx.block_detection import state_dependencies
    from gotranx.linearization import diagonal_jacobian

    from .structure import hh_gates

    return hh_gates(model.ode, state_dependencies(model.ode), diagonal_jacobian(model.ode))


class CStepper:
    def __init__(self, ode, config: SchemeConfig, build_dir: Path, cc: str | None = None) -> None:
        from gotranx.codegen import CCodeGenerator
        from gotranx.codegen.c import Format
        from gotranx.schemes import get_scheme

        cc = cc or compiler()
        codegen = CCodeGenerator(ode, format=Format.none)
        source = "\n".join(
            [
                codegen.imports(),
                codegen.parameter_index(),
                codegen.state_index(),
                codegen.initial_parameter_values(),
                codegen.initial_state_values(),
                codegen.scheme(get_scheme(config.gotranx_scheme), **config.kwargs),
                f"#define SCHEME {config.gotranx_scheme}",
                DRIVER,
            ]
        )
        build_dir = Path(build_dir)
        build_dir.mkdir(parents=True, exist_ok=True)
        c_file = build_dir / f"{config.key}.c"
        library = c_file.with_suffix(".so")
        if not library.exists() or not c_file.exists() or c_file.read_text() != source:
            tmp_c = build_dir / f"{config.key}.{os.getpid()}.c"
            tmp_lib = tmp_c.with_suffix(".so")
            tmp_c.write_text(source)
            _compile(cc, tmp_c, tmp_lib)
            os.replace(tmp_lib, library)
            os.replace(tmp_c, c_file)
        lib = ctypes.CDLL(str(library))
        void_p, double, long = ctypes.c_void_p, ctypes.c_double, ctypes.c_long
        lib.run_beat.restype = long
        lib.run_beat.argtypes = [ctypes.c_int, void_p, void_p, double, long, long, void_p]
        lib.time_beats.restype = double
        lib.time_beats.argtypes = [ctypes.c_int, void_p, void_p, double, long, long, void_p]
        lib.init_state_values.argtypes = [void_p]
        lib.init_parameter_values.argtypes = [void_p]
        lib.state_index.argtypes = [ctypes.c_char_p]
        lib.state_index.restype = ctypes.c_int
        lib.parameter_index.argtypes = [ctypes.c_char_p]
        lib.parameter_index.restype = ctypes.c_int
        self.lib = lib
        self.config = config
        self.num_states = ode.num_states
        self.u0 = np.zeros(ode.num_states)
        self.parameters = np.zeros(ode.num_parameters)
        lib.init_state_values(_pointer(self.u0))
        lib.init_parameter_values(_pointer(self.parameters))
        self.state_names = sorted((s.name for s in ode.states), key=self.state_index)
        self.c_file = c_file

    def state_index(self, name: str) -> int:
        index = self.lib.state_index(name.encode())
        if index < 0:
            raise KeyError(f"{name!r} is not a state of the model")
        return index

    def run(self, dt: float, n_steps: int, save_every: int, parameters=None, u0=None):
        p = self.parameters if parameters is None else np.ascontiguousarray(parameters, dtype=float)
        y0 = self.u0 if u0 is None else np.ascontiguousarray(u0, dtype=float)
        num_saves = n_steps // save_every + 1
        out = np.full((num_saves, self.num_states), np.nan)
        fail = self.lib.run_beat(self.num_states, _pointer(y0), _pointer(p), dt, n_steps, save_every, _pointer(out))
        return np.arange(num_saves) * save_every * dt, out, int(fail)

    def time(self, dt: float, n_steps: int, repeats: int = 5, parameters=None, u0=None) -> float:
        p = self.parameters if parameters is None else np.ascontiguousarray(parameters, dtype=float)
        y0 = self.u0 if u0 is None else np.ascontiguousarray(u0, dtype=float)
        last = np.zeros(self.num_states)
        return min(self.lib.time_beats(self.num_states, _pointer(y0), _pointer(p), dt, n_steps, 1, _pointer(last)) for _ in range(repeats))
