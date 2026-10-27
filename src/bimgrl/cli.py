"""``bimgrl <command> <model_dir> [options]``: the deterministic steps of the skill"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _ws(args):
    from .workspace import Workspace

    return Workspace(Path(args.model_dir))


def cmd_setup(args):
    from . import setup

    setup.run(Path(args.skill_dir) if args.skill_dir else None)


def cmd_intake(args):
    from . import intake

    ws = _ws(args)
    meta = intake.run(ws, voltage=args.voltage, force=args.force)
    print((ws.out / "intake.md").read_text())
    print(f"wrote {ws.model_dir / 'meta.toml'} (status {meta.status})")


def cmd_structure(args):
    from . import structure

    ws = _ws(args)
    structure.run(ws)
    print((ws.out / "structure.md").read_text())


def cmd_invariants(args):
    ws = _ws(args)
    if args.adapter:
        from . import native
        from .adapters import open_adapter

        groups = [g.split(",") for g in args.group] if args.group else None
        native.invariants_numeric(ws, open_adapter(args.adapter, ws), Path(args.trace), groups)
        return
    from . import invariants

    invariants.run(ws, numeric_check=not args.no_numeric)
    print((ws.out / "invariants.md").read_text())


def cmd_reference(args):
    from . import reference

    ws = _ws(args)
    reference.run(ws, protocols_=args.protocol.split(","), beats=args.beats, rtol=args.rtol, atol=args.atol)


def cmd_stiffness(args):
    ws = _ws(args)
    if args.adapter:
        from . import native
        from .adapters import open_adapter

        native.stiffness(ws, open_adapter(args.adapter, ws), Path(args.trace), args.threshold)
        return
    from . import stiffness

    stiffness.run(ws, protocol=args.protocol, beats=args.beats, every_ms=args.every, threshold=args.threshold)
    print((ws.out / "stiffness.md").read_text())


def cmd_convert(args):
    from . import convert
    from .meta import Stim

    ws = _ws(args)
    stim = Stim(args.stim_start, args.stim_duration, args.stim_period, args.stim_amplitude) if args.stim_variable else None
    convert.convert(ws, Path(args.cellml), args.voltage, args.stim_variable, stim, args.id)


def cmd_fidelity(args):
    from . import fidelity

    fidelity.run(_ws(args), Path(args.samples) if args.samples else None, Path(args.params) if args.params else None,
                 Path(args.trace) if args.trace else None, Path(args.map) if args.map else None, args.rhs_tol)


def cmd_sweep(args):
    from . import sweep

    ws = _ws(args)
    dts = [float(x) for x in args.dt.split(",")] if args.dt else None
    labels = args.scheme.split(",")
    sweep.sweep(ws, labels, args.protocol, args.beats, dts, args.budget, args.blocks, args.suggested,
                stability=args.stability, save_traces=args.save_traces, force=args.force)
    if args.time:
        sweep.time_rows(ws, labels, args.protocol, blocks_file=args.blocks)


def cmd_time(args):
    from . import sweep

    sweep.time_rows(_ws(args), args.scheme.split(","), args.protocol, beats_timed=args.beats, repeats=args.repeats, blocks_file=args.blocks)


def cmd_propose(args):
    from . import propose

    ws = _ws(args)
    propose.run(ws, dt=args.dt, stiff_threshold=args.threshold, nest_v_markov=args.nest_v_markov, write=not args.no_write)
    print((ws.out / "propose.md").read_text())


def cmd_codegen(args):
    from . import codegen

    codegen.run(_ws(args), args.backend, args.scheme, args.out, args.blocks, args.extra or [])


def cmd_validate(args):
    from . import validate_cli

    validate_cli.run(_ws(args), Path(args.test), Path(args.reference) if args.reference else None, args.level, args.name)


def cmd_log(args):
    from . import log

    ws = _ws(args)
    if args.section:
        print(f"wrote {log.add_section_line(ws, args.section, args.line)}")
        return
    path = log.append_round(ws, args.round, args.change, args.result, args.decision, args.blocks)
    if args.freeze:
        log.set_round_header(ws, args.round, args.blocks)
    print(f"wrote {path}")


def cmd_report(args):
    from . import report

    report.run(_ws(args), protocol=args.protocol, baseline=args.baseline)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="bimgrl", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("setup", help="install gotranx into the skill environment from config.toml")
    p.add_argument("--skill-dir", default=None)
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("intake", help="draft meta.toml from model.ode and report what to confirm")
    p.add_argument("model_dir")
    p.add_argument("--voltage", default=None, help="name of the membrane potential state")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_intake)

    p = sub.add_parser("structure", help="dependency graph, SCCs, update types, linear subsystems")
    p.add_argument("model_dir")
    p.set_defaults(func=cmd_structure)

    p = sub.add_parser("invariants", help="linear invariants, conservation-law checks, ion bookkeeping")
    p.add_argument("model_dir")
    p.add_argument("--no-numeric", action="store_true", help="skip the numeric cross-check")
    p.add_argument("--adapter", default=None, help="native path: ode | cmd:<command> | so:<lib>:<meta.json>")
    p.add_argument("--trace", default=None, help="native path: a trace of the model for the checkpoints")
    p.add_argument("--group", action="append", help="native path: comma list of states to test together (repeatable)")
    p.set_defaults(func=cmd_invariants)

    p = sub.add_parser("reference", help="LSODA reference traces for protocols A and B")
    p.add_argument("model_dir")
    p.add_argument("--protocol", default="A,B")
    p.add_argument("--beats", type=int, default=10)
    p.add_argument("--rtol", type=float, default=1e-10)
    p.add_argument("--atol", type=float, default=1e-12)
    p.set_defaults(func=cmd_reference)

    p = sub.add_parser("stiffness", help="eigenmodes, participation factors, owners, FE and GRL limits")
    p.add_argument("model_dir")
    p.add_argument("--protocol", default="B")
    p.add_argument("--beats", type=int, default=2)
    p.add_argument("--every", type=float, default=1.0, help="checkpoint spacing in ms")
    p.add_argument("--threshold", type=float, default=-20.0, help="Re lambda below which a mode is stiff (/ms)")
    p.add_argument("--adapter", default=None, help="native path: ode | cmd:<command> | so:<lib>:<meta.json>")
    p.add_argument("--trace", default=None, help="native path: a trace of the model for the checkpoints")
    p.set_defaults(func=cmd_stiffness)

    p = sub.add_parser("convert", help="CellML -> model.ode (ms, mV, half-open stimulus) and a draft meta.toml")
    p.add_argument("model_dir")
    p.add_argument("--cellml", required=True)
    p.add_argument("--voltage", required=True, help="qualified myokit name of V, e.g. membrane.V")
    p.add_argument("--stim-variable", default=None, help="qualified name of the stimulus current, e.g. membrane.i_Stim")
    p.add_argument("--stim-start", type=float, default=0.0)
    p.add_argument("--stim-duration", type=float, default=1.0)
    p.add_argument("--stim-period", type=float, default=1000.0)
    p.add_argument("--stim-amplitude", type=float, default=-52.0, help="value of the current when on (model units)")
    p.add_argument("--id", default=None)
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("fidelity", help="prove a translated model.ode against the original code")
    p.add_argument("model_dir")
    p.add_argument("--samples", default=None, help="rhs_samples.csv: t, y_<state>..., dy_<state>... from the original code")
    p.add_argument("--params", default=None, help="params.json from the original code")
    p.add_argument("--trace", default=None, help="a trace of the original code (trace contract)")
    p.add_argument("--map", default=None, help="translation_map.json")
    p.add_argument("--rhs-tol", type=float, default=1e-9)
    p.set_defaults(func=cmd_fidelity)

    p = sub.add_parser("sweep", help="run schemes over the dt grid against the reference, coarsest first")
    p.add_argument("model_dir")
    p.add_argument("--scheme", default="GRL,BIMGRL", help="comma list of FE,RL,GRL,BIGRL,BIMGRL,BIMGRL-auto")
    p.add_argument("--protocol", default="B")
    p.add_argument("--beats", type=int, default=10)
    p.add_argument("--dt", default=None, help="comma list of time steps (default: the 1-2-5 grid, descending until two pass tight)")
    p.add_argument("--budget", type=float, default=600.0, help="seconds per run above which a run is skipped")
    p.add_argument("--blocks", default=None, help="blocks.toml (default <model_dir>/blocks.toml)")
    p.add_argument("--suggested", default=None, help="blocks file for BIMGRL-auto")
    p.add_argument("--stability", action="store_true", help="bisect the stability limit of the explicit schemes")
    p.add_argument("--save-traces", action="store_true")
    p.add_argument("--time", action="store_true", help="also time every ok run (serial, best of 5)")
    p.add_argument("--force", action="store_true", help="rerun rows that exist")
    p.set_defaults(func=cmd_sweep)

    p = sub.add_parser("time", help="time the ok runs of the sweep (one beat, best of N)")
    p.add_argument("model_dir")
    p.add_argument("--scheme", default="GRL,BIMGRL")
    p.add_argument("--protocol", default="B")
    p.add_argument("--beats", type=int, default=1)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--blocks", default=None)
    p.set_defaults(func=cmd_time)

    p = sub.add_parser("propose", help="propose blocks.toml by the paper's rules from structure, invariants, stiffness")
    p.add_argument("model_dir")
    p.add_argument("--dt", type=float, default=0.1, help="target time step for the gotranx suggest-blocks comparison")
    p.add_argument("--threshold", type=float, default=-20.0)
    p.add_argument("--nest-v-markov", action="store_true", help="nest voltage-dependent Markov blocks in the V block")
    p.add_argument("--no-write", action="store_true", help="do not write <model_dir>/blocks.toml")
    p.set_defaults(func=cmd_propose)

    p = sub.add_parser("codegen", help="generate scheme code with gotranx (c, python, julia, cytozoo)")
    p.add_argument("model_dir")
    p.add_argument("--backend", default="c", choices=["c", "python", "julia", "cytozoo"])
    p.add_argument("--scheme", default="BIMGRL", choices=["FE", "RL", "GRL", "BIGRL", "BIMGRL"])
    p.add_argument("--blocks", default=None)
    p.add_argument("-o", "--out", default=None)
    p.add_argument("extra", nargs="*", help="extra gotranx flags after --")
    p.set_defaults(func=cmd_codegen)

    p = sub.add_parser("validate", help="validate a trace file (any language) against a reference trace")
    p.add_argument("model_dir")
    p.add_argument("--test", required=True, help="trace .csv with its .json run meta")
    p.add_argument("--reference", default=None, help="reference trace .csv (default bimgrl/reference_<protocol>.csv)")
    p.add_argument("--level", default="moderate", choices=["moderate", "tight"])
    p.add_argument("--name", default=None)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("log", help="append a tuning round (or a section line) to NOTES.md")
    p.add_argument("model_dir")
    p.add_argument("--round", type=int, default=1)
    p.add_argument("--change", default="", help="what changed in the blocks (membership / tolerances)")
    p.add_argument("--result", default="", help="the numbers: status and errors at the time steps tried")
    p.add_argument("--decision", default="", help="adopted / rejected, and why")
    p.add_argument("--blocks", default="blocks.toml")
    p.add_argument("--freeze", action="store_true", help="also stamp the round into the blocks.toml header")
    p.add_argument("--section", default=None, help="instead: add --line to this NOTES.md section")
    p.add_argument("--line", default="")
    p.set_defaults(func=cmd_log)

    p = sub.add_parser("report", help="cheapest passing runs, speed-ups, Table-S1 row, error-vs-cost figure")
    p.add_argument("model_dir")
    p.add_argument("--protocol", default="B")
    p.add_argument("--baseline", default="GRL")
    p.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    try:
        args.func(args)
    except SystemExit:
        raise
    except KeyboardInterrupt:
        sys.exit(130)


if __name__ == "__main__":
    main()
