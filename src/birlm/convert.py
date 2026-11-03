"""CellML -> model.ode through myokit and gotranx, in ms and mV, with the half-open stimulus
(vendored from the survey's convert.py @ 58af805, without the per-model hook machinery)"""

from __future__ import annotations

from pathlib import Path

from .meta import ModelMeta, Stim, dumps_meta, save_meta
from .vendored.stimulus import rewrite_stimulus
from .workspace import Workspace


def load_cellml(path: Path):
    import myokit.formats.cellml

    return myokit.formats.cellml.CellMLImporter().model(str(path))


def merge_duplicate_connections(text: str) -> tuple[str, list[str]]:
    from lxml import etree

    root = etree.fromstring(text.encode())
    namespace = root.nsmap.get(None) or root.nsmap.get("cellml")
    tag = (lambda name: f"{{{namespace}}}{name}") if namespace else (lambda name: name)
    first: dict = {}
    edits = []
    for connection in root.findall(tag("connection")):
        components = connection.find(tag("map_components"))
        c1, c2 = components.get("component_1"), components.get("component_2")
        pair = frozenset((c1, c2))
        if pair not in first:
            first[pair] = (connection, c1)
            continue
        target, target_c1 = first[pair]
        for mapping in connection.findall(tag("map_variables")):
            if c1 != target_c1:
                v1, v2 = mapping.get("variable_1"), mapping.get("variable_2")
                mapping.set("variable_1", v2)
                mapping.set("variable_2", v1)
            target.append(mapping)
        connection.getparent().remove(connection)
        edits.append(f"XML repair: a second <connection> of {c1} and {c2} merged into the first (CellML 1.0, 3.4.5.4)")
    if not edits:
        return text, []
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8").decode(), edits


def _find(model, name: str):
    if "." in name:
        return model.get(name)
    model.create_unique_names()
    matches = [v for v in model.variables(deep=True) if v.uname() == name]
    if len(matches) != 1:
        raise KeyError(f"{name!r} is not a unique variable name of the model")
    return matches[0]


def normalize_units(model, voltage: str) -> list[str]:
    import myokit

    edits = []
    for var, unit, what in ((model.time(), "ms", "time"), (_find(model, voltage), "mV", "membrane potential")):
        new = myokit.Unit.parse_simple(unit)
        old = var.unit()
        if old is None:
            raise ValueError(f"The {what} {var.qname()} has no unit, convert it by hand")
        if old == new or old == myokit.units.dimensionless:
            continue
        var.convert_unit(new)
        edits.append(f"{what} {var.qname()} converted from {old} to {unit} (myokit convert_unit)")
    return edits


def inline_derivative_references(model) -> list[str]:
    import myokit

    edits = []
    for var in model.variables(deep=True):
        if var.is_state() or var.rhs() is None:
            continue
        derivatives = [ref for ref in var.rhs().references() if isinstance(ref, myokit.Derivative)]
        if not derivatives:
            continue
        old = var.rhs().code()
        var.set_rhs(var.rhs().clone(subst={d: d.var().rhs() for d in derivatives}))
        edits.append(f"{var.qname()} = {old}: dot() inlined (gotranx cannot export it)")
    return edits


def convert(ws: Workspace, cellml: Path, voltage: str, stim_variable: str | None, stim: Stim | None, model_id: str | None = None) -> tuple[Path, list[str]]:
    import gotranx
    from gotranx.myokit import myokit_to_gotran

    text = Path(cellml).read_text()
    text, edits = merge_duplicate_connections(text)
    src = ws.out / "import" / Path(cellml).name
    src.parent.mkdir(parents=True, exist_ok=True)
    src.write_text(text)
    model = load_cellml(src)
    edits += normalize_units(model, voltage)
    edits += inline_derivative_references(model)
    if stim is not None and stim_variable:
        edits += rewrite_stimulus(model, stim_variable, stim)
    path = ws.model_dir / "model.ode"
    myokit_to_gotran(model).save(path)
    ode = gotranx.load_ode(path)
    v_names = [s.name for s in ode.states if s.name == voltage.split(".")[-1]]
    meta = ModelMeta(id=model_id or ws.model_dir.name, language="cellml", source=str(cellml), voltage=v_names[0] if v_names else voltage,
                     stim_variable=stim_variable, stim=stim, bcl=stim.period if stim else 1000.0, edits=edits, status="converted")
    if not (ws.model_dir / "meta.toml").exists():
        save_meta(meta, ws.model_dir)
    (ws.out / "convert.md").write_text("# Conversion\n\n" + "\n".join(f"- {e}" for e in edits) + f"\n\n{ode.num_states} states, {ode.num_parameters} parameters -> `{path}`\n")
    print(f"wrote {path} ({ode.num_states} states); edits:\n" + "\n".join(f"- {e}" for e in edits))
    return path, edits
