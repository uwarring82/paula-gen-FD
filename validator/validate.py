#!/usr/bin/env python3
"""
iontrap-reference validator
============================

Two layers of checking, run over the whole repository:

  1. Field-level validity via JSON Schema (schema/*.schema.json) — required
     fields, enums, types, cross-field conditionals.
  2. Graph-level invariants JSON Schema cannot express:
       - reference resolution  (source.ref -> sources; generation ->
         generations; configuration slot keys -> configuration_slots;
         derived_from -> records; supersedes -> sibling key)
       - the WALL test          (invariants 1, 2, 4, 5: kinds well-formed;
         differential AC Stark shift is never `input`; benchmark-by-default
         and input-by-default lints)
       - the INHERITANCE test   (invariant 3: no `kind: input` record whose
         transitive derived_from closure touches a `benchmark`)
       - acyclic derived_from
       - benchmarks carry a valid `measured_on`; all dates are real ISO dates

Usage:
    python validator/validate.py [REPO_ROOT]

Exit codes: 0 = clean, 1 = one or more ERRORs, 2 = usage / IO failure.

Dependencies: pyyaml, jsonschema  (see pyproject.toml).
"""
from __future__ import annotations

import datetime as _dt
import json
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.stderr.write("missing dependency: pyyaml. Install with `pip install 'pyyaml>=6'`\n")
    raise SystemExit(2)
try:
    from jsonschema import Draft202012Validator, FormatChecker
except ImportError:  # pragma: no cover
    import importlib.metadata as _md

    try:
        _ver = _md.version("jsonschema")
    except Exception:
        _ver = "not installed"
    sys.stderr.write(
        f"jsonschema >= 4.18 is required (Draft 2020-12 support); found {_ver}. "
        'This project needs Python >= 3.10. Install with `pip install -e ".[dev]"` '
        "in a virtualenv, or `pip install 'jsonschema>=4.18'`.\n"
    )
    raise SystemExit(2)


# --------------------------------------------------------------------------- #
# YAML loader: parse unsigned-exponent scientific notation (2.213e6) as float, #
# and keep ISO dates as plain strings so JSON-Schema string validation applies #
# --------------------------------------------------------------------------- #
class RecordLoader(yaml.SafeLoader):
    pass


RecordLoader.yaml_implicit_resolvers = {
    k: list(v) for k, v in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
# Drop the timestamp resolver -> "2026-06-18" stays a string, not a date object.
for _ch, _mappers in list(RecordLoader.yaml_implicit_resolvers.items()):
    RecordLoader.yaml_implicit_resolvers[_ch] = [
        (tag, rx) for (tag, rx) in _mappers if tag != "tag:yaml.org,2002:timestamp"
    ]
# Add a YAML-1.2-style float resolver (handles 2.213e6, 1e3, .5e-2, etc.).
_FLOAT_RX = re.compile(
    r"""^[-+]?(?:
            [0-9][0-9_]*\.[0-9_]*(?:[eE][-+]?[0-9]+)?
          | \.[0-9_]+(?:[eE][-+]?[0-9]+)?
          | [0-9][0-9_]*[eE][-+]?[0-9]+
        )$
        | ^[-+]?\.(?:inf|Inf|INF)$
        | ^\.(?:nan|NaN|NAN)$""",
    re.X,
)
RecordLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float", _FLOAT_RX, list("-+0123456789.")
)


# --------------------------------------------------------------------------- #
# Reporting                                                                    #
# --------------------------------------------------------------------------- #
class Report:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def error(self, where: str, msg: str) -> None:
        self.errors.append(f"{where}: {msg}")

    def warn(self, where: str, msg: str) -> None:
        self.warnings.append(f"{where}: {msg}")

    def ok(self) -> bool:
        return not self.errors


_DATE_RX = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def valid_date(value) -> bool:
    if not isinstance(value, str) or not _DATE_RX.match(value):
        return False
    try:
        _dt.date.fromisoformat(value)
        return True
    except ValueError:
        return False


def load_yaml(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return yaml.load(fh, Loader=RecordLoader)


def load_registry(path: Path, rep: Report):
    """Load a registry YAML, reporting (not raising) if it is missing/malformed."""
    if not path.exists():
        rep.error(path.name, "registry file missing")
        return {}
    doc = load_yaml(path)
    if doc is None:
        return {}
    if not isinstance(doc, dict):
        rep.error(path.name, "registry must be a mapping")
        return {}
    return doc


# --------------------------------------------------------------------------- #
# Field-level: JSON Schema                                                     #
# --------------------------------------------------------------------------- #
def schema_check(root: Path, rep: Report, records: list[dict]) -> None:
    schema_dir = root / "schema"
    fmt = FormatChecker()

    def validator_for(name: str):
        return Draft202012Validator(
            json.loads((schema_dir / name).read_text(encoding="utf-8")),
            format_checker=fmt,
        )

    # Records (each list item against the record schema)
    rec_v = validator_for("record.schema.json")
    for i, rec in enumerate(records):
        if not isinstance(rec, dict):
            rep.error(f"records[{i}]", f"record is not a mapping: {rec!r}")
            continue
        label = rec.get("_where", "record")
        payload = {k: v for k, v in rec.items() if k != "_where"}
        for err in sorted(rec_v.iter_errors(payload), key=str):
            loc = "/".join(str(p) for p in err.absolute_path) or "<root>"
            rep.error(label, f"schema: {loc}: {err.message}")

    # Registries (whole-document)
    for fname, schema in (
        ("sources.yaml", "sources.schema.json"),
        ("generations.yaml", "generations.schema.json"),
        ("configuration_slots.yaml", "configuration_slots.schema.json"),
    ):
        path = root / "registries" / fname
        if not path.exists():
            rep.error(fname, "registry file missing")
            continue
        doc = load_yaml(path) or {}
        v = validator_for(schema)
        for err in sorted(v.iter_errors(doc), key=str):
            loc = "/".join(str(p) for p in err.absolute_path) or "<root>"
            rep.error(fname, f"schema: {loc}: {err.message}")


# --------------------------------------------------------------------------- #
# Graph-level invariants                                                       #
# --------------------------------------------------------------------------- #
# Benchmark-by-default / input-by-default lint patterns (invariants 4 & 5).
_BENCHMARK_BY_DEFAULT = [
    re.compile(p)
    for p in (
        r"decoher",
        r"dephas",
        r"heating_rate",
        r"residual",
    )
]
# Differential AC Stark shift: benchmark *never* input (Out-of-scope + inv. 4).
_STARK_RX = re.compile(r"(differential.*stark|stark.*shift)")
_INPUT_BY_DEFAULT = [
    re.compile(p)
    for p in (r"_power$", r"waist", r"detuning", r"k_vector", r"polaris", r"polariz")
]


def graph_check(root: Path, rep: Report, records: list[dict]) -> None:
    sources = load_registry(root / "registries" / "sources.yaml", rep)
    generations = load_registry(root / "registries" / "generations.yaml", rep)
    slots = load_registry(root / "registries" / "configuration_slots.yaml", rep)

    by_name: dict[str, dict] = {}
    for rec in records:
        if not isinstance(rec, dict):
            continue
        name = rec.get("name")
        where = rec.get("_where", "record")
        if not isinstance(name, str):
            continue
        if name in by_name:
            rep.error(where, f"duplicate record name '{name}'")
            continue
        by_name[name] = rec

    # --- supersedes resolution in the registries -------------------------- #
    for gkey, gval in (generations or {}).items():
        sup = (gval or {}).get("supersedes")
        if sup is not None and sup not in generations:
            rep.error("generations.yaml", f"'{gkey}'.supersedes '{sup}' does not resolve")
    for slot_name, slot in (slots or {}).items():
        for skey, sval in (slot or {}).items():
            sup = (sval or {}).get("supersedes")
            if sup is not None and sup not in slot:
                rep.error(
                    "configuration_slots.yaml",
                    f"{slot_name}.{skey}.supersedes '{sup}' does not resolve within slot",
                )

    # --- per-record reference resolution + dates + kind lints ------------- #
    for name, rec in by_name.items():
        where = rec.get("_where", name)

        src = rec.get("source")
        if not isinstance(src, dict):
            src = {}
        ref = src.get("ref")
        if ref is not None and ref not in sources:
            rep.error(where, f"source.ref '{ref}' does not resolve in sources registry")
        if not valid_date(src.get("extracted_on")):
            rep.error(where, f"source.extracted_on '{src.get('extracted_on')}' is not a YYYY-MM-DD date")

        gen = rec.get("generation")
        if gen is not None and gen not in generations:
            rep.error(where, f"generation '{gen}' does not resolve in generations registry")

        cfg = rec.get("configuration")
        if isinstance(cfg, dict):
            if cfg and all(v is None for v in cfg.values()):
                rep.error(
                    where,
                    "configuration is an object but every slot is null — it carries no "
                    "resolvable key. Use configuration: null for apparatus-independent "
                    "quantities, or set at least one slot (trap/beams/b_field).",
                )
            for slot_name, slot_key in cfg.items():
                if slot_key is None:
                    continue
                if slot_name not in slots or slot_key not in (slots.get(slot_name) or {}):
                    rep.error(
                        where,
                        f"configuration.{slot_name} '{slot_key}' does not resolve in configuration_slots",
                    )

        if rec.get("kind") == "benchmark":
            mo = rec.get("measured_on")
            if mo is None:
                rep.error(where, "benchmark record is missing measured_on [req-b]")
            elif not valid_date(mo):
                rep.error(where, f"measured_on '{mo}' is not a YYYY-MM-DD date")
        elif "measured_on" in rec and not valid_date(rec["measured_on"]):
            rep.error(where, f"measured_on '{rec['measured_on']}' is not a YYYY-MM-DD date")

        for dep in rec.get("derived_from", []) or []:
            if dep not in by_name:
                rep.error(where, f"derived_from '{dep}' does not resolve to a known record")

    # --- the WALL, soft side: AC-Stark hard rule + default-kind lints ----- #
    _check_kind_wall(rep, by_name)

    # --- source resolvability of REFERENCED sources (warn, don't fail) ---- #
    def _ref_of(r):
        s = r.get("source")
        return s.get("ref") if isinstance(s, dict) else None

    referenced = {_ref_of(r) for r in by_name.values() if _ref_of(r) in sources}
    for ref in sorted(x for x in referenced if x):
        entry = sources.get(ref) or {}
        if not any(entry.get(k) for k in ("link", "doi", "urn")):
            rep.warn(
                f"sources.yaml[{ref}]",
                "referenced source has no resolvable identifier (link/doi/urn) — "
                "traceability is degraded (FAIR Accessible)",
            )
        elif entry.get("verified") is not True:
            rep.warn(
                f"sources.yaml[{ref}]",
                "referenced source identifier is not marked verified — confirm it resolves",
            )

    # --- acyclic derived_from + transitive benchmark inheritance ---------- #
    _check_inheritance_and_cycles(rep, by_name)


def _check_kind_wall(rep: Report, by_name: dict[str, dict]) -> None:
    """The WALL, soft side (invariants 2, 4, 5). A differential AC Stark shift
    named `input` is a hard error (it is both calibrated and predictable, so
    feeding it in closes the loop); benchmark-by-default and input-by-default
    names carrying the opposite `kind` raise warnings to flag likely mistakes."""
    for name, rec in by_name.items():
        where = rec.get("_where", name)
        lname = name.lower()
        kind = rec.get("kind")
        if _STARK_RX.search(lname) and kind == "input":
            rep.error(
                where,
                "differential AC Stark shift must be kind:benchmark, never input "
                "(it is both calibrated and predictable -> feeding it in closes the loop)",
            )
        elif kind == "input" and any(p.search(lname) for p in _BENCHMARK_BY_DEFAULT):
            rep.warn(
                where,
                f"'{name}' looks benchmark-by-default (decoherence/heating/residual) "
                "but is marked input — confirm this is a deliberate hold-out",
            )
        elif kind == "benchmark" and any(p.search(lname) for p in _INPUT_BY_DEFAULT):
            rep.warn(
                where,
                f"'{name}' looks input-by-default (beam power/waist/detuning/…) "
                "but is marked benchmark — confirm this is a deliberate hold-out",
            )


def _check_inheritance_and_cycles(rep: Report, by_name: dict[str, dict]) -> None:
    WHITE, GREY, BLACK = 0, 1, 2
    colour = {n: WHITE for n in by_name}

    def closure_has_benchmark(start: str) -> tuple[bool, list[str]]:
        """Return (touches_benchmark, the benchmark names found) over the
        transitive derived_from closure of `start` (excluding `start`)."""
        seen: set[str] = set()
        found: list[str] = []
        stack = list(by_name[start].get("derived_from", []) or [])
        while stack:
            cur = stack.pop()
            if cur in seen or cur not in by_name:
                continue
            seen.add(cur)
            if by_name[cur].get("kind") == "benchmark":
                found.append(cur)
            stack.extend(by_name[cur].get("derived_from", []) or [])
        return (bool(found), found)

    def visit(n: str, path: list[str]) -> None:
        colour[n] = GREY
        for dep in by_name[n].get("derived_from", []) or []:
            if dep not in by_name:
                continue
            if colour[dep] == GREY:
                cyc = " -> ".join(path + [n, dep])
                rep.error(by_name[n].get("_where", n), f"derived_from cycle: {cyc}")
            elif colour[dep] == WHITE:
                visit(dep, path + [n])
        colour[n] = BLACK

    for n in by_name:
        if colour[n] == WHITE:
            visit(n, [])

    # Inheritance test (invariant 3)
    for name, rec in by_name.items():
        if rec.get("kind") != "input":
            continue
        touches, found = closure_has_benchmark(name)
        if touches:
            rep.error(
                rec.get("_where", name),
                "INHERITANCE VIOLATION: kind:input but its transitive derived_from "
                f"closure contains benchmark record(s) {found}. Such a quantity "
                "inherits benchmark status unless re-derived exclusively from inputs.",
            )


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
def collect_records(root: Path, rep: Report) -> list[dict]:
    records: list[dict] = []
    rec_dir = root / "records"
    if not rec_dir.is_dir():
        rep.error("records/", "records directory missing")
        return records
    for path in sorted(rec_dir.glob("*.yaml")):
        doc = load_yaml(path)
        if doc is None:
            continue
        if not isinstance(doc, list):
            rep.error(path.name, "records file must be a YAML list of records")
            continue
        for i, rec in enumerate(doc):
            if isinstance(rec, dict):
                rec = dict(rec)
                rec["_where"] = f"{path.name}[{i}] ({rec.get('name', '?')})"
            records.append(rec)
    return records


def main(argv: list[str]) -> int:
    root = Path(argv[1]).resolve() if len(argv) > 1 else Path(__file__).resolve().parent.parent
    rep = Report()

    records = collect_records(root, rep)
    schema_check(root, rep, records)
    # Graph-level checks assume schema-valid, well-typed data; running them on
    # malformed records would be noisy at best and crash at worst. Defer them
    # until field-level errors are fixed.
    deferred = not rep.ok()
    if not deferred:
        graph_check(root, rep, records)

    n_rec = sum(1 for r in records if isinstance(r, dict))
    print(f"iontrap-reference validator — {n_rec} record(s) checked under {root}")
    for w in rep.warnings:
        print(f"  WARN  {w}")
    for e in rep.errors:
        print(f"  ERROR {e}")

    if rep.ok():
        print(f"OK — no errors ({len(rep.warnings)} warning(s)).")
        return 0
    if deferred:
        print("  (graph-level checks deferred until the field-level errors above are fixed)")
    print(f"FAILED — {len(rep.errors)} error(s), {len(rep.warnings)} warning(s).")
    return 1


def main_cli() -> None:
    """Zero-argument console-script entry point (see pyproject.toml)."""
    raise SystemExit(main(sys.argv))


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
