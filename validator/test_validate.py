"""
Unit tests for the iontrap-reference validator.

These prove that the two graph-level invariants actually FIRE — i.e. that the
validator would reject a real violation, not just pass clean data. They also
pin the custom YAML loader behaviour (scientific-notation floats; ISO dates as
strings) that the JSON-Schema layer depends on.

Run:  pytest validator/
"""
import io

import yaml

from validator.validate import (
    RecordLoader,
    Report,
    _check_inheritance_and_cycles,
    valid_date,
)


# --------------------------------------------------------------------------- #
# YAML loader contract                                                         #
# --------------------------------------------------------------------------- #
def _load(text):
    return yaml.load(io.StringIO(text), Loader=RecordLoader)


def test_unsigned_exponent_parses_as_float():
    doc = _load("a: 2.213e6\nb: 1e3\nc: 5.0e-6\n")
    assert doc["a"] == 2.213e6 and isinstance(doc["a"], float)
    assert doc["b"] == 1000.0 and isinstance(doc["b"], float)
    assert doc["c"] == 5.0e-6 and isinstance(doc["c"], float)


def test_iso_date_stays_string():
    doc = _load("d: 2026-06-18\n")
    assert doc["d"] == "2026-06-18" and isinstance(doc["d"], str)


def test_plain_int_stays_int():
    doc = _load("n: 2213000\n")
    assert doc["n"] == 2213000 and isinstance(doc["n"], int)


def test_valid_date():
    assert valid_date("2026-06-18")
    assert not valid_date("2026-13-01")   # bad month
    assert not valid_date("2026-6-1")     # not zero-padded
    assert not valid_date("not-a-date")
    assert not valid_date(20260618)       # not a string


# --------------------------------------------------------------------------- #
# Graph invariant: transitive benchmark inheritance (invariant 3)             #
# --------------------------------------------------------------------------- #
def _rec(name, kind, derived_from=()):
    return {"name": name, "kind": kind, "derived_from": list(derived_from), "_where": name}


def test_inheritance_rejects_input_derived_from_benchmark():
    by_name = {
        "bench": _rec("bench", "benchmark"),
        "mid": _rec("mid", "input", ["bench"]),       # transitively touches benchmark
        "leaf": _rec("leaf", "input", ["mid"]),       # also touches it, two hops away
    }
    rep = Report()
    _check_inheritance_and_cycles(rep, by_name)
    assert not rep.ok()
    msgs = " ".join(rep.errors)
    assert "INHERITANCE VIOLATION" in msgs
    assert "mid" in msgs and "leaf" in msgs


def test_inheritance_allows_input_closure_all_input():
    by_name = {
        "a": _rec("a", "input"),
        "b": _rec("b", "input", ["a"]),
        "c": _rec("c", "input", ["a", "b"]),
    }
    rep = Report()
    _check_inheritance_and_cycles(rep, by_name)
    assert rep.ok(), rep.errors


def test_benchmark_derived_from_benchmark_is_fine():
    # A benchmark whose ancestor is a benchmark is allowed (it inherits, and it
    # already IS benchmark) — only kind:input ancestors-touching-benchmark fail.
    by_name = {
        "bench": _rec("bench", "benchmark"),
        "down": _rec("down", "benchmark", ["bench"]),
    }
    rep = Report()
    _check_inheritance_and_cycles(rep, by_name)
    assert rep.ok(), rep.errors


# --------------------------------------------------------------------------- #
# Graph invariant: acyclic derived_from                                       #
# --------------------------------------------------------------------------- #
def test_cycle_detected():
    by_name = {
        "x": _rec("x", "input", ["y"]),
        "y": _rec("y", "input", ["x"]),
    }
    rep = Report()
    _check_inheritance_and_cycles(rep, by_name)
    assert not rep.ok()
    assert any("cycle" in e for e in rep.errors)
