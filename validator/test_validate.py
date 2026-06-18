"""
Unit tests for the iontrap-reference validator.

These prove that the two graph-level invariants actually FIRE — i.e. that the
validator would reject a real violation, not just pass clean data. They also
pin the custom YAML loader behaviour (scientific-notation floats; ISO dates as
strings) that the JSON-Schema layer depends on.

Run:  pytest validator/
"""
import io
from pathlib import Path

import yaml

from validator.validate import (
    RecordLoader,
    Report,
    _check_inheritance_and_cycles,
    _check_kind_wall,
    graph_check,
    schema_check,
    valid_date,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


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


# --------------------------------------------------------------------------- #
# Graph invariant: the WALL (soft side) — AC-Stark hard rule + default lints   #
# --------------------------------------------------------------------------- #
def _by_name(*recs):
    return {r["name"]: r for r in recs}


def test_wall_rejects_differential_ac_stark_as_input():
    rep = Report()
    _check_kind_wall(rep, _by_name(_rec("differential_ac_stark_shift_qubit", "input")))
    assert not rep.ok()
    assert any("stark" in e.lower() for e in rep.errors)


def test_wall_allows_differential_ac_stark_as_benchmark():
    rep = Report()
    _check_kind_wall(rep, _by_name(_rec("differential_ac_stark_shift_qubit", "benchmark")))
    assert rep.ok() and not rep.warnings


def test_wall_warns_benchmark_by_default_marked_input():
    rep = Report()
    _check_kind_wall(rep, _by_name(_rec("motional_decoherence_rate", "input")))
    assert rep.ok()  # a warning, not an error
    assert any("benchmark-by-default" in w for w in rep.warnings)


def test_wall_warns_input_by_default_marked_benchmark():
    rep = Report()
    _check_kind_wall(rep, _by_name(_rec("raman_beam_power", "benchmark")))
    assert rep.ok()
    assert any("input-by-default" in w for w in rep.warnings)


def test_wall_clean_record_no_findings():
    rep = Report()
    _check_kind_wall(rep, _by_name(_rec("omega_z_axial_com", "input")))
    assert rep.ok() and not rep.warnings


# --------------------------------------------------------------------------- #
# Robustness: malformed records must be reported, never crash the validator    #
# --------------------------------------------------------------------------- #
def _full_record(**overrides):
    rec = {
        "name": "x",
        "value": 1.0,
        "units": "Hz",
        "kind": "input",
        "uncertainty": {"value": 0, "type": "exact"},
        "source": {"ref": "clos2017", "loc": "p.1", "extracted_by": "t", "extracted_on": "2026-06-18"},
        "scope": {"isotope": "apparatus", "subsystem": "fields"},
        "observation_type": "direct",
        "derived_from": [],
        "generation": "legacy",
        "status": "provisional",
        "_where": "test",
    }
    rec.update(overrides)
    return rec


def test_schema_check_survives_non_dict_record():
    rep = Report()
    schema_check(REPO_ROOT, rep, ["i am a string", _full_record()])  # must not raise
    assert any("not a mapping" in e for e in rep.errors)


def test_graph_check_survives_malformed_nested_fields():
    rep = Report()
    bad = [
        "not a mapping",
        {"name": "r1", "source": "oops-a-string", "configuration": "oops", "_where": "r1"},
    ]
    graph_check(REPO_ROOT, rep, bad)  # the bug: .get() on a str -> AttributeError
    assert isinstance(rep.errors, list)  # reached here == no crash


# --------------------------------------------------------------------------- #
# Schema hole: a confirmed record must not pass with an EMPTY configuration    #
# --------------------------------------------------------------------------- #
def test_configuration_empty_object_rejected():
    rep = Report()
    schema_check(REPO_ROOT, rep, [_full_record(status="confirmed", configuration={})])
    assert any("configuration" in e for e in rep.errors)


def test_configuration_with_one_slot_accepted():
    rep = Report()
    schema_check(REPO_ROOT, rep, [_full_record(status="confirmed", configuration={"b_field": "nominal"})])
    assert not any("configuration" in e for e in rep.errors)


def test_configuration_null_accepted_for_confirmed():
    rep = Report()
    schema_check(REPO_ROOT, rep, [_full_record(status="confirmed", configuration=None)])
    assert not any("configuration" in e for e in rep.errors)
