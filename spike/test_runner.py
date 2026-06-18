"""
Tests for the composition root (spike/runner.py): the ValidationResult math, the
real engine validations being consistent, and that the runner DETECTS tension,
SURFACES a wall violation as an error, and FLAGS uncovered benchmarks.

Run:  pytest spike/
"""
import math

import pytest

from spike.ledger import Ledger
from spike.runner import (
    REGISTRY,
    Validation,
    ValidationResult,
    _validate_clock,
    render_table,
    run_all,
    uncovered_benchmarks,
)


# --- ValidationResult math --------------------------------------------------
def test_result_residual_sigma_nsigma():
    r = ValidationResult(benchmark="x", engine="e", predicted=1.002e9, measured=1.0e9,
                         units="Hz", sigma_pred=1e6, sigma_meas=1e6)
    assert r.residual == pytest.approx(2e6)
    assert r.sigma == pytest.approx(math.hypot(1e6, 1e6))
    assert r.nsigma == pytest.approx(2e6 / math.hypot(1e6, 1e6))
    assert r.ok  # ~1.4 sigma


def test_result_exact_match_zero_sigma_is_ok():
    r = ValidationResult(benchmark="x", engine="e", predicted=1.0, measured=1.0,
                         sigma_pred=0.0, sigma_meas=0.0)
    assert r.nsigma == 0.0 and r.ok


def test_result_error_is_not_ok():
    r = ValidationResult(benchmark="x", engine="e", error="boom")
    assert r.status == "ERROR" and not r.ok and math.isnan(r.nsigma)


# --- real validations -------------------------------------------------------
def test_real_validations_consistent():
    results = run_all(Ledger.load())
    assert len(results) == len(REGISTRY)
    assert all(r.ok for r in results)
    names = {r.benchmark for r in results}
    assert names == {"clock_transition_25mg", "clock_transition_weber_25mg",
                     "omega_z_axial_stretch_2ion_25mg", "omega_radial_rocking_2ion_25mg"}


def test_render_table_has_rows():
    table = render_table(run_all(Ledger.load()))
    assert "clock_transition_25mg" in table and "omega_z_axial_stretch_2ion_25mg" in table
    assert "n_sigma" in table and "status" in table


# --- the runner catches problems --------------------------------------------
def test_runner_detects_tension():
    def _tension(_ledger):
        # predicted 1% off a tightly-measured benchmark -> many sigma
        return ValidationResult(benchmark="fake", engine="fake", units="Hz",
                                predicted=2.0e6, measured=1.0e6,
                                sigma_pred=0.0, sigma_meas=1.0e3)
    results = run_all(Ledger.load(), registry=[Validation("fake", "fake", _tension)])
    assert results[0].status == "TENSION" and not results[0].ok


def test_runner_surfaces_wall_violation_as_error():
    def _bad(ledger):
        ledger.input_quantity("clock_transition_25mg")   # consume a benchmark -> ValueError
        return ValidationResult(benchmark="clock_transition_25mg", engine="fake")
    results = run_all(Ledger.load(), registry=[Validation("clock_transition_25mg", "fake", _bad)])
    assert results[0].status == "ERROR" and not results[0].ok
    assert "must be kind:input" in results[0].error


def test_uncovered_benchmarks_flags_missing_engine():
    ledger = Ledger.load()
    # a registry that only covers the clock leaves the stretch uncovered
    reg = [Validation("clock_transition_25mg", "levels", _validate_clock)]
    unc = uncovered_benchmarks(ledger, registry=reg)
    assert "omega_z_axial_stretch_2ion_25mg" in unc
    # derived benchmarks (derived_from != []) are bookkeeping, not direct targets
    assert "clock_transition_residual_25mg" not in unc
    # the full registry leaves nothing uncovered
    assert uncovered_benchmarks(ledger) == []
