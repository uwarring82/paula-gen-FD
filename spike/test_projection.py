"""
Tests for the mode-projection engine — Raman (TPSR) combination -> motional mode.

The geometry is fixed by the four Delta_k directions and the radial-mode tilt; the
engine must reproduce Doerr 2024's documented addressing and the 45 deg AC-axial /
30 deg radial direction cosines.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.projection import ModeProjection, _unit, mode_axes
from spike.ledger import Ledger

COS30, SIN30, COS45 = math.cos(math.radians(30)), math.sin(math.radians(30)), math.cos(math.radians(45))

# canonical normalised Delta_k directions (trap x,y,z)
DK = {
    "CC": (0.0, 0.0, 0.0),
    "OC": (0.0, 0.0, 1.0),
    "AC": (-1 / math.sqrt(2), 0.0, -1 / math.sqrt(2)),
    "ROC": (1.0, 0.0, 0.0),
}
DOERR = {"CC": (), "OC": ("lf",), "AC": ("lf", "mf", "hf"), "ROC": ("mf", "hf")}


def _eng(tilt=30.0):
    return ModeProjection(tilt_deg=tilt, dk_directions=DK)


def test_mode_axes_orthonormal_and_tilted():
    ax = mode_axes(30.0)
    # lf is axial
    assert ax["lf"] == (0.0, 0.0, 1.0)
    # radial modes are in the x-y plane and mutually orthogonal
    assert ax["mf"][2] == 0.0 and ax["hf"][2] == 0.0
    assert sum(a * b for a, b in zip(ax["mf"], ax["hf"])) == pytest.approx(0.0, abs=1e-12)
    # mf is tilted 30 deg from +x
    assert ax["mf"][0] == pytest.approx(COS30) and ax["mf"][1] == pytest.approx(SIN30)


def test_addressed_modes_reproduce_doerr():
    eng = _eng()
    for comb in ("CC", "OC", "AC", "ROC"):
        assert eng.addressed_modes(comb) == DOERR[comb]


def test_cc_is_carrier_only():
    eng = _eng()
    assert eng.addressed_modes("CC") == ()
    assert all(v == pytest.approx(0.0) for v in eng.projections("CC").values())


def test_oc_addresses_axial_only():
    eng = _eng()
    p = eng.projections("OC")
    assert p["lf"] == pytest.approx(1.0)
    assert p["mf"] == pytest.approx(0.0, abs=1e-12) and p["hf"] == pytest.approx(0.0, abs=1e-12)


def test_ac_addresses_all_three_with_45deg_axial():
    eng = _eng()
    p = eng.projections("AC")
    assert p["lf"] == pytest.approx(COS45)            # 45 deg to the lf axis (Doerr)
    assert p["mf"] == pytest.approx(COS30 / math.sqrt(2))  # two-step: x-comp (1/sqrt2) x cos30
    assert p["hf"] == pytest.approx(SIN30 / math.sqrt(2))  # x-comp (1/sqrt2) x sin30


def test_roc_addresses_radials_with_tilt_cosines():
    eng = _eng()
    p = eng.projections("ROC")
    assert p["lf"] == pytest.approx(0.0, abs=1e-12)
    assert p["mf"] == pytest.approx(COS30) and p["hf"] == pytest.approx(SIN30)


def test_tilt_zero_collapses_one_radial():
    # with no tilt, x-directed ROC fully projects onto mf (along x) and nothing onto hf
    eng = _eng(tilt=0.0)
    p = eng.projections("ROC")
    assert p["mf"] == pytest.approx(1.0) and p["hf"] == pytest.approx(0.0, abs=1e-12)


def test_from_ledger_consumes_vectors_and_reproduces_doerr():
    eng = ModeProjection.from_ledger(Ledger.load())
    assert eng.tilt_deg == pytest.approx(30.0)
    for comb in ("CC", "OC", "AC", "ROC"):
        assert eng.addressed_modes(comb) == DOERR[comb]


def test_from_ledger_refuses_benchmark_as_tilt():
    with pytest.raises(ValueError):
        ModeProjection.from_ledger(Ledger.load(), tilt_name="doppler_cooling_limit_25mg")


def test_unit_normalises_and_preserves_zero():
    assert _unit((2.0, 0.0, 0.0)) == pytest.approx((1.0, 0.0, 0.0))
    assert _unit((0.0, 0.0, 0.0)) == (0.0, 0.0, 0.0)   # CC carrier: zero stays zero


def test_non_unit_dk_gives_same_direction_cosine():
    # an un-normalised dk (5*z) must yield the same projection as the unit z
    eng = ModeProjection(30.0, {"X": (0.0, 0.0, 5.0)})
    assert eng.projection("X", "lf") == pytest.approx(1.0)


def test_addressed_threshold_drops_fp_residual_at_degeneracy():
    # ROC (+x) onto radial modes: at 89 deg both addressed; at exactly 90 deg the
    # mf projection is FP residual (~6e-17) and must be dropped, hf stays.
    lo = ModeProjection(89.0, {"ROC": (1.0, 0.0, 0.0)})
    assert "mf" in lo.addressed_modes("ROC")          # cos89 = 0.017 > 1e-9
    deg = ModeProjection(90.0, {"ROC": (1.0, 0.0, 0.0)})
    assert deg.addressed_modes("ROC") == ("hf",)       # mf=cos90~0 dropped, hf=1 kept


def test_malformed_dk_vector_rejected():
    with pytest.raises(ValueError):
        ModeProjection(30.0, {"BAD": (0.0, 1.0)})      # 2-vector, not (x,y,z)
