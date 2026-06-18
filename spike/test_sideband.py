"""
Tests for the sideband engine — absolute Lamb-Dicke, sideband Rabi rates, and the
Raman differential AC-Stark factor.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.sideband import Sideband, raman_differential_stark_factor
from spike.ledger import Ledger

# magnitude-bearing Delta_k/k vectors (trap x,y,z)
SQRT2 = math.sqrt(2.0)
DK = {
    "CC": (0.0, 0.0, 0.0),
    "OC": (0.0, 0.0, SQRT2),
    "AC": (-SQRT2, 0.0, -SQRT2),
    "ROC": (SQRT2, 0.0, 0.0),
}


def _sb(ref_eta=0.32, tilt=30.0):
    return Sideband(DK, tilt_deg=tilt, ref_eta=ref_eta)


def test_anchor_reproduced():
    # eta(OC, lf, omega_ref) must equal the anchor eta_ref by construction
    assert _sb().lamb_dicke("OC", "lf", 1.92e6) == pytest.approx(0.32)


def test_eta_scales_as_inverse_sqrt_omega():
    sb = _sb()
    e1 = sb.lamb_dicke("OC", "lf", 1.0e6)
    e4 = sb.lamb_dicke("OC", "lf", 4.0e6)
    assert e1 / e4 == pytest.approx(2.0)   # eta ~ omega^-1/2


def test_ac_and_oc_share_axial_eta():
    # AC's z-component equals OC's full |Delta_k| (both sqrt2 onto lf)
    sb = _sb()
    assert sb.lamb_dicke("AC", "lf", 1.3e6) == pytest.approx(sb.lamb_dicke("OC", "lf", 1.3e6))


def test_cc_has_zero_lamb_dicke():
    sb = _sb()
    assert sb.lamb_dicke("CC", "lf", 1.3e6) == pytest.approx(0.0)
    assert sb.lamb_dicke("CC", "mf", 3.0e6) == pytest.approx(0.0)


def test_sideband_rabi_blue_red():
    sb = _sb()
    eta = sb.lamb_dicke("OC", "lf", 1.3e6)
    assert sb.sideband_rabi("OC", "lf", 1.3e6, 1.0e5, n=0, order="blue") == pytest.approx(eta * 1.0e5)
    assert sb.sideband_rabi("OC", "lf", 1.3e6, 1.0e5, n=4, order="blue") == pytest.approx(eta * math.sqrt(5) * 1.0e5)
    assert sb.sideband_rabi("OC", "lf", 1.3e6, 1.0e5, n=0, order="red") == pytest.approx(0.0)


def test_raman_differential_stark_factor():
    # omega_HF / Delta_R ~ 1.789 GHz / 20 GHz ~ 0.089
    assert raman_differential_stark_factor(1.789e9, 20e9) == pytest.approx(0.0894, abs=1e-3)


def test_from_ledger_anchor_and_addressed_etas():
    ledger = Ledger.load()
    sb = Sideband.from_ledger(ledger)
    assert sb.lamb_dicke("OC", "lf", 1.92e6) == pytest.approx(0.32)
    # AC addresses all three with decreasing eta (higher-freq modes -> smaller extent)
    lf = sb.lamb_dicke("AC", "lf", ledger.value("omega_z_axial_com_25mg"))
    mf = sb.lamb_dicke("AC", "mf", ledger.value("omega_radial_mf_25mg"))
    hf = sb.lamb_dicke("AC", "hf", ledger.value("omega_radial_hf_25mg"))
    assert lf > mf > hf > 0


def test_from_ledger_refuses_benchmark():
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        Sideband.from_ledger(ledger, tilt_name="doppler_cooled_occupation_25mg")


def test_ledger_oc_vector_magnitude_is_sqrt2():
    # the anchor (eta=0.32 for OC) relies on |Delta_k_OC/k| = sqrt2; pin the record
    oc = Ledger.load().value("raman_oc_combination_25mg")
    assert math.hypot(*oc) == pytest.approx(SQRT2)


def test_lamb_dicke_rejects_zero_frequency():
    with pytest.raises(ValueError):
        _sb().lamb_dicke("OC", "lf", 0.0)
