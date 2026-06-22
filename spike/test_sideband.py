"""
Tests for the sideband engine — absolute Lamb-Dicke, sideband Rabi rates, and the
Raman differential AC-Stark factor.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.sideband import (
    Sideband,
    _laguerre,
    carrier_rabi_factor,
    raman_differential_stark_factor,
    thermal_carrier_flip,
    thermal_dephasing_rate,
    thermal_pn,
)
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


# --- carrier Debye-Waller thermal dephasing --------------------------------
def test_laguerre_low_orders():
    x = 0.151
    assert _laguerre(0, x) == pytest.approx(1.0)
    assert _laguerre(1, x) == pytest.approx(1.0 - x)
    assert _laguerre(2, x) == pytest.approx(1.0 - 2.0 * x + x * x / 2.0)
    assert _laguerre(3, x) == pytest.approx(1.0 - 3 * x + 1.5 * x * x - x ** 3 / 6.0)


def test_thermal_pn_is_a_distribution():
    nbar = 0.7
    ps = [thermal_pn(nbar, n) for n in range(200)]
    assert sum(ps) == pytest.approx(1.0, abs=1e-6)
    assert thermal_pn(nbar, 0) == pytest.approx(1.0 / (nbar + 1.0))
    assert thermal_pn(0.0, 0) == 1.0 and thermal_pn(0.0, 1) == 0.0
    # geometric ratio P_{n+1}/P_n = nbar/(nbar+1)
    assert ps[5] / ps[4] == pytest.approx(nbar / (nbar + 1.0))


def test_carrier_rabi_factor_limits():
    eta = 0.389
    assert carrier_rabi_factor(eta, 0) == pytest.approx(math.exp(-eta * eta / 2.0))
    assert carrier_rabi_factor(0.0, 5) == pytest.approx(1.0)            # no coupling spread
    # n=1 is slower than n=0 in the Lamb-Dicke regime (L_1 = 1-eta^2 < 1)
    assert carrier_rabi_factor(eta, 1) < carrier_rabi_factor(eta, 0)


def test_thermal_flip_nbar0_is_coherent_debye_waller_flop():
    eta, rabi = 0.389, 2.0e5
    om0 = rabi * math.exp(-eta * eta / 2.0)        # Debye-Waller-reduced carrier rate
    t_pi = 1.0 / (2.0 * om0)
    assert thermal_carrier_flip(0.0, rabi, eta, 0.0) == pytest.approx(0.0)
    assert thermal_carrier_flip(t_pi, rabi, eta, 0.0) == pytest.approx(1.0, abs=1e-6)


def test_thermal_flip_decays_with_nbar():
    eta, rabi = 0.389, 2.0e5
    om0 = rabi * math.exp(-eta * eta / 2.0)
    t2pi = 1.0 / om0                                # one full period: coherent returns to ~0
    cold = thermal_carrier_flip(t2pi, rabi, eta, 0.0)
    warm = thermal_carrier_flip(t2pi, rabi, eta, 1.0)
    hot = thermal_carrier_flip(t2pi, rabi, eta, 5.0)
    assert cold < warm < hot                        # dephasing fills in the trough toward 1/2
    assert hot == pytest.approx(0.5, abs=0.1)        # large nbar -> mixed within a period


def test_thermal_flip_detuning_caps_amplitude():
    eta, rabi = 0.389, 2.0e5
    om0 = rabi * math.exp(-eta * eta / 2.0)
    # with nbar=0 and a detuning, the peak is capped below 1
    peak_on = thermal_carrier_flip(1.0 / (2 * om0), rabi, eta, 0.0, detuning_hz=0.0)
    eff = math.hypot(om0, 0.5 * rabi)
    peak_off = thermal_carrier_flip(1.0 / (2 * eff), rabi, eta, 0.0, detuning_hz=0.5 * rabi)
    assert peak_on == pytest.approx(1.0, abs=1e-6)
    assert peak_off < 0.95


def test_thermal_dephasing_rate_scaling():
    eta, rabi = 0.389, 2.0e5
    assert thermal_dephasing_rate(rabi, eta, 0.0) == 0.0
    base = thermal_dephasing_rate(rabi, eta, 0.5)
    assert thermal_dephasing_rate(2 * rabi, eta, 0.5) == pytest.approx(2 * base)      # ~ Omega_0
    assert thermal_dephasing_rate(rabi, 2 * eta, 0.5) == pytest.approx(4 * base)      # ~ eta^2
    # ~ sqrt(nbar(nbar+1))
    r = thermal_dephasing_rate(rabi, eta, 2.0) / base
    assert r == pytest.approx(math.sqrt(2 * 3) / math.sqrt(0.5 * 1.5))


def test_carrier_thermal_flip_method_matches_function():
    sb = _sb()
    eta = sb.lamb_dicke("OC", "lf", 1.30e6)
    assert sb.carrier_thermal_flip("OC", "lf", 1.30e6, 2.0e5, 0.3, 3e-6) == pytest.approx(
        thermal_carrier_flip(3e-6, 2.0e5, eta, 0.3))


# --- thermal sideband flop (RSB/BSB thermometry) ---------------------------
def test_sideband_rabi_factor_red_blue():
    from spike.engines.sideband import sideband_rabi_factor
    assert sideband_rabi_factor(0.39, 0, "red") == 0.0           # ground state can't subtract
    assert sideband_rabi_factor(0.39, 0, "blue") == pytest.approx(0.39)   # n=0 -> n=1 ok
    assert sideband_rabi_factor(0.39, 4, "red") == pytest.approx(0.39 * 2.0)   # sqrt(4)=2
    assert sideband_rabi_factor(0.39, 3, "blue") == pytest.approx(0.39 * 2.0)  # sqrt(3+1)=2


def test_thermal_sideband_red_suppressed_when_cold():
    from spike.engines.sideband import thermal_sideband_flip
    ts = [0.5e-6 * k for k in range(40)]
    red_cold = max(thermal_sideband_flip(t, 80e3, 0.39, 0.05, "red") for t in ts)
    blue_cold = max(thermal_sideband_flip(t, 80e3, 0.39, 0.05, "blue") for t in ts)
    assert red_cold < 0.1                       # cold RSB barely flops (n=0 stuck)
    assert blue_cold > 0.8                       # cold BSB flops nearly fully
    # red amplitude grows with nbar (thermometer)
    red_warm = max(thermal_sideband_flip(t, 80e3, 0.39, 1.0, "red") for t in ts)
    assert red_warm > red_cold


def test_thermal_sideband_red_amplitude_tracks_1_minus_P0():
    from spike.engines.sideband import thermal_sideband_flip, thermal_pn
    ts = [0.3e-6 * k for k in range(80)]
    for nbar in (0.1, 0.5, 2.0):
        red_peak = max(thermal_sideband_flip(t, 80e3, 0.39, nbar, "red") for t in ts)
        # the red flop saturates near (1 - P_0)/2 <= 1 - P_0 (only n>=1 can flop)
        assert red_peak <= (1.0 - thermal_pn(nbar, 0)) + 1e-6
