"""
Tests for the Raman-beam (relative-phase) dephasing engine: the two contrast
envelopes (Lorentzian -> exponential, Gaussian -> quasi-static), the rate <->
mutual-linewidth inversion used to read out the OC-flop residual, the T2 <-> sigma_nu
relations, and the path-imbalance Wiener phase variance.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.raman_dephasing import (
    C_LIGHT,
    coherence_gaussian,
    coherence_lorentzian,
    coherence_time_from_rate,
    dephasing_rate,
    mutual_linewidth_from_rate,
    path_imbalance_phase_variance,
    sigma_nu_from_t2,
    static_contrast,
    t2_from_sigma_nu,
)


def test_lorentzian_envelope():
    dnu = 10e3
    assert coherence_lorentzian(0.0, dnu) == pytest.approx(1.0)
    assert coherence_lorentzian(5e-6, dnu) == pytest.approx(math.exp(-math.pi * dnu * 5e-6))
    assert coherence_lorentzian(1e-6, dnu) < 1.0          # decays
    assert coherence_lorentzian(1e-6, 0.0) == 1.0         # no linewidth -> no decay


def test_gaussian_envelope():
    t2 = 10e-6
    assert coherence_gaussian(0.0, t2) == pytest.approx(1.0)
    assert coherence_gaussian(t2, t2) == pytest.approx(math.exp(-1.0))
    assert coherence_gaussian(5e-6, float("inf")) == 1.0  # infinite T2 -> no decay


def test_rate_linewidth_roundtrip():
    dnu = 28e3
    g = dephasing_rate(dnu)
    assert g == pytest.approx(math.pi * dnu)
    assert mutual_linewidth_from_rate(g) == pytest.approx(dnu)     # inverse
    assert coherence_time_from_rate(g) == pytest.approx(1.0 / g)
    # the exponential envelope at t = T_phi has dropped to 1/e? NO: T_phi = 1/Gamma,
    # and C(t)=exp(-Gamma t), so C(T_phi)=1/e by construction
    assert coherence_lorentzian(coherence_time_from_rate(g), dnu) == pytest.approx(math.exp(-1.0))


def test_t2_sigma_nu_roundtrip():
    sigma_nu = 5e3
    t2 = t2_from_sigma_nu(sigma_nu)
    assert t2 == pytest.approx(math.sqrt(2.0) / (2 * math.pi * sigma_nu))
    assert sigma_nu_from_t2(t2) == pytest.approx(sigma_nu)
    assert t2_from_sigma_nu(0.0) == float("inf")


def test_path_imbalance_phase_variance():
    # <dphi^2> = 2 pi dnu_laser tau, tau = dL/c (Wiener); static contrast = exp(-var/2)
    dnu_laser, dL = 1e3, 0.3       # 1 kHz laser, 30 cm path imbalance
    var = path_imbalance_phase_variance(dnu_laser, dL)
    assert var == pytest.approx(2 * math.pi * dnu_laser * dL / C_LIGHT)
    assert 0.0 < static_contrast(var) <= 1.0
    assert path_imbalance_phase_variance(0.0, dL) == 0.0          # perfect laser -> no jitter
    # scales linearly in both laser linewidth and path difference
    assert path_imbalance_phase_variance(2 * dnu_laser, dL) == pytest.approx(2 * var)
    assert path_imbalance_phase_variance(dnu_laser, 2 * dL) == pytest.approx(2 * var)


def test_residual_readout_scale():
    # a ~9e4 /s residual (the OC flop) -> a ~28 kHz mutual linewidth, T_phi ~ 11 us
    g_resid = 8.8e4
    assert mutual_linewidth_from_rate(g_resid) == pytest.approx(2.8e4, rel=0.1)
    assert coherence_time_from_rate(g_resid) * 1e6 == pytest.approx(11.4, rel=0.1)


def test_path_variation_for_phase():
    from spike.engines.raman_dephasing import path_variation_for_phase
    lam = 279.6e-9
    assert path_variation_for_phase(2 * math.pi, lam) == pytest.approx(lam)        # full 2pi = lambda
    assert path_variation_for_phase(1.0, lam) * 1e9 == pytest.approx(44.5, rel=1e-2)  # 1 rad ~ 45 nm
    assert path_variation_for_phase(0.0, lam) == 0.0


def test_path_jitter_velocity():
    from spike.engines.raman_dephasing import path_jitter_velocity
    lam = 279.6e-9
    # v = lambda * df ; a 21 kHz mutual linewidth -> ~5.9 mm/s
    assert path_jitter_velocity(21e3, lam) * 1e3 == pytest.approx(5.87, rel=1e-2)
    assert path_jitter_velocity(2 * 21e3, lam) == pytest.approx(2 * path_jitter_velocity(21e3, lam))
