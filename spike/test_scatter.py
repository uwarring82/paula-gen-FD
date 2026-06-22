"""
Tests for the Raman scatter / differential-AC-Stark engine: the off-resonant
scattering rate, the detuning-only SE-per-pi figure of merit, the balanced-beam
collapse Gamma_sc = (Gamma/Delta_R) Omega, the 1/Delta_R scaling, and the forward
flop (AC-Stark cap + scattering envelope).

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.scatter import (
    CONTRAST_DECAY_FACTOR,
    RamanScatter,
    contrast_decay_rate,
    differential_stark_shift,
    flip_probability,
    scatter_rate,
    scatter_rate_from_rabi,
    se_probability_per_pi,
    two_photon_rabi,
)
from spike.ledger import Ledger

GAMMA = 41.8e6      # 3P_3/2 linewidth (mg_p32_natural_linewidth)
DELTA_R = 20.0e9    # raman_detuning_from_p32
OMEGA_HF = 1.79e9   # hyperfine_splitting_25mg_f2_f3


def test_two_photon_rabi_formula():
    assert two_photon_rabi(1e9, 1e9, DELTA_R) == pytest.approx(1e18 / (2 * DELTA_R))


def test_scatter_rate_formula():
    # true rate in 1/s: Gamma_decay = 2pi*gamma turns the linewidth into a rate
    ob, orr = 1.2e9, 0.8e9
    assert scatter_rate(ob, orr, DELTA_R, GAMMA) == pytest.approx(
        2 * math.pi * GAMMA * (ob ** 2 + orr ** 2) / (4 * DELTA_R ** 2))


def test_balanced_collapse_gamma_sc_equals_gamma_over_delta_times_omega():
    # For balanced beams: Gamma_sc = (Gamma/Delta_R) * Omega (Omega in rad/s).
    omega_single = 1.0e9                       # single-beam Rabi/2pi [Hz]
    rabi = two_photon_rabi(omega_single, omega_single, DELTA_R)   # Omega/2pi [Hz]
    direct = scatter_rate(omega_single, omega_single, DELTA_R, GAMMA)
    via_rabi = scatter_rate_from_rabi(rabi, DELTA_R, GAMMA, balance=1.0)
    closed = (GAMMA / DELTA_R) * (2 * math.pi * rabi)
    assert via_rabi == pytest.approx(direct)
    assert via_rabi == pytest.approx(closed)


def test_se_per_pi_is_detuning_only_balanced():
    # P_SE(pi) = pi * Gamma / Delta_R, independent of the Rabi rate (balanced).
    assert se_probability_per_pi(DELTA_R, GAMMA) == pytest.approx(math.pi * GAMMA / DELTA_R)
    # ~0.66% at 20 GHz / 41.8 MHz -- the large-detuning SE floor
    assert se_probability_per_pi(DELTA_R, GAMMA) == pytest.approx(6.57e-3, rel=1e-2)


def test_se_per_pi_equals_gamma_sc_times_t_pi():
    rabi = 170e3
    t_pi = 1.0 / (2.0 * rabi)
    g_sc = scatter_rate_from_rabi(rabi, DELTA_R, GAMMA)
    assert g_sc * t_pi == pytest.approx(se_probability_per_pi(DELTA_R, GAMMA))


def test_scaling_with_detuning():
    # At fixed MEASURED Rabi, Gamma_sc and P_SE/pi both fall as 1/Delta_R.
    rabi = 170e3
    g1 = scatter_rate_from_rabi(rabi, DELTA_R, GAMMA)
    g2 = scatter_rate_from_rabi(rabi, 2 * DELTA_R, GAMMA)
    assert g2 == pytest.approx(g1 / 2)
    p1 = se_probability_per_pi(DELTA_R, GAMMA)
    p2 = se_probability_per_pi(2 * DELTA_R, GAMMA)
    assert p2 == pytest.approx(p1 / 2)


def test_imbalance_increases_scatter_min_at_balanced():
    rabi = 170e3
    balanced = scatter_rate_from_rabi(rabi, DELTA_R, GAMMA, balance=1.0)
    for r in (0.5, 2.0, 0.55 / 0.998):       # the data's IR/green power ratio region
        assert scatter_rate_from_rabi(rabi, DELTA_R, GAMMA, balance=r) > balanced


def test_contrast_decay_rate():
    rabi = 170e3
    assert contrast_decay_rate(rabi, DELTA_R, GAMMA) == pytest.approx(
        CONTRAST_DECAY_FACTOR * scatter_rate_from_rabi(rabi, DELTA_R, GAMMA))


def test_differential_stark_shift_leading_order():
    rabi = 170e3
    assert differential_stark_shift(rabi, OMEGA_HF, DELTA_R) == pytest.approx(
        (OMEGA_HF / DELTA_R) * rabi)


def test_flip_probability_resonant_no_scatter_is_sin2():
    rabi = 170e3
    t_pi = 1.0 / (2.0 * rabi)
    assert flip_probability(0.0, rabi) == pytest.approx(0.0)
    assert flip_probability(t_pi, rabi) == pytest.approx(1.0)        # full pi flop
    assert flip_probability(2 * t_pi, rabi) == pytest.approx(0.0, abs=1e-9)  # 2pi back


def test_flip_probability_scatter_settles_to_half():
    rabi = 170e3
    g = 1.0e5                                  # heavy scatter for a clear settle
    late = flip_probability(50e-6, rabi, gamma_sc_hz=g)
    assert late == pytest.approx(0.5, abs=0.05)


def test_flip_probability_stark_caps_amplitude():
    rabi = 170e3
    delta = 0.5 * rabi                         # large detuning for a visible cap
    t_pi = 1.0 / (2.0 * rabi)
    eff = math.hypot(rabi, delta)
    # peak of the capped flop is Omega^2/Omega_eff^2 < 1, reached at t = 1/(2 eff)
    peak = flip_probability(1.0 / (2.0 * eff), rabi, stark_detuning_hz=delta)
    assert peak == pytest.approx(rabi ** 2 / eff ** 2, rel=1e-6)
    assert peak < 1.0
    assert flip_probability(t_pi, rabi, stark_detuning_hz=delta) < flip_probability(t_pi, rabi)


def test_from_ledger_consumes_inputs_and_numbers_sane():
    ledger = Ledger.load()
    rs = RamanScatter.from_ledger(ledger)
    assert abs(rs.delta_r) == pytest.approx(20.0e9)   # SIGNED red (-20 GHz); rates use |.|
    assert rs.delta_r < 0
    assert rs.gamma == pytest.approx(41.8e6)
    # SE/pi ~ 0.66%; at 170 kHz the scattering coherence time is sub-ms (hundreds of us)
    assert rs.se_per_pi() == pytest.approx(6.57e-3, rel=2e-2)
    g_sc = rs.scatter_rate(170e3)
    assert 1.0e3 < g_sc < 5.0e3
    assert 1e-4 < 1.0 / g_sc < 1e-2            # 0.1-10 ms coherence floor
    # differential AC-Stark ~ (1.79/20) * 170 kHz ~ 15 kHz; SIGN follows the red detuning (<0)
    assert rs.stark_detuning(170e3) == pytest.approx(-15.2e3, rel=5e-2)


def test_from_ledger_object_flip_matches_free_functions():
    ledger = Ledger.load()
    rs = RamanScatter.from_ledger(ledger)
    rabi, t = 170e3, 3e-6
    g = rs.scatter_rate(rabi)
    d = rs.stark_detuning(rabi)
    assert rs.flip_probability(t, rabi) == pytest.approx(
        flip_probability(t, rabi, gamma_sc_hz=g, stark_detuning_hz=d))
