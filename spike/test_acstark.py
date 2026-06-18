"""
Tests for the AC-Stark engine — the far-detuned single-beam light shift, and the
BDD benchmark vs Hasse 2025.

Run:  pytest spike/
"""
import pytest

from spike.engines.acstark import ac_stark_shift, is_far_detuned
from spike.ledger import Ledger

GAMMA = 41.8e6


def test_formula_matches_omega2_over_4delta():
    # delta_AC = s Gamma^2/(8 delta) = Omega^2/(4 delta) with Omega^2 = s Gamma^2/2
    s, delta = 3.0, 500e6
    omega_sq = s * GAMMA ** 2 / 2.0
    assert ac_stark_shift(s, GAMMA, delta) == pytest.approx(omega_sq / (4.0 * delta))


def test_bdd_shift_matches_hasse_10MHz():
    # BDD: s=20, delta=-10*Gamma -> ~ -10.4 MHz (Hasse: ~2pi x 10 MHz)
    shift = ac_stark_shift(20.0, GAMMA, -10.0 * GAMMA)
    assert abs(shift) == pytest.approx(10.45e6, rel=1e-2)
    assert shift < 0   # red detuning -> negative


def test_shift_sign_follows_detuning():
    assert ac_stark_shift(1.0, GAMMA, +100e6) > 0
    assert ac_stark_shift(1.0, GAMMA, -100e6) < 0


def test_far_detuned_predicate():
    assert is_far_detuned(GAMMA, -10 * GAMMA)        # BDD: far
    assert not is_far_detuned(GAMMA, -0.5 * GAMMA)   # BD: near resonance
    assert not is_far_detuned(GAMMA, -0.1 * GAMMA)   # BDX: near resonance
    assert is_far_detuned(GAMMA, 20e9)               # Raman: very far


def test_bdd_from_ledger_reproduces_benchmark():
    ledger = Ledger.load()
    s = ledger.value("bdd_far_cooling_saturation")
    delta = ledger.value("bdd_far_cooling_detuning")
    gamma = ledger.value("mg_p32_natural_linewidth")
    bench = ledger.benchmark_quantity("bdd_ac_stark_shift_25mg")
    pred = abs(ac_stark_shift(s, gamma, delta))
    assert abs(pred - bench.value) < bench.sigma   # within the (power-dependent) error bar


def test_raman_beams_are_far_detuned():
    # the far-detuned formula applies to the Raman beams (Delta_R = 20 GHz)
    ledger = Ledger.load()
    d_r = ledger.value("raman_detuning_from_p32")
    assert is_far_detuned(41.8e6, d_r)
