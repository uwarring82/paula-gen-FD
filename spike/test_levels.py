"""
Tests for the levels engine. They pin the physics against independent limits:
zero-field splitting = 3|A|, the quadratic Zeeman coefficient, the measured
clock benchmark, stretched-state linearity, the inverted (A<0) manifold, and
that the engine refuses to consume a benchmark (the wall).

Run:  pytest spike/
"""
import math

import pytest

from spike import constants as C
from spike.engines.levels import GroundStateZeeman
from spike.ledger import Ledger

A_25MG = -5.96254376e8       # Hz, Itano & Wineland
I_25MG = 2.5
NU0 = 3 * abs(A_25MG)        # 1788.763128 MHz


def eng():
    return GroundStateZeeman(A_hz=A_25MG, I=I_25MG)


def test_zero_field_splitting_is_3A():
    assert eng().zero_field_splitting() == pytest.approx(NU0, abs=1.0)


def test_clock_at_zero_field_equals_3A():
    assert eng().clock_transition(0.0) == pytest.approx(NU0, abs=1.0)


def test_quadratic_zeeman_coefficient():
    # K ~ 2195 Hz/G^2 ; check in Hz/G^2 to ~1%
    K_per_G2 = eng().quadratic_zeeman_coeff() / 1e8
    assert K_per_G2 == pytest.approx(2195.0, rel=0.01)


def test_clock_shift_is_quadratic_at_low_field():
    e = eng()
    B = 5.5 * C.GAUSS
    shift = e.clock_transition(B) - e.clock_transition(0.0)
    K = e.quadratic_zeeman_coeff()
    assert shift == pytest.approx(K * B * B, rel=2e-3)   # quartic correction is tiny


def test_clock_at_5p5_gauss():
    # exact Breit-Rabi value (closed form, incl. the g_I term), ~1788.829549 MHz
    val = eng().clock_transition(5.5 * C.GAUSS)
    assert val == pytest.approx(1788829549.0, abs=20.0)


def test_level_structure_counts():
    lv = eng().levels(5.5 * C.GAUSS)
    assert len(lv) == 12
    f3 = sorted(mF for (F, mF) in lv if F == 3.0)
    f2 = sorted(mF for (F, mF) in lv if F == 2.0)
    assert f3 == [-3, -2, -1, 0, 1, 2, 3]
    assert f2 == [-2, -1, 0, 1, 2]


def test_inverted_manifold_F3_below_F2():
    # A < 0 -> the F=I+1/2 (=3) manifold lies BELOW F=I-1/2 (=2) at low field
    e = eng()
    assert e.energy(3.0, 0.0, 1e-5) < e.energy(2.0, 0.0, 1e-5)


def test_stretched_state_is_linear_in_field():
    e = eng()
    B = 3.0 * C.GAUSS
    e0 = e.energy(3.0, 3.0, 0.0)
    s1 = e.energy(3.0, 3.0, B) - e0
    s2 = e.energy(3.0, 3.0, 2 * B) - e0
    assert s2 == pytest.approx(2 * s1, rel=1e-9)


def test_from_ledger_reproduces_clock_benchmark():
    ledger = Ledger.load()
    e = GroundStateZeeman.from_ledger(ledger)
    B = ledger.value("b_field_quantization_freddy")
    predicted = e.clock_transition(B)
    bench = ledger.quantity("clock_transition_25mg")
    # residual should be a few kHz (the field is only stated as ~5.5 G)
    assert abs(predicted - bench.value) < 4.0e3


def test_from_ledger_refuses_a_benchmark_input():
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        GroundStateZeeman.from_ledger(ledger, a_name="clock_transition_25mg")


def test_rejects_zero_nuclear_spin():
    # 24Mg/26Mg are I=0 inputs in the ledger: must fail clearly, not KeyError
    with pytest.raises(ValueError):
        GroundStateZeeman(A_hz=A_25MG, I=0.0)
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        GroundStateZeeman.from_ledger(ledger, i_name="nuclear_spin_24mg")


def test_clock_requires_half_integer_spin():
    # integer I builds (valid level structure) but has no m_F=0 clock sublevel
    e = GroundStateZeeman(A_hz=A_25MG, I=1.0)
    with pytest.raises(ValueError):
        e.clock_transition(5.5 * C.GAUSS)


def test_rejects_nonfinite_A():
    with pytest.raises(ValueError):
        GroundStateZeeman(A_hz=float("nan"), I=2.5)


# --- ground-state hyperfine / Zeeman transition spectrum --------------------
def test_clock_is_the_00_hyperfine_transition():
    e = eng()
    B = 5.5 * C.GAUSS
    hf = e.hyperfine_transitions(B)
    assert hf[(0.0, 0.0)] == pytest.approx(e.clock_transition(B))
    assert e.clock_transition_from_spectrum(B) == pytest.approx(e.clock_transition(B))


def test_hyperfine_transition_count():
    # F=3 (mF -3..3) <-> F=2 (mF -2..2) with |dmF| <= 1  ->  15 transitions
    assert len(eng().hyperfine_transitions(5.5 * C.GAUSS)) == 15


def test_zeeman_splitting_matches_lande_g_factor():
    e = eng()
    B = 5.5 * C.GAUSS
    dz = e.zeeman_splitting(3.0, 0.0, B)
    gF3 = C.G_J_2S12 * (3 * 4 + 0.5 * 1.5 - 2.5 * 3.5) / (2 * 3 * 4)   # Landé g_F (electronic)
    assert dz == pytest.approx(gF3 * C.MU_B_OVER_H * B, rel=0.02)


def test_hyperfine_transitions_require_half_integer_I():
    with pytest.raises(ValueError):
        GroundStateZeeman(A_hz=A_25MG, I=1.0).hyperfine_transitions(5.5 * C.GAUSS)


def test_input_quantity_refuses_benchmark():
    ledger = Ledger.load()
    # the wall, at the ledger boundary: a benchmark cannot be consumed as an input
    with pytest.raises(ValueError):
        ledger.input_quantity("clock_transition_25mg")
    assert ledger.input_quantity("b_field_quantization_freddy").kind == "input"


def test_benchmark_quantity_refuses_input():
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        ledger.benchmark_quantity("hyperfine_a_constant_25mg")
    assert ledger.benchmark_quantity("clock_transition_25mg").kind == "benchmark"


def test_from_ledger_sources_g_factors_from_ledger():
    # the wall now covers the constants too: g_J/g_I come from input records
    ledger = Ledger.load()
    e = GroundStateZeeman.from_ledger(ledger)
    assert e.g_J == ledger.input_quantity("g_factor_electron_2s12").value
    assert e.g_I == ledger.input_quantity("g_factor_nuclear_25mg").value
    # and the prediction is unchanged from the constants.py values
    assert e.clock_transition(5.5 * C.GAUSS) == pytest.approx(1788829549.0, abs=20.0)
