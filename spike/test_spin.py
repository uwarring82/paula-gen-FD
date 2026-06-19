"""
Tests for the spin state + coherent operations.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.spin import free, p_up, pi_pulse_time, prepare, pulse

R = 50e3


def test_prepare_population():
    assert p_up(prepare(0.0)) == pytest.approx(0.0)       # pure |down>
    assert p_up(prepare(0.1)) == pytest.approx(0.1)       # 10% prep error


def test_pi_pi2_2pi_pulses():
    assert p_up(pulse(prepare(0), R, 0.0, pi_pulse_time(R))) == pytest.approx(1.0, abs=1e-9)
    assert p_up(pulse(prepare(0), R, 0.0, 1 / (4 * R))) == pytest.approx(0.5, abs=1e-9)
    assert p_up(pulse(prepare(0), R, 0.0, 1 / R)) == pytest.approx(0.0, abs=1e-9)   # 2pi


def test_detuned_pulse_caps_contrast():
    # detuning = Rabi -> max flip = Omega^2/(Omega^2+delta^2) = 1/2
    eff = math.hypot(R, R)
    assert p_up(pulse(prepare(0), R, R, 1 / (2 * eff))) == pytest.approx(0.5)


def test_free_precession_preserves_population_rotates_coherence():
    s = pulse(prepare(0), R, 0.0, 1 / (4 * R))             # pi/2 -> equatorial
    s2 = free(s, 1000.0, 100e-6)
    assert p_up(s2) == pytest.approx(p_up(s))              # z (population) unchanged
    assert (s2.x, s2.y) != pytest.approx((s.x, s.y))       # x,y (coherence) rotate
