"""
Tests for the Rabi-rate engine — the damped-cosine fit and the Rabi relations.

Run:  pytest spike/
"""
import math

import pytest

from spike.datfile import DatFile
from spike.engines.rabi import fit_rabi, pi_time_from_rabi, rabi_from_pi_time

_DUR = "sources/data/microwave/13_28_39_15_06_2026.dat"


def test_rabi_pi_time_relations():
    assert rabi_from_pi_time(1e-5) == pytest.approx(5e4)        # t_pi 10 us -> 50 kHz
    assert pi_time_from_rabi(5e4) == pytest.approx(1e-5)
    assert pi_time_from_rabi(rabi_from_pi_time(7e-6)) == pytest.approx(7e-6)


def test_fit_recovers_synthetic_flop():
    f_true = 0.055   # MHz (55 kHz)
    t = [0.5 * k for k in range(40)]                            # 0..19.5 us
    y = [3.0 + 2.5 * math.cos(2 * math.pi * f_true * ti) for ti in t]
    fit = fit_rabi(t, y, [0.05] * len(t))
    assert fit["freq_hz"] == pytest.approx(55e3, rel=2e-2)
    assert fit["t_pi_us"] == pytest.approx(1 / (2 * f_true), rel=2e-2)


def test_fit_recovers_damped_flop():
    f_true, g_true = 0.06, 0.02
    t = [0.5 * k for k in range(60)]
    y = [3.0 + 2.5 * math.exp(-g_true * ti) * math.cos(2 * math.pi * f_true * ti) for ti in t]
    fit = fit_rabi(t, y, [0.05] * len(t))
    assert fit["freq_hz"] == pytest.approx(60e3, rel=2e-2)
    assert fit["decay_per_s"] == pytest.approx(g_true * 1e6, rel=0.5)   # decay is loosely constrained


def test_fit_kalis_duration_scan():
    dur = DatFile(_DUR)
    t, y, s = dur.signal()
    fit = fit_rabi(t, y, s)
    # ~53 kHz from the ~18.8 us min-to-min period; ~10% below the Doerr 59.45 kHz
    assert fit["freq_hz"] == pytest.approx(53.3e3, rel=3e-2)
    assert 8.5 < fit["t_pi_us"] < 10.0
