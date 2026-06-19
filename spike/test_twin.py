"""
Tests for the integrated experiment-cycle twin (prepare -> drive -> detect).

Run:  pytest spike/
"""
import math
import random

import pytest

from spike.twin import (
    MWModel,
    detection_levels,
    ensemble_p_up,
    fit_ramsey,
    make_seq_ramsey_freq,
    seq_rabi_freq,
    seq_rabi_time,
    seq_ramsey_time,
    simulate_counts,
)


def test_sigma_delta_from_t2():
    m = MWModel(rabi_hz=50e3, t2star_s=800e-6)
    assert m.sigma_delta_hz == pytest.approx(math.sqrt(2) / (2 * math.pi * 800e-6))
    assert MWModel(rabi_hz=50e3).sigma_delta_hz == 0.0      # T2 = inf -> no dephasing


def test_rabi_flop_on_resonance():
    m = MWModel(rabi_hz=50e3, delta_set_hz=0.0, acz_hz=0.0)
    assert ensemble_p_up(seq_rabi_time, 10.0, m) == pytest.approx(1.0, abs=1e-6)   # pi at 10 us
    assert ensemble_p_up(seq_rabi_time, 0.0, m) == pytest.approx(0.0, abs=1e-6)


def test_ramsey_fringe_oscillates_at_ac_zeeman():
    m = MWModel(rabi_hz=50e3, acz_hz=2000.0, delta_set_hz=2000.0)   # drive on the Rabi resonance
    assert ensemble_p_up(seq_ramsey_time, 0.0, m) == pytest.approx(1.0, abs=1e-6)
    assert ensemble_p_up(seq_ramsey_time, 250.0, m) == pytest.approx(0.0, abs=1e-3)  # half period (2 kHz)


def test_dephasing_decays_the_fringe():
    m = MWModel(rabi_hz=50e3, acz_hz=2000.0, delta_set_hz=2000.0, t2star_s=300e-6)
    peak = ensemble_p_up(seq_ramsey_time, 500.0, m)        # full-period peak (would be 1 without dephasing)
    assert peak == pytest.approx(0.5 * (1 + math.exp(-(500 / 300) ** 2)), abs=0.02)
    assert peak < 0.6


def test_fit_ramsey_recovers_injected_acz_and_t2():
    m = MWModel(rabi_hz=50e3, acz_hz=3000.0, delta_set_hz=3000.0, t2star_s=800e-6)
    tau = [25.0 * k for k in range(41)]
    fit = fit_ramsey(tau, [ensemble_p_up(seq_ramsey_time, t, m) for t in tau])
    assert fit["freq_hz"] == pytest.approx(3000.0, rel=0.05)
    assert fit["t2_us"] == pytest.approx(800.0, rel=0.2)


def test_ramsey_freq_comb_on_bare_resonance():
    # Ramsey frequency scan: comb teeth (P_up max) at f0_bare + n/tau_eff, nulls halfway;
    # the comb is centred on the BARE resonance (free precession, MW off), not f0 + acz.
    f0 = 1775.6e6
    m = MWModel(rabi_hz=20e3, f0_bare_hz=f0, acz_hz=2000.0, eps_prep=0.0)
    seq = make_seq_ramsey_freq(300.0)
    sp = 1.0 / (300e-6 + 1.0 / (math.pi * 20e3))           # 1/tau_eff (finite-pulse corrected)
    assert ensemble_p_up(seq, f0 / 1e6, m) > 0.95          # tooth at the bare resonance
    assert ensemble_p_up(seq, (f0 + sp) / 1e6, m) > 0.95   # next tooth one spacing away
    assert ensemble_p_up(seq, (f0 + 0.5 * sp) / 1e6, m) < 0.05   # null halfway between


def test_rabi_freq_dip_pulled_to_ac_zeeman():
    # Rabi (MW-on) spectroscopy line is pulled to f0_bare + acz, NOT the bare f0.
    f0 = 1775.6e6
    m = MWModel(rabi_hz=20e3, f0_bare_hz=f0, acz_hz=2000.0, pulse_us=25.0, eps_prep=0.0)
    on = ensemble_p_up(seq_rabi_freq, (f0 + 2000.0) / 1e6, m)    # MW-on resonance
    off = ensemble_p_up(seq_rabi_freq, (f0 - 15000.0) / 1e6, m)  # far detuned
    assert on > 0.98 and off < on - 0.3


def test_simulate_counts_cloud():
    m = MWModel(rabi_hz=50e3, mu_bright=6.0, mu_dark=0.1, n_shots=200)
    counts = simulate_counts(seq_rabi_time, 0.0, m, random.Random(0))   # t=0 -> |down> bright
    assert len(counts) == 200 and all(c >= 0 for c in counts)
    assert 5.0 < sum(counts) / len(counts) < 7.0                        # ~ Poisson(6)


def test_depump_bright_adds_low_count_tail():
    # A pure bright state (t=0 -> P_up=0) with depumping has more zero/few-count
    # shots than the pure-Poisson case (the Thomm bright-histogram tail).
    base = dict(rabi_hz=50e3, mu_bright=6.0, mu_dark=0.1, eps_prep=0.0, n_shots=4000)
    c0 = simulate_counts(seq_rabi_time, 0.0, MWModel(depump_bright=0.0, **base), random.Random(0))
    c1 = simulate_counts(seq_rabi_time, 0.0, MWModel(depump_bright=0.4, **base), random.Random(0))
    low0 = sum(1 for x in c0 if x <= 1) / len(c0)
    low1 = sum(1 for x in c1 if x <= 1) / len(c1)
    assert low1 > low0 + 0.05
    assert sum(c1) / len(c1) < sum(c0) / len(c0)         # depumping also lowers the bright mean


def test_detection_levels_from_efficiency():
    mu_b, mu_d = detection_levels(scatter_rate_hz=2.0e7, detection_eff=5.6e-3, t_det_s=5.0e-5)
    assert mu_b == pytest.approx(2.0e7 * 5.6e-3 * 5.0e-5)               # = 5.6 counts
    assert mu_d == pytest.approx(0.01 * mu_b)
