"""
Tests for the detection engine — Poisson statistics and state discrimination.

Run:  pytest spike/
"""
import pytest

from spike.datfile import DatFile
import math

import random

from spike.engines.detection import (
    detection_fidelity,
    empirical_fidelity,
    expected_bright_counts,
    mandel_q,
    ml_estimate_p_down,
    optimal_threshold,
    poisson_cdf,
    poisson_pmf,
    qpn,
    transition_count_pmf,
)
from spike.twin import _detect_count


def test_transition_pmf_normalizes_and_recovers_poisson():
    for k in range(8):                                   # decay=0 -> pure Poisson
        assert transition_count_pmf(k, 6.0, 0.1, 0.0) == pytest.approx(poisson_pmf(k, 6.0))
    assert sum(transition_count_pmf(k, 6.0, 0.1, 0.3) for k in range(80)) == pytest.approx(1.0, abs=1e-6)


def test_depumping_enhances_low_count_tail():
    tail_dp = sum(transition_count_pmf(k, 6.0, 0.1, 0.3) for k in range(2))
    tail_po = sum(poisson_pmf(k, 6.0) for k in range(2))
    assert tail_dp > tail_po + 0.03                      # bright->dark depumping -> more zero/few photons


def test_dark_leak_enhances_high_count_tail():
    hi_leak = sum(transition_count_pmf(k, 0.1, 6.0, 0.3) for k in range(3, 40))
    hi_po = sum(poisson_pmf(k, 0.1) for k in range(3, 40))
    assert hi_leak > hi_po + 0.05                        # dark->bright leak into the cycle -> high-count tail


def test_ml_readout_recovers_p_down_and_beats_naive():
    rng = random.Random(3)
    lam_b, lam_d, dp, p_true = 6.0, 0.1, 0.3, 0.65
    counts = []
    for _ in range(3000):
        if rng.random() < p_true:
            counts.append(_detect_count(lam_b, lam_d, dp, rng))      # bright (may depump)
        else:
            counts.append(_detect_count(lam_d, lam_b, 0.0, rng))     # dark
    p_ml = ml_estimate_p_down(counts, lam_b, lam_d, depump_bright=dp)
    p_naive = ml_estimate_p_down(counts, lam_b, lam_d, depump_bright=0.0)   # ignores depump
    assert p_ml == pytest.approx(p_true, abs=0.05)
    assert p_naive < p_ml                                # ignoring the tail underestimates P_down

_DUR = "sources/data/microwave/13_28_39_15_06_2026.dat"


def test_poisson_pmf_and_cdf():
    assert poisson_pmf(0, 0.0) == 1.0
    assert poisson_pmf(3, 0.0) == 0.0
    assert sum(poisson_pmf(n, 4.0) for n in range(60)) == pytest.approx(1.0, abs=1e-9)
    assert poisson_cdf(59, 4.0) == pytest.approx(1.0, abs=1e-9)
    assert poisson_pmf(0, 5.0) == pytest.approx(2.718281828 ** -5)


def test_optimal_threshold_and_fidelity_well_separated():
    nc, err = optimal_threshold(6.0, 0.01)
    assert nc == 1                            # bright if >=1 count
    assert detection_fidelity(6.0, 0.01) == pytest.approx(1.0 - err)
    assert detection_fidelity(6.0, 0.01) > 0.99


def test_fidelity_degrades_when_means_overlap():
    assert detection_fidelity(6.0, 0.01) > detection_fidelity(2.0, 1.0)


def test_mandel_q_sign():
    assert mandel_q(5.0, 5.0) == pytest.approx(0.0)     # Poissonian
    assert mandel_q(5.0, 9.0) > 0                       # super-Poissonian
    assert mandel_q(5.0, 3.0) < 0                       # sub-Poissonian


def test_empirical_fidelity_from_histograms():
    bright = {0: 1, 5: 40, 6: 34}     # 1/75 misread as dark at threshold 1
    dark = {0: 74, 1: 1}              # 1/75 misread as bright
    f = empirical_fidelity(bright, dark, threshold=1)
    assert f == pytest.approx(1.0 - 0.5 * (1 / 75 + 1 / 75))


def test_expected_bright_counts():
    assert expected_bright_counts(2.0e7, 1.0e-3, 3.0e-4) == pytest.approx(6.0)


def test_qpn_binomial_standard_error():
    assert qpn(0.5, 100) == pytest.approx((0.25 / 100) ** 0.5)   # 0.05, the worst case
    assert qpn(0.5, 75) == pytest.approx((0.25 / 75) ** 0.5)
    assert qpn(0.0, 75) == 0.0 and qpn(1.0, 75) == 0.0           # no projection noise at the poles
    assert math.isnan(qpn(0.5, 0))                                # no shots -> undefined


def test_kalis_detection_fidelity_high_and_super_poissonian():
    dur = DatFile(_DUR)
    hists = dur.histograms()
    bright = max((h for h in hists
                  if mandel_q(DatFile.hist_mean(h), DatFile.hist_variance(h)) < 0.5),
                 key=DatFile.hist_mean)
    dark = min(hists, key=DatFile.hist_mean)
    mb, md = DatFile.hist_mean(bright), DatFile.hist_mean(dark)
    nc, _ = optimal_threshold(mb, md)
    assert nc == 1
    assert detection_fidelity(mb, md) > 0.98
    assert empirical_fidelity(bright, dark, nc) > 0.97
