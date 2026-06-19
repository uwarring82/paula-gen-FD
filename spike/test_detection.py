"""
Tests for the detection engine — Poisson statistics and state discrimination.

Run:  pytest spike/
"""
import pytest

from spike.datfile import DatFile
from spike.engines.detection import (
    detection_fidelity,
    empirical_fidelity,
    expected_bright_counts,
    mandel_q,
    optimal_threshold,
    poisson_cdf,
    poisson_pmf,
)

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
