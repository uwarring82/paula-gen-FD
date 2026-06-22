"""
Tests for the bootstrap uncertainty module + the fitters' opt-in error bars: robust
summary stats, parametric (Gaussian) and non-parametric (shot) resamplers, recovery
of injected parameters within the bootstrapped uncertainty, and the key PROPERTY that
the error bar GROWS with the noise.

Run:  pytest spike/
"""
import math
import random

import pytest

from spike.bootstrap import (
    gaussian_bootstrap,
    robust_summary,
    shot_bootstrap,
    summarize,
)
from spike.engines.rabi import fit_rabi
from spike.engines.tickle import fit_tickle


# --- robust summary ---------------------------------------------------------
def test_robust_summary_basics():
    s = robust_summary([1, 2, 3, 4, 5, 6, 7, 8, 9])
    assert s["median"] == pytest.approx(5.0)
    assert s["n"] == 9
    assert s["sigma"] > 0
    # robust to a wild outlier (median/percentile, not mean/std)
    s2 = robust_summary([1, 2, 3, 4, 5, 6, 7, 8, 9, 1e6])
    assert s2["median"] == pytest.approx(5.5)
    assert s2["sigma"] < 10            # 16-84 half-width ignores the outlier


def test_robust_summary_drops_nan_inf():
    s = robust_summary([1.0, 2.0, float("nan"), float("inf"), 3.0])
    assert s["n"] == 3 and s["median"] == pytest.approx(2.0)


def test_robust_summary_empty():
    s = robust_summary([])
    assert s["n"] == 0 and math.isnan(s["median"])


def test_summarize_drops_failed_runs():
    runs = [{"a": 1.0}, None, {"a": 2.0}, {"a": None}, {"a": 3.0}]
    out = summarize(runs, ("a",))
    assert out["a"]["n"] == 3 and out["a"]["median"] == pytest.approx(2.0)


# --- synthetic flop helpers -------------------------------------------------
def _flop(t, f=0.15, g=0.05, amp=0.5, off=0.5):
    return [off + math.exp(-g * ti) * amp * math.cos(2 * math.pi * f * ti) for ti in t]


def _noisy(y0, noise, seed):
    rng = random.Random(seed)
    return [yi + rng.gauss(0, noise) for yi in y0]


def _mean_fit(x, y, sigma):
    """A cheap weighted-mean 'fit' (estimates the DC offset) — exercises the bootstrap
    MACHINERY without the cost of the fit_rabi grid (which the two integration tests
    below cover). err of the mean ~ noise/sqrt(N)."""
    w = [1.0 / (s * s) if s > 0 else 1.0 for s in sigma]
    return {"offset": sum(wi * yi for wi, yi in zip(w, y)) / sum(w)}


# --- gaussian bootstrap recovers + scales with noise (cheap machinery test) --
def test_gaussian_bootstrap_recovers_within_uncertainty():
    x = list(range(40))
    y = _noisy([5.0] * 40, 0.5, seed=3)
    s = summarize(gaussian_bootstrap(_mean_fit, x, y, [0.5] * 40, n_boot=300, seed=1),
                  ("offset",))["offset"]
    assert s["sigma"] > 0
    assert abs(s["median"] - 5.0) < 4 * s["sigma"]        # truth within the error bar
    assert s["sigma"] == pytest.approx(0.5 / math.sqrt(40), rel=0.4)   # ~ noise/sqrt(N)


def test_bootstrap_error_grows_with_noise():
    x = list(range(40))
    lo = summarize(gaussian_bootstrap(_mean_fit, x, [5.0] * 40, [0.2] * 40, n_boot=300, seed=2),
                   ("offset",))["offset"]["sigma"]
    hi = summarize(gaussian_bootstrap(_mean_fit, x, [5.0] * 40, [0.8] * 40, n_boot=300, seed=2),
                   ("offset",))["offset"]["sigma"]
    assert hi > 2.5 * lo                                  # 4x the noise -> ~4x the error bar


def test_shot_bootstrap_recovers_mean():
    rng = random.Random(0)
    # two scan points, each 200 shots with means 2 and 8
    points = [(0.0, [rng.gauss(2, 1) for _ in range(200)]),
              (1.0, [rng.gauss(8, 1) for _ in range(200)])]
    runs = shot_bootstrap(lambda x, y, s: {"y0": y[0], "y1": y[1]}, points, n_boot=50, seed=5)
    s = summarize(runs, ("y0", "y1"))
    assert s["y0"]["median"] == pytest.approx(2.0, abs=0.3)
    assert s["y1"]["median"] == pytest.approx(8.0, abs=0.3)
    assert s["y0"]["sigma"] > 0                           # resampling spread ~ SEM ~ 1/sqrt(200)


# --- fitter *_err keys (integration with the real grid fitters; small n_boot) --
def test_fit_rabi_n_boot_adds_err_keys():
    t = [10 * k / 16 for k in range(17)]
    y = _noisy(_flop(t), 0.05, seed=7)
    fit = fit_rabi(t, y, [0.05] * len(t), f_lo_khz=80, f_hi_khz=260, n_boot=30, seed=0)
    for k in ("freq_hz_err", "decay_per_s_err", "t_pi_us_err", "amplitude_err", "offset_err"):
        assert k in fit and fit[k] > 0
    # no bootstrap by default -> no err keys (fast path unchanged)
    assert "freq_hz_err" not in fit_rabi(t, y, [0.05] * len(t), f_lo_khz=80, f_hi_khz=260)


def test_fit_tickle_n_boot_adds_f0_err():
    from spike.engines.tickle import _amp_shape
    texc, f0 = 200e-6, 1.30e6
    freqs = [f0 - 20e3 + 40e3 * k / 24 for k in range(25)]
    counts = _noisy([5.0 - 3.0 * _amp_shape(f, f0, texc) ** 2 for f in freqs], 0.2, seed=9)
    fit = fit_tickle(freqs, counts, texc, sigma=[0.2] * len(freqs), n_boot=30, seed=0)
    assert "f0_hz_err" in fit and fit["f0_hz_err"] > 0
    assert abs(fit["f0_hz"] - f0) < 4 * fit["f0_hz_err"]
