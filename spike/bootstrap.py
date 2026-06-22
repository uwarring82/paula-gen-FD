"""
Bootstrap uncertainties for the grid-search fitters (rabi, tickle) and any pipeline
built on them — so a fitted Omega, decay rate, secular frequency, or the twin's
inverted nbar_eff carries a real error bar, not just a chi-square.

Two resamplers, both re-running the SAME fit on perturbed data and summarising the
spread of the parameters (no analytic covariance, so it works for the non-linear
grid fits and propagates cleanly through downstream inversions):

  * gaussian_bootstrap — PARAMETRIC: perturb each point y_i -> y_i + N(0, sigma_i)
    using the per-point errors the .dat already carries. Use when you have (x, y,
    sigma) but not the raw shots.
  * shot_bootstrap — NON-PARAMETRIC: resample the per-shot counts at each scan point
    with replacement and recompute (mean, standard error). Use when the raw per-shot
    histograms are available (DatFile.point_shots) — it captures the true (possibly
    non-Gaussian, low-count Poisson) sampling distribution.

Summaries are ROBUST (median + the 16-84 percentile half-width), matching the
heavy-tail-aware convention already used in twin_freqscan: a few catastrophic
grid-fit outliers (e.g. a mis-identified period) then do not blow up the error bar.

Pure Python; seeded stdlib RNG (reproducible).
"""
from __future__ import annotations

import math
import random


def _percentile(sorted_vals, p: float) -> float:
    """Linear-interpolated percentile p in [0,1] of an already-sorted list."""
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    if n == 1:
        return sorted_vals[0]
    i = p * (n - 1)
    lo = int(i)
    frac = i - lo
    return sorted_vals[lo] if lo + 1 >= n else sorted_vals[lo] * (1 - frac) + sorted_vals[lo + 1] * frac


def robust_summary(samples) -> dict:
    """Robust location/spread of a sample list: median, the 16-84 percentile half-
    width `sigma` (= (p84 - p16)/2, the 1-sigma-equivalent immune to outliers), the
    plain mean/std, and n. NaN/inf samples are dropped first."""
    vals = sorted(v for v in samples if v == v and abs(v) != float("inf"))
    n = len(vals)
    if n == 0:
        return {"median": float("nan"), "sigma": float("nan"), "mean": float("nan"),
                "std": float("nan"), "p16": float("nan"), "p84": float("nan"), "n": 0}
    median = _percentile(vals, 0.5)
    p16, p84 = _percentile(vals, 0.16), _percentile(vals, 0.84)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0
    return {"median": median, "sigma": 0.5 * (p84 - p16), "mean": mean,
            "std": math.sqrt(var), "p16": p16, "p84": p84, "n": n}


def summarize(runs, keys) -> dict:
    """Summarise a list of fit-result dicts into {key: robust_summary} for `keys`.
    Missing/None values are dropped per key (a fit that failed contributes nothing)."""
    return {k: robust_summary([r[k] for r in runs if r is not None and r.get(k) is not None])
            for k in keys}


def _mean_sem(shots):
    """(mean, standard error of the mean) of a per-point shot list."""
    n = len(shots)
    if n == 0:
        return 0.0, 1.0
    m = sum(shots) / n
    if n == 1:
        return m, 1.0
    var = sum((s - m) ** 2 for s in shots) / (n - 1)
    return m, math.sqrt(var / n)


def gaussian_bootstrap(fit_fn, x, y, sigma, n_boot: int = 200, seed: int = 0):
    """Re-fit `n_boot` times with each y_i resampled as N(y_i, sigma_i); return the
    list of result dicts (a failed/None fit becomes None). `fit_fn(x, y, sigma)` ->
    dict. Use `summarize(runs, keys)` to reduce."""
    rng = random.Random(seed)
    runs = []
    for _ in range(n_boot):
        yb = [yi + rng.gauss(0.0, si if si > 0 else 0.0) for yi, si in zip(y, sigma)]
        try:
            runs.append(fit_fn(x, yb, sigma))
        except Exception:
            runs.append(None)
    return runs


def shot_bootstrap(fit_fn, points, n_boot: int = 200, seed: int = 0):
    """Re-fit `n_boot` times, resampling the per-shot counts at each scan point with
    replacement and recomputing (mean, SEM). `points` = [(x, [shot, ...]), ...];
    `fit_fn(x, y, sigma)` -> dict. Returns the list of result dicts (None on failure)."""
    rng = random.Random(seed)
    xs = [p[0] for p in points]
    shotlists = [p[1] for p in points]
    runs = []
    for _ in range(n_boot):
        yb, sb = [], []
        for shots in shotlists:
            n = len(shots)
            if n == 0:
                yb.append(0.0)
                sb.append(1.0)
                continue
            resampled = [shots[rng.randrange(n)] for _ in range(n)]
            m, sem = _mean_sem(resampled)
            yb.append(m)
            sb.append(sem if sem > 0 else 1.0)
        try:
            runs.append(fit_fn(xs, yb, sb))
        except Exception:
            runs.append(None)
    return runs
