"""
Detection engine — state-dependent fluorescence readout: photon-count statistics
and bright/dark discrimination.

On the detection beam (BDX, near resonance) the bright qubit state |down> scatters
photons while |up> stays dark, giving per-shot count distributions that are
(ideally) Poissonian with means mu_bright, mu_dark. A count threshold n_c assigns
each shot; the readout fidelity is

    F = 1 - 1/2 [ P(count < n_c | bright) + P(count >= n_c | dark) ],

maximised over n_c. This engine predicts the optimal threshold and fidelity from
the two means (Poisson model) and measures them from the empirical per-shot
histograms (DatFile.histograms()). For the kalis2017 example the Poisson model
(0.992) overestimates the empirical fidelity (~0.97): the gap is one-sided, almost
entirely the BRIGHT histogram's excess of LOW-count shots (a few shots at 0-2
counts), i.e. bright-state loss / depumping during the detection window — NOT a
dark-channel effect (the dark/reference channel matches Poisson) and NOT captured
by the Mandel Q = var/mean - 1 (which measures overall variance, here near-zero).
The bright rate itself follows from the cooling engine's BDX scatter rate times the
(apparatus) photon-collection efficiency and detection time.
"""
from __future__ import annotations

import math


def poisson_pmf(n: int, mu: float) -> float:
    """P(count = n) for a Poisson distribution of mean mu."""
    if mu <= 0.0:
        return 1.0 if n == 0 else 0.0
    return math.exp(-mu + n * math.log(mu) - math.lgamma(n + 1))


def poisson_cdf(n: int, mu: float) -> float:
    """P(count <= n)."""
    if n < 0:
        return 0.0
    return sum(poisson_pmf(k, mu) for k in range(0, n + 1))


def optimal_threshold(mu_bright: float, mu_dark: float):
    """Threshold n_c (declare 'bright' if count >= n_c) minimising the mean
    misassignment error, for two Poissonians. Returns (n_c, error)."""
    nmax = int(mu_bright + 8.0 * math.sqrt(mu_bright + 1.0)) + 3
    best = None
    for nc in range(0, nmax + 1):
        e_bright = poisson_cdf(nc - 1, mu_bright)        # bright read as dark
        e_dark = 1.0 - poisson_cdf(nc - 1, mu_dark)      # dark read as bright
        err = 0.5 * (e_bright + e_dark)
        if best is None or err < best[1]:
            best = (nc, err)
    return best


def detection_fidelity(mu_bright: float, mu_dark: float, threshold: int | None = None) -> float:
    """Readout fidelity F = 1 - 1/2[P(<n_c|bright) + P(>=n_c|dark)] (Poisson).
    If threshold is None, the optimal threshold is used."""
    if threshold is None:
        return 1.0 - optimal_threshold(mu_bright, mu_dark)[1]
    e_bright = poisson_cdf(threshold - 1, mu_bright)
    e_dark = 1.0 - poisson_cdf(threshold - 1, mu_dark)
    return 1.0 - 0.5 * (e_bright + e_dark)


def mandel_q(mean: float, variance: float) -> float:
    """Mandel Q = var/mean - 1. Q=0 Poissonian, Q>0 super-Poissonian (broadened),
    Q<0 sub-Poissonian."""
    return variance / mean - 1.0 if mean > 0 else float("nan")


def empirical_fidelity(bright_hist: dict, dark_hist: dict, threshold: int) -> float:
    """Readout fidelity computed directly from two per-shot count histograms
    {count: n_shots}, with the same threshold convention (bright if count>=n_c)."""
    nb = sum(bright_hist.values())
    nd = sum(dark_hist.values())
    if nb == 0 or nd == 0:
        return float("nan")
    e_bright = sum(o for c, o in bright_hist.items() if c < threshold) / nb
    e_dark = sum(o for c, o in dark_hist.items() if c >= threshold) / nd
    return 1.0 - 0.5 * (e_bright + e_dark)


def qpn(p: float, n: int) -> float:
    """Quantum projection noise: the binomial standard error sqrt(p(1-p)/n) on a
    spin-flip probability p estimated from n projective shots — the fundamental
    measurement uncertainty for sample size n."""
    if n <= 0 or p < 0.0 or p > 1.0:
        return float("nan")
    return (p * (1.0 - p) / n) ** 0.5


def expected_bright_counts(scatter_rate_hz: float, collection_eff: float, t_det_s: float) -> float:
    """Mean detected counts of the bright state: R_scatter * eta_collection * t_det.
    collection_eff folds solid angle, optics transmission and detector QE (apparatus)."""
    return scatter_rate_hz * collection_eff * t_det_s
