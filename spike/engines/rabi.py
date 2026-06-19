"""
Rabi-rate engine — extract the Rabi frequency from a duration (Rabi-flop) scan,
and the standard Rabi relations.

A resonant drive flops the bright-state population as
    y(t) = c + e^{-gamma t} (a cos(2pi f t) + b sin(2pi f t)),
with f = Omega/(2pi) the Rabi frequency, gamma a phenomenological decay, and the
amplitude sqrt(a^2+b^2). `fit_rabi` finds (f, gamma) by a grid scan with an
exact weighted linear least-squares solve for (c, a, b) at each grid point
(pure Python, via linalg.solve) — no SciPy. The pi-time is t_pi = 1/(2f).

This is the MEASURED rate from raw data; the *relative* rates across transitions
come from the drive engine's Clebsch-Gordan couplings, and the absolute magnetic-
dipole rate Omega = (mu_B/hbar) g B_mw |CG| needs the (apparatus-dependent) drive
field B_mw — provided as a calibration helper, not a from-first-principles
prediction.
"""
from __future__ import annotations

import math

from ..linalg import solve


def _ls_solve(t, y, sigma, f, gamma):
    """Weighted linear least squares for (c, a, b) in
    y = c + e^{-gamma t}(a cos 2pi f t + b sin 2pi f t); returns (c, a, b, chi2).
    t in us, f in MHz, gamma in 1/us."""
    ata = [[0.0] * 3 for _ in range(3)]
    aty = [0.0, 0.0, 0.0]
    basis = []
    for ti, yi, si in zip(t, y, sigma):
        env = math.exp(-gamma * ti)
        ph = 2.0 * math.pi * f * ti
        r = (1.0, env * math.cos(ph), env * math.sin(ph))
        w = 1.0 / (si * si) if si > 0 else 1.0
        basis.append((r, yi, w))
        for i in range(3):
            aty[i] += w * r[i] * yi
            for j in range(3):
                ata[i][j] += w * r[i] * r[j]
    c, a, b = solve(ata, aty)
    chi2 = 0.0
    for r, yi, w in basis:
        pred = c * r[0] + a * r[1] + b * r[2]
        chi2 += w * (yi - pred) ** 2
    return c, a, b, chi2


def fit_rabi(t, y, sigma=None, f_lo_khz=None, f_hi_khz=None):
    """Fit a damped Rabi flop. t in microseconds. When the frequency bounds are not
    given they are auto-derived from the time span (so a slow 7 kHz flop and a fast
    60 kHz flop both work): f_lo ~ 0.4/T_max, f_hi ~ 0.95*Nyquist (capped). Returns a
    dict with the Rabi frequency (Hz), Omega (rad/s), decay rate (1/s), pi-time (us),
    amplitude/offset, the reduced chi-square, and a grid-edge flag."""
    if sigma is None:
        sigma = [1.0] * len(t)
    tmax = max(t) if t else 1.0
    dt = tmax / (len(t) - 1) if len(t) > 1 and tmax > 0 else (tmax or 1.0)
    if f_lo_khz is None:
        f_lo_khz = max(0.5, 0.4 * 1000.0 / tmax)               # >= a few-tenths of a period over T_max
    if f_hi_khz is None:
        f_hi_khz = min(0.95 * 1000.0 / (2.0 * dt), 150.0)      # <= ~Nyquist, capped for MW rates
        if f_hi_khz <= f_lo_khz * 1.5:
            f_hi_khz = f_lo_khz * 6.0
    f_lo, f_hi = f_lo_khz * 1e-3, f_hi_khz * 1e-3              # MHz
    nf = 240
    df = (f_hi - f_lo) / nf
    g_grid = [0.005 * k for k in range(21)]                    # 0 .. 0.1 /us
    best = None
    for k in range(nf + 1):
        f = f_lo + df * k
        for g in g_grid:
            c, a, b, chi2 = _ls_solve(t, y, sigma, f, g)
            if best is None or chi2 < best[0]:
                best = (chi2, f, g, c, a, b)
    # local refinement (1/5 of a grid cell in f; 0.001/us in gamma)
    _, f0, g0, *_ = best
    for kf in range(-6, 7):
        f = f0 + df * kf / 5.0
        if f <= 0:
            continue
        for g in [max(0.0, g0 + 0.001 * kg) for kg in range(-6, 7)]:
            c, a, b, chi2 = _ls_solve(t, y, sigma, f, g)
            if chi2 < best[0]:
                best = (chi2, f, g, c, a, b)
    chi2, f, g, c, a, b = best
    dof = max(1, len(t) - 5)            # 5 free params: c, a, b, f, gamma
    edge = f <= f_lo + df or f >= f_hi - df
    return {
        "freq_hz": f * 1e6,
        "omega": 2.0 * math.pi * f * 1e6,
        "decay_per_s": g * 1e6,
        "offset": c,
        "amplitude": math.hypot(a, b),
        "t_pi_us": 1.0 / (2.0 * f),
        "chi2_reduced": chi2 / dof,      # >1 here: residual scatter exceeds shot noise
        "ndata": len(t),
        "grid_edge": edge,               # True -> optimum hit the [f_lo, f_hi] edge (widen the range)
    }


def generalized_rabi(t_s: float, detuning_hz: float, rabi_hz: float) -> float:
    """Spin-flip probability after a square pulse of duration t_s for a drive of
    Rabi frequency rabi_hz detuned by detuning_hz (the twin's forward prediction):

        P = (Omega^2 / Omega_eff^2) sin^2(Omega_eff t / 2),  Omega_eff = sqrt(Omega^2 + delta^2).

    Inputs in ordinary Hz: the 2pi cancels in the prefactor and the sin argument is
    pi * Omega_eff_hz * t (= (2pi Omega_eff) t / 2)."""
    eff = math.hypot(rabi_hz, detuning_hz)
    if eff == 0.0:
        return 0.0
    return (rabi_hz * rabi_hz) / (eff * eff) * math.sin(math.pi * eff * t_s) ** 2


def rabi_from_pi_time(t_pi_s: float) -> float:
    """Rabi frequency Omega/(2pi) [Hz] from a pi-time [s]: 1/(2 t_pi)."""
    return 1.0 / (2.0 * t_pi_s)


def pi_time_from_rabi(freq_hz: float) -> float:
    """pi-time [s] from the Rabi frequency Omega/(2pi) [Hz]: 1/(2 f)."""
    return 1.0 / (2.0 * freq_hz)
