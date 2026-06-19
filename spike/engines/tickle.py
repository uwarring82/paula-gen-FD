"""
Tickle engine — secular-frequency spectroscopy by motional excitation, following
the group's own method (Kalis et al., PRA 94, 023401 (2016); thesis kalis2017).

A finite resonant excitation pulse of duration texc drives a motional mode as a
CLASSICAL driven harmonic oscillator. The final motional amplitude is (Kalis 2016
Eq. 3)

    A_i(texc) ∝ sin([w_exc - w_i] texc / 2) / (w_exc^2 - w_i^2),

a SINC in the excitation frequency (FWHM ≈ 1/texc, Fourier-limited) — NOT a
Lorentzian. The oscillating ion Doppler-modulates the detection (BD) scattering with
modulation index beta_i = |<u_i, k_BD>| A_i, transferring carrier population to
motional sidebands at +-v*w_i with Bessel weights J_v(beta_i)^2, so the normalized
fluorescence is (Eq. 2)

    F = prod_i  sum_{v=-inf..inf}  J_v(beta_i)^2 (Gw/2)^2 / ((Delta_BD + v w_i)^2 + (Gw/2)^2),

truncated at |v| = 15. Scanning w_exc across one mode gives a sinc-driven Bessel DIP
centred on the secular frequency; its width is 1/texc (broadened by Gaussian
mode-frequency noise, kalis2017 Eq. 4.10, for texc > ~1 ms). This engine builds that
lineshape and fits a tickle .dat frequency scan for the mode frequency. Pure Python.
"""
from __future__ import annotations

import math

from ..constants import H_PLANCK
from ..linalg import solve

HBAR = H_PLANCK / (2.0 * math.pi)


def besselj(v: int, x: float, terms: int = 100) -> float:
    """J_|v|(x) via its power series (stable for |x| <~ 20, |v| <= 15). Returns the
    ABSOLUTE order |v| (so J_{-v} loses the (-1)^v sign) — harmless here as it is only
    ever used squared (J_v(beta)^2 = J_{-v}(beta)^2)."""
    v = abs(int(v))
    half = 0.5 * x
    if half == 0.0:
        return 1.0 if v == 0 else 0.0
    term = half ** v / math.factorial(v)         # k = 0
    s = term
    for k in range(1, terms):
        term *= -half * half / (k * (k + v))
        s += term
        if abs(term) <= 1e-17 * abs(s):
            break
    return s


def excitation_amplitude(f_exc_hz: float, f0_hz: float, texc_s: float, scale: float = 1.0) -> float:
    """Final motional amplitude after a texc resonant pulse (Kalis 2016 Eq. 3):
    scale * sin([w_exc - w0] texc/2) / (w_exc^2 - w0^2), with the resonant L'Hopital
    limit scale*texc/(4 w0). scale = (Q Uexc/m)|<u,El>| (units of length if SI)."""
    we, w0 = 2.0 * math.pi * f_exc_hz, 2.0 * math.pi * f0_hz
    denom = we * we - w0 * w0
    if abs(denom) < 1e-9 * (w0 * w0 + 1.0):
        return scale * texc_s / (4.0 * w0) if w0 else math.inf
    return scale * math.sin((we - w0) * texc_s / 2.0) / denom


def _amp_shape(f_exc_hz: float, f0_hz: float, texc_s: float) -> float:
    """Excitation amplitude relative to resonance (=1 at f0): the sinc shape."""
    a0 = excitation_amplitude(f0_hz, f0_hz, texc_s)
    return excitation_amplitude(f_exc_hz, f0_hz, texc_s) / a0 if a0 else 0.0


def coherent_occupation(amplitude_m: float, omega_hz: float, mass_kg: float) -> float:
    """Mean occupation of the coherent state at amplitude A (kalis2017 Eq. 2.63):
    <n> = 1/2 (m*omega*A^2/hbar - 1)."""
    return 0.5 * (mass_kg * 2.0 * math.pi * omega_hz * amplitude_m ** 2 / HBAR - 1.0)


def fluorescence(beta: float, detuning_bd_hz: float, omega_hz: float, gamma_w_hz: float,
                 v_max: int = 15) -> float:
    """Normalised single-mode fluorescence (Kalis 2016 Eq. 2, one factor of the
    product), F(beta)/F(0): sum_v J_v(beta)^2 g(Delta_BD + v omega) / g(Delta_BD),
    g(x) = (Gw/2)^2/(x^2+(Gw/2)^2). =1 at beta=0, DIPS as the modulation grows ONLY
    near Delta_BD=0 (on-resonance detection, where any sideband spread leaves the line
    peak); a red/blue Delta_BD on a slope can instead give a peak."""
    hw = 0.5 * gamma_w_hz

    def g(x):
        return hw * hw / (x * x + hw * hw)

    num = sum(besselj(v, beta) ** 2 * g(detuning_bd_hz + v * omega_hz) for v in range(-v_max, v_max + 1))
    den = g(detuning_bd_hz)
    return num / den if den > 0.0 else math.nan


def tickle_lineshape(f_exc_hz: float, f0_hz: float, texc_s: float, beta_res: float,
                     detuning_bd_hz: float, gamma_w_hz: float, baseline: float,
                     v_max: int = 15) -> float:
    """Fluorescence vs excitation frequency: baseline * F(beta_res*|sinc shape|)."""
    beta = beta_res * abs(_amp_shape(f_exc_hz, f0_hz, texc_s))
    return baseline * fluorescence(beta, detuning_bd_hz, f0_hz, gamma_w_hz, v_max)


def fit_tickle(freqs_hz, counts, texc_s, sigma=None, n_f: int = 400):
    """Fit a tickle frequency scan for the mode frequency f0, with the Kalis-2016
    finite-pulse lineshape: counts = baseline - depth * |A(f; f0, texc)|^2, where
    A is the SINC excitation amplitude (Eq. 3; FWHM ~ 1/texc, NOT a Lorentzian). The
    dip profile is the leading-order modulation (fluorescence loss ~ motional energy ~
    A^2); the full Bessel saturation (fluorescence()) only sets the absolute depth,
    which depends on the detection sensitivity (Delta_BD, Gw, beam waist) and is left
    free as `depth`. Grid over f0 with an EXACT linear solve for (baseline, depth),
    then refine. Returns {f0_hz, depth, baseline, fwhm_hz, chi2}. f0 (the secular
    frequency) is the robust deliverable (< 1 kHz scatter)."""
    n = len(freqs_hz)
    if sigma is None:
        sigma = [1.0] * n
    fmin, fmax = min(freqs_hz), max(freqs_hz)
    span = (fmax - fmin) or 1.0

    def _eval(f0):
        a00 = a01 = a11 = b0 = b1 = 0.0
        rows = []
        for f, c, s in zip(freqs_hz, counts, sigma):
            s2 = _amp_shape(f, f0, texc_s) ** 2                 # |sinc amplitude|^2
            w = 1.0 / (s * s) if s > 0 else 1.0
            rows.append((s2, c, w))
            a00 += w
            a01 += -w * s2
            a11 += w * s2 * s2
            b0 += w * c
            b1 += -w * c * s2
        try:
            base, depth = solve([[a00, a01], [a01, a11]], [b0, b1])
        except Exception:
            return None
        if depth < 0.0:                                        # a tickle is fluorescence LOSS:
            base, depth = b0 / a00, 0.0                        # constrain depth>=0 (no spurious peak)
        chi2 = sum(w * (c - (base - depth * s2)) ** 2 for s2, c, w in rows)
        return chi2, base, depth

    best = None
    for i in range(n_f + 1):
        f0 = fmin + span * i / n_f
        ev = _eval(f0)
        if ev and (best is None or ev[0] < best[0]):
            best = (ev[0], f0, ev[1], ev[2])
    for i in range(-8, 9):                                     # refine f0
        f0 = best[1] + (span / n_f) * i / 8.0
        ev = _eval(f0)
        if ev and ev[0] < best[0]:
            best = (ev[0], f0, ev[1], ev[2])
    chi2, f0, base, depth = best
    # significance: F-test of the 3-parameter dip vs a flat line (is a dip resolved?)
    ws = [1.0 / (s * s) if s > 0 else 1.0 for s in sigma]
    ybar = sum(wi * c for wi, c in zip(ws, counts)) / sum(ws)
    chi2_flat = sum(wi * (c - ybar) ** 2 for wi, c in zip(ws, counts))
    dof = max(1, n - 3)
    f_stat = ((chi2_flat - chi2) / 2.0) / (chi2 / dof) if chi2 > 0 else 0.0
    step = span / (n - 1) if n > 1 else span
    edge_pinned = (f0 - fmin) < 1.5 * step or (fmax - f0) < 1.5 * step   # resonance not bracketed
    return {"f0_hz": f0, "depth": depth, "baseline": base, "fwhm_hz": 1.0 / texc_s,
            "chi2": chi2, "f_stat": f_stat, "edge_pinned": edge_pinned,
            "resolved": f_stat > 4.0 and depth > 0.0 and not edge_pinned}
