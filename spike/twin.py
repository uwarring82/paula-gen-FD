"""
Integrated experiment-cycle twin: prepare -> drive -> detect, on one spin STATE,
simulated per run so the count cloud emerges, plus the analytic ensemble average.

The microwave drive carries an AC-Zeeman shift and (quasi-static) dephasing. KEY:
the AC-Zeeman shift exists only while the MW is ON (the pulses) — during a Ramsey
FREE gap (MW off) the spin precesses at the BARE detuning. So the Rabi (MW-on)
resonance includes the AC-Zeeman shift, while the Ramsey fringe (MW-off free
precession) oscillates at the bare detuning delta_set; when the drive is tuned to
the Rabi (MW-on) resonance this means delta_set = acz, so the fringe frequency
THEN equals the AC-Zeeman shift, and the fringe's contrast decay gives the
dephasing T2*. That asymmetry is exactly what a Rabi-vs-Ramsey comparison exploits.

Per run a quasi-static detuning offset delta_noise ~ N(0, sigma_delta) is drawn
(sigma_delta = sqrt(2)/(2 pi T2*)); the spin is evolved, projectively measured
(QPN), and a detection count is drawn from Poisson(mu_bright) for |down> or
Poisson(mu_dark) for |up>. Detection levels come from the cooling scatter rate x
the collection efficiency (Friedenauer 2010: 5.6 per-mille) x the detection time.

Pure Python (no numpy/scipy); randomness via a seeded stdlib RNG (reproducible).
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .engines.detection import expected_bright_counts
from .engines.spin import Bloch, free, p_up, prepare, pulse
from .linalg import solve


@dataclass
class MWModel:
    """Microwave qubit model. Frequencies in Hz, times in s unless noted."""
    rabi_hz: float                       # Omega/2pi
    f0_bare_hz: float = 0.0              # bare transition (levels prediction); scans are relative to it
    acz_hz: float = 0.0                 # AC-Zeeman shift (present only with MW ON)
    delta_set_hz: float = 0.0           # f_drive - f0_bare (the SET detuning; = acz if tuned to the Rabi resonance)
    t2star_s: float = float("inf")      # dephasing time (quasi-static Gaussian)
    eps_prep: float = 0.0               # residual |up> after preparation
    mu_bright: float = 6.0              # mean bright (|down>) counts
    mu_dark: float = 0.05               # mean dark (|up>) counts
    n_shots: int = 75                   # runs per scan point
    pulse_us: float = 50.0              # fixed pulse length for frequency spectroscopy (seq_rabi_freq)
    ramsey_phase: float = 0.0           # phase of the 2nd Ramsey pi/2 pulse [rad]
    depump_bright: float = 0.0          # bright->dark depumping during detection (Gamma*t_det); low-count tail
    leak_dark: float = 0.0              # dark->bright leak into the cycle (Gamma*t_det); high-count tail

    @property
    def sigma_delta_hz(self) -> float:
        return math.sqrt(2.0) / (2.0 * math.pi * self.t2star_s) if math.isfinite(self.t2star_s) else 0.0


# --- sequences: (spin, scan_value, model, delta_noise) -> spin -------------
def seq_rabi_time(s: Bloch, t_us: float, m: MWModel, dn: float) -> Bloch:
    """Rabi flop vs pulse duration (MW on throughout)."""
    det = m.delta_set_hz - m.acz_hz - dn          # detuning from the MW-on resonance
    return pulse(s, m.rabi_hz, det, t_us * 1e-6, 0.0)


def seq_rabi_freq(s: Bloch, f_mhz: float, m: MWModel, dn: float) -> Bloch:
    """Rabi spectroscopy vs drive frequency (MHz), fixed pulse length m.pulse_us.
    Detuning from the MW-on resonance f0_bare + acz."""
    det = (f_mhz * 1e6 - m.f0_bare_hz) - m.acz_hz - dn
    return pulse(s, m.rabi_hz, det, m.pulse_us * 1e-6, 0.0)


def seq_ramsey_time(s: Bloch, tau_us: float, m: MWModel, dn: float) -> Bloch:
    """Ramsey: pi/2 -- free(tau) -- pi/2. AC-Zeeman acts only during the pulses, so
    the free precession (hence the fringe) runs at delta_set; it equals the
    AC-Zeeman shift only when the drive is on the Rabi resonance (delta_set = acz).
    Two equal-phase pi/2 pulses add to a pi, so P_up(tau=0) = 1 (fringe = (1+cos)/2)."""
    det_pulse = m.delta_set_hz - m.acz_hz - dn    # MW on
    det_free = m.delta_set_hz - dn                # MW off -> no AC-Zeeman
    t_pi2 = 1.0 / (4.0 * m.rabi_hz)
    s = pulse(s, m.rabi_hz, det_pulse, t_pi2, 0.0)
    s = free(s, det_free, tau_us * 1e-6)
    s = pulse(s, m.rabi_hz, det_pulse, t_pi2, m.ramsey_phase)
    return s


def make_seq_ramsey_freq(tau_us: float):
    """Factory: a Ramsey FREQUENCY scan at fixed free time tau_us. Scans the drive
    frequency f (MHz). During the pi/2 pulses (MW on) the detuning is from the
    AC-Zeeman-shifted resonance f0_bare + acz; during the free gap (MW off) it is
    from the BARE resonance f0_bare. So the fringe comb (free-precession phase
    2 pi (f - f0_bare) tau) is centred on f0_bare with spacing 1/tau, while the
    pulse envelope peaks at f0_bare + acz."""
    def seq(s: Bloch, f_mhz: float, m: MWModel, dn: float) -> Bloch:
        det_pulse = (f_mhz * 1e6 - m.f0_bare_hz) - m.acz_hz - dn
        det_free = (f_mhz * 1e6 - m.f0_bare_hz) - dn
        t_pi2 = 1.0 / (4.0 * m.rabi_hz)
        s = pulse(s, m.rabi_hz, det_pulse, t_pi2, 0.0)
        s = free(s, det_free, tau_us * 1e-6)
        s = pulse(s, m.rabi_hz, det_pulse, t_pi2, m.ramsey_phase)
        return s
    return seq


# --- ensemble (analytic) + per-run Monte Carlo -----------------------------
def ensemble_p_up(seq_fn, scan_val, m: MWModel, n_grid: int = 41) -> float:
    """Ensemble-averaged P(|up>): integrate the per-run P_up over the quasi-static
    Gaussian detuning noise (the dephasing). n_grid points over +-3 sigma; the grid
    must resolve cos(2 pi * delta_noise * tau) at the largest tau (default 41 is
    ample for sigma_delta * tau_max <~ 1)."""
    sig = m.sigma_delta_hz
    if sig <= 0.0:
        return p_up(seq_fn(prepare(m.eps_prep), scan_val, m, 0.0))
    num = den = 0.0
    for i in range(n_grid):
        dn = sig * (-3.0 + 6.0 * i / (n_grid - 1))
        w = math.exp(-dn * dn / (2.0 * sig * sig))
        num += w * p_up(seq_fn(prepare(m.eps_prep), scan_val, m, dn))
        den += w
    return num / den


def _poisson(mu: float, rng: random.Random) -> int:
    """Knuth's Poisson sampler."""
    if mu <= 0.0:
        return 0
    el, k, p = math.exp(-mu), 0, 1.0
    while True:
        k += 1
        p *= rng.random()
        if p <= el:
            return k - 1


def _detect_count(lam_start: float, lam_other: float, decay: float, rng: random.Random) -> int:
    """Sample one detection count: Poisson(lam_start), but if the state SWITCHES to
    lam_other at fraction u~Exp(decay) of the window (u<1), Poisson(lam_start*u +
    lam_other*(1-u)). decay = Gamma*t_det; decay<=0 -> pure Poisson(lam_start).
    (Matches engines.detection.transition_count_pmf.)"""
    if decay <= 0.0:
        return _poisson(lam_start, rng)
    u = -math.log(rng.random()) / decay
    if u >= 1.0:
        return _poisson(lam_start, rng)                # survived the window
    return _poisson(lam_start * u + lam_other * (1.0 - u), rng)


def simulate_counts(seq_fn, scan_val, m: MWModel, rng: random.Random):
    """One scan point -> the list of per-run detection counts (the count cloud).
    The detection step is realistic: a bright |down> ion may DEPUMP to dark
    (m.depump_bright -> low-count tail) and a dark |up> ion may LEAK into the cycling
    transition (m.leak_dark -> high-count tail); both default to 0 (pure Poisson)."""
    out = []
    for _ in range(m.n_shots):
        dn = rng.gauss(0.0, m.sigma_delta_hz) if m.sigma_delta_hz > 0 else 0.0
        s = seq_fn(prepare(m.eps_prep), scan_val, m, dn)
        if rng.random() < p_up(s):                     # projective measurement (QPN) -> dark
            out.append(_detect_count(m.mu_dark, m.mu_bright, m.leak_dark, rng))
        else:                                          # bright
            out.append(_detect_count(m.mu_bright, m.mu_dark, m.depump_bright, rng))
    return out


def simulate_scan(seq_fn, scan_vals, m: MWModel, seed: int = 0):
    """-> [(scan_val, [per-run counts])] mirroring a .dat scan."""
    rng = random.Random(seed)
    return [(x, simulate_counts(seq_fn, x, m, rng)) for x in scan_vals]


def detection_levels(scatter_rate_hz: float, detection_eff: float, t_det_s: float,
                     dark_fraction: float = 0.01):
    """(mu_bright, mu_dark) from first principles: bright = R_scatter * eta * t_det,
    where eta is the TOTAL photon detection efficiency (solid angle x filter/mirror
    transmission x PMT QE; Friedenauer 5.6e-3, mg_detection_efficiency_25mg), not
    the bare solid-angle collection. dark ~ a small background fraction of bright.
    Delegates the bright-count formula to engines.detection.expected_bright_counts."""
    mu_b = expected_bright_counts(scatter_rate_hz, detection_eff, t_det_s)
    return mu_b, dark_fraction * mu_b


# --- inference: fit a Ramsey fringe (cosine x Gaussian envelope) ------------
def fit_ramsey(tau_us, y, sigma=None):
    """Fit y = c + e^{-(t/T2)^2}(a cos 2pi f t + b sin 2pi f t) by a grid over
    (f, T2) with an exact weighted linear solve for (c, a, b). Returns the fringe
    frequency (Hz; this is the free-precession detuning delta_set, which equals the
    AC-Zeeman shift when the drive is on the Rabi resonance), T2* (us), contrast and
    offset. t in us."""
    n = len(tau_us)
    if sigma is None:
        sigma = [1.0] * n
    tmax = max(tau_us) or 1.0
    dt = tmax / (n - 1) if n > 1 else tmax
    f_lo, f_hi = max(2e-4, 0.25 / tmax), 0.95 / (2.0 * dt)      # MHz (cycles/us)

    def _eval(f, t2):
        ata = [[0.0] * 3 for _ in range(3)]
        aty = [0.0, 0.0, 0.0]
        basis = []
        for ti, yi, si in zip(tau_us, y, sigma):
            env = math.exp(-(ti / t2) ** 2)
            ph = 2.0 * math.pi * f * ti
            r = (1.0, env * math.cos(ph), env * math.sin(ph))
            w = 1.0 / (si * si) if si > 0 else 1.0
            basis.append((r, yi, w))
            for i in range(3):
                aty[i] += w * r[i] * yi
                for j in range(3):
                    ata[i][j] += w * r[i] * r[j]
        try:
            c, a, b = solve(ata, aty)
        except Exception:
            return None
        chi2 = sum(w * (yi - (c * r[0] + a * r[1] + b * r[2])) ** 2 for r, yi, w in basis)
        return chi2, c, a, b

    best = None
    nf, df = 200, (f_hi - f_lo) / 200
    t2_step = 0.25 * tmax
    for kf in range(nf + 1):
        f = f_lo + df * kf
        for kt in range(25):
            t2 = (0.25 + 0.25 * kt) * tmax
            ev = _eval(f, t2)
            if ev and (best is None or ev[0] < best[0]):
                best = (ev[0], f, t2, ev[1], ev[2], ev[3])
    _, f0, t20, *_ = best                                       # refine in (f, T2)
    for kf in range(-5, 6):
        f = f0 + df * kf / 5.0
        for kt in range(-5, 6):
            t2 = max(0.05 * tmax, t20 + t2_step * kt / 5.0)
            ev = _eval(f, t2)
            if ev and ev[0] < best[0]:
                best = (ev[0], f, t2, ev[1], ev[2], ev[3])
    chi2, f, t2, c, a, b = best
    return {"freq_hz": f * 1e6, "t2_us": t2, "contrast": 2.0 * math.hypot(a, b),
            "offset": c, "chi2": chi2}
