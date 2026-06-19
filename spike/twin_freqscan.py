"""
Frequency-domain test case: Rabi vs Ramsey SPECTROSCOPY on the integrated twin.

A microwave-on Rabi frequency scan (fixed pulse) gives a single fluorescence dip
PULLED to the AC-Zeeman-shifted resonance f0_bare + acz. A Ramsey frequency scan
at fixed free time tau gives a fringe COMB whose teeth sit at f0_bare + n/tau_eff
(the free precession, MW off, is referenced to the BARE resonance) under an
envelope centred on f0_bare + acz. Three free times tau = 100, 300, 600 us show
the fringes sharpening as ~1/tau (with the finite-pulse correction
tau_eff = tau + 1/(pi*Omega)) and the dephasing T2* shrinking the contrast at long
tau. The offset between the Rabi dip and the Ramsey comb centre is the AC-Zeeman
shift; the multi-tau ladder (coarse tau disambiguates, fine tau refines) pins the
bare resonance — a robust-phase-estimation-style hierarchy.

    python -m spike.twin_freqscan  ->  docs/figures/twin_freqscan_rabi_ramsey.png
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .linalg import solve  # noqa: E402
from .twin import (  # noqa: E402
    MWModel,
    detection_levels,
    ensemble_p_up,
    make_seq_ramsey_freq,
    seq_rabi_freq,
    simulate_scan,
)

_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_freqscan_rabi_ramsey.png"
_BLUE, _RED, _GRAY = "#1f77b4", "#d62728", "#888888"
_TAUS = [100.0, 300.0, 600.0]            # us
_F0 = 1775.6e6                            # bare |3,+3>-|2,+2> resonance (levels engine)
_WIN_RABI, _NPTS_RABI = 25.0, 101         # Rabi: +-25 kHz (~1.25 x Omega), 101 points
_K_FRINGES, _NPTS_RM = 2.5, 73            # Ramsey window = +-2.5 fringe spacings (resolution-matched)
_N_SEED = 100                            # Monte-Carlo replicas (heavy-tailed centre -> need many + robust stats)


def _model():
    mu_b, mu_d = detection_levels(scatter_rate_hz=4.5e7, detection_eff=5.6e-3, t_det_s=25e-6)
    return MWModel(rabi_hz=20e3, f0_bare_hz=_F0, acz_hz=3000.0, delta_set_hz=0.0,
                   t2star_s=800e-6, eps_prep=0.02, mu_bright=mu_b, mu_dark=max(mu_d, 0.1),
                   n_shots=60, pulse_us=25.0)


def _dets(win, npts):
    return [-win + 2 * win * k / (npts - 1) for k in range(npts)]


def _f_mhz(det_khz):
    return (_F0 + det_khz * 1e3) / 1e6


def _means_in_det(scan):
    return [(f * 1e6 - _F0) / 1e3 for f, _ in scan], [sum(c) / len(c) for _, c in scan]


def _parab_min(xs, ys, i):
    """Sub-grid minimum near index i via a 3-point parabola (equal spacing)."""
    if i <= 0 or i >= len(xs) - 1:
        return xs[i]
    denom = ys[i - 1] - 2 * ys[i] + ys[i + 1]
    if denom <= 0:
        return xs[i]
    return xs[i] + 0.5 * (xs[i] - xs[i - 1]) * (ys[i - 1] - ys[i + 1]) / denom


def _quad_vertex(xs, ys, x_centre, half_width):
    """Vertex of a least-squares quadratic fit to the points within half_width of
    x_centre (robust to noise for a broad dip; falls back to x_centre)."""
    pts = [(x, y) for x, y in zip(xs, ys) if abs(x - x_centre) <= half_width]
    if len(pts) < 3:
        return x_centre
    s4 = s3 = s2 = s1 = s0 = t2 = t1 = t0 = 0.0
    for x, y in pts:
        x2 = x * x
        s4 += x2 * x2; s3 += x2 * x; s2 += x2; s1 += x; s0 += 1.0
        t2 += x2 * y; t1 += x * y; t0 += y
    try:
        a, b, _ = solve([[s4, s3, s2], [s3, s2, s1], [s2, s1, s0]], [t2, t1, t0])
    except Exception:
        return x_centre
    return -b / (2.0 * a) if a > 0 else x_centre


def _rabi_dip(xs, ys, half_width=5.0):
    """Centre of the (broad) Rabi fluorescence dip by a windowed quadratic fit."""
    i = min(range(len(ys)), key=lambda k: ys[k])
    return _quad_vertex(xs, ys, xs[i], half_width)


def _ramsey_centre(xs, ys, half_width):
    """Centre of the n=0 Ramsey fringe (sharp dip nearest 0), windowed quadratic."""
    mins = [i for i in range(1, len(ys) - 1) if ys[i] < ys[i - 1] and ys[i] <= ys[i + 1]]
    if not mins:
        return 0.0
    i0 = min(mins, key=lambda k: abs(xs[k]))
    return _quad_vertex(xs, ys, xs[i0], half_width)


def _curve_fringe_spacing(seqfn, m, lo, hi, n=801):
    """Fringe spacing (kHz) of the NOISELESS ensemble curve over [lo, hi] kHz."""
    xs = [lo + (hi - lo) * k / (n - 1) for k in range(n)]
    ys = [ensemble_p_up(seqfn, _f_mhz(d), m) for d in xs]      # maxima of P_up = fringe teeth
    peaks = [xs[i] for i in range(1, n - 1) if ys[i] > ys[i - 1] and ys[i] >= ys[i + 1]]
    if len(peaks) < 2:
        return 0.0
    return sum(peaks[i + 1] - peaks[i] for i in range(len(peaks) - 1)) / (len(peaks) - 1)


def _robust(v):
    """(median, robust sigma = half the 16-84 percentile range). The Ramsey-centre
    estimator is heavy-tailed (rare catastrophic fringe mis-identification), so a
    plain std over few seeds is unstable; the 16-84 spread ignores the ~1% tail."""
    s = sorted(v)
    n = len(s)

    def q(p):
        i = p * (n - 1)
        lo = int(i)
        return s[lo] if lo + 1 >= n else s[lo] * (1 - (i - lo)) + s[lo + 1] * (i - lo)

    return q(0.5), 0.5 * (q(0.84) - q(0.16))


def _model_pull(seqfn, m, win, n=801):
    """The envelope-pull bias (kHz): on the NOISELESS ensemble curve the bare
    resonance is 0, so the comb tooth (P_up maximum) nearest 0 is the pure
    systematic shift the AC-Zeeman envelope imprints on the fringe comb."""
    xs = [-win + 2 * win * k / (n - 1) for k in range(n)]
    ys = [ensemble_p_up(seqfn, _f_mhz(d), m) for d in xs]
    teeth = [xs[i] for i in range(1, n - 1) if ys[i] > ys[i - 1] and ys[i] >= ys[i + 1]]
    return min(teeth, key=abs) if teeth else 0.0


def main(argv=None) -> int:
    m = _model()
    tau_eff = {t: t + 1e6 / (math.pi * m.rabi_hz) for t in _TAUS}    # us
    rm_win = {t: _K_FRINGES * 1000.0 / tau_eff[t] for t in _TAUS}    # +-window (kHz), ~2.5 fringe spacings
    fvals_rabi = [_f_mhz(d) for d in _dets(_WIN_RABI, _NPTS_RABI)]
    fvals_rm = {t: [_f_mhz(d) for d in _dets(rm_win[t], _NPTS_RM)] for t in _TAUS}

    # ---- inference over replicas (Rabi dip [wide scan] + Ramsey comb centre) ----
    # The AC-Zeeman envelope shifts the whole Ramsey comb by a known per-tau bias
    # (_model_pull, from the noiseless comb x envelope); subtract it (for real data
    # one fits the comb x envelope, e.g. fit_ramsey in the frequency domain). Stats
    # are robust (median, 16-84%) because the comb-centre estimator is heavy-tailed.
    pulls = {t: _model_pull(make_seq_ramsey_freq(t), m, rm_win[t]) for t in _TAUS}
    rabi_dips, centres = [], {t: [] for t in _TAUS}
    for sd in range(_N_SEED):
        xr, yr = _means_in_det(simulate_scan(seq_rabi_freq, fvals_rabi, m, seed=10 + sd))
        rabi_dips.append(_rabi_dip(xr, yr, half_width=8.0))
        for t in _TAUS:
            xs, ys = _means_in_det(simulate_scan(make_seq_ramsey_freq(t), fvals_rm[t], m, seed=50 + sd))
            centres[t].append(_ramsey_centre(xs, ys, 0.4 * 1000.0 / tau_eff[t]) - pulls[t])
    rd_mu, rd_sd = _robust(rabi_dips)
    bare_mu, bare_sd = _robust(centres[max(_TAUS)])                  # finest tau -> bare resonance
    acz_mu = rd_mu - bare_mu
    acz_sd = (rd_sd ** 2 + bare_sd ** 2) ** 0.5
    sharper = rd_sd / max(bare_sd, 1e-9)

    print(f"Rabi-vs-Ramsey FREQUENCY scans (Omega/2pi = {m.rabi_hz / 1e3:.0f} kHz, "
          f"injected AC-Zeeman = {m.acz_hz / 1e3:.2f} kHz; N={_N_SEED} replicas, robust median/16-84%):")
    print(f"  Rabi dip (MW-on, broad ±{_WIN_RABI:.0f} kHz)      : {rd_mu:+.3f} +/- {rd_sd:.3f} kHz  (= f0_bare + acz)")
    print(f"  Ramsey comb centre, tau={int(max(_TAUS))}us (bare)  : {bare_mu:+.3f} +/- {bare_sd:.3f} kHz "
          f"(envelope-pull {pulls[max(_TAUS)]:+.3f} removed; {sharper:.0f}x sharper than Rabi)")
    print(f"  => AC-Zeeman = Rabi - Ramsey        : {acz_mu:.3f} +/- {acz_sd:.3f} kHz")
    for t in _TAUS:
        sp = _curve_fringe_spacing(make_seq_ramsey_freq(t), m, -rm_win[t], rm_win[t])
        print(f"  tau={int(t):3d}us: window ±{rm_win[t]:4.1f} kHz, model fringe spacing {sp:.2f} kHz "
              f"(1/tau={1000 / t:.2f}, 1/tau_eff={1000 / tau_eff[t]:.2f}); envelope pull {pulls[t]:+.3f}")

    # ---- plot: Rabi (wide) + three Ramsey (resolution-matched zoom) --------
    panels = [("Rabi", seq_rabi_freq, fvals_rabi, _WIN_RABI)] + \
             [(f"Ramsey {int(t)} us", make_seq_ramsey_freq(t), fvals_rm[t], rm_win[t]) for t in _TAUS]
    fig, axes = plt.subplots(len(panels), 1, figsize=(10.5, 11.5))

    def _to_counts(p):
        return m.mu_bright * (1.0 - p) + m.mu_dark * p

    for ax, (name, seqfn, fvals, win) in zip(axes, panels):
        scan = simulate_scan(seqfn, fvals, m, seed=0)
        for f, counts in scan:
            d = (f * 1e6 - _F0) / 1e3
            ax.scatter([d] * len(counts), counts, s=5, color=_BLUE, alpha=0.09, lw=0, zorder=2)
        xs_, means_ = _means_in_det(scan)
        ax.scatter(xs_, means_, s=20, color=_BLUE, edgecolor="white", lw=0.5, zorder=4,
                   label=f"per-point mean (N={m.n_shots})")
        grid_d = [-win + 2 * win * k / 600 for k in range(601)]
        ax.plot(grid_d, [_to_counts(ensemble_p_up(seqfn, _f_mhz(d), m)) for d in grid_d],
                "-", color=_RED, lw=1.5, zorder=3, label=r"$\langle$twin$\rangle$ ensemble")
        ax.axvline(0.0, color=_GRAY, ls=":", lw=1.2, zorder=1)
        ax.axvline(m.acz_hz / 1e3, color=_RED, ls="--", lw=1.0, alpha=0.7, zorder=1)
        ax.set_xlim(-win, win)
        ax.set_ylim(-0.4, 10.5)
        ax.set_ylabel("counts")
        if name == "Rabi":
            ax.set_title(rf"Rabi spectroscopy (fixed {m.pulse_us:.0f} $\mu$s pulse, ±{win:.0f} kHz) — broad "
                         rf"dip at $f_0+\delta_{{\rm ACZ}}$ = {rd_mu:+.2f} ± {rd_sd:.2f} kHz", fontsize=10)
            ax.legend(fontsize=7.5, loc="lower right", ncol=2)
        else:
            t = float(name.split()[1])
            ax.set_title(rf"Ramsey $\tau$={int(t)} $\mu$s (±{win:.0f} kHz) — fringe comb on bare $f_0$, "
                         rf"spacing $1/\tau_{{\rm eff}}$={1000 / tau_eff[t]:.2f} kHz", fontsize=10)

    axes[-1].set_xlabel(r"drive detuning from bare resonance $f_0$ (kHz)")
    axes[0].axvline(-999, color=_GRAY, ls=":", label=r"bare $f_0$ (Ramsey comb centre)")
    axes[0].axvline(-999, color=_RED, ls="--", alpha=0.7, label=r"$f_0+\delta_{\rm ACZ}$ (Rabi dip)")
    axes[0].legend(fontsize=7.5, loc="lower right", ncol=2)
    fig.suptitle(rf"Integrated twin — Rabi vs Ramsey frequency scans: AC-Zeeman = {acz_mu:.2f} ± "
                 rf"{acz_sd:.2f} kHz (injected {m.acz_hz / 1e3:.2f}); Ramsey ~{sharper:.0f}× "
                 rf"sharper, $T_2^*$={m.t2star_s * 1e6:.0f} $\mu$s", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT, dpi=130)
    plt.close(fig)
    print("wrote", _OUT)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
