"""
Integrated-twin test case: a Rabi-vs-Ramsey comparison INFERS the AC-Zeeman shift
and the dephasing T2* that the twin was given.

The twin runs the full cycle (prepare -> drive -> detect) per shot, so the count
cloud emerges. The Rabi flop (MW on) pins the Rabi rate Omega and is essentially
blind to the small AC-Zeeman shift and the slow dephasing; the Ramsey fringe
(free precession, MW off) oscillates at the AC-Zeeman shift and its Gaussian
contrast decay gives T2*. Fitting both recovers (Omega, AC-Zeeman, T2*).

    python -m spike.twin_demo   ->  docs/figures/twin_rabi_ramsey_inference.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .engines.rabi import fit_rabi  # noqa: E402
from .twin import (  # noqa: E402
    MWModel,
    detection_levels,
    ensemble_p_up,
    fit_ramsey,
    seq_rabi_time,
    seq_ramsey_time,
    simulate_scan,
)

_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_rabi_ramsey_inference.png"
_BLUE, _RED, _GREEN = "#1f77b4", "#d62728", "#2ca02c"
_N_SEED = 16                       # Monte-Carlo replicas -> inference uncertainty


def _to_counts(p_up, mu_b, mu_d):
    return mu_b * (1.0 - p_up) + mu_d * p_up


def _means(scan):
    return [x for x, _ in scan], [sum(c) / len(c) for _, c in scan]


def _stats(xs):
    mu = sum(xs) / len(xs)
    sd = (sum((x - mu) ** 2 for x in xs) / max(1, len(xs) - 1)) ** 0.5
    return mu, sd


def _infer_one(t_rabi, tau, m, seed):
    """One Monte-Carlo replica: simulate Rabi + Ramsey, return (Omega, acz, T2*)."""
    xr, yr = _means(simulate_scan(seq_rabi_time, t_rabi, m, seed=100 + seed))
    br, dr = max(yr), min(yr)
    omega = fit_rabi(xr, [(br - v) / (br - dr) for v in yr])["freq_hz"]
    xs, ys = _means(simulate_scan(seq_ramsey_time, tau, m, seed=200 + seed))
    bs, ds = max(ys), min(ys)
    fm = fit_ramsey(xs, [(bs - v) / (bs - ds) for v in ys])
    return omega, fm["freq_hz"], fm["t2_us"]


def main(argv=None) -> int:
    # ---- the "true" apparatus the twin stands in for --------------------
    # Detection levels are NOT hand-set: bright counts follow from the cooling
    # scatter rate x the Friedenauer detection efficiency x the detection window.
    mu_b, mu_d = detection_levels(scatter_rate_hz=4.5e7, detection_eff=5.6e-3, t_det_s=25e-6)
    m = MWModel(rabi_hz=50e3, acz_hz=3000.0, delta_set_hz=3000.0, t2star_s=800e-6,
                eps_prep=0.02, mu_bright=mu_b, mu_dark=max(mu_d, 0.1), n_shots=75)

    t_rabi = [2.0 * k for k in range(21)]                   # 0..40 us
    tau = [25.0 * k for k in range(41)]                     # 0..1000 us

    # ---- inference over N Monte-Carlo replicas (so we quote an uncertainty) ----
    reps = [_infer_one(t_rabi, tau, m, sd) for sd in range(_N_SEED)]
    om_mu, om_sd = _stats([r[0] for r in reps])
    acz_mu, acz_sd = _stats([r[1] for r in reps])
    t2_mu, t2_sd = _stats([r[2] for r in reps])

    print(f"Rabi-vs-Ramsey inference (integrated twin, N={_N_SEED} Monte-Carlo replicas):")
    print(f"  injected : Omega/2pi = {m.rabi_hz / 1e3:6.2f} kHz | AC-Zeeman = {m.acz_hz / 1e3:5.2f} kHz | "
          f"T2* = {m.t2star_s * 1e6:5.0f} us")
    print(f"  inferred : Omega/2pi = {om_mu / 1e3:5.2f} +/- {om_sd / 1e3:.2f} kHz (Rabi) | "
          f"AC-Zeeman = {acz_mu / 1e3:.2f} +/- {acz_sd / 1e3:.2f} kHz | "
          f"T2* = {t2_mu:.0f} +/- {t2_sd:.0f} us (Ramsey)")
    print(f"  detection: mu_bright = {m.mu_bright:.2f}, mu_dark = {m.mu_dark:.2f} counts "
          f"(R_scatter 45 MHz x eta 5.6e-3 x t_det 25 us; Friedenauer)")

    # one representative replica for the figure (clouds + means)
    rabi = simulate_scan(seq_rabi_time, t_rabi, m, seed=100)
    rams = simulate_scan(seq_ramsey_time, tau, m, seed=200)

    # ---- plot -----------------------------------------------------------
    fig, (ax_r, ax_m) = plt.subplots(1, 2, figsize=(12, 4.8))

    def _panel(ax, scan, gx, gp, xlabel, twin_label):
        for x, counts in scan:
            ax.scatter([x] * len(counts), counts, s=8, color=_BLUE, alpha=0.16, lw=0, zorder=2)
        xs_, means_ = _means(scan)
        ax.scatter(xs_, means_, s=30, color=_BLUE, edgecolor="white", lw=0.8, zorder=4,
                   label=f"per-point mean (N={m.n_shots})")
        ax.plot(gx, [_to_counts(p, m.mu_bright, m.mu_dark) for p in gp], "-", color=_RED, lw=1.8,
                label=twin_label, zorder=3)
        ax.scatter([], [], s=10, color=_BLUE, alpha=0.5, label=f"individual runs ({m.n_shots}/pt)")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("photon counts")
        ax.legend(fontsize=7.5, loc="upper right")

    grid_t = [40.0 * k / 300 for k in range(301)]
    _panel(ax_r, rabi, grid_t, [ensemble_p_up(seq_rabi_time, t, m) for t in grid_t],
           r"MW pulse duration ($\mu$s)", r"twin $\langle$Rabi$\rangle$")
    ax_r.set_title(rf"Rabi flop  →  $\Omega/2\pi$ = {om_mu / 1e3:.2f} ± {om_sd / 1e3:.2f} kHz  "
                   r"(AC-Zeeman/dephasing hidden)", fontsize=10)

    grid_tau = [1000.0 * k / 300 for k in range(301)]
    _panel(ax_m, rams, grid_tau, [ensemble_p_up(seq_ramsey_time, t, m) for t in grid_tau],
           r"Ramsey free time $\tau$ ($\mu$s)", r"twin $\langle$Ramsey$\rangle$")
    ax_m.set_title(rf"Ramsey fringe  →  AC-Zeeman = {acz_mu / 1e3:.2f} ± {acz_sd / 1e3:.2f} kHz,  "
                   rf"$T_2^*$ = {t2_mu:.0f} ± {t2_sd:.0f} $\mu$s", fontsize=10)

    fig.suptitle(rf"Integrated twin (prepare→drive→detect, per shot): inferring AC-Zeeman + dephasing "
                 rf"from Rabi vs Ramsey  ($N$={_N_SEED} replicas)", fontsize=12)
    fig.text(0.5, 0.005, rf"injected: $\Omega/2\pi$ {m.rabi_hz / 1e3:.0f} kHz, AC-Zeeman {m.acz_hz / 1e3:.2f} "
             rf"kHz, $T_2^*$ {m.t2star_s * 1e6:.0f} $\mu$s; detection $\mu_b$={m.mu_bright:.1f}/$\mu_d$="
             rf"{m.mu_dark:.1f} from Friedenauer $\eta$=5.6×10⁻³", ha="center", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT, dpi=130)
    plt.close(fig)
    print("wrote", _OUT)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
