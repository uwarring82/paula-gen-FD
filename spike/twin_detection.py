"""
Realistic detection: the depumping tail + the OPTIONAL maximum-likelihood readout.

During detection the BD beam pumps a fraction of BRIGHT ions out of the cycling
transition, so the bright photon-count histogram is a Poisson core PLUS a tail of
zero/few-photon events (observed on PAULA; cf. Thomm 2021). The raw counts stay the
primary observable; when a state probability is wanted, a maximum-likelihood readout
that MODELS the depumping recovers P_down accurately, whereas a fixed threshold (or
mean-count normalisation) is biased low because the bright tail masquerades as dark.

This demo uses Thomm's measured single-ion levels lambda_down = 2.682, lambda_up =
0.036 counts per 30 us and a representative depumping (Gamma*t_det = 0.3):

    python -m spike.twin_detection  ->  docs/figures/twin_detection_depumping.png
"""
from __future__ import annotations

import random
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .engines.detection import (  # noqa: E402
    ml_estimate_p_down,
    optimal_threshold,
    poisson_pmf,
    transition_count_pmf,
)
from .twin import (  # noqa: E402
    MWModel,
    _detect_count,
    ensemble_p_up,
    seq_rabi_time,
    simulate_scan,
)

_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_detection_depumping.png"
_LAM_B, _LAM_D, _DEPUMP = 2.682, 0.036, 0.3      # Thomm 2021 lambda_down/up; representative depumping
_BLUE, _RED, _GREEN, _GRAY = "#1f77b4", "#d62728", "#2ca02c", "#888888"


def main(argv=None) -> int:
    # ---- Panel A: bright-state count histogram with depumping ----------------
    rng = random.Random(0)
    n_hist = 2500                                  # Thomm's shot count
    bright = [_detect_count(_LAM_B, _LAM_D, _DEPUMP, rng) for _ in range(n_hist)]
    hist = Counter(bright)
    kmax = max(bright)
    ks = list(range(kmax + 1))
    frac = [hist.get(k, 0) / n_hist for k in ks]
    pure = [poisson_pmf(k, _LAM_B) for k in ks]
    model = [transition_count_pmf(k, _LAM_B, _LAM_D, _DEPUMP) for k in ks]
    tail_obs = sum(frac[k] for k in range(2))
    tail_pure = sum(pure[k] for k in range(2))

    # ---- Panel B: a Rabi flop read out three ways ----------------------------
    m = MWModel(rabi_hz=50e3, mu_bright=_LAM_B, mu_dark=_LAM_D, depump_bright=_DEPUMP,
                eps_prep=0.0, n_shots=300)
    ts = [1.0 * k for k in range(41)]              # 0..40 us
    scan = simulate_scan(seq_rabi_time, ts, m, seed=1)
    nc = optimal_threshold(_LAM_B, _LAM_D)[0]
    p_true = [1.0 - ensemble_p_up(seq_rabi_time, t, m) for t in ts]       # spin P_down (no detection error)
    p_ml = [ml_estimate_p_down(c, _LAM_B, _LAM_D, depump_bright=_DEPUMP) for _, c in scan]
    p_thr = [sum(1 for x in c if x >= nc) / len(c) for _, c in scan]      # threshold (ignores depump)

    def _rmse(a, b):
        return (sum((x - y) ** 2 for x, y in zip(a, b)) / len(a)) ** 0.5

    print("Realistic detection (Thomm levels lam_down=2.682, lam_up=0.036, depump Gamma*t_det=0.3):")
    print(f"  bright low-count tail P(k<=1): observed {tail_obs:.3f} vs pure Poisson {tail_pure:.3f}")
    print(f"  Rabi-flop P_down readout RMSE vs truth: ML {_rmse(p_ml, p_true):.3f}, "
          f"threshold(n_c={nc}) {_rmse(p_thr, p_true):.3f}")

    fig, (axh, axf) = plt.subplots(1, 2, figsize=(12, 4.6))

    axh.bar(ks, frac, width=0.9, color=_BLUE, alpha=0.45, label=f"bright shots (N={n_hist})")
    axh.plot(ks, pure, "o--", color=_GRAY, lw=1.4, ms=4, label=r"pure Poisson($\lambda_\downarrow$)")
    axh.plot(ks, model, "s-", color=_RED, lw=1.8, ms=4, label="depumping model")
    axh.set_xlabel("photon counts $k$")
    axh.set_ylabel("probability")
    axh.set_title(rf"Bright histogram: depumping tail  ($P(k{{\leq}}1)$={tail_obs:.2f} vs "
                  rf"Poisson {tail_pure:.2f})", fontsize=10)
    axh.legend(fontsize=8)

    axf.plot(ts, p_true, "-", color="black", lw=1.6, label=r"true $P_\downarrow$ (spin)")
    axf.plot(ts, p_ml, "o", color=_GREEN, ms=4, label=f"ML readout (RMSE {_rmse(p_ml, p_true):.3f})")
    axf.plot(ts, p_thr, "x", color=_RED, ms=5, label=f"threshold $n_c$={nc} (RMSE {_rmse(p_thr, p_true):.3f})")
    axf.set_xlabel(r"MW pulse duration ($\mu$s)")
    axf.set_ylabel(r"$P_\downarrow$")
    axf.set_ylim(-0.05, 1.08)
    axf.set_title("Optional ML readout recovers $P_\\downarrow$; threshold biased low by the tail", fontsize=10)
    axf.legend(fontsize=8, loc="upper right")

    fig.suptitle("Realistic detection — depumping during BD readout + optional maximum-likelihood "
                 "state readout (Thomm 2021)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT, dpi=130)
    plt.close(fig)
    print("wrote", _OUT)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
