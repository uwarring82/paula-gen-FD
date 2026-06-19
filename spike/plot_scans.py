"""
Plot the kalis2017 frequency + duration scans against the digital-twin
(generalized-Rabi) prediction, with the quantum-projection-noise (QPN) band.

EVERY twin parameter is read from the .dat ion properties: the sample size
(exp_point), the pi-time t_mw_3p3_2p2 -> Rabi frequency Omega/2pi = 1/(2 t_pi), the
operating resonance fr_mw_3p3_2p2, and the freq-scan pulse duration. Counts are
mapped to the spin-flip probability P(|up>) via the bright/dark levels from the
duration-scan fit; QPN = sqrt(P(1-P)/N).

    python -m spike.plot_scans   ->  docs/figures/kalis_twin_vs_data.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .datfile import DatFile  # noqa: E402
from .engines.detection import qpn  # noqa: E402
from .engines.levels import GroundStateZeeman  # noqa: E402
from .engines.rabi import fit_rabi, generalized_rabi  # noqa: E402
from .ledger import Ledger  # noqa: E402

_DD = Path(__file__).resolve().parent.parent / "sources" / "data" / "microwave"
_FREQ = _DD / "13_28_34_15_06_2026.dat"
_DUR = _DD / "13_28_39_15_06_2026.dat"
_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "kalis_twin_vs_data.png"

_BLUE, _RED = "#1f77b4", "#d62728"


def _grid(lo, hi, n=401):
    return [lo + (hi - lo) * k / n for k in range(n + 1)]


def main(argv=None) -> int:
    if not _DUR.exists():
        print("no kalis2017 example data found")
        return 0
    freq, dur = DatFile(_FREQ), DatFile(_DUR)

    # --- twin parameters, ALL from the ion properties ----------------------
    n_shots = dur.scan["shots"]                             # exp_point
    t_pi_us = dur.settings["t_mw_3p3_2p2"]                  # pi-time
    rabi_hz = 1.0 / (2.0 * t_pi_us * 1e-6)                  # Omega/2pi
    f0_dur = dur.settings["fr_mw_3p3_2p2"]                  # MHz
    f0_frq = freq.settings["fr_mw_3p3_2p2"]                 # MHz
    t_pulse_s = freq.settings["t_mw_3p3_2p2"] * 1e-6        # freq-scan pulse

    # counts -> P(|up>) from the duration-fit bright/dark asymptotes
    xt, yt, st = dur.signal()
    fit = fit_rabi(xt, yt, st)
    s_br, s_dk = fit["offset"] + fit["amplitude"], fit["offset"] - fit["amplitude"]
    to_p = lambda s: (s_br - s) / (s_br - s_dk)             # noqa: E731

    # levels-engine prediction of the transition (independent field cross-check)
    f_levels = None
    try:
        ledger = Ledger.load()
        eng = GroundStateZeeman.from_ledger(ledger)
        f_levels = eng.hyperfine_transitions(ledger.value("b_field_zeeman_weber_25mg"))[(3.0, 2.0)] / 1e6
    except Exception:  # pragma: no cover
        pass

    fig, (ax_f, ax_d) = plt.subplots(1, 2, figsize=(12, 4.8))

    # --- frequency scan -----------------------------------------------------
    fx, fy, _ = freq.signal()
    p_f = [to_p(y) for y in fy]
    ax_f.errorbar(fx, p_f, yerr=[qpn(p, n_shots) for p in p_f], fmt="o", ms=5, capsize=3,
                  color=_BLUE, label=f"experiment (QPN, N={n_shots})", zorder=3)
    gx = _grid(min(fx), max(fx))
    gp = [generalized_rabi(t_pulse_s, (f - f0_frq) * 1e6, rabi_hz) for f in gx]
    ax_f.plot(gx, gp, "-", color=_RED, lw=1.8, label="digital twin (gen. Rabi)", zorder=2)
    ax_f.fill_between(gx, [p - qpn(p, n_shots) for p in gp], [p + qpn(p, n_shots) for p in gp],
                      color=_RED, alpha=0.15, lw=0, label="twin QPN band", zorder=1)
    ax_f.axvline(f0_frq, color="gray", ls=":", lw=1, label=f"twin f0 = {f0_frq:.3f} (ion props)")
    if f_levels:
        ax_f.axvline(f_levels, color="green", ls="--", lw=1, label=f"levels f0 = {f_levels:.3f} (Weber B)")
    ax_f.set_xlabel("MW frequency (MHz)")
    ax_f.set_ylabel(r"spin-flip probability $P(|{\uparrow}\rangle)$")
    ax_f.set_title(f"Frequency scan  (pulse {t_pulse_s * 1e6:.1f} $\\mu$s)")
    ax_f.set_ylim(-0.12, 1.12)
    ax_f.legend(fontsize=7.5, loc="upper right")

    # --- duration scan ------------------------------------------------------
    dx, dy, _ = dur.signal()
    p_d = [to_p(y) for y in dy]
    ax_d.errorbar(dx, p_d, yerr=[qpn(p, n_shots) for p in p_d], fmt="o", ms=5, capsize=3,
                  color=_BLUE, label=f"experiment (QPN, N={n_shots})", zorder=3)
    gt = _grid(0.0, max(dx))
    gpd = [generalized_rabi(t * 1e-6, (f0_dur - f0_frq) * 1e6, rabi_hz) for t in gt]
    ax_d.plot(gt, gpd, "-", color=_RED, lw=1.8, label=r"digital twin $\sin^2(\Omega t/2)$", zorder=2)
    ax_d.fill_between(gt, [p - qpn(p, n_shots) for p in gpd], [p + qpn(p, n_shots) for p in gpd],
                      color=_RED, alpha=0.15, lw=0, label="twin QPN band", zorder=1)
    ax_d.set_xlabel(r"MW pulse duration ($\mu$s)")
    ax_d.set_ylabel(r"spin-flip probability $P(|{\uparrow}\rangle)$")
    ax_d.set_title(rf"Duration scan  ($\Omega/2\pi$={rabi_hz / 1e3:.1f} kHz from $t_\pi$={t_pi_us:.1f} $\mu$s)")
    ax_d.set_ylim(-0.12, 1.12)
    ax_d.legend(fontsize=7.5, loc="upper right")

    fig.suptitle(r"kalis2017  $|3,+3\rangle \leftrightarrow |2,+2\rangle$  microwave: experiment vs digital twin",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT, dpi=130)
    print("wrote", _OUT)
    print(f"twin (from ion properties): N={n_shots}, t_pi={t_pi_us:.2f} us -> Omega/2pi={rabi_hz / 1e3:.2f} kHz, "
          f"f0={f0_frq:.4f} MHz, pulse={t_pulse_s * 1e6:.2f} us; bright/dark = {s_br:.2f}/{s_dk:.2f} counts")
    if f_levels:
        print(f"levels f0 = {f_levels:.4f} MHz (Weber B) vs operating {f0_frq:.4f} MHz "
              f"-> {(f0_frq - f_levels) * 1e3:+.0f} kHz")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
