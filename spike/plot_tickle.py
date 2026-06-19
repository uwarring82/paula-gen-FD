"""
Tickle (secular-excitation) spectroscopy of the three motional modes, analysed with
the Kalis-2016 finite-pulse model (sinc excitation amplitude, FWHM ~ 1/texc; engine
spike.engines.tickle). Each PDQ_*_FScan directory holds several scans (some narrow
calibration scans that MISS the mode); we fit ALL of them, keep only files with a
SIGNIFICANT, bracketed dip (engine F-test + edge guard), and report the consensus
secular frequency vs the ledger's nominal mode-frequency input.

    python -m spike.plot_tickle  ->  docs/figures/tickle_modes.png
"""
from __future__ import annotations

import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .datfile import DatFile  # noqa: E402
from .engines.tickle import _amp_shape, fit_tickle  # noqa: E402

_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "tickle_modes.png"
_BLUE, _RED, _GRAY = "#1f77b4", "#d62728", "#888888"
# (label, data dir, texc [s], ledger nominal [MHz])
_MODES = [("LF axial", "PDQ_LF_FScan", 200e-6, 1.30),
          ("MF radial", "PDQ_MF_FScan", 200e-6, 3.00),
          ("HF radial", "PDQ_HF_FScan", 100e-6, 4.50)]


def _load(fn):
    pts = DatFile(fn).point_shots()
    order = sorted(range(len(pts)), key=lambda i: pts[i][0])
    xs = [pts[i][0] for i in order]
    clouds = [pts[i][1] for i in order]
    ys = [sum(c) / len(c) for c in clouds]
    return xs, clouds, ys


def main(argv=None) -> int:
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    print("Tickle motional-mode spectroscopy (Kalis 2016 sinc model; significant-dip files only):")
    for ax, (lab, d, texc, nominal) in zip(axes, _MODES):
        files = sorted(glob.glob(str(Path("sources/data/Tickle") / d / "*.dat")))
        fits = []
        for fn in files:
            xs, clouds, ys = _load(fn)
            fit = fit_tickle([x * 1e6 for x in xs], ys, texc_s=texc)
            fits.append((fit, xs, clouds, ys))
        good = [t for t in fits if t[0]["resolved"]]

        if good:
            f0s = sorted(t[0]["f0_hz"] / 1e6 for t in good)
            mid = len(f0s) // 2
            f0 = f0s[mid] if len(f0s) % 2 else 0.5 * (f0s[mid - 1] + f0s[mid])   # median (robust)
            spread = (sum((x - sum(f0s) / len(f0s)) ** 2 for x in f0s) / len(f0s)) ** 0.5
            off = (f0 - nominal) / nominal * 100.0
            # display the resolved file closest to the median (so the curve matches f0),
            # breaking ties toward the widest scan for context
            best = min(good, key=lambda t: (abs(t[0]["f0_hz"] / 1e6 - f0), -(t[1][-1] - t[1][0])))
            note = (f"f0 = {f0:.5f} +/- {spread * 1e3:.1f} kHz ({off:+.1f}% vs {nominal:.2f}); "
                    f"{len(good)}/{len(files)} files resolve")
            title = rf"{lab}: $f_0$={f0:.4f} MHz ({off:+.1f}% vs {nominal:.2f}); {len(good)}/{len(files)} files"
        else:
            best = max(fits, key=lambda t: t[1][-1] - t[1][0])
            note = f"UNRESOLVED — no significant dip in {len(files)} files"
            title = f"{lab}: UNRESOLVED (no significant dip)"
        print(f"  {lab:9s}: {note}")

        fit, xs, clouds, ys = best
        for x, c in zip(xs, clouds):
            ax.scatter([x] * len(c), c, s=7, color=_BLUE, alpha=0.13, lw=0, zorder=2)
        ax.scatter(xs, ys, s=22, color=_BLUE, edgecolor="white", lw=0.6, zorder=4, label="per-point mean")
        gx = [xs[0] + (xs[-1] - xs[0]) * k / 400 for k in range(401)]
        gy = [fit["baseline"] - fit["depth"] * _amp_shape(g * 1e6, fit["f0_hz"], texc) ** 2 for g in gx]
        ax.plot(gx, gy, "-", color=_RED, lw=1.7, zorder=3, label="Kalis sinc fit")
        if good:
            ax.axvline(f0, color=_RED, ls="--", lw=1.0, alpha=0.7)                # the reported median
        if xs[0] <= nominal <= xs[-1]:
            ax.axvline(nominal, color=_GRAY, ls=":", lw=1.2, label="ledger nominal")
        ax.set_xlim(xs[0], xs[-1])
        ax.set_xlabel("tickle frequency (MHz)")
        ax.set_ylabel("photon counts")
        ax.set_title(title, fontsize=9.5)
        ax.legend(fontsize=7.5, loc="lower left")

    fig.suptitle("Tickle spectroscopy of the motional modes — Kalis 2016 sinc lineshape (FWHM ~ 1/t_exc), "
                 "significant-dip files only", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT, dpi=130)
    plt.close(fig)
    print("wrote", _OUT)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
