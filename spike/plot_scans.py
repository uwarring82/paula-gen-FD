"""
Plot PAULA microwave frequency + duration scans (the group's own data) against the
digital-twin (generalized-Rabi) prediction, with the quantum-projection-noise (QPN)
band. Works for ANY |F,mF> <-> |F',mF'> transition: the transition is parsed from
the scan parameter (e.g. fr_mw_3p1_2p2 / t_mw_3p1_2p2), and every twin parameter is
read from the .dat ion properties — sample size (exp_point), the pi-time t_mw_<tr>
-> Rabi rate, the resonance fr_mw_<tr>, the pulse duration. The .dat FORMAT is
documented in kalis2017; the data are the group's own measurements.

    python -m spike.plot_scans   ->  docs/figures/mw_<transition>_<date>_twin_vs_data.png
"""
from __future__ import annotations

import glob
import os
import random
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .datfile import DatFile  # noqa: E402
from .engines.detection import qpn  # noqa: E402
from .engines.levels import GroundStateZeeman  # noqa: E402
from .engines.rabi import fit_rabi, generalized_rabi  # noqa: E402
from .ledger import Ledger  # noqa: E402

_ROOT = Path(__file__).resolve().parent.parent
_FIGDIR = _ROOT / "docs" / "figures"
_DATADIRS = [_ROOT / "sources" / "data" / "microwave", _ROOT / "sources" / "data" / "MW_Control"]
_BLUE, _RED = "#1f77b4", "#d62728"
_MF = re.compile(r"^(\d)([pm])(\d+)$")


# --- transition parsing -----------------------------------------------------
def _base_transition(scan_name: str) -> str:
    """'fr_mw_3p1_2p2_s' -> '3p1_2p2' (the two |F,mF> tokens, dropping variant tags)."""
    body = scan_name.replace("fr_mw_", "").replace("t_mw_", "")
    mf = [t for t in body.split("_") if _MF.match(t)][:2]
    return "_".join(mf) if len(mf) == 2 else body


def _mf(tok: str):
    m = _MF.match(tok)
    return int(m.group(1)), (-1 if m.group(2) == "m" else 1) * int(m.group(3))


def _levels_key(base: str):
    """'3p1_2p2' -> (mF of the F=3 state, mF of the F=2 state) for hyperfine_transitions."""
    a, b = base.split("_")
    fa, ma = _mf(a)
    fb, mb = _mf(b)
    return (float(ma if fa == 3 else mb), float(mb if fb == 2 else ma))


def _label(base: str) -> str:
    a, b = base.split("_")

    def ket(tok):
        f, m = _mf(tok)
        return rf"|{f},{'+' if m >= 0 else '-'}{abs(m)}\rangle"

    return rf"${ket(a)} \leftrightarrow {ket(b)}$"


def _secs(path) -> int:
    h, m, s = os.path.basename(str(path)).split("_")[:3]
    return int(h) * 3600 + int(m) * 60 + int(s)


def _grid(lo, hi, n=401):
    return [lo + (hi - lo) * k / n for k in range(n + 1)]


def _discover(data_dir):
    """-> [(base, freq_path, dur_path)]: the most-complete duration scan per transition,
    paired with the nearest-in-time frequency scan (same calibration session)."""
    freqs, durs = {}, {}
    for fn in glob.glob(os.path.join(str(data_dir), "**", "*.dat"), recursive=True):
        sc = DatFile(fn).scan
        if not sc:
            continue
        base = _base_transition(sc["name"])
        (freqs if sc["name"].startswith("fr_") else durs).setdefault(base, []).append((fn, sc["points"]))
    out = []
    for base in sorted(set(freqs) & set(durs)):
        dur = max(durs[base], key=lambda x: (x[1], _secs(x[0])))[0]
        freq = min(freqs[base], key=lambda x: abs(_secs(x[0]) - _secs(dur)))[0]
        out.append((base, freq, dur))
    return out


# --- the twin-vs-data plot for one transition pair --------------------------
def plot_pair(base, freq_path, dur_path, out_path, ledger=None):
    freq, dur = DatFile(freq_path), DatFile(dur_path)
    n_shots = dur.scan["shots"]
    t_pi_us = dur.settings.get(f"t_mw_{base}")
    f0 = dur.settings.get(f"fr_mw_{base}") or freq.settings.get(f"fr_mw_{base}")
    t_pulse_s = (freq.settings.get(f"t_mw_{base}") or t_pi_us) * 1e-6

    if t_pi_us:
        rabi_hz = 1.0 / (2.0 * t_pi_us * 1e-6)             # from the configured pi-time (ion props)
    else:
        rabi_hz = fit_rabi(*dur.signal())["freq_hz"]       # fall back to the duration-scan fit
        t_pi_us = 1e6 / (2.0 * rabi_hz)

    def _panel(ax, pts, gx, gp, twin_label):
        """Raw photon counts, EVERY experimental run as a point (the count variation),
        plus the per-point mean, the twin curve and the QPN band — the last mapped into
        count space via the per-point bright/dark means (only to PLACE the twin; the
        data are untouched). QPN in counts = span * sqrt(P(1-P)/N)."""
        xs = [x for x, _ in pts]
        means = [sum(s) / len(s) for _, s in pts]
        s_br, s_dk = max(means), min(means)
        span = (s_br - s_dk) or 1.0
        xr = (max(xs) - min(xs)) or 1.0
        width = 0.55 * xr / max(1, len(xs) - 1)
        rng = random.Random(0)
        for x, shots in pts:                                  # each run = a point (jittered in x)
            ax.scatter([x + width * (rng.random() - 0.5) for _ in shots], shots,
                       s=8, color=_BLUE, alpha=0.22, lw=0, zorder=2)
        ax.scatter(xs, means, s=30, color=_BLUE, edgecolor="white", lw=0.8, zorder=4,
                   label=f"per-point mean (N={n_shots})")
        tw = [s_br - p * span for p in gp]
        band = [span * qpn(p, n_shots) for p in gp]
        ax.plot(gx, tw, "-", color=_RED, lw=1.8, label=twin_label, zorder=3)
        ax.fill_between(gx, [t - b for t, b in zip(tw, band)], [t + b for t, b in zip(tw, band)],
                        color=_RED, alpha=0.12, lw=0, label="twin QPN band", zorder=1)
        ax.scatter([], [], s=10, color=_BLUE, alpha=0.5, label=f"individual runs ({n_shots}/pt)")

    f_levels = None
    if ledger is not None:
        try:
            eng = GroundStateZeeman.from_ledger(ledger)
            f_levels = eng.hyperfine_transitions(ledger.value("b_field_zeeman_weber_25mg"))[_levels_key(base)] / 1e6
        except Exception:  # pragma: no cover
            pass

    fig, (ax_f, ax_d) = plt.subplots(1, 2, figsize=(12, 4.8))

    # frequency scan (raw counts, individual runs)
    pts_f = freq.point_shots()
    fx = [x for x, _ in pts_f]
    gx = _grid(min(fx), max(fx))
    gp = [generalized_rabi(t_pulse_s, (f - f0) * 1e6, rabi_hz) for f in gx]
    _panel(ax_f, pts_f, gx, gp, "digital twin (gen. Rabi)")
    ax_f.axvline(f0, color="gray", ls=":", lw=1, label=f"twin f0 = {f0:.3f} (ion props)")
    if f_levels:
        ax_f.axvline(f_levels, color="green", ls="--", lw=1, label=f"levels f0 = {f_levels:.3f} (Weber B)")
    _ft = (freq.timestamp or "").split()
    ax_f.set_xlabel("MW frequency (MHz)")
    ax_f.set_ylabel("photon counts")
    ax_f.set_title(f"Frequency scan  {_ft[1] if len(_ft) > 1 else ''}  (pulse {t_pulse_s * 1e6:.1f} $\\mu$s)")
    ax_f.legend(fontsize=7.5, loc="upper right")

    # duration scan (raw counts, individual runs)
    pts_d = dur.point_shots()
    dx = [x for x, _ in pts_d]
    gt = _grid(0.0, max(dx))
    gpd = [generalized_rabi(t * 1e-6, 0.0, rabi_hz) for t in gt]
    _panel(ax_d, pts_d, gt, gpd, r"digital twin $\sin^2(\Omega t/2)$")
    _dt = (dur.timestamp or "").split()
    ax_d.set_xlabel(r"MW pulse duration ($\mu$s)")
    ax_d.set_ylabel("photon counts")
    ax_d.set_title(rf"Duration scan  {_dt[1] if len(_dt) > 1 else ''}  "
                   rf"($\Omega/2\pi$={rabi_hz / 1e3:.1f} kHz, $t_\pi$={t_pi_us:.1f} $\mu$s)")
    ax_d.legend(fontsize=7.5, loc="upper right")

    date = _dt[0] if _dt else ""
    fig.suptitle(rf"$^{{25}}$Mg$^+$  {_label(base)}  microwave — PAULA, {date}   (experiment vs digital twin)",
                 fontsize=12)
    fig.text(0.5, 0.005, ".dat DAQ format: kalis2017 (data are the group's own measurement)",
             ha="center", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path.name, rabi_hz, f0, f_levels


def main(argv=None) -> int:
    ledger = None
    try:
        ledger = Ledger.load()
    except Exception:  # pragma: no cover
        pass
    made = []
    for dd in _DATADIRS:
        if not dd.exists():
            continue
        for base, fpath, dpath in _discover(dd):
            date = (DatFile(dpath).timestamp or "????-??-??").split()[0]
            out = _FIGDIR / f"mw_{base}_{date}_twin_vs_data.png"
            name, rabi_hz, f0, f_lev = plot_pair(base, fpath, dpath, out, ledger)
            off = f" ; levels f0 {f_lev:.4f} ({(f0 - f_lev) * 1e3:+.0f} kHz)" if f_lev else ""
            made.append(f"{name}: Omega/2pi={rabi_hz / 1e3:.2f} kHz, f0={f0:.4f} MHz{off}")
    print("wrote:")
    for m in made:
        print("  " + m)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
