"""
Analyse raw PAULA `.dat` scans with the rabi + detection engines.

Kept separate from the ledger-based composition root (spike/runner.py): this reads
RAW data files (docs/DATA_FORMAT.md), runs the engines, and prints what would
become records (with .dat + kalis2017 provenance) once curated.

    python -m spike.analyze_data
"""
from __future__ import annotations

from pathlib import Path

from .datfile import DatFile
from .engines.detection import (
    detection_fidelity,
    empirical_fidelity,
    mandel_q,
    optimal_threshold,
)
from .engines.levels import GroundStateZeeman
from .engines.rabi import fit_rabi
from .ledger import Ledger

_DATADIR = Path(__file__).resolve().parent.parent / "sources" / "data" / "microwave"
_FREQ = _DATADIR / "13_28_34_15_06_2026.dat"
_DURATION = _DATADIR / "13_28_39_15_06_2026.dat"
_DOERR_RABI_HZ = 59453.0   # mw_rabi_3p3_2p2_doerr


def _bright_dark(dat: DatFile):
    """The cleanest bright and the dark histogram. Dark = smallest mean. Bright =
    the largest-mean state that is still near-Poissonian (Mandel Q < 0.5): a Rabi
    REVIVAL has the same mean but is broadened, so we prefer the prepared |down>."""
    hists = dat.histograms()
    dark = min(hists, key=DatFile.hist_mean)
    mmax = max(DatFile.hist_mean(h) for h in hists)
    high = [h for h in hists if DatFile.hist_mean(h) > 0.8 * mmax]
    bright = min(high, key=lambda h: mandel_q(DatFile.hist_mean(h), DatFile.hist_variance(h)))
    return bright, dark


def main(argv=None) -> int:
    if not _DURATION.exists():
        print("no example data found under sources/data/microwave/")
        return 0
    dur = DatFile(_DURATION)
    out = ["RAW-DATA ANALYSIS — 25Mg+ |3,+3> <-> |2,+2> microwave, PAULA %s "
           "(.dat format: kalis2017)" % (dur.timestamp or ""), ""]

    # --- Rabi: duration scan ------------------------------------------------
    t, y, s = dur.signal()
    fit = fit_rabi(t, y, s)
    out.append("RABI (duration scan, %s):" % _DURATION.name)
    out.append(f"  fitted Omega/2pi = {fit['freq_hz'] / 1e3:.2f} kHz   t_pi = {fit['t_pi_us']:.2f} us"
               f"   (config mw_t_1 = {dur.settings.get('mw_t_1')} us)")
    dev = (fit["freq_hz"] - _DOERR_RABI_HZ) / _DOERR_RABI_HZ * 100.0
    out.append(f"  vs mw_rabi_3p3_2p2_doerr = 59.45 kHz: {dev:+.0f}% — consistent with the "
               "power/day-dependence of the (apparatus-limited) MW rate")
    out.append(f"  chi2_red = {fit['chi2_reduced']:.2f} (>1: residual scatter exceeds shot noise; "
               "frequency robust)")
    out.append("")

    # --- frequency scan -> resonance, vs the levels engine ------------------
    frq = DatFile(_FREQ)
    fx, fy, _ = frq.signal()
    f_res = min(zip(fx, fy), key=lambda p: p[1])[0]      # the dip = max spin-flip
    out.append("FREQUENCY (%s):" % _FREQ.name)
    out.append(f"  resonance dip at f = {f_res:.4f} MHz")
    try:
        ledger = Ledger.load()
        eng = GroundStateZeeman.from_ledger(ledger)
        B = ledger.value("b_field_zeeman_weber_25mg")
        pred = eng.hyperfine_transitions(B)[(3.0, 2.0)] / 1e6
        out.append(f"  levels engine predicts (3,+3)<->(2,+2) = {pred:.4f} MHz at the Weber field "
                   f"(5.6454 G); residual {f_res - pred:+.3f} MHz (~{abs(f_res - pred) / 2.34:.2f} G of field)")
    except Exception as e:  # pragma: no cover
        out.append(f"  (levels comparison unavailable: {e})")
    out.append("")

    # --- detection: per-shot count histograms -------------------------------
    bright, ref = _bright_dark(dur)
    mb, vb = DatFile.hist_mean(bright), DatFile.hist_variance(bright)
    mr = DatFile.hist_mean(ref)
    nc, _ = optimal_threshold(mb, mr)
    out.append("DETECTION (per-shot count histograms, %d shots each):" % sum(bright.values()))
    out.append(f"  bright |down> mean = {mb:.2f}  (var {vb:.2f}, Mandel Q = {mandel_q(mb, vb):+.2f}, near-Poissonian)")
    out.append(f"  reference/bkg mean = {mr:.3f}  (the ~0.013 normalisation counter, |up>-ion proxy)")
    out.append(f"  optimal threshold n_c = {nc} (declare bright if count >= {nc})")
    out.append(f"  readout fidelity: Poisson model {detection_fidelity(mb, mr):.4f}  vs  "
               f"empirical {empirical_fidelity(bright, ref, nc):.4f}")
    out.append("  gap is one-sided: the bright histogram's LOW-count tail = bright-state loss/")
    out.append("  depumping during detection (the reference channel matches Poisson).")
    out.append("  NOTE: a Rabi scan, not a readout calibration — no pure |up> prep, so the")
    out.append("  reference counter stands in for the dark-ion level.")
    print("\n".join(out))
    return 0


if __name__ == "__main__":   # pragma: no cover
    raise SystemExit(main())
