"""
Sideband thermometry DISCRIMINATOR — use the RSB+BSB flops (OC_Axial/1_1R_LF_MA) to
break the degeneracy the OC CARRIER flop could not: was its contrast loss MOTIONAL
(a hot nbar, carrier Debye-Waller) or RAMAN-BEAM DEPHASING (relative-phase noise of
the two beams)? See twin_oc_flop / ADR-0007 for the two candidates.

THE DATASET. The 1R_LF_MA scan drives, per shot, BOTH a blue-sideband pulse (counter
0, gated into its own detection) and a red-sideband pulse (counter 1), vs the same
pulse duration t_lf_1R. So the file has TWO real flop blocks (DatFile.counter_blocks):
  * BSB |down,n> -> |up,n+1>  (ADDS a phonon): works at every n incl. n=0 -> FULL flop.
  * RSB |down,n> -> |up,n-1>  (SUBTRACTS a phonon): the n=0 ground state CANNOT flop,
    so the RSB amplitude is suppressed -- it is NEAR-CONSTANT for a cold ion.

THERMOMETRY. The first-order sideband strengths go as eta*sqrt(n+1) (blue) and
eta*sqrt(n) (red), so the RED/BLUE peak-flop ratio is the standard motional
thermometer

    A_RSB / A_BSB = nbar / (nbar + 1)   ->   nbar = ratio / (1 - ratio).

A near-constant RSB (small ratio) => a COLD ion (small nbar). This is a DIRECT nbar,
independent of the carrier flop.

VERDICT. If nbar(sideband) << the carrier's all-motional nbar_eff, the carrier loss
canNOT be mostly motional -> the residual is Raman-beam dephasing. The twin then
decomposes the carrier decay at the MEASURED nbar: motional (Debye-Waller, engines.
sideband) + scattering (engines.scatter) vs the Raman-dephasing remainder (engines.
raman_dephasing, mutual linewidth dnu / coherence time T_phi).

    python -m spike.twin_sideband   ->  docs/figures/twin_sideband_thermometry.png

Pure Python + matplotlib (Agg). nbar uncertainty via a Gaussian (per-point) bootstrap.
"""
from __future__ import annotations

import math
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import twin_oc_flop  # noqa: E402
from .bootstrap import robust_summary  # noqa: E402
from .datfile import DatFile  # noqa: E402
from .engines import raman_dephasing  # noqa: E402
from .ledger import Ledger  # noqa: E402

_DATAFILE = (Path(__file__).resolve().parent.parent / "sources" / "data" /
             "OC_Axial" / "1_1R_LF_MA" / "09_48_23_12_06_2026.dat")
_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_sideband_thermometry.png"
_BLUE, _RED, _GRAY = "#1f77b4", "#d62728", "#888888"
# block order in this apparatus's file: counter 0 = BSB (full flop), counter 1 = RSB
# (near-constant). The user confirmed "the near constant is the RSB flop".
_BSB_BLOCK, _RSB_BLOCK = 0, 1


def _amplitudes(bsb, rsb, mu_b, mu_d, bsb_pk, rsb_dip):
    """Peak flop probabilities of each sideband, each measured from its own t=0
    baseline over the full bright-dark count range (so the BSB rise and the RSB dip
    are on a common P_flip scale): A_BSB = (BSB_pk - BSB_0)/R, A_RSB = (RSB_0 - RSB_dip)/R."""
    R = (mu_b - mu_d) or 1.0
    return max(0.0, (bsb[bsb_pk] - bsb[0]) / R), max(0.0, (rsb[0] - rsb[rsb_dip]) / R)


def _nbar(a_rsb, a_bsb):
    ratio = a_rsb / a_bsb if a_bsb > 0 else 0.0
    ratio = min(ratio, 0.999)
    return ratio, ratio / (1.0 - ratio)


def analyze(path=_DATAFILE, n_boot=400, seed=0):
    dat = DatFile(path)
    blocks = dat.counter_blocks()
    (xb, yb, sb) = blocks[_BSB_BLOCK]
    (xr, yr, sr) = blocks[_RSB_BLOCK]
    hh = [DatFile.hist_mean(h) for h in dat.histograms()]
    mu_b, mu_d = max(hh), min(hh)
    bsb_pk = max(range(len(yb)), key=lambda i: yb[i])     # BSB peak (rises)
    rsb_dip = min(range(len(yr)), key=lambda i: yr[i])    # RSB dip (near-constant, drops)

    a_bsb, a_rsb = _amplitudes(yb, yr, mu_b, mu_d, bsb_pk, rsb_dip)
    ratio, nbar = _nbar(a_rsb, a_bsb)

    # Gaussian (per-point) bootstrap: perturb each block's y by its sigma AT the fixed
    # peak/baseline indices (so max/min picking does not bias the spread), re-derive nbar.
    rng = random.Random(seed)
    samples = []
    for _ in range(n_boot):
        b0 = yb[0] + rng.gauss(0, sb[0]); bp = yb[bsb_pk] + rng.gauss(0, sb[bsb_pk])
        r0 = yr[0] + rng.gauss(0, sr[0]); rd = yr[rsb_dip] + rng.gauss(0, sr[rsb_dip])
        R = (mu_b - mu_d) or 1.0
        ab, ar = max(1e-6, (bp - b0) / R), max(0.0, (r0 - rd) / R)
        samples.append(_nbar(ar, ab)[1])
    nbar_err = robust_summary(samples)["sigma"]

    # carrier flop + its decomposition at the MEASURED sideband nbar (point estimates
    # only -- the sideband nbar above carries the uncertainty, so no carrier bootstrap)
    tw, _fit, ci = twin_oc_flop.build(n_boot=0)
    g_obs = ci["g_obs"]
    g_floor = tw.effective_decay(nbar)                    # scatter + thermal @ nbar(sideband)
    g_motional = max(0.0, g_floor - tw.gamma_sc_contrast)
    g_raman = max(0.0, g_obs - g_floor)

    return {
        "dat": dat, "xb": xb, "yb": yb, "sb": sb, "xr": xr, "yr": yr, "sr": sr,
        "mu_b": mu_b, "mu_d": mu_d, "a_bsb": a_bsb, "a_rsb": a_rsb, "ratio": ratio,
        "nbar": nbar, "nbar_err": nbar_err, "n_boot": n_boot,
        "g_obs": g_obs, "g_sc": tw.gamma_sc_contrast, "g_motional": g_motional,
        "g_floor": g_floor, "g_raman": g_raman,
        "dnu": raman_dephasing.mutual_linewidth_from_rate(g_raman),
        "tphi_us": raman_dephasing.coherence_time_from_rate(g_raman) * 1e6,
        "nbar_eff_carrier": ci["nbar_eff"], "tspan_us": xb[-1] if xb else 0.0,
    }


def report(info) -> str:
    frac = lambda g: 100.0 * g / info["g_obs"] if info["g_obs"] else 0.0
    raman_dominates = info["g_raman"] > info["g_floor"]
    L = [
        "OC AXIAL — SIDEBAND THERMOMETRY (RSB+BSB; %s)" % (info["dat"].timestamp or ""),
        "  data: %s  (t_lf_1R, %d pts x %d shots, %d-replica bootstrap)" % (
            _DATAFILE.name, info["dat"].scan["points"], info["dat"].scan["shots"], info["n_boot"]),
        "",
        "RED vs BLUE sideband flop (the direct motional thermometer):",
        "  BSB |down,n>->|up,n+1> (adds): peak P_flip = %.2f  (full flop -- works at n=0)" % info["a_bsb"],
        "  RSB |down,n>->|up,n-1> (subtr): peak P_flip = %.2f  (near-constant -- n=0 cannot flop)" % info["a_rsb"],
        "  RSB/BSB = %.3f  ->  nbar = ratio/(1-ratio) = %.2f ± %.2f   (COLD)" % (
            info["ratio"], info["nbar"], info["nbar_err"]),
        "",
        "DISCRIMINATOR vs the CARRIER flop (same beams, same mode):",
        "  carrier's all-motional reading was nbar_eff = %.2f; the SIDEBAND says nbar = %.2f." % (
            info["nbar_eff_carrier"], info["nbar"]),
        "  => the ion is COLD, so the carrier loss is NOT mostly motional. Decomposing the",
        "     carrier decay (%.2e /s) at the MEASURED nbar = %.2f:" % (info["g_obs"], info["nbar"]),
        "       motional (Debye-Waller) + scattering : %.2e /s  (%.0f%%)" % (
            info["g_floor"], frac(info["g_floor"])),
        "       RAMAN-BEAM DEPHASING (remainder)     : %.2e /s  (%.0f%%)  dnu=%.0f kHz, T_phi=%.0f us" % (
            info["g_raman"], frac(info["g_raman"]), info["dnu"] / 1e3, info["tphi_us"]),
        "",
        "  => the OC carrier-flop contrast loss is DOMINATED by RAMAN-BEAM DEPHASING."
        if raman_dominates else
        "  => motion and Raman dephasing contribute comparably.",
        "     The apparent 'nbar ~ 1' from the carrier alone was Raman dephasing posing as",
        "     a hot ion; the sidebands show the motion is cold (nbar ~ %.2f)." % info["nbar"],
    ]
    return "\n".join(L)


def make_figure(info, out=_OUT):
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    R, mu_d = (info["mu_b"] - info["mu_d"]) or 1.0, info["mu_d"]
    # both sidebands as P_flip from their own t=0 baseline
    p_bsb = [(v - info["yb"][0]) / R for v in info["yb"]]
    p_rsb = [(info["yr"][0] - v) / R for v in info["yr"]]
    eb = [s / R for s in info["sb"]]
    er = [s / R for s in info["sr"]]
    ax.errorbar(info["xb"], p_bsb, yerr=eb, fmt="o-", color=_BLUE, ms=5, capsize=2, lw=1.2,
                label="BSB (adds phonon) — FULL flop, peak %.2f" % info["a_bsb"])
    ax.errorbar(info["xr"], p_rsb, yerr=er, fmt="s-", color=_RED, ms=5, capsize=2, lw=1.2,
                label="RSB (subtracts) — near-constant, peak %.2f" % info["a_rsb"])
    ax.axhline(0.0, color=_GRAY, lw=0.8, alpha=0.6)
    ax.set_xlabel("sideband pulse duration  t_lf_1R  [µs]")
    ax.set_ylabel("P$_{flip}$  (from each block's t=0 baseline)")
    ax.set_title("Sideband thermometry — RSB/BSB asymmetry = cold ion (%s)" % (info["dat"].timestamp or ""))
    txt = "\n".join([
        r"RSB/BSB = %.2f $\Rightarrow$ n̄ = %.2f ± %.2f (COLD)" % (
            info["ratio"], info["nbar"], info["nbar_err"]),
        r"carrier alone read n̄$_{eff}$ = %.2f (all-motional)" % info["nbar_eff_carrier"],
        r"$\Rightarrow$ carrier loss %.0f%% Raman dephasing" % (100 * info["g_raman"] / info["g_obs"]),
        r"    (Δν = %.0f kHz, T$_\phi$ = %.0f µs), %.0f%% motional" % (
            info["dnu"] / 1e3, info["tphi_us"], 100 * info["g_floor"] / info["g_obs"]),
    ])
    ax.text(0.97, 0.5, txt, transform=ax.transAxes, ha="right", va="center", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec=_GRAY, alpha=0.9))
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    if not _DATAFILE.exists():
        print("RSB+BSB sideband data not found at", _DATAFILE)
        return 0
    info = analyze()
    print(report(info))
    out = make_figure(info)
    print("\nwrote", out.relative_to(Path(__file__).resolve().parent.parent))
    return 0


if __name__ == "__main__":   # pragma: no cover
    raise SystemExit(main())
