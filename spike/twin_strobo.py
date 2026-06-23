"""
Twin of the stroboscopic OC carrier flop -- the displ=0 baseline of the "active phase
grating" (Strobo2.0/1_FlopN_3p3_2p2_PDQ_displ_strobo).

THE SEQUENCE (decoded from the .dat script): cool -> sideband-cool to |n~0> -> MW
pi-pulse (spin prep) -> tickle displacement of the lf axial mode (BUT u_displ = 0 here
-> NO displacement, motion stays in |0>) -> active_phase_grating_Laser: N = 50
stroboscopic OC pulses (B1 continuous, R2 pulsed for `delta_t`), one per period
DELTA_t = 0.769 us = 1/(1.30 MHz) = the lf MOTIONAL PERIOD, at fr_oc_strobo ~ the
qubit carrier (-40 kHz) -> detect spin. The scan is `delta_t` (the per-cycle pulse
width), 0..0.1 us.

PHYSICS (UW-confirmed). The stroboscopic OC drive is a CARRIER flop locked to the
motion: its k-vector gradient makes the flop rate sample the instantaneous motional
position -> the "phase grating". With u_displ = 0 the mean position is zero, so this
run is the GRATING'S n=0 CALIBRATION: a clean stroboscopic carrier flop on |0>,

    P_up(delta_t) = (1/2)(1 - cos(2 pi * N * Omega_strobo * delta_t)),

with N = 50 amplifying the per-cycle Rabi Omega_strobo (the total flop angle is
N * Omega_strobo * delta_t). The ground state gives a single Debye-Waller-reduced
rate Omega_strobo = Omega_0 * e^{-eta^2/2} (no collapse/revival -- that appears only
once the motion is displaced). The twin fits Omega_strobo and the contrast.

KEY: the whole train lasts a FIXED N * DELTA_t = 38.5 us (independent of delta_t), yet
the flop keeps a sizeable contrast -- much more than the CONTINUOUS carrier flop's
T_phi ~ 15 us (engines.raman_dephasing, from twin_sideband) would leave over 38.5 us
(e^{-38.5/15} ~ 8%). The stroboscopic structure DECOUPLES the slow Raman-beam
dephasing (a dynamical-decoupling benefit of the phase-grating technique).

    python -m spike.twin_strobo   ->  docs/figures/twin_strobo_flop.png

Pure Python + matplotlib (Agg). Uncertainties via the rabi-fit bootstrap.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .datfile import DatFile  # noqa: E402
from .engines.rabi import _ls_solve, fit_rabi  # noqa: E402
from .engines.sideband import Sideband  # noqa: E402
from .engines.raman_optical import RamanBeam, RamanOptics  # noqa: E402
from .engines.scatter import differential_stark_shift  # noqa: E402
from .engines.strobo_sim import strobo_detuning_scan  # noqa: E402
from .ledger import Ledger  # noqa: E402

_DATAFILE = (Path(__file__).resolve().parent.parent / "sources" / "data" / "Strobo2.0" /
             "1_FlopN_3p3_2p2_PDQ_displ_strobo" / "10_22_22_12_06_2026.dat")
_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_strobo_flop.png"
_OUT_SCAN = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_strobo_detuning_scan.png"
_OUT_ACS = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_strobo_acstark_vs_N.png"
_BLUE, _RED, _GRAY = "#1f77b4", "#d62728", "#888888"
_F_LO_KHZ, _F_HI_KHZ = 10_000.0, 60_000.0      # fast flop: ~25 cycles/us in delta_t
_TPHI_CONT_US = 15.0                            # continuous-carrier-flop T_phi (twin_sideband)


def analyze(path=_DATAFILE, n_boot=300, seed=0):
    dat = DatFile(path)
    x, y, s = dat.counter_blocks()[0]              # the strobo flop (counter 0)
    N = dat.settings.get("N_strobo_OC_lf_PDQ", 50)
    DELTA_t = dat.settings.get("DELTA_t", 0.769)   # us, = motional period
    det_khz = (dat.settings.get("fr_oc_strobo", 0) - dat.settings.get("fr_mw_3p3_2p2", 0)) * 1e3

    fit = fit_rabi(x, y, s, f_lo_khz=_F_LO_KHZ, f_hi_khz=_F_HI_KHZ, n_boot=n_boot, seed=seed)
    f_cyc = fit["freq_hz"] / 1e6                    # cycles per us (in delta_t)
    omega_strobo = fit["freq_hz"] / N              # per-cycle Rabi Omega/2pi [Hz]
    # full damped-cosine (c, a, b) at the fitted frequency, for the model curve + phase
    c, a, b, _ = _ls_solve(x, y, s, f_cyc, 0.0)

    hh = [DatFile.hist_mean(h) for h in dat.histograms()]
    mu_b, mu_d = max(hh), min(hh)
    contrast = 2.0 * math.hypot(a, b) / (mu_b - mu_d) if mu_b > mu_d else 0.0

    eta = Sideband.from_ledger(Ledger.load()).lamb_dicke("OC", "lf",
            Ledger.load().input_quantity("omega_z_axial_com_25mg").value)
    dw = math.exp(-eta * eta / 2.0)                 # ground-state Debye-Waller
    omega0 = omega_strobo / dw                      # bare per-cycle Rabi (before Debye-Waller)

    total_us = N * DELTA_t
    # dephasing: a FIXED-time contrast factor. What the CONTINUOUS flop (T_phi~15us)
    # would leave over total_us, vs what the strobo keeps -> decoupling factor.
    cont_contrast = math.exp(-total_us / _TPHI_CONT_US)
    tphi_eff_us = -total_us / math.log(contrast) if 0 < contrast < 1 else float("inf")

    return {
        "dat": dat, "x": x, "y": y, "s": s, "c": c, "a": a, "b": b, "f_cyc": f_cyc,
        "N": N, "DELTA_t": DELTA_t, "det_khz": det_khz, "total_us": total_us,
        "omega_strobo": omega_strobo, "omega_strobo_err": fit.get("freq_hz_err", 0.0) / N,
        "omega0": omega0, "eta": eta, "dw": dw, "contrast": contrast,
        "mu_b": mu_b, "mu_d": mu_d, "turns": f_cyc * (x[-1] if x else 0.0),
        "cont_contrast": cont_contrast, "tphi_eff_us": tphi_eff_us,
        "tspan_us": x[-1] if x else 0.0, "n_boot": n_boot,
    }


def report(info) -> str:
    decoupled = info["contrast"] > 2.0 * info["cont_contrast"]
    L = [
        "STROBOSCOPIC OC CARRIER FLOP — phase-grating n=0 baseline (displ off; %s)" % (
            info["dat"].timestamp or ""),
        "  data: %s  (delta_t 0..%.2f us, %d pts x %d shots, %d-replica bootstrap)" % (
            _DATAFILE.name, info["tspan_us"], info["dat"].scan["points"],
            info["dat"].scan["shots"], info["n_boot"]),
        "  strobo: N=%g OC pulses, one per DELTA_t=%.3f us (= lf motional period); "
        "fr_oc_strobo at carrier %+.0f kHz; total %.1f us." % (
            info["N"], info["DELTA_t"], info["det_khz"], info["total_us"]),
        "",
        "FLOP (vs the per-cycle pulse width delta_t):",
        "  %.1f cycles/us -> %.1f flop turns over the scan; N x Omega amplification" % (
            info["f_cyc"], info["turns"]),
        "  per-cycle Omega_strobo/2pi = freq/N = %.0f ± %.0f kHz" % (
            info["omega_strobo"] / 1e3, info["omega_strobo_err"] / 1e3),
        "  Debye-Waller (|0>, eta=%.3f): e^{-eta^2/2}=%.3f -> bare Omega_0/2pi = %.0f kHz" % (
            info["eta"], info["dw"], info["omega0"] / 1e3),
        "",
        "CONTRAST = %.2f (of the bright-dark range). The train lasts a FIXED %.1f us:" % (
            info["contrast"], info["total_us"]),
        "  a CONTINUOUS carrier flop at T_phi=%.0f us (twin_sideband) would keep only "
        "~%.0f%% over %.1f us;" % (_TPHI_CONT_US, 100 * info["cont_contrast"], info["total_us"]),
        "  the strobo keeps ~%.0f%% (effective T_phi ~ %.0f us) =>" % (
            100 * info["contrast"], info["tphi_eff_us"]),
        "  the STROBOSCOPIC structure DECOUPLES the slow Raman-beam dephasing." if decoupled
        else "  contrast consistent with the continuous-flop dephasing.",
        "",
        "  This is the GRATING's n=0 calibration (u_displ=0): a clean single-rate flop. With",
        "  a displacement the flop rate would be MODULATED by the motional position (the",
        "  position-sensitive phase grating) -> collapse/revival; that is the displ!=0 twin.",
    ]
    return "\n".join(L)


def make_figure(info, out=_OUT):
    x, y, s = info["x"], info["y"], info["s"]
    ts = [info["tspan_us"] * k / 300 for k in range(301)]
    model = [info["c"] + info["a"] * math.cos(2 * math.pi * info["f_cyc"] * t)
             + info["b"] * math.sin(2 * math.pi * info["f_cyc"] * t) for t in ts]

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.errorbar(x, y, yerr=s, fmt="o", color=_BLUE, ms=5, capsize=2, zorder=3,
                label="data (mean ± s.e.)")
    ax.plot(ts, model, "-", color=_RED, lw=1.8, zorder=4,
            label=r"strobo flop: $\Omega_{strobo}/2\pi$=%.0f kHz/cycle (×N=%g)" % (
                info["omega_strobo"] / 1e3, info["N"]))
    ax.set_xlabel("per-cycle OC pulse width  delta_t  [µs]")
    ax.set_ylabel("fluorescence counts  (∝ P$_{flip}$)")
    ax.set_title("Stroboscopic OC carrier flop — phase-grating n=0 baseline (%s)" % (
        info["dat"].timestamp or ""))
    txt = "\n".join([
        r"N=%g pulses @ DELTA_t=%.2f µs (motional period)" % (info["N"], info["DELTA_t"]),
        r"$\Omega_{strobo}$=%.0f kHz/cycle, %.1f turns over the scan" % (
            info["omega_strobo"] / 1e3, info["turns"]),
        r"contrast %.2f over %.0f µs (cont. flop: ~%.0f%%)" % (
            info["contrast"], info["total_us"], 100 * info["cont_contrast"]),
        r"$\Rightarrow$ stroboscopic dephasing-decoupling (T$_\phi^{eff}$~%.0f µs)" % info["tphi_eff_us"],
    ])
    ax.text(0.97, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec=_GRAY, alpha=0.9))
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


# --- forward SIMULATION: detuning scan of the pulse train (delta_t fixed) ----
def simulate_detuning_scan(omega_strobo_hz=4.99e5, delta_t_us=0.02, span_flf=1.6,
                           n_pts=641, F=10):
    """SIMULATE P_flip vs the drive detuning of the stroboscopic pulse train, delta_t
    FIXED at delta_t_us, all other params from the apparatus (engines.strobo_sim). The
    scan covers +- span_flf * f_lf, i.e. the carrier (delta=0) and the FIRST SIDEBANDS
    (+- the motional frequency). Because DELTA_t = the motional period, the motional
    sidebands ALIAS onto the carrier comb -> equal narrow teeth at delta = k * f_lf."""
    ledger = Ledger.load()
    eta = Sideband.from_ledger(ledger).lamb_dicke(
        "OC", "lf", ledger.input_quantity("omega_z_axial_com_25mg").value)
    dat = DatFile(_DATAFILE) if _DATAFILE.exists() else None
    f_lf = (dat.settings.get("fr_lf", 1.3001) * 1e6) if dat else \
        ledger.input_quantity("omega_z_axial_com_25mg").value
    DELTA_t = dat.settings.get("DELTA_t", 0.769172) if dat else 0.769172
    N = int(dat.settings.get("N_strobo_OC_lf_PDQ", 50)) if dat else 50

    dets = [(-span_flf + 2 * span_flf * k / (n_pts - 1)) * f_lf for k in range(n_pts)]
    P = strobo_detuning_scan(eta, omega_strobo_hz, delta_t_us, DELTA_t, N, f_lf, dets, F=F)
    peaks = [(dets[k], P[k]) for k in range(1, n_pts - 1)
             if P[k] > 0.5 and P[k] >= P[k - 1] and P[k] >= P[k + 1]]
    return {"dets": dets, "P": P, "f_lf": f_lf, "DELTA_t": DELTA_t, "N": N, "eta": eta,
            "omega_strobo": omega_strobo_hz, "delta_t_us": delta_t_us, "peaks": peaks,
            "fwhm_khz": 1e3 / (N * DELTA_t)}      # tooth width ~ 1/(N*DELTA_t)


def make_detuning_figure(sim, out=_OUT_SCAN):
    f_lf_mhz = sim["f_lf"] / 1e6
    dmhz = [d / 1e6 for d in sim["dets"]]
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.plot(dmhz, sim["P"], "-", color=_RED, lw=1.4)
    for d, _h in sim["peaks"]:
        k = round(d / sim["f_lf"])
        lab = "carrier" if k == 0 else ("%+d × f$_{lf}$" % k)
        ax.annotate(lab, (d / 1e6, 1.02), ha="center", va="bottom", fontsize=8, color=_GRAY)
    ax.set_xlabel("drive detuning from the carrier  δ  [MHz]")
    ax.set_ylabel("P$_{flip}$ (simulated)")
    ax.set_ylim(-0.03, 1.15)
    ax.set_title("Stroboscopic pulse train — simulated detuning scan (delta_t=%.2f µs, N=%g)" % (
        sim["delta_t_us"], sim["N"]))
    txt = "\n".join([
        r"strobo comb: teeth at δ = k·f$_{lf}$ (f$_{lf}$=%.3f MHz)" % f_lf_mhz,
        r"carrier (k=0) + FIRST SIDEBANDS (k=±1, ±%.2f MHz)" % f_lf_mhz,
        r"tooth width ~ 1/(N·DELTA_t) = %.0f kHz" % sim["fwhm_khz"],
        r"DELTA_t = motional period ⇒ motional sidebands alias onto the comb",
    ])
    ax.text(0.5, 0.5, txt, transform=ax.transAxes, ha="center", va="center", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec=_GRAY, alpha=0.9))
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


# --- AC-Stark systematics vs N (B1 continuous, only R2 pulsed for a pi pulse) -
def ac_stark_vs_N(Ns=(2, 5, 10, 20, 50, 100, 200, 500), omega_strobo_hz=4.99e5,
                  eta=None):
    """Estimate the (NEGATIVE) differential AC-Stark shift of the OC pi pulse vs N (the
    number of stroboscopic cycles). R2 is pulsed, its width set for a pi pulse,
    delta_t = 1/(2 N Omega_strobo); B1 is ON CONTINUOUSLY. Hence:

      * B1 acts for the WHOLE train t_seq = N*DELTA_t -> phase grows ~ N (its shift is a
        constant detuning floor, present even between pulses);
      * R2's on-time t_R2 = N*delta_t = 1/(2 Omega_strobo) is FIXED by the pi condition
        -> its AC-Stark phase is N-INDEPENDENT, but its duty delta_t/DELTA_t (its weight
        in the time-averaged shift) falls ~ 1/N.

    ABSOLUTE scale: the per-beam shifts are NEGATIVE (|2,2> is closer to P3/2 -> the
    qubit resonance moves DOWN -- the sign of the -40 kHz fr_oc_strobo offset and of our
    earlier scatter estimate). The NOMINAL powers/waists (s_B1~100, s_R2~300) over-
    predict the OBSERVED two-photon flop rate (the true measure of the intensity at the
    ion), so we RE-ANCHOR the effective intensity to it: kappa = Omega_0(observed) /
    Omega_2gamma(nominal). The Rabi-anchored shifts are reported alongside the nominal
    upper bound and the scalar (omega_HF/Delta) cross-check (a known factor-~2.7 spread
    between the scalar and the polarization-resolved differential -- an open item)."""
    L = Ledger.load()
    ro = RamanOptics.from_ledger(L)
    Isat = L.input_quantity("mg_saturation_intensity").value
    sB1 = ro.saturation(L.input_quantity("raman_b1_power_25mg").value,
                        L.input_quantity("raman_b1_waist_25mg").value, Isat)
    sR2 = ro.saturation(L.input_quantity("raman_r2_power_25mg").value,
                        L.input_quantity("raman_r2_waist_25mg").value, Isat)
    polB1 = tuple(L.input_quantity("raman_b1_polarization_25mg").value)   # (0,1,0)=pi
    polR2 = tuple(L.input_quantity("raman_r2_polarization_25mg").value)   # (.5,0,.5)=lin-perp
    B1n, R2n = RamanBeam(sB1, polB1), RamanBeam(sR2, polR2)
    dB1_nom = ro.differential_stark_hz(B1n)                           # nominal upper bound
    dR2_nom = ro.differential_stark_hz(R2n)
    # re-anchor the effective intensity to the observed bare carrier Rabi Omega_0
    if eta is None:
        eta = Sideband.from_ledger(L).lamb_dicke(
            "OC", "lf", L.input_quantity("omega_z_axial_com_25mg").value)
    om0_obs = omega_strobo_hz / math.exp(-eta * eta / 2.0)            # un-Debye-Waller bare Rabi
    om2g_nom = ro.two_photon_rabi_hz(B1n, R2n)                        # nominal-intensity prediction
    kappa = om0_obs / om2g_nom if om2g_nom else 1.0                   # effective/nominal intensity
    dB1, dR2 = dB1_nom * kappa, dR2_nom * kappa                       # RABI-ANCHORED (used below)
    d_scalar = differential_stark_shift(om0_obs, ro.omega_hf, ro.delta_p32)  # scalar cross-check
    dat = DatFile(_DATAFILE) if _DATAFILE.exists() else None
    DELTA_t = (dat.settings.get("DELTA_t", 0.769172) if dat else 0.769172) * 1e-6   # s
    Om = omega_strobo_hz
    t_R2 = 1.0 / (2.0 * Om)                           # fixed pi-pulse R2 on-time [s]
    phi_R2 = 2.0 * math.pi * dR2 * t_R2              # R2 AC-Stark phase [rad], N-independent
    rows = []
    for N in Ns:
        delta_t = 1.0 / (2.0 * N * Om)              # per-cycle R2 width for a pi pulse [s]
        duty = delta_t / DELTA_t                     # R2 fraction of each period
        t_seq = N * DELTA_t                           # full train duration [s]
        d_eff = dB1 + dR2 * duty                     # time-averaged shift over the train
        phi_B1 = 2.0 * math.pi * dB1 * t_seq         # B1 AC-Stark phase [rad] (~N)
        rows.append({"N": N, "delta_t_us": delta_t * 1e6, "duty": duty,
                     "t_seq_us": t_seq * 1e6, "d_eff_khz": d_eff / 1e3,
                     "phi_B1_cyc": phi_B1 / (2.0 * math.pi),
                     "phi_tot_cyc": (phi_B1 + phi_R2) / (2.0 * math.pi)})
    return {"dB1": dB1, "dR2": dR2, "dB1_nom": dB1_nom, "dR2_nom": dR2_nom,
            "d_scalar": d_scalar, "sB1": sB1, "sR2": sR2, "kappa": kappa,
            "om0_obs": om0_obs, "om2g_nom": om2g_nom, "t_R2_us": t_R2 * 1e6,
            "phi_R2_cyc": phi_R2 / (2.0 * math.pi), "DELTA_t_us": DELTA_t * 1e6, "rows": rows}


def report_acstark(acs) -> str:
    L = ["", "OC pi-pulse AC-Stark systematics vs N  (B1 continuous, only R2 pulsed)",
         "  differential shift delta_AC = LS(|2,2>)-LS(|3,3>) is NEGATIVE (resonance",
         "  moves down; same sign as the -40 kHz fr_oc_strobo offset & our scatter est.):",
         "  intensity anchor: nominal s_B1=%.0f, s_R2=%.0f predict Omega_2g=%.2f MHz, but the" % (
             acs["sB1"], acs["sR2"], acs["om2g_nom"] / 1e6),
         "    OBSERVED bare Rabi Omega_0=%.2f MHz -> re-anchor kappa=%.2f (ion below nominal peak)" % (
             acs["om0_obs"] / 1e6, acs["kappa"]),
         "                              nominal(upper)   Rabi-anchored",
         "    B1 (pi,  CONTINUOUS):       %+7.1f kHz     %+7.1f kHz   [constant floor]" % (
             acs["dB1_nom"] / 1e3, acs["dB1"] / 1e3),
         "    R2 (lin, PULSED):           %+7.1f kHz     %+7.1f kHz   [t_R2=%.2f us fixed]" % (
             acs["dR2_nom"] / 1e3, acs["dR2"] / 1e3, acs["t_R2_us"]),
         "    sum (both on, during pulse):%+7.1f kHz     %+7.1f kHz" % (
             (acs["dB1_nom"] + acs["dR2_nom"]) / 1e3, (acs["dB1"] + acs["dR2"]) / 1e3),
         "    scalar (omega_HF/Delta)*Omega_0 cross-check: %+.1f kHz  (~ the -40 kHz offset)" % (
             acs["d_scalar"] / 1e3),
         "    R2 pi-pulse phase (N-independent): %+.2f cycles" % acs["phi_R2_cyc"],
         "    N    delta_t/us  R2-duty  t_seq/us  <delta_AC>/kHz  phi_B1/cyc  phi_tot/cyc"]
    for r in acs["rows"]:
        L.append("  %4d   %8.4f   %5.2f   %7.1f   %+11.1f   %+8.2f   %+8.2f" % (
            r["N"], r["delta_t_us"], r["duty"], r["t_seq_us"], r["d_eff_khz"],
            r["phi_B1_cyc"], r["phi_tot_cyc"]))
    L.append("  => the SHAPE is robust: B1's continuous shift dominates and its phase")
    L.append("     grows ~N (few cycles already at N=50); R2's pi-pulse phase is fixed.")
    L.append("     Small N keeps the AC-Stark phase small; the grating wants many -> trade-off.")
    return "\n".join(L)


def make_acstark_figure(acs, out=_OUT_ACS):
    Ns = [r["N"] for r in acs["rows"]]
    deff = [r["d_eff_khz"] for r in acs["rows"]]
    phiB1 = [abs(r["phi_B1_cyc"]) for r in acs["rows"]]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))
    ax1.semilogx(Ns, deff, "o-", color=_RED, label="Rabi-anchored ⟨δ$_{AC}$⟩")
    ax1.axhline(acs["dB1"] / 1e3, ls="--", color=_BLUE,
                label="B1 floor (continuous): %+.0f kHz" % (acs["dB1"] / 1e3))
    ax1.axhline(acs["d_scalar"] / 1e3, ls=":", color=_GRAY,
                label="scalar (ω$_{HF}$/Δ)·Ω$_0$: %+.0f kHz" % (acs["d_scalar"] / 1e3))
    ax1.set_xlabel("N  (stroboscopic cycles)")
    ax1.set_ylabel("time-averaged differential AC-Stark shift  [kHz]")
    ax1.set_title("Effective shift over the train  (δ$_{B1}$ + δ$_{R2}$·duty)")
    ax1.grid(alpha=0.3, which="both"); ax1.legend(fontsize=8)
    ax2.loglog(Ns, phiB1, "o-", color=_RED, label="B1 (continuous) ∝ N")
    ax2.axhline(abs(acs["phi_R2_cyc"]), ls="--", color=_BLUE,
                label="R2 (pulsed π) fixed: %.2f cyc" % abs(acs["phi_R2_cyc"]))
    ax2.set_xlabel("N  (stroboscopic cycles)")
    ax2.set_ylabel("|accumulated AC-Stark phase|  [cycles]")
    ax2.set_title("AC-Stark phase of the π pulse")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=8)
    fig.suptitle("OC π-pulse AC-Stark systematics vs N  —  B1 continuous, only R2 pulsed "
                 "(Rabi-anchored, κ=%.2f)" % acs["kappa"], fontsize=10.5)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    if not _DATAFILE.exists():
        print("strobo data not found at", _DATAFILE)
        return 0
    info = analyze()
    print(report(info))
    out = make_figure(info)
    print("\nwrote", out.relative_to(Path(__file__).resolve().parent.parent))
    # forward simulation: detuning scan at delta_t = 0.02 us (covers the first sidebands)
    sim = simulate_detuning_scan(omega_strobo_hz=info["omega_strobo"], delta_t_us=0.02)
    print("\nSIMULATED DETUNING SCAN (delta_t=0.02 us, N=%g): comb teeth at" % sim["N"],
          ", ".join("%+.0f kHz" % (d / 1e3) for d, _ in sim["peaks"]),
          "(width ~%.0f kHz)" % sim["fwhm_khz"])
    out2 = make_detuning_figure(sim)
    print("wrote", out2.relative_to(Path(__file__).resolve().parent.parent))
    # AC-Stark systematics vs N (B1 continuous, only R2 pulsed for a pi pulse)
    acs = ac_stark_vs_N(omega_strobo_hz=info["omega_strobo"])
    print(report_acstark(acs))
    out3 = make_acstark_figure(acs)
    print("\nwrote", out3.relative_to(Path(__file__).resolve().parent.parent))
    return 0


if __name__ == "__main__":   # pragma: no cover
    raise SystemExit(main())
