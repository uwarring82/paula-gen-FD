"""
Digital twin of the OC (Orthogonal-Carrier) axial Raman carrier flop, PAULA
sources/data/OC_Axial/0_Car_Flop/09_46_23_12_06_2026.dat.

The OC combination is a two-photon stimulated-Raman (TPSR) drive: B1 (pi) + R2
(sigma), both ~20 GHz red of 3P_3/2, with Delta_k || z (axial). `0_Car_Flop` is a
CARRIER flop (Delta n = 0) vs the pulse duration t_oc (0..9.75 us, 21 points, 40
shots), so the addressed |down>=|3,3> <-> |up>=|2,2> population Rabi-oscillates and
its contrast decays.

This twin composes FOUR mechanisms, each labelled by what anchors it:

  COHERENT       carrier Rabi Omega = Omega_B Omega_R / (2 Delta_R); here taken
                 from the measured flop (engines.rabi.fit_rabi -> Omega, t_pi).

  AC-STARK       the differential light shift delta_AC = (omega_HF/Delta_R) Omega
                 between |down> and |up> (engines.scatter). LEDGER-DERIVED from
                 Delta_R + omega_HF. Speeds Omega_eff and caps the amplitude.

  SCATTERING     off-resonant spontaneous photon scattering Gamma_sc = (Gamma/
                 Delta_R) Omega (engines.scatter), envelope e^{-(3/4)Gamma_sc t}.
                 LEDGER-DERIVED from Delta_R + Gamma.

  MOTIONAL       the carrier Debye-Waller dephasing over the THERMAL phonon
                 distribution: Omega_{n,n} = Omega_0 e^{-eta^2/2} L_n(eta^2), so a
                 thermal state (mean nbar) is a spread of Rabi frequencies and the
                 flop dephases (engines.sideband.thermal_carrier_flip). LEDGER-
                 DERIVED from eta (sideband, OC->lf at omega_lf) + nbar (the RSB-
                 cooled benchmark mg_rsb_cooled_nbar_axial_lf_25mg).

This replaces the earlier *empirical* exponential residual (ADR-0007) with the
ledger-anchored MOTIONAL channel, and adds the diagnostic that makes the result
honest: the **effective nbar inversion**. At the RSB-cooled nbar=0.07 the motional
dephasing is still small (~8% of the observed decay; ~4x the 2% scattering floor),
so the ledger floor under-predicts the data. Inverting the observed decay through
the same thermal model gives the nbar that WOULD reproduce it: here nbar_eff ~ 1,
i.e. ~14x the cooled benchmark. So this OC flop is consistent with a near-unity-nbar
motional state -- sideband-cooling underperformance / heating in this run and/or
technical dephasing (Raman intensity-phase or B-field noise) -- NOT with the
RSB-cooled 0.07. The twin SEPARATES and QUANTIFIES; it does not fabricate a fit.

    python -m spike.twin_oc_flop   ->  docs/figures/twin_oc_axial_carrier_flop.png

Pure Python + matplotlib (Agg); randomness via a seeded stdlib RNG.
"""
from __future__ import annotations

import math
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .datfile import DatFile  # noqa: E402
from .engines.rabi import fit_rabi  # noqa: E402
from .engines.raman_optical import RamanOptics, beams_from_ledger  # noqa: E402
from .engines.scatter import CONTRAST_DECAY_FACTOR, RamanScatter  # noqa: E402
from .engines.sideband import Sideband, thermal_coherence  # noqa: E402
from .ledger import Ledger  # noqa: E402

_DATAFILE = (Path(__file__).resolve().parent.parent / "sources" / "data" /
             "OC_Axial" / "0_Car_Flop" / "09_46_23_12_06_2026.dat")
_OUT = Path(__file__).resolve().parent.parent / "docs" / "figures" / "twin_oc_axial_carrier_flop.png"
_BLUE, _RED, _GRAY, _GREEN = "#1f77b4", "#d62728", "#888888", "#2ca02c"
_F_LO_KHZ, _F_HI_KHZ = 80.0, 260.0      # flop-fit window (the peak sits near t_pi ~ 2.5 us)
_TSPAN_US = 9.75                        # scan span; the contrast-decay reference window
_OMEGA_LF_NAME = "omega_z_axial_com_25mg"           # input: freddy axial lf mode (1.30 MHz)
_NBAR_NAME = "mg_rsb_cooled_nbar_axial_lf_25mg"     # benchmark: RSB-cooled nbar (0.07)
_NBAR_GRID = [0.0, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]


# --- twin model -------------------------------------------------------------
class OCFlopTwin:
    """Forward model of the OC carrier flop: coherent Rabi x AC-Stark x scattering
    x thermal motional dephasing, plus per-shot detection counts. The motional
    channel uses the carrier Debye-Waller thermal sum; the scattering channel an
    exponential envelope toward the mixed state."""

    def __init__(self, rabi_hz: float, scatter: RamanScatter, sideband: Sideband,
                 omega_lf_hz: float, nbar: float, mu_bright: float, mu_dark: float,
                 balance: float = 1.0, n_shots: int = 40,
                 gamma_sc_hz: float = None, delta_ac_hz: float = None):
        self.rabi = rabi_hz
        self.sc = scatter
        self.sb = sideband
        self.omega_lf = omega_lf_hz
        self.nbar = nbar                    # the RSB-cooled benchmark nbar
        self.mu_bright = mu_bright          # mean counts when FLIPPED (|up>, bright here)
        self.mu_dark = mu_dark              # mean counts when UN-flipped (|down>)
        self.balance = balance
        self.n_shots = n_shots
        self.eta = sideband.lamb_dicke("OC", "lf", omega_lf_hz)
        # scattering + differential AC-Stark: the polarization+power-resolved engine
        # (raman_optical) when supplied, else the scalar/leading-order scatter engine.
        self.gamma_sc = scatter.scatter_rate(rabi_hz, balance) if gamma_sc_hz is None else gamma_sc_hz
        self.delta_ac = scatter.stark_detuning(rabi_hz) if delta_ac_hz is None else delta_ac_hz
        self.gamma_sc_contrast = CONTRAST_DECAY_FACTOR * self.gamma_sc
        self.omega_eff = math.hypot(rabi_hz, self.delta_ac)
        self.amp_cap = (rabi_hz * rabi_hz) / (self.omega_eff * self.omega_eff)

    def flip_prob(self, t_us: float, nbar: float | None = None,
                  with_scatter: bool = True, with_motional: bool = True) -> float:
        """P(flipped) at duration t_us. MOTIONAL: the carrier Debye-Waller thermal
        sum (engines.sideband) with the AC-Stark detuning folded in per component.
        SCATTERING: an exponential envelope damping the oscillation toward 1/2."""
        t = t_us * 1e-6
        nb = self.nbar if nbar is None else nbar
        if with_motional:
            p = self.sb.carrier_thermal_flip("OC", "lf", self.omega_lf, self.rabi,
                                             nb, t, detuning_hz=self.delta_ac)
        else:
            eff, amp = self.omega_eff, self.amp_cap
            p = 0.5 * amp * (1.0 - math.cos(2.0 * math.pi * eff * t))
        if with_scatter and self.gamma_sc_contrast > 0.0:
            p = 0.5 + (p - 0.5) * math.exp(-self.gamma_sc_contrast * t)
        return p

    def curve(self, ts_us, **kw):
        return [self.flip_prob(t, **kw) for t in ts_us]

    def signal_curve(self, ts_us, **kw):
        """Expected COUNT signal = mu_dark + P_flip (mu_bright - mu_dark)."""
        return [self.mu_dark + p * (self.mu_bright - self.mu_dark)
                for p in self.curve(ts_us, **kw)]

    def effective_decay(self, nbar: float, t_ref_us: float = _TSPAN_US) -> float:
        """Effective contrast-decay rate [1/s] of the full predicted flop at a given
        nbar, over the window [0, t_ref]: the scattering rate (3/4)Gamma_sc PLUS the
        thermal dephasing -ln|C(t_ref)|/t_ref (engines.sideband.thermal_coherence).
        The two add because the model contrast at t_ref is |C|*e^{-(3/4)Gamma_sc t},
        and the data's fit_rabi exponential gives e^{-gamma_obs t} -- so matching this
        to gamma_obs matches the contrast at the end of the scan (same estimator,
        analytically, no per-nbar curve fit)."""
        t = t_ref_us * 1e-6
        c = thermal_coherence(t, self.rabi, self.eta, nbar)
        gamma_motional = -math.log(max(c, 1e-12)) / t
        return self.gamma_sc_contrast + gamma_motional

    def invert_nbar(self, gamma_target: float):
        """The nbar whose predicted flop decays at gamma_target (the observed rate),
        by monotone interpolation on _NBAR_GRID. Returns (nbar_eff, saturated): if
        the decay never reaches gamma_target on the grid (fast-dephasing limit, the
        flop barely completes a period), returns the grid max flagged saturated."""
        grid = [(nb, self.effective_decay(nb)) for nb in _NBAR_GRID]
        if gamma_target <= grid[0][1]:
            return 0.0, False
        for (n0, g0), (n1, g1) in zip(grid, grid[1:]):
            if g0 <= gamma_target <= g1 and g1 > g0:
                return n0 + (n1 - n0) * (gamma_target - g0) / (g1 - g0), False
        return grid[-1][0], True            # never reached target -> saturated

    def simulate_counts(self, t_us: float, rng: random.Random, nbar: float | None = None):
        """Per-shot detection cloud at one duration: QPN projection at P_flip, then a
        Poisson count from the bright (flipped) or dark (un-flipped) level."""
        p = self.flip_prob(t_us, nbar=nbar)
        out = []
        for _ in range(self.n_shots):
            mu = self.mu_bright if rng.random() < p else self.mu_dark
            out.append(_poisson(mu, rng))
        return out


def _poisson(mu: float, rng: random.Random) -> int:
    """Knuth's Poisson sampler (same as twin._poisson)."""
    if mu <= 0.0:
        return 0
    el, k, p = math.exp(-mu), 0, 1.0
    while True:
        k += 1
        p *= rng.random()
        if p <= el:
            return k - 1


# --- build from data + ledger ----------------------------------------------
def build(path=_DATAFILE, n_boot=0, seed=0):
    """Fit the measured flop, read the four loss channels from the ledger, and
    decompose the observed decay: the ledger floor (scattering + thermal at the
    RSB-cooled nbar) vs the effective-nbar inversion. Returns (twin, fit, info).

    With n_boot > 0, a non-parametric shot bootstrap (resample the per-shot counts at
    each scan point, re-fit, re-invert nbar_eff) adds robust 1-sigma uncertainties
    `rabi_hz_err`, `g_obs_err`, `nbar_eff_err` to info -- a real error bar on nbar_eff,
    not just a chi-square."""
    dat = DatFile(path)
    t, y, s = dat.signal()
    fit = fit_rabi(t, y, s, f_lo_khz=_F_LO_KHZ, f_hi_khz=_F_HI_KHZ)
    rabi = fit["freq_hz"]
    g_obs = fit["decay_per_s"]

    ledger = Ledger.load()
    scatter = RamanScatter.from_ledger(ledger)
    sideband = Sideband.from_ledger(ledger)
    omega_lf = ledger.input_quantity(_OMEGA_LF_NAME).value
    nbar_rsb = ledger.benchmark_quantity(_NBAR_NAME).value

    hists = [DatFile.hist_mean(h) for h in dat.histograms()]
    mu_bright, mu_dark = max(hists), min(hists)

    # POLARIZATION + POWER-RESOLVED differential AC-Stark + scattering (raman_optical):
    # OC = B1 (pi) + R2 (linear), per-beam powers from the .dat (pwr_b1/pwr_r2). The
    # dimensionless ratios are anchored to the measured Rabi (delta_AC[Hz] = ratio*rabi;
    # Gamma_sc[1/s] = ratio*2pi*rabi -- matching the scalar engine's conventions).
    optics = RamanOptics.from_ledger(ledger)
    powers = {"b1": dat.settings.get("pwr_b1", 1.0), "r2": dat.settings.get("pwr_r2", 1.0)}
    bB, bR = beams_from_ledger(ledger, "OC", powers=powers)
    beams = [bB, bR]
    stark_ratio = optics.differential_stark_per_rabi(beams, bB, bR)        # delta_AC = ratio*rabi
    scatter_ratio = optics.scatter_per_rabi(beams, bB, bR, state=(3.0, 3.0))  # Gamma_sc = ratio*2pi*rabi
    delta_ac = stark_ratio * rabi
    gamma_sc = scatter_ratio * 2.0 * math.pi * rabi
    # leading-order (scalar) values, for the old-vs-new comparison in the report
    delta_ac_scalar = scatter.stark_detuning(rabi)
    gamma_sc_scalar = scatter.scatter_rate(rabi)

    def _make_twin(rabi_hz):
        """A twin at a given Rabi, with the (rabi-scaled) resolved loss channels."""
        return OCFlopTwin(rabi_hz=rabi_hz, scatter=scatter, sideband=sideband,
                          omega_lf_hz=omega_lf, nbar=nbar_rsb, mu_bright=mu_bright,
                          mu_dark=mu_dark, n_shots=dat.scan["shots"],
                          gamma_sc_hz=scatter_ratio * 2.0 * math.pi * rabi_hz,
                          delta_ac_hz=stark_ratio * rabi_hz)

    twin = _make_twin(rabi)
    g_floor = twin.effective_decay(nbar_rsb)                      # scatter + thermal@nbar_rsb
    g_thermal = max(0.0, g_floor - twin.gamma_sc_contrast)        # thermal-only contribution
    nbar_eff, saturated = twin.invert_nbar(g_obs)

    info = {
        "t": t, "y": y, "s": s, "dat": dat,
        "rabi_hz": rabi, "t_pi_us": fit["t_pi_us"], "chi2r": fit["chi2_reduced"],
        "g_obs": g_obs, "g_sc": twin.gamma_sc_contrast, "g_thermal": g_thermal,
        "g_floor": g_floor, "eta": twin.eta, "omega_lf": omega_lf,
        "nbar_rsb": nbar_rsb, "nbar_eff": nbar_eff, "nbar_saturated": saturated,
        "delta_ac": twin.delta_ac, "amp_cap": twin.amp_cap, "se_per_pi": scatter.se_per_pi(),
        "gamma_sc_rate": twin.gamma_sc, "mu_bright": mu_bright, "mu_dark": mu_dark,
        "tspan_us": t[-1] if t else 0.0,
        # polarization+power-resolved vs leading-order (scalar) comparison
        "delta_ac_scalar": delta_ac_scalar, "gamma_sc_scalar": gamma_sc_scalar,
        "pwr_b1": powers["b1"], "pwr_r2": powers["r2"],
        "circ_b1": optics.circular_fraction(bB), "circ_r2": optics.circular_fraction(bR),
    }

    if n_boot > 0:
        from .bootstrap import shot_bootstrap, summarize

        def _replica(xs, yb, sb):
            f = fit_rabi(xs, yb, sb, f_lo_khz=_F_LO_KHZ, f_hi_khz=_F_HI_KHZ)
            tw = _make_twin(f["freq_hz"])
            nb, _sat = tw.invert_nbar(f["decay_per_s"])
            return {"rabi_hz": f["freq_hz"], "g_obs": f["decay_per_s"], "nbar_eff": nb}

        runs = shot_bootstrap(_replica, dat.point_shots(), n_boot=n_boot, seed=seed)
        summ = summarize(runs, ("rabi_hz", "g_obs", "nbar_eff"))
        info["rabi_hz_err"] = summ["rabi_hz"]["sigma"]
        info["g_obs_err"] = summ["g_obs"]["sigma"]
        info["nbar_eff_err"] = summ["nbar_eff"]["sigma"]
        info["n_boot"] = summ["nbar_eff"]["n"]

    return twin, fit, info


# --- report -----------------------------------------------------------------
def _loss(rate, tspan_us):
    return 100.0 * (1.0 - math.exp(-rate * tspan_us * 1e-6))


def report(info) -> str:
    tspan = info["tspan_us"]
    frac = lambda g: 100.0 * g / info["g_obs"] if info["g_obs"] else 0.0

    def pm(key, scale=1.0, fmt="%.1f"):
        """' ± <err>' for a bootstrapped quantity, '' if no bootstrap was run."""
        e = info.get(key)
        return "" if e is None else " ± " + (fmt % (e / scale))

    boot = (" (%d-replica shot bootstrap)" % info["n_boot"]) if info.get("n_boot") else ""
    nbar_eff = ((">= %.1f" % info["nbar_eff"]) if info["nbar_saturated"]
                else ("%.2f%s" % (info["nbar_eff"], pm("nbar_eff_err", 1.0, "%.2f"))))
    L = [
        "OC AXIAL CARRIER FLOP — digital twin (TPSR, Delta_k || z; %s)" % (info["dat"].timestamp or ""),
        "  data: %s  (t_oc 0..%.2f us, %d pts x %d shots)%s" % (
            _DATAFILE.name, tspan, info["dat"].scan["points"], info["dat"].scan["shots"], boot),
        "",
        "COHERENT (from the measured flop):",
        "  Omega/2pi = %.1f%s kHz   t_pi = %.2f us   (fit chi2_red = %.1f)" % (
            info["rabi_hz"] / 1e3, pm("rabi_hz_err", 1e3, "%.1f"), info["t_pi_us"], info["chi2r"]),
        "",
        "AC-STARK — polarization+power resolved (raman_optical, full |J',mJ'> sum):",
        "  OC = B1(pi, pwr %.3f, C=%+.2f) + R2(linear, pwr %.3f, C=%+.2f)" % (
            info["pwr_b1"], info["circ_b1"], info["pwr_r2"], info["circ_r2"]),
        "  delta_AC = %+.1f kHz -> amplitude cap %.4f (%.2f%% loss)   "
        "[scalar omega_HF/Delta_R: %+.1f kHz]" % (
            info["delta_ac"] / 1e3, info["amp_cap"], 100.0 * (1.0 - info["amp_cap"]),
            info["delta_ac_scalar"] / 1e3),
        "",
        "LEDGER-ANCHORED DECAY CHANNELS (effective rate; %% of the observed decay):",
        "  scattering  (full mJ sum):           %.2e /s  (%.0f%%)  [scalar Gamma/Delta_R: %.2e /s]" % (
            info["g_sc"], frac(info["g_sc"]), CONTRAST_DECAY_FACTOR * info["gamma_sc_scalar"]),
        "  motional    (eta=%.3f, nbar=%.2f):  %.2e /s  (%.0f%%)  carrier Debye-Waller" % (
            info["eta"], info["nbar_rsb"], info["g_thermal"], frac(info["g_thermal"])),
        "  -> ledger floor (scatter+motional):  %.2e /s  (%.0f%%)" % (
            info["g_floor"], frac(info["g_floor"])),
        "",
        "OBSERVED vs FLOOR:",
        "  observed decay = %.2e%s /s  (tau = %.1f us, %.0f%% loss over the scan)" % (
            info["g_obs"], pm("g_obs_err", 1.0, "%.1e"), 1e6 / info["g_obs"], _loss(info["g_obs"], tspan)),
        "  the ledger floor explains only ~%.0f%% of it at the RSB-cooled nbar=%.2f." % (
            frac(info["g_floor"]), info["nbar_rsb"]),
        "",
        "EFFECTIVE-nbar INVERSION (attribute the decay to motion through the same model):",
        "  nbar_eff = %s   vs RSB-cooled nbar = %.2f  (~%.0fx hotter)" % (
            nbar_eff, info["nbar_rsb"],
            (info["nbar_eff"] / info["nbar_rsb"]) if info["nbar_rsb"] else float("inf")),
        "  => consistent with a near-unity-nbar motional state: sideband-cooling",
        "     underperformance / heating in this run and/or technical dephasing",
        "     (Raman intensity-phase or B-field noise) — NOT the cooled 0.07. Flagged.",
        "",
        "DETECTION: bright(flipped) mu = %.2f counts, dark(un-flipped) mu = %.3f counts." % (
            info["mu_bright"], info["mu_dark"]),
    ]
    return "\n".join(L)


# --- figure -----------------------------------------------------------------
def make_figure(twin, info, out=_OUT, seed: int = 7):
    t, y, s = info["t"], info["y"], info["s"]
    rng = random.Random(seed)
    ts = [0.02 * k * info["tspan_us"] for k in range(51)]
    nbar_eff = info["nbar_eff"]

    floor = twin.signal_curve(ts, nbar=info["nbar_rsb"])         # ledger floor: scatter + thermal@0.07
    best = twin.signal_curve(ts, nbar=nbar_eff)                  # decay attributed to motion: nbar_eff

    cloud_x, cloud_y = [], []
    for ti in t:
        for c in twin.simulate_counts(ti, rng, nbar=nbar_eff):
            cloud_x.append(ti)
            cloud_y.append(c)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    ax.scatter(cloud_x, cloud_y, s=8, color=_GRAY, alpha=0.18, edgecolors="none",
               zorder=1, label="twin per-shot counts (QPN+Poisson, n̄$_{eff}$)")
    ax.errorbar(t, y, yerr=s, fmt="o", color=_BLUE, ms=5, capsize=2, zorder=3,
                label="data (mean ± s.e.)")
    ax.plot(ts, best, "-", color=_RED, lw=2.0, zorder=4,
            label="twin @ n̄$_{eff}$ = %.2f (decay attributed to motion)" % nbar_eff)
    ax.plot(ts, floor, "--", color=_GREEN, lw=1.6, zorder=2,
            label="ledger floor: scatter + thermal @ n̄ = %.2f" % info["nbar_rsb"])

    ax.set_xlabel("OC pulse duration  t_oc  [µs]")
    ax.set_ylabel("fluorescence counts  (∝ P$_{flip}$)")
    ax.set_title("OC axial carrier flop — twin vs PAULA data (%s)" % (info["dat"].timestamp or ""))
    txt = "\n".join([
        r"$\Omega/2\pi$ = %.0f kHz   $\eta$ = %.3f" % (info["rabi_hz"] / 1e3, info["eta"]),
        r"$\delta_{AC}$ = %.1f kHz   $P_{SE}/\pi$ = %.2f%%" % (
            info["delta_ac"] / 1e3, 100 * info["se_per_pi"]),
        r"scatter %.0f%% + motional@%.2f %.0f%% of decay" % (
            100 * info["g_sc"] / info["g_obs"], info["nbar_rsb"],
            100 * info["g_thermal"] / info["g_obs"]),
        r"observed $\tau$ = %.1f µs  $\Rightarrow$  n̄$_{eff}$ = %.2f%s" % (
            1e6 / info["g_obs"], nbar_eff,
            (r" $\pm$ %.2f" % info["nbar_eff_err"]) if info.get("nbar_eff_err") else ""),
        r"(RSB-cooled n̄ = %.2f $\rightarrow$ ~%.0f× hotter)" % (
            info["nbar_rsb"], nbar_eff / info["nbar_rsb"] if info["nbar_rsb"] else 0),
    ])
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec=_GRAY, alpha=0.85))
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main(argv=None) -> int:
    if not _DATAFILE.exists():
        print("OC carrier-flop data not found at", _DATAFILE)
        return 0
    twin, fit, info = build(n_boot=300)        # error bars on Omega, tau_obs, nbar_eff
    print(report(info))
    out = make_figure(twin, info)
    print("\nwrote", out.relative_to(Path(__file__).resolve().parent.parent))
    return 0


if __name__ == "__main__":   # pragma: no cover
    raise SystemExit(main())
