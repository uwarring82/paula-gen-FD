"""
Generate the figures embedded in
docs/notes/aom_finite_sound_velocity_rabi.md:
  aom_single_vs_double_pass.png   -- rotation-area & peak efficiency and delivered
                                     pulse shapes for n=1 (R2 single-pass, current),
                                     n=2, n=3 (R2 single + B1 double, alternative)
  aom_flop_single_vs_double_pass.png -- predicted carrier flop P_down vs delta_t for the
                                     current (n=1) and alternative (n=3) schemes at EQUAL
                                     drive power, with the measured 0-100 ns flop overlaid.

Device constants from sources/spec sheets/intraact_aom220uv_specs.pdf (IntraAction
ASM-2202B3): acoustic velocity V = 5.95 mm/us, fused silica.  Run:
    python docs/figures/make_aom_rise_figs.py
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_REPO, "spike"))
from engines import aom_rise as aom  # noqa: E402  (canonical pure-Python model; single source of truth)

# All plotting helpers below work in NANOSECONDS (tf via aom.tau_f[us] -> ns); the canonical
# model lives in spike/engines/aom_rise.py (pure Python, no numpy.trapz -> portable across
# numpy 1.x/2.x). These thin wrappers just delegate, so the note tables and tests agree.
_C = {1: "#1f77b4", 2: "#ff7f0e", 3: "#d62728"}
_LBL = {1: "n=1  current: R2 single-pass switched (B1 on)",
        2: "n=2  e.g. R2 single + B1 single, both switched",
        3: "n=3  alternative: R2 single + B1 double, both switched"}


def tau_f(D_mm):
    """field 1/e rise time constant [ns]; w = D/2 (D = 1/e^2 intensity diameter)."""
    return aom.tau_f(D_mm) * 1e3


_af = np.vectorize(aom.field_envelope)


def a_field(t, dt, D_mm):
    """single-pass diffracted FIELD envelope for a nominal gate of width dt [ns]."""
    return _af(t, dt, tau_f(D_mm))


def area_eff(dt, D_mm, n):
    return aom.area_ratio(dt, tau_f(D_mm), n)


def peak_eff(dt, D_mm, n):
    return aom.peak(dt, tau_f(D_mm), n)


def w_equiv(dt, D_mm, n):
    return aom.w_equiv(dt, tau_f(D_mm), n)


def fwhm(dt, D_mm, n):
    return aom.fwhm(dt, tau_f(D_mm), n)


# --------------------------------------------------------------------------- #
def fig_efficiency(path, D_mm=1.0):
    tf = tau_f(D_mm)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.5))

    # (a) efficiency vs pulse width ------------------------------------------
    dt = np.logspace(np.log10(5), np.log10(4000), 260)
    axL.axvspan(5, 100, color="0.85", alpha=.7, lw=0)
    axL.text(22, 0.06, "scan\n0-100 ns", fontsize=8, ha="center", color="0.3")
    for n in (1, 2, 3):
        ra = np.array([area_eff(x, D_mm, n) for x in dt])
        axL.plot(dt, ra, _C[n], lw=2, label=_LBL[n])
        rp = np.array([peak_eff(x, D_mm, n) for x in dt])
        axL.plot(dt, rp, _C[n], lw=1, ls=":", alpha=.7)
    axL.plot([], [], "k-", lw=2, label="rotation area  $R_{\\rm area}$ (solid)")
    axL.plot([], [], "k:", lw=1, label="peak amplitude $R_{\\rm peak}$ (dotted)")
    axL.axvline(60, color="0.4", lw=1, ls="--")
    axL.text(64, 0.9, "60 ns\n(set)", fontsize=8, color="0.3")
    axL.set_xscale("log"); axL.set_xlim(5, 4000); axL.set_ylim(0, 1.1)
    axL.set_xlabel("nominal pulse width  $\\delta t$  [ns]")
    axL.set_ylabel("efficiency  /  ideal rectangular pulse")
    axL.set_title(f"(a) effective Rabi vs pulse width  (D={D_mm} mm, $\\tau_f$={tf:.0f} ns)")
    axL.legend(fontsize=7.4, loc="center right"); axL.grid(alpha=.3)

    # (b) delivered pulse shapes at dt = 60 ns -------------------------------
    dt0 = 60.0
    t = np.linspace(-2.5 * tf, dt0 + 2.5 * tf, 1200)
    a = a_field(t, dt0, D_mm)
    for n in (1, 2, 3):
        axR.plot(t, a ** n, _C[n], lw=2,
                 label=f"n={n}: peak {peak_eff(dt0, D_mm, n):.2f}, "
                       f"area {area_eff(dt0, D_mm, n):.2f}")
    axR.axvspan(0, dt0, color="0.9", lw=0)
    axR.text(dt0 / 2, -0.07, "nominal gate", fontsize=8, ha="center", color="0.4")
    axR.set_xlabel("time [ns]"); axR.set_ylabel("Raman coupling $\\Omega(t)/\\Omega_0 \\propto a^n$")
    axR.set_title("(b) delivered pulse at $\\delta t$=60 ns: amplitude lost, edges sharpened")
    axR.legend(fontsize=8); axR.grid(alpha=.3); axR.set_xlim(-160, 220); axR.set_ylim(-0.12, 1.0)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", path)


# --------------------------------------------------------------------------- #
def _measured_flop():
    """measured 0-100 ns flop (file 10_22_22); returns (delta_t[ns], P_norm in [0,1])."""
    from datfile import DatFile
    f = os.path.join(_REPO, "sources/data/Strobo2.0/"
                     "1_FlopN_3p3_2p2_PDQ_displ_strobo/10_22_22_12_06_2026.dat")
    df = DatFile(f)
    x, y, _ = df.signal()
    x = np.array(x) * 1e3                       # us -> ns
    y = np.array(y)
    return x, (y - y.min()) / (y.max() - y.min())


def fig_flop(path, D_mm=1.0, N=50, T_period=38.0):
    """carrier flop  P_down = sin^2(Theta/2),  Theta(dt) = N*Omega0*dt*R_area_n(dt),
    Omega0 fixed by the measured n=1 fringe spacing (~T_period ns).  Same Omega0 for n=3."""
    Omega0 = 2 * np.pi / (N * T_period)         # rad/ns, from the data fringe spacing
    dt = np.linspace(0, 150, 600)
    fig, ax = plt.subplots(figsize=(8.6, 4.6))

    try:
        mx, my = _measured_flop()
        # coarse phase alignment of the model to the displaced-state start
        ax.plot(mx, my, "o", ms=5, color="0.45", label="measured (norm. fluorescence)")
        phi0 = -0.17  # set so n=1 fringe maxima land on the measured peaks (~20/58/98 ns)
    except Exception as e:                       # data optional
        print("data overlay skipped:", e); phi0 = -0.17

    for n, ls in ((1, "-"), (3, "-")):
        ra = np.array([area_eff(x, D_mm, n) if x > 0 else 0.0 for x in dt])
        Theta = N * Omega0 * dt * ra
        P = np.sin((Theta + phi0) / 2) ** 2
        ax.plot(dt, P, ls, color=_C[n], lw=2.2, label=_LBL[n].split("  ", 1)[1])

    ax.axvline(60, color="0.4", lw=1, ls="--"); ax.text(61, 0.02, "60 ns", fontsize=8, color="0.3")
    ax.set_xlabel("nominal pulse width  $\\delta t$  [ns]")
    ax.set_ylabel("$P_\\downarrow$  (flip probability)")
    ax.set_title("Predicted carrier flop at EQUAL drive power — current (n=1) vs alternative (n=3)")
    ax.set_xlim(0, 150); ax.set_ylim(-0.03, 1.03); ax.grid(alpha=.3); ax.legend(fontsize=8.5)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", path)


def fig_effective_length(path, D_mm=1.0, f_lf=1.3001e-3):
    """The ion sees a LONGER pulse than the electronic gate, and its width SATURATES at a
    floor set by the acoustic transit (it cannot be made arbitrarily short)."""
    tf = tau_f(D_mm)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.6))

    # (a) effective width vs nominal gate ------------------------------------
    dt = np.linspace(1, 420, 240)
    we1 = np.array([w_equiv(x, D_mm, 1) for x in dt])
    fw1 = np.array([fwhm(x, D_mm, 1) for x in dt])
    we3 = np.array([w_equiv(x, D_mm, 3) for x in dt])
    floor1 = np.sqrt(np.pi) * tf                 # equiv-width floor, single pass (n=1)
    floor3 = np.sqrt(np.pi / 3) * tf             # equiv-width floor, R2single+B1double (n=3)
    axL.plot(dt, dt, "0.6", lw=1.2, ls="--", label="ideal: ion pulse = gate $\\delta t$")
    axL.plot(dt, we1, _C[1], lw=2.4, label="$w_{\\rm eff}$ current (n=1, single-pass)")
    axL.plot(dt, fw1, _C[1], lw=1.2, ls=":", label="FWHM current (n=1)")
    axL.plot(dt, we3, _C[3], lw=2.4, label="$w_{\\rm eff}$ alternative (n=3, +B1 double)")
    axL.axhline(floor1, color=_C[1], lw=1, alpha=.5)
    axL.axhline(floor3, color=_C[3], lw=1, alpha=.5)
    axL.text(300, floor1 + 4, f"floor $\\sqrt{{\\pi}}\\,\\tau_f$ = {floor1:.0f} ns", fontsize=8, color=_C[1])
    axL.text(300, floor3 + 4, f"floor $\\sqrt{{\\pi/3}}\\,\\tau_f$ = {floor3:.0f} ns", fontsize=8, color=_C[3])
    axL.axvspan(0, 100, color="0.85", alpha=.7, lw=0)
    axL.text(50, 360, "scan 0-100 ns:\nion pulse pinned at the floor", fontsize=8, ha="center", color="0.3")
    axL.set_xlim(0, 420); axL.set_ylim(0, 420)
    axL.set_xlabel("nominal electronic gate  $\\delta t$  [ns]")
    axL.set_ylabel("effective pulse length seen by the ion  [ns]")
    axL.set_title(f"(a) effective pulse length saturates  (D={D_mm} mm, $\\tau_f$={tf:.0f} ns)")
    axL.legend(fontsize=8, loc="lower right"); axL.grid(alpha=.3)
    ax2 = axL.twinx()                            # motional phase smeared per pulse
    ax2.set_ylim(0, 420 * 2 * np.pi * f_lf)
    ax2.set_ylabel("motional phase per pulse  $\\omega_{\\rm lf} w_{\\rm eff}$  [rad]", color="0.4")
    ax2.tick_params(colors="0.4")

    # (b) normalised pulses collapse onto the floor at small dt --------------
    for dt0, c in [(10, "#7b3294"), (30, "#1f77b4"), (60, "#2ca02c"),
                   (120, "#ff7f0e"), (240, "#d62728")]:
        t = np.linspace(-3 * tf, dt0 + 3 * tf, 1400)
        p = a_field(t, dt0, D_mm)
        axR.plot(t - dt0 / 2, p / p.max(), c, lw=1.9,
                 label=f"$\\delta t$={dt0} ns  ($w_{{\\rm eff}}$={w_equiv(dt0, D_mm, 1):.0f} ns)")
    axR.set_xlabel("time (centred) [ns]"); axR.set_ylabel("normalised $\\Omega(t)/\\Omega_{\\rm peak}$  (n=1)")
    axR.set_title("(b) short gates deliver the SAME minimum-width pulse")
    axR.legend(fontsize=8); axR.grid(alpha=.3); axR.set_xlim(-220, 220)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", path)


def fig_acoustic_delay(path, D_mm=1.0, d_offset_mm=5.0):
    """Schematic: the fixed propagation delay tau_d = d/V (transducer -> beam) shifts the
    whole optical pulse; the rise tau_f shapes it. The apparatus compensates tau_d in the
    pulse timing (delay_Raman_r2 = 0.46 us, delay_Raman_b1 = 0.50 us)."""
    tf = tau_f(D_mm)
    td_r2 = 460.0                                 # ns, delay_Raman_r2 from the .dat (=> d ~2.7 mm)
    td_b1 = 500.0                                 # ns, delay_Raman_b1 from the .dat (=> d ~3.0 mm)
    td_user = d_offset_mm / (aom.V_FUSED_SILICA * 1e-3)   # ns, user's "~5 mm" => ~840 ns
    dt0 = 200.0
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(9.6, 5.4), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 1]))

    # (a) one beam: RF gate vs delivered light, delayed by tau_d and shaped by tau_f
    axA.axvspan(0, dt0, ymin=0.55, ymax=0.95, color="#888", alpha=.7, lw=0)
    axA.text(dt0 / 2, 0.99, "RF gate (TTL) — width $\\delta t$", ha="center", fontsize=9)
    t = np.linspace(-100, td_r2 + dt0 + 6 * tf, 2400)
    light = 0.42 * a_field(t - td_r2, dt0, D_mm)
    axA.plot(t, light, _C[1], lw=2.2)
    axA.fill_between(t, 0, light, color=_C[1], alpha=.18)
    axA.text(td_r2 + dt0 / 2, 0.46, "light at the ion\n(area preserved)", ha="center", fontsize=8.5, color=_C[1])
    axA.annotate("", xy=(td_r2, 0.5), xytext=(0, 0.5), arrowprops=dict(arrowstyle="<->", color="k", lw=1.3))
    axA.text(td_r2 / 2, 0.52, f"$\\tau_d=d/V$  ($d$=2.7 mm $\\to$ {td_r2:.0f} ns;  $d$=5 mm $\\to$ {td_user:.0f} ns)",
             ha="center", fontsize=8.5)
    axA.set_ylim(0, 1.12); axA.set_yticks([]); axA.set_ylabel("R2 (single-pass)")
    axA.set_title("(a) finite sound speed = fixed DELAY $\\tau_d$ (transducer$\\to$beam) + rise $\\tau_f$ (beam size)")

    # (b) two beams must OVERLAP: differential delay tau_d(B1)-tau_d(R2) compensated by timing
    lr2 = 0.42 * a_field(t - td_r2, dt0, D_mm)
    lb1 = 0.42 * a_field(t - td_b1, dt0, D_mm)
    axB.plot(t, lr2, _C[1], lw=2.2, label=f"R2 light (delay {td_r2:.0f} ns)")
    axB.plot(t, lb1, _C[3], lw=2.2, label=f"B1 light if also switched (delay {td_b1:.0f} ns)")
    axB.fill_between(t, 0, np.minimum(lr2, lb1), color="0.5", alpha=.3)
    axB.annotate("", xy=(td_b1, 0.30), xytext=(td_r2, 0.30), arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
    axB.text((td_r2 + td_b1) / 2, 0.33, f"$\\Delta\\tau_d$={td_b1 - td_r2:.0f} ns\n(must be matched)",
             ha="center", fontsize=8)
    axB.text(td_r2 + dt0 + 40, 0.30, "Raman coupling $\\propto$ OVERLAP\n(only where both beams are on)",
             fontsize=8.5, color="0.3")
    axB.set_ylim(0, 0.55); axB.set_yticks([]); axB.set_ylabel("two-beam overlap")
    axB.set_xlabel("time [ns]"); axB.set_xlim(-100, td_b1 + dt0 + 320); axB.legend(fontsize=8, loc="upper right")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", path)


def fig_detuning_comb(path, D_mm=1.0, N=50, f_lf=1.3001, Dt=0.769172, dt0=60e-3):
    """Frequency-domain fingerprint. A train of N pulses (period Dt) gives a Floquet comb:
    P(delta) = |array factor|^2 * |single-pulse spectrum|^2.  The array factor (teeth at
    delta = k f_lf, width ~1/(N Dt)) is the SAME with/without the AOM effect; the single-pulse
    envelope (aom.comb_envelope, with/without the Gaussian roll-off) differs. All times in us,
    frequencies in MHz."""
    tf = tau_f(D_mm) * 1e-3                        # ns -> us
    _env = np.vectorize(aom.comb_envelope)

    def array_factor(d):
        num = np.sin(N * np.pi * d * Dt)
        den = N * np.sin(np.pi * d * Dt)
        return np.where(np.abs(den) < 1e-9, 1.0, num / den) ** 2

    def env_ideal(d, dt):
        return _env(d, dt, tf, False)

    def env_aom(d, dt):
        return _env(d, dt, tf, True)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.8, 4.6))

    # (a) the measured lineshape at dt = 60 ns -------------------------------
    d = np.linspace(-6.2 * f_lf, 6.2 * f_lf, 20000)
    af = array_factor(d)
    axL.plot(d / f_lf, af * env_ideal(d, dt0), "0.6", lw=1.0,
             label="without AOM effect (rect $\\delta t$)")
    axL.plot(d / f_lf, af * env_aom(d, dt0), _C[1], lw=1.0,
             label="with AOM transit (rounded)")
    axL.plot(d / f_lf, env_ideal(d, dt0), "0.6", lw=1.6, ls="--", alpha=.9)
    axL.plot(d / f_lf, env_aom(d, dt0), _C[1], lw=1.8, ls="--", alpha=.9)
    axL.text(3.4, 0.55, "sinc envelope\n(no effect)", fontsize=8, color="0.4")
    axL.text(0.25, 0.30, "Gaussian-narrowed\nenvelope (AOM)", fontsize=8, color=_C[1])
    axL.set_xlabel("drive detuning  $\\delta / f_{\\rm lf}$"); axL.set_ylabel("$P_\\downarrow$ (norm.)")
    axL.set_title(f"(a) detuning-scan comb at $\\delta t$=60 ns ($\\tau_f$={tf*1e3:.0f} ns, N={N})")
    axL.legend(fontsize=8, loc="upper right"); axL.grid(alpha=.3); axL.set_xlim(-6, 6); axL.set_ylim(0, 1.05)

    # (b) comb spectral width vs dt: 1/e half-width of the envelope ----------
    def width_1e(dt, which):
        return aom.comb_halfwidth(dt, tf, aom=(which == "aom"), fmax_units=60 * f_lf) / f_lf

    dts = np.logspace(np.log10(5e-3), np.log10(2.0), 90)   # us
    wi = np.array([width_1e(x, "ideal") for x in dts])
    wa = np.array([width_1e(x, "aom") for x in dts])
    axR.axvspan(5, 100, color="0.85", alpha=.7, lw=0)
    axR.text(22, 1.2, "scan\n0-100 ns", fontsize=8, ha="center", color="0.3")
    axR.plot(dts * 1e3, wi, "0.5", lw=2.2, label="without effect $\\propto 1/\\delta t$ (keeps broadening)")
    axR.plot(dts * 1e3, wa, _C[1], lw=2.4, label="with AOM transit (saturates)")
    floor = 1.0 / (np.pi * tf * np.sqrt(2) * f_lf)
    axR.axhline(floor, color=_C[1], lw=1, alpha=.6)
    axR.text(300, floor + 0.4, f"floor $\\approx${floor:.1f} $f_{{\\rm lf}}$", fontsize=8, color=_C[1])
    axR.set_xscale("log"); axR.set_yscale("log"); axR.set_xlim(5, 2000); axR.set_ylim(0.8, 60)
    axR.set_xlabel("nominal pulse width  $\\delta t$  [ns]")
    axR.set_ylabel("comb half-width  [units of $f_{\\rm lf}$]")
    axR.set_title("(b) comb width: broadens (ideal) vs saturates (AOM)")
    axR.legend(fontsize=8, loc="upper right"); axR.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", path)


def fig_strobe_sequence(path, D_mm=1.0, DELTA_t=769.0, dt0=60.0, td=460.0, N=50, n_show=3):
    """Annotated timing diagram of a steady-state section of the active phase grating, mapping
    the script parameters onto B1 (always on), the R2 TTL gate, the delivered R2 light, and the
    Rabi rate the ion actually sees. Times in ns; the x-axis is shown in us.

    The fixed propagation delay tau_d (= delay_Raman_r2) is a uniform, compensated shift, so for
    readability the light/Rabi lanes are drawn ALIGNED under their gates (the delay is annotated,
    not offset; see aom_acoustic_delay.png for the delay itself). Only the rise tau_f (which
    rounds and lengthens each pulse) is rendered, since that is the per-pulse shaping effect."""
    tf = tau_f(D_mm)                                  # ns
    gates = [k * DELTA_t for k in range(n_show)]
    t = np.linspace(-0.45 * DELTA_t, (n_show - 1) * DELTA_t + td + dt0 + 7 * tf, 5000)
    tu = t / 1e3                                      # us for plotting
    colB1, colR2, colRabi = "#2ca02c", _C[1], "#7b3294"
    we = w_equiv(dt0, D_mm, 1)

    fig, (axB1, axR2, axL, axR) = plt.subplots(
        4, 1, figsize=(12.0, 6.6), sharex=True,
        gridspec_kw=dict(height_ratios=[0.5, 1.0, 1.0, 1.25], hspace=0.16))
    # dotted vertical guides at the TTL gate times -> the light/Rabi pulses sit one tau_d to the right
    for ax in (axB1, axR2, axL, axR):
        for g in gates:
            ax.axvline(g / 1e3, color="0.8", lw=0.8, ls=":", zorder=0)

    # --- Lane 1: B1 always on ----------------------------------------------
    axB1.axhspan(0.0, 1.0, xmin=0.0, xmax=1.0, color=colB1, alpha=0.45)
    axB1.set_ylim(0, 1.3); axB1.set_yticks([]); axB1.set_ylabel("B1\n(blue Raman)", fontsize=9)
    axB1.text(tu[len(tu) // 2], 1.12,
              "continuous for the whole train:  turn(b1,1) … turn(b1,0)   (double-pass AOM, NOT gated)",
              ha="center", fontsize=8.5, color=colB1)

    # --- Lane 2: R2 TTL gate -----------------------------------------------
    for g in gates:
        axR2.fill_between([g / 1e3, (g + dt0) / 1e3], 0, 1, color=colR2, alpha=0.9, lw=0)
    axR2.set_ylim(0, 1.5); axR2.set_yticks([]); axR2.set_ylabel("R2 TTL\ngate", fontsize=9)
    # delta_t (on)
    axR2.annotate("", xy=(gates[0] / 1e3, 1.12), xytext=((gates[0] + dt0) / 1e3, 1.12),
                  arrowprops=dict(arrowstyle="<->", color="k", lw=1))
    axR2.text((gates[0] + dt0) / 1e3 + 0.02, 1.12, "delta_t  (60 ns; scanned 0–100 ns)", fontsize=8.5, va="center")
    # DELTA_t (period)
    axR2.annotate("", xy=(gates[0] / 1e3, 1.34), xytext=(gates[1] / 1e3, 1.34),
                  arrowprops=dict(arrowstyle="<->", color="0.3", lw=1))
    axR2.text((gates[0] + gates[1]) / 2 / 1e3, 1.38, "DELTA_t = 0.769 µs  (= 1/f_lf)",
              ha="center", fontsize=8.5, color="0.3")
    # gap = DELTA_t - delta_t
    axR2.annotate("", xy=((gates[0] + dt0) / 1e3, 0.5), xytext=(gates[1] / 1e3, 0.5),
                  arrowprops=dict(arrowstyle="<->", color="0.45", lw=0.9))
    axR2.text((gates[0] + dt0 + gates[1]) / 2 / 1e3, 0.6, "sleep(DELTA_t − delta_t)", ha="center",
              fontsize=7.6, color="0.45")
    axR2.text(gates[-1] / 1e3 + 0.13, 0.5, "…  × N_strobo = 50", fontsize=9, va="center", color=colR2)

    # --- Lane 3: R2 light at the ion (delayed by tau_d, rounded by tau_f) ----
    for g in gates:
        light = a_field(t - (g + td), dt0, D_mm)
        axL.plot(tu, light, colR2, lw=1.8)
        axL.fill_between(tu, 0, light, color=colR2, alpha=0.18)
    axL.set_ylim(0, 1.3); axL.set_yticks([]); axL.set_ylabel("R2 light\nat ion", fontsize=9)
    axL.annotate("", xy=((gates[0] + td) / 1e3, 0.52), xytext=(gates[0] / 1e3, 0.52),
                 arrowprops=dict(arrowstyle="->", color="0.25", lw=1.3))
    axL.text(gates[0] / 1e3 + 0.01, 0.62,
             "τ_d = delay_Raman_r2 ≈ 0.46 µs\n(propagation delay, compensated)", fontsize=7.6, color="0.25")
    axL.text((gates[1] + td) / 1e3 + 0.02, 0.55,
             f"rounded by rise τ_f={tf:.0f} ns →\neffective width w_eff≈{we:.0f} ns", fontsize=8, color=colR2)

    # --- Lane 4: Rabi rate on the ion (also delayed by tau_d) ---------------
    for g in gates:
        rabi = a_field(t - (g + td), dt0, D_mm)
        axR.plot(tu, rabi, colRabi, lw=2.0)
        axR.fill_between(tu, 0, rabi, color=colRabi, alpha=0.2)
        axR.fill_between([(g + td) / 1e3, (g + td + dt0) / 1e3], 0, 1, color="0.5", alpha=0.18, lw=0)  # ideal box
        axR.plot([(g + td) / 1e3, (g + td) / 1e3, (g + td + dt0) / 1e3, (g + td + dt0) / 1e3],
                 [0, 1, 1, 0], "0.5", lw=0.9, ls="--")
    axR.set_ylim(0, 1.3); axR.set_yticks([]); axR.set_ylabel("Rabi rate\nΩ(t)/Ω₀ on ion", fontsize=9)
    axR.set_xlabel("time  [µs]")
    axR.text(gates[1] / 1e3, 0.66,
             "Ω(t) ∝ E_B1·E_R2 = (B1 const)·a(t)\npeak = Ω₀·erf(δt/2τ_f) ≈ 0.39 Ω₀\n"
             "area = Ω₀·δt  (PRESERVED, = grey box)", fontsize=8, color=colRabi)
    axR.text(gates[0] / 1e3 - 0.1, 1.12, "dashed grey = ideal\ninstant-switch pulse\n(same area)",
             fontsize=7.4, color="0.4")

    axB1.set_title("Active phase grating (active_phase_grating_Laser): one section of the N_strobo train  "
                   "—  DDS fr_oc_strobo=1775.46 MHz, pwr_strobo=0.46", fontsize=10)
    axR.set_xlim(tu[0], tu[-1])
    fig.subplots_adjust(top=0.93, bottom=0.09, left=0.08, right=0.98)
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", path)


def write_table(path, D_mm=1.0):
    """Dump the note's numeric tables so they are reproducible/verifiable (not just baked
    into the markdown). One row per nominal gate width, for D = D_mm."""
    gates = [10, 30, 60, 100, 200, 400, 1000]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        # standard CSV header on row 1 (no leading '#', so strict parsers read it cleanly);
        # provenance carried in dedicated columns repeated per row.
        w.writerow(["delta_t_ns", "peak_n1", "R_area_n1", "R_area_n2", "R_area_n3",
                    "w_eff_n1_ns", "fwhm_n1_ns", "w_eff_n3_ns", "D_mm", "tau_f_ns"])
        tfn = round(tau_f(D_mm), 2)
        for g in gates:
            w.writerow([g,
                        round(peak_eff(g, D_mm, 1), 4),
                        round(area_eff(g, D_mm, 1), 4),
                        round(area_eff(g, D_mm, 2), 4),
                        round(area_eff(g, D_mm, 3), 4),
                        round(w_equiv(g, D_mm, 1), 1),
                        round(fwhm(g, D_mm, 1), 1),
                        round(w_equiv(g, D_mm, 3), 1),
                        D_mm, tfn])
    print("wrote", path)


if __name__ == "__main__":
    fig_efficiency(os.path.join(_HERE, "aom_single_vs_double_pass.png"))
    fig_effective_length(os.path.join(_HERE, "aom_effective_pulse_length.png"))
    fig_flop(os.path.join(_HERE, "aom_flop_single_vs_double_pass.png"))
    fig_acoustic_delay(os.path.join(_HERE, "aom_acoustic_delay.png"))
    fig_detuning_comb(os.path.join(_HERE, "aom_detuning_comb.png"))
    fig_strobe_sequence(os.path.join(_HERE, "aom_strobe_sequence.png"))
    write_table(os.path.join(_HERE, "aom_rise_table.csv"))
