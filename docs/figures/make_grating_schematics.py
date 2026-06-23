"""
Generate the pedagogical schematics embedded in
docs/notes/strobo_grating_transfer_function.md:
  grating_pulse_sequence.png  -- the stroboscopic grating train locked to the motion
  grating_ramsey_phasespace.png -- the two-pulse Ramsey sequence + the phase-space disk
These are illustrations (not data); run:  python docs/figures/make_grating_schematics.py
"""
import cmath
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_RED, _BLUE, _GRAY, _GREEN = "#d62728", "#1f77b4", "#888888", "#2ca02c"


def pulse_sequence(path):
    N, Dt, dt = 6, 1.0, 0.11                      # illustrative: 6 pulses, period 1, narrow pulse
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.4, 4.1), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 1.4]))
    # --- top: the drive pulse train ---
    for k in range(N):
        ax1.axvspan(k * Dt, k * Dt + dt, ymin=0.0, ymax=0.85, color=_RED, alpha=0.85)
    ax1.set_ylim(0, 1.25); ax1.set_yticks([]); ax1.set_ylabel("OC drive\n(qubit flip + kick)")
    # annotate one pulse width (delta t) and the period (Delta t)
    ax1.annotate("", xy=(0, 1.02), xytext=(dt, 1.02),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=1))
    ax1.text(dt + 0.04, 1.04, r"$\delta t$ (pulse width)", fontsize=8, va="center")
    ax1.annotate("", xy=(1.0, 1.13), xytext=(2.0, 1.13),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=1))
    ax1.text(1.5, 1.16, r"$\Delta_t$ (strobe period $= 1/f_{\rm lf}$)", fontsize=8, ha="center")
    ax1.set_title(r"Stroboscopic grating: $N$ pulses, one per motional period "
                  r"(the strobe is locked to the motion)", fontsize=10)
    # --- bottom: the axial motion, sampled at the pulse times ---
    t = np.linspace(0, N - 1 + 0.001, 1000)
    ax2.plot(t, np.cos(2 * np.pi * t), color=_BLUE, lw=1.4)
    ax2.plot([k for k in range(N)], [np.cos(2 * np.pi * k) for k in range(N)], "o",
             color=_RED, ms=8, zorder=5, label="pulse times → same motional phase")
    for k in range(N):
        ax2.axvline(k, color=_GRAY, ls=":", lw=0.7)
    ax2.set_ylim(-1.5, 1.5); ax2.set_ylabel("axial motion  $x(t)$")
    ax2.set_xlabel(r"time (units of the motional period $\Delta_t = 1/f_{\rm lf}$)")
    ax2.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def ramsey_phasespace(path, eta=0.45):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2),
                                   gridspec_kw=dict(width_ratios=[1.25, 1]))
    # --- left: the two-pi/2-pulse Ramsey timeline ---
    axL.set_xlim(0, 10); axL.set_ylim(0, 3); axL.axis("off")
    axL.plot([0.3, 9.7], [1, 1], color="k", lw=1)                     # time axis
    def box(x, w, label, col):
        axL.add_patch(plt.Rectangle((x, 0.7), w, 0.6, fc=col, ec="k", alpha=0.85))
        axL.text(x + w / 2, 1.0, label, ha="center", va="center", fontsize=9, color="white")
    box(1.0, 1.4, r"$\pi/2$ ref" + "\n" + r"$\beta_r$", _BLUE)
    box(5.0, 1.7, r"$\pi/2$ grating" + "\n" + r"$\beta_g,\ \varphi$", _RED)
    axL.add_patch(plt.Rectangle((8.0, 0.7), 1.4, 0.6, fc=_GREEN, ec="k", alpha=0.85))
    axL.text(8.7, 1.0, "measure\n$P_\\downarrow$", ha="center", va="center", fontsize=8.5, color="white")
    axL.annotate("", xy=(5.0, 1.9), xytext=(2.4, 1.9),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=1))
    axL.text(3.7, 2.1, r"$T_R$ (Ramsey separation)", ha="center", fontsize=8.5)
    axL.text(5.0, 0.3, r"scan $\varphi$: $P_\downarrow(\varphi)=\frac{1}{2}[1+\mathrm{Re}(e^{i\varphi'}\chi(\Delta\beta))]$",
             ha="center", fontsize=8.5)
    axL.set_title("Two-pulse Ramsey: reference + phase-coherent grating", fontsize=10)
    # --- right: phase space (the ring and the reachable disk) ---
    th = np.linspace(0, 2 * np.pi, 200)
    axR.add_patch(plt.Circle((0, 0), 2 * eta, fc="0.92", ec="none", zorder=0))   # disk |Db|<=2eta
    axR.plot(2 * eta * np.cos(th), 2 * eta * np.sin(th), color="0.55", lw=1.0)
    axR.plot(eta * np.cos(th), eta * np.sin(th), "--", color=_GRAY, lw=1.1)       # ring |b|=eta
    br = eta * np.exp(1j * 2.3); bg = eta * np.exp(1j * 0.5)                       # two kicks on the ring
    for b, c, lab in ((br, _BLUE, r"$\beta_r$"), (bg, _RED, r"$\beta_g$")):
        axR.annotate("", xy=(b.real, b.imag), xytext=(0, 0),
                     arrowprops=dict(arrowstyle="->", color=c, lw=1.8))
        axR.text(b.real * 1.18, b.imag * 1.18, lab, color=c, fontsize=11, ha="center")
    db = bg - br
    axR.annotate("", xy=(db.real, db.imag), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color="k", lw=2.0))
    axR.text(db.real * 1.1 + 0.03, db.imag * 1.1 - 0.06, r"$\Delta\beta=\beta_g-\beta_r$", fontsize=10)
    axR.text(0, -2 * eta - 0.12, r"reachable disk $|\Delta\beta|\leq 2\eta$", ha="center", fontsize=9, color="0.4")
    axR.text(eta * 0.72, eta * 0.72, r"ring $|\beta|=\eta$", fontsize=8, color=_GRAY, rotation=-45)
    L = 2 * eta + 0.25
    axR.set_xlim(-L, L); axR.set_ylim(-L, L); axR.set_aspect("equal")
    axR.set_xlabel(r"$\mathrm{Re}\,\beta$"); axR.set_ylabel(r"$\mathrm{Im}\,\beta$")
    axR.set_title("Phase space: kicks on the ring → disk of differences", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def tomography_sequence(path, gamma=1.3 * cmath.exp(1j * 0.785), eta=0.45):
    """The experimental measurement sequence for ONE probe point of the Wigner-tomography
    scan (state prep -> qubit pi/2 -> conditional displacement D(beta) -> analysis pi/2 ->
    detect), plus the phase-space picture of the scanned probe."""
    fig = plt.figure(figsize=(11.5, 4.4))
    axS = fig.add_axes([0.04, 0.05, 0.63, 0.9]); axS.set_xlim(0, 12); axS.set_ylim(0, 3); axS.axis("off")
    axS.plot([0.2, 11.8], [1.25, 1.25], color="k", lw=1)               # time axis
    def box(x, w, label, col, fs=8.5):
        axS.add_patch(plt.Rectangle((x, 0.82), w, 0.86, fc=col, ec="k", alpha=0.88))
        axS.text(x + w / 2, 1.25, label, ha="center", va="center", fontsize=fs, color="white")
    box(0.4, 2.0, "cool +\nprepare\n" + r"$|\gamma\rangle$", _GRAY)
    box(2.9, 0.9, r"$\pi/2$", _BLUE)
    box(4.2, 2.2, "conditional\ndisplacement\n" + r"$D(\beta)$", _RED)
    box(7.1, 1.3, r"$\pi/2$" + "\n" + r"$(\varphi)$", _BLUE)
    box(9.0, 1.6, "detect\n" + r"$\to P_\downarrow$", _GREEN)
    axS.annotate("", xy=(2.9, 0.55), xytext=(10.6, 0.55),
                 arrowprops=dict(arrowstyle="<->", color="0.4", lw=1))
    axS.text(6.7, 0.30, r"repeat $M$ shots $\to P_\downarrow$;   scan probe "
             r"$\beta=\mathrm{mag}\cdot e^{i\phi_g}$ and readout $\varphi\in\{0,\pi/2\}$",
             ha="center", fontsize=8.5, color="0.25")
    axS.text(6.0, 2.62, "Wigner-tomography measurement sequence (one probe point)",
             ha="center", fontsize=10.5)
    # phase-space inset: the prepared state and the scanned probe displacement
    axP = fig.add_axes([0.73, 0.16, 0.25, 0.66])
    th = np.linspace(0, 2 * np.pi, 200)
    axP.add_patch(plt.Circle((gamma.real, gamma.imag), 0.7, fc=_GRAY, ec="none", alpha=0.35))
    axP.text(gamma.real, gamma.imag, r"$|\gamma\rangle$", ha="center", va="center", fontsize=10)
    beta = 1.6 * cmath.exp(1j * 2.4)
    axP.annotate("", xy=(beta.real, beta.imag), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color=_RED, lw=2))
    axP.text(beta.real * 1.15, beta.imag * 1.15, r"probe $\beta$", color=_RED, fontsize=9, ha="center")
    axP.plot(0, 0, "k+", ms=8)
    axP.set_xlim(-2.6, 2.6); axP.set_ylim(-2.6, 2.6); axP.set_aspect("equal")
    axP.set_xlabel(r"Re", fontsize=8); axP.set_ylabel(r"Im", fontsize=8)
    axP.set_title(r"scan $\beta$ over phase space", fontsize=9)
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    pulse_sequence(os.path.join(_HERE, "grating_pulse_sequence.png"))
    ramsey_phasespace(os.path.join(_HERE, "grating_ramsey_phasespace.png"))
    tomography_sequence(os.path.join(_HERE, "wigner_tomography_sequence.png"))
    print("wrote grating_pulse_sequence.png, grating_ramsey_phasespace.png, "
          "wigner_tomography_sequence.png")
