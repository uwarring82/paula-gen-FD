"""
End-to-end DIGITAL-TWIN demonstration of motional-state Wigner tomography for a DISPLACED
(coherent) state of fixed amplitude and phase, via the characteristic-function / Ramsey
route (docs/notes/strobo_grating_transfer_function.md, §6-7; engine: grating_tomography).

PIPELINE (run: python -m spike.twin_wigner_tomography):
  1. PREPARE a known coherent state |gamma>, gamma = r e^{i phi} (fixed amplitude+phase).
  2. MEASURE (twin-simulated): the spin-dependent-force / Ramsey characteristic-function
     readout. For each probe displacement beta on a phase-space grid, reading the qubit in
     two bases gives  P_x = 1/2(1 + C Re chi(beta)),  P_y = 1/2(1 + C Im chi(beta)), with
     contrast C; each is sampled with SHOT NOISE (M repetitions -> binomial counts). The
     raw counts are written to docs/examples/wigner_tomography/.
  3. ANALYSE: reconstruct chi_hat(beta) = [(2 n_x/M - 1) + i(2 n_y/M - 1)] / C from the
     counts, then the Wigner function W(alpha) = (1/pi^2) FT[chi_hat] (engine transform).
  4. VALIDATE: compare W to the analytic W of |gamma>; recover the displacement from the
     reconstruction (centroid <alpha> = gamma) -> amplitude and phase.

IMPORTANT (FAIR): the raw data here is SIMULATED by the digital twin, NOT a real
measurement. It is written under docs/examples/ (separate from sources/data/, the real
DAQ path) and every file is labelled twin-simulated.

SCOPE -- what this DOES and does NOT use (read this). This validates the RECONSTRUCTION
PIPELINE (populations -> chi -> Wigner) and the SDF/conditional-displacement MEASUREMENT
MODEL. It does **NOT run the stroboscopic phase grating** (engines.strobo_sim): the
characteristic function chi(beta) is the ANALYTIC chi of the prepared coherent state, and
the readout is the idealized model P_down = 1/2(1 + C * Re/Im chi(beta)) + shot noise. So
the actual grating physics -- the Floquet comb, the per-cycle kicks beta_k, the finite reach
|Delta beta| <= 2 eta ~ 0.78, the multipulse phase-programming criterion, and the
finite-pulse corrections (all in docs/notes/strobo_grating_transfer_function.md) -- is
BYPASSED. A grating-faithful version would instead generate the chi samples from the
propagator running in the SDF / conditional-displacement mode (an engine extension not yet
built), inheriting those reach/phase/pulse limits. As written, this demonstrates the
tomography MATHEMATICS and the measurement MODEL, not the grating hardware.

Pure Python + matplotlib (Agg). Reproducible (seeded shot noise).
"""
from __future__ import annotations

import cmath
import math
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .engines import grating_tomography as gt  # noqa: E402

_ROOT = Path(__file__).resolve().parent.parent
_DATADIR = _ROOT / "docs" / "examples" / "wigner_tomography"
_FIGDIR = _ROOT / "docs" / "figures"
_RAW = _DATADIR / "coherent_tomography_raw.dat"

# --- the prepared state and the (simulated) measurement settings --------------
GAMMA_R, GAMMA_PHI = 1.3, math.pi / 4.0       # displaced coherent state: amplitude, phase
CONTRAST = 0.90                                # Ramsey fringe contrast (imperfection)
SHOTS = 500                                    # repetitions per (beta, quadrature)
BETA_HALF, BETA_N = 4.0, 25                    # chi-sampling grid: +-4, 25x25 (disk-masked)
ALPHA_HALF, ALPHA_N = 2.5, 51                  # Wigner reconstruction grid
SEED = 20260623
_BLUE, _RED, _GRAY = "#1f77b4", "#d62728", "#888888"


def _binomial_count(p: float, M: int, rng: random.Random) -> int:
    """Shot-noise count out of M (Normal approximation to Binomial(M, p), clamped)."""
    n = round(M * p + math.sqrt(max(M * p * (1.0 - p), 1e-12)) * rng.gauss(0.0, 1.0))
    return min(M, max(0, n))


# --- step 1+2: prepare |gamma>, simulate the SDF/Ramsey chi measurement -------
def simulate(gamma: complex, contrast=CONTRAST, shots=SHOTS,
             beta_half=BETA_HALF, beta_n=BETA_N, seed=SEED):
    """Return raw 'data': list of (beta, n_x, n_y) shot counts on the disk |beta|<=beta_half,
    plus the uniform grid spacing h. chi from the engine; SDF/Ramsey readout + shot noise."""
    rng = random.Random(seed)
    pts, h = gt.square_grid(beta_half, beta_n)
    raw = []
    for b in pts:
        if abs(b) > beta_half:                 # sample a disk, not the square corners
            continue
        chi = gt.chi_coherent(b, gamma)
        px = 0.5 * (1.0 + contrast * chi.real)
        py = 0.5 * (1.0 + contrast * chi.imag)
        raw.append((b, _binomial_count(px, shots, rng), _binomial_count(py, shots, rng)))
    return raw, h


# --- raw data file I/O (self-documenting; twin-simulated provenance) ----------
def write_raw(raw, gamma, h, path=_RAW, contrast=CONTRAST, shots=SHOTS):
    """Write the raw scan log: one row per experimental SETTING (the parameters varied --
    the probe-displacement magnitude and phase, and the readout-quadrature phase) and the
    measured spin-DOWN POPULATION (the dependent, y-axis quantity)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("# Wigner-tomography RAW DATA -- SIMULATED BY THE DIGITAL TWIN (NOT a real "
                "measurement).\n")
        f.write("# Generated by spike/twin_wigner_tomography.py (engine: grating_tomography).\n")
        f.write("# Prepared state: coherent |gamma>, gamma = %.6f %+.6fj "
                "(|gamma|=%.4f, arg=%.4f rad).\n" % (gamma.real, gamma.imag, abs(gamma), cmath.phase(gamma)))
        f.write("# Each row is ONE experimental run. SCANNED PARAMETERS (independent vars):\n")
        f.write("#   disp_mag    = |beta|, the conditional-displacement magnitude (set by the drive area)\n")
        f.write("#   disp_phase  = arg(beta), the grating optical phase phi_g [rad]  (so beta = disp_mag * e^{i disp_phase})\n")
        f.write("#   readout_phase = the analysis-pulse phase varphi [rad]: 0 -> Re-chi quadrature, pi/2 -> Im-chi quadrature\n")
        f.write("# MEASURED (dependent var, y-axis): P_down = spin-down probability = counts_down / M_shots.\n")
        f.write("#   model: P_down = 1/2 (1 + C * [Re or Im] chi(beta)),  C=%.3f, M=%d shots/run.\n" % (contrast, shots))
        f.write("#   beta sampled on the disk |beta|<=%.2f, uniform grid spacing h=%.5f.\n" % (BETA_HALF, h))
        f.write("# columns:  disp_mag  disp_phase  readout_phase  M_shots  counts_down  P_down\n")
        for (b, nx, ny) in raw:
            r, ph = abs(b), cmath.phase(b)
            f.write("%.5f %+.5f %.5f %d %d %.5f\n" % (r, ph, 0.0, shots, nx, nx / shots))            # Re-chi quad
            f.write("%.5f %+.5f %.5f %d %d %.5f\n" % (r, ph, math.pi / 2, shots, ny, ny / shots))     # Im-chi quad
    return path


def read_raw(path=_RAW):
    """Read the scan log back, recombining the two readout quadratures per displacement
    setting into (beta, M, n_x, n_y) rows (+ the grid spacing h from the header)."""
    h, groups = None, {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "grid spacing h=" in line:
                    h = float(line.split("h=")[1].strip().rstrip("."))
                continue
            p = line.split()
            r, ph, vphi, M, n = float(p[0]), float(p[1]), float(p[2]), int(p[3]), int(p[4])
            key = (round(r, 5), round(ph, 5))
            g = groups.setdefault(key, {"beta": cmath.rect(r, ph), "M": M})
            g["nx" if abs(vphi) < 1e-6 else "ny"] = n
    rows = [(g["beta"], g["M"], g["nx"], g["ny"]) for g in groups.values()]
    return rows, h


# --- step 3: reconstruct chi_hat, then the Wigner function --------------------
def reconstruct(rows, contrast, alpha_half=ALPHA_HALF, alpha_n=ALPHA_N, h=None):
    """chi_hat(beta) from the counts, then W(alpha) via the engine's 2-D transform."""
    beta = [r[0] for r in rows]
    chi_hat = [((2.0 * nx / M - 1.0) + 1j * (2.0 * ny / M - 1.0)) / contrast
               for (_b, M, nx, ny) in rows]
    step = 2.0 * alpha_half / (alpha_n - 1)
    ax = [-alpha_half + i * step for i in range(alpha_n)]
    alpha_pts = [complex(x, y) for y in ax for x in ax]
    W = gt.wigner_from_samples(beta, chi_hat, alpha_pts, h * h)
    return chi_hat, beta, W, alpha_pts, ax, step


def recovered_gamma(W, alpha_pts, frac=0.1):
    """Displacement = centroid <alpha> over the W>frac*max region (= gamma for a coherent
    state). Returns the recovered complex amplitude."""
    wmax = max(W)
    num, den = 0j, 0.0
    for w, a in zip(W, alpha_pts):
        if w > frac * wmax:
            num += w * a
            den += w
    return num / den if den else 0j


# --- step 4: figures + report -------------------------------------------------
def make_raw_data_figure(rows, gamma, contrast=CONTRAST, out=_FIGDIR / "twin_wigner_raw_data.png"):
    """The RAW DATA: the measured spin-down POPULATION (y-axis) vs the scanned probe
    displacement, for both readout quadratures."""
    bx = [r[0].real for r in rows]
    by = [r[0].imag for r in rows]
    Px = [nx / M for (_b, M, nx, _ny) in rows]
    Py = [ny / M for (_b, M, _nx, ny) in rows]
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))
    for a, vals, t in [(ax[0], Px, r"$P_\downarrow$, readout $\varphi$=0 (Re-$\chi$ quad.)"),
                       (ax[1], Py, r"$P_\downarrow$, readout $\varphi=\pi/2$ (Im-$\chi$ quad.)")]:
        sc = a.scatter(bx, by, c=vals, cmap="RdBu_r", vmin=0.5 - 0.5 * contrast,
                       vmax=0.5 + 0.5 * contrast, s=22)
        a.set_aspect("equal"); a.set_title(t, fontsize=9.5)
        a.set_xlabel(r"Re $\beta$ (= disp. mag.$\times\cos\phi_g$)"); a.set_ylabel(r"Im $\beta$")
        fig.colorbar(sc, ax=a, fraction=0.046, label=r"$P_\downarrow$")
    # 1-D cut along the real-beta axis: measured population (+ shot error bars) vs the model
    cut = sorted([(b.real, nx / M) for (b, M, nx, _ny) in rows if abs(b.imag) < 1e-6])
    xs = [x for x, _ in cut]
    ys = [y for _, y in cut]
    sig = [math.sqrt(p * (1 - p) / SHOTS) for p in ys]            # binomial shot error on P
    ax[2].errorbar(xs, ys, yerr=sig, fmt="o", color=_RED, ms=4, capsize=2, label="measured $P_\\downarrow$")
    fine = [(-4 + 8 * i / 300) for i in range(301)]
    ax[2].plot(fine, [0.5 * (1 + contrast * gt.chi_coherent(complex(x, 0), gamma).real) for x in fine],
               _BLUE, lw=1.4, label="ideal $P_\\downarrow$")
    ax[2].axhline(0.5, color=_GRAY, lw=0.7, ls=":")
    ax[2].set_xlabel(r"disp. magnitude along Re $\beta$ (Im $\beta$=0)"); ax[2].set_ylabel(r"$P_\downarrow$")
    ax[2].set_ylim(0, 1); ax[2].set_title("raw populations (1-D cut)")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)
    fig.suptitle(r"Step 2: raw measured populations $P_\downarrow$ vs the scanned probe "
                 r"displacement (twin-simulated, shot noise)", y=1.02, fontsize=11)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    return out


def make_wigner_figure(W, alpha_pts, ax_lin, gamma, g_rec,
                       out=_FIGDIR / "twin_wigner_reconstruction.png"):
    na = len(ax_lin)
    def grid(vals):
        return [[vals[j * na + i] for i in range(na)] for j in range(na)]
    W_ana = [gt.wigner_coherent(a, gamma) for a in alpha_pts]
    ext = [ax_lin[0], ax_lin[-1], ax_lin[0], ax_lin[-1]]
    vmax = max(W_ana)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    for a, dat, t in [(ax[0], grid(W_ana), "analytic $W$ (input $|\\gamma\\rangle$)"),
                      (ax[1], grid(W), "reconstructed $W$ (from raw data)"),
                      (ax[2], grid([r - b for r, b in zip(W, W_ana)]), "residual")]:
        im = a.imshow(dat, origin="lower", extent=ext, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        a.set_xlabel("Re $\\alpha$"); a.set_ylabel("Im $\\alpha$"); a.set_aspect("equal")
        fig.colorbar(im, ax=a, fraction=0.046)
        a.set_title(t)
    for a in ax[:2]:
        a.plot(gamma.real, gamma.imag, "k+", ms=11, mew=2)
        a.plot(g_rec.real, g_rec.imag, "x", color="lime", ms=9, mew=2)
    fig.suptitle("Step 4: reconstructed Wigner function vs the prepared state  "
                 "(+ = input $\\gamma$, × = recovered)", y=1.02, fontsize=11)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    return out


def report(gamma, g_rec, W, alpha_pts, n_pts):
    W_ana = [gt.wigner_coherent(a, gamma) for a in alpha_pts]
    err = max(abs(r - b) for r, b in zip(W, W_ana))
    L = ["", "Digital-twin Wigner tomography of a displaced coherent state",
         "  Step 1  prepared state : |gamma>, gamma = %.3f e^{i %.3f} = %.3f %+.3fj"
         % (GAMMA_R, GAMMA_PHI, gamma.real, gamma.imag),
         "  Step 2  measured       : SDF/Ramsey chi-readout, %d probe points, C=%.2f, M=%d "
         "shots/quadrature -> raw counts written" % (n_pts, CONTRAST, SHOTS),
         "  Step 3  reconstructed  : chi_hat = (2 n/M - 1)/C ; W = (1/pi^2) FT[chi_hat]",
         "  Step 4  recovered      : <alpha> = %.3f %+.3fj  -> |gamma|=%.3f (input %.3f), "
         "arg=%.3f rad (input %.3f)" % (g_rec.real, g_rec.imag, abs(g_rec), GAMMA_R,
                                        cmath.phase(g_rec), GAMMA_PHI),
         "          max |W_recon - W_analytic| = %.4f  (peak W = %.3f)" % (err, max(W_ana))]
    return "\n".join(L)


def main(argv=None) -> int:
    gamma = GAMMA_R * cmath.exp(1j * GAMMA_PHI)
    raw, h = simulate(gamma)
    write_raw(raw, gamma, h)
    rows, h_read = read_raw()                  # round-trip through the raw file
    chi_hat, beta, W, alpha_pts, ax_lin, _step = reconstruct(rows, CONTRAST, h=h_read)
    g_rec = recovered_gamma(W, alpha_pts)
    print(report(gamma, g_rec, W, alpha_pts, len(rows)))
    f1 = make_raw_data_figure(rows, gamma)
    f2 = make_wigner_figure(W, alpha_pts, ax_lin, gamma, g_rec)
    print("\nwrote", _RAW.relative_to(_ROOT))
    for f in (f1, f2):
        print("wrote", f.relative_to(_ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
