"""
Raman-beam dephasing — decoherence of a two-photon Raman carrier flop from RELATIVE
optical phase/frequency noise between the two beams.

A TPSR carrier is driven by the BEAT NOTE of the two beams, so the spin phase tracks
the relative optical phase phi(t) = phi_B(t) - phi_R(t). Noise in THIS relative phase
decoheres the flop -- a channel INDEPENDENT of the carrier Debye-Waller motional
dephasing (engines/sideband), the off-resonant scattering (engines/scatter,
raman_optical), and the AC-Stark shift. It is the leading candidate for the
"technical" residual that the OC-flop twin otherwise has to absorb into an
effective nbar (engines/../twin_oc_flop, ADR-0007).

WHERE IT COMES FROM. For beams split from ONE laser and offset by AOMs/DDS, the
common-mode laser phase CANCELS in the difference, so the two-photon FREQUENCY is set
by the (RF-stable) AOM offset, not the optical laser. What survives:
  * PATH IMBALANCE: the two beams reach the ion via different path lengths, so the
    differential samples the laser phase at a delay tau = dL/c. A Lorentzian-linewidth
    laser (phase random walk, <[phi(t)-phi(0)]^2> = 2 pi dnu_laser |t|) then imprints a
    relative-phase variance ~ 2 pi dnu_laser tau (path_imbalance_phase_variance).
  * RF (AOM/DDS) phase noise, fibre noise, beam-pointing -> intensity/phase jitter.

CONTRAST ENVELOPE C(t) (the factor multiplying the coherent oscillation), two limits:
  * white relative-frequency noise / LORENTZIAN mutual line of FWHM dnu:
        C(t) = exp(-pi dnu t)            -> EXPONENTIAL, rate Gamma_phi = pi dnu.
  * quasi-static / slow (1/f) noise, GAUSSIAN relative-frequency spread sigma_nu:
        C(t) = exp(-(t/T2)^2),  T2 = sqrt(2)/(2 pi sigma_nu).
Real Raman noise sits between; the exponential (Lorentzian) limit matches the
phenomenological exponential the rabi fit already uses, so it is the default.

CAPABILITY + DIAGNOSTIC: there is NO measured Raman mutual linewidth in the theses to
anchor against, so this engine is not a sigma-validation. Its job is to convert an
OBSERVED residual decay rate into the mutual linewidth dnu / coherence time T_phi it
implies (mutual_linewidth_from_rate / coherence_time_from_rate), as the alternative
to a hot motional state -- the two are degenerate in a single flop (both reproduce
the same envelope); an independent probe (vary the path imbalance, a spin echo, or a
direct beat-note linewidth) breaks the degeneracy.

Pure Python; frequencies in Hz, times in s.
"""
from __future__ import annotations

import math

C_LIGHT = 2.99792458e8     # m/s


def coherence_lorentzian(t_s: float, mutual_linewidth_hz: float) -> float:
    """Contrast envelope C(t) = exp(-pi * dnu * t) for a Lorentzian mutual (two-photon)
    line of FWHM `mutual_linewidth_hz` (white relative-frequency noise). Exponential."""
    if mutual_linewidth_hz <= 0.0:
        return 1.0
    return math.exp(-math.pi * mutual_linewidth_hz * t_s)


def coherence_gaussian(t_s: float, t2_s: float) -> float:
    """Contrast envelope C(t) = exp(-(t/T2)^2) for quasi-static / slow relative-
    frequency noise (Gaussian). T2 relates to the relative-frequency spread via
    t2_from_sigma_nu."""
    if not math.isfinite(t2_s) or t2_s <= 0.0:
        return 1.0
    return math.exp(-(t_s / t2_s) ** 2)


def dephasing_rate(mutual_linewidth_hz: float) -> float:
    """Exponential contrast-decay rate Gamma_phi = pi * dnu [1/s] (the 1/e rate of the
    Lorentzian envelope)."""
    return math.pi * abs(mutual_linewidth_hz)


def mutual_linewidth_from_rate(gamma_phi_hz: float) -> float:
    """Invert: the mutual (two-photon) linewidth dnu = Gamma_phi / pi [Hz] implied by
    an observed exponential contrast-decay rate. Use this to read out the Raman
    dephasing a measured flop residual corresponds to."""
    return abs(gamma_phi_hz) / math.pi


def coherence_time_from_rate(gamma_phi_hz: float) -> float:
    """Mutual coherence time T_phi = 1/Gamma_phi [s] for an exponential envelope."""
    return 1.0 / gamma_phi_hz if gamma_phi_hz > 0.0 else float("inf")


def t2_from_sigma_nu(sigma_nu_hz: float) -> float:
    """Gaussian-dephasing T2 = sqrt(2)/(2 pi sigma_nu) [s] from the rms relative-
    frequency spread sigma_nu (the quasi-static limit)."""
    return math.sqrt(2.0) / (2.0 * math.pi * sigma_nu_hz) if sigma_nu_hz > 0.0 else float("inf")


def sigma_nu_from_t2(t2_s: float) -> float:
    """Inverse of t2_from_sigma_nu: rms relative-frequency spread from a Gaussian T2."""
    return math.sqrt(2.0) / (2.0 * math.pi * t2_s) if (math.isfinite(t2_s) and t2_s > 0.0) else 0.0


def path_imbalance_phase_variance(laser_linewidth_hz: float, path_diff_m: float) -> float:
    """Relative-phase variance <dphi^2> = 2 pi * dnu_laser * tau imprinted by a
    path-length difference dL (tau = dL/c) for a Lorentzian-linewidth laser (Wiener
    phase diffusion). This is the STATIC, shot-to-shot phase jitter from the common
    laser's path imbalance -> a constant contrast factor exp(-<dphi^2>/2); the
    time-dependent decay during the pulse needs the laser frequency-noise PSD and is
    lumped into the mutual linewidth `dnu` instead. ORDER-OF-MAGNITUDE scale estimate."""
    tau = abs(path_diff_m) / C_LIGHT
    return 2.0 * math.pi * abs(laser_linewidth_hz) * tau


def static_contrast(phase_variance: float) -> float:
    """Constant contrast factor exp(-<dphi^2>/2) from a static (per-shot) relative-
    phase jitter of variance `phase_variance`."""
    return math.exp(-0.5 * phase_variance)


# --- map the dephasing to a beam-path-length variation ---------------------
# The two-photon (carrier) phase = phi_B - phi_R tracks the relative OPTICAL phase of
# the two beams, so a DIFFERENTIAL path-length change dL between the two arms shifts it
# by dphi = (2 pi / lambda) dL  (the FULL optical k, lambda ~ 280 nm; a COMMON-mode
# path change before the beam split cancels in phi_B - phi_R and would couple only
# through the much smaller Delta_k). NOTE: a STATIC dL does NOT dephase a Rabi flop
# (the drive phase only sets the rotation axis; P_flip is phase-independent). What
# decoheres the flop is dL CHANGING during the pulse -> a relative-frequency (two-
# photon detuning) offset df. So these convert the dephasing's frequency/phase budget
# into the length/velocity scale of the DYNAMIC path jitter that would cause it.

def path_variation_for_phase(phase_rad: float, wavelength_m: float) -> float:
    """Differential beam-path-length change dL = (lambda/2pi) * phase imprinting a
    relative-phase slip `phase_rad`. lambda/2pi (~45 nm at 280 nm) is one radian."""
    return wavelength_m * phase_rad / (2.0 * math.pi)


def path_jitter_velocity(freq_hz: float, wavelength_m: float) -> float:
    """Relative path-length drift speed v = lambda * df that produces a two-photon
    frequency offset df (a moving path Doppler-shifts the beat note) [m/s]. The mutual
    linewidth / relative-frequency noise df maps to this path-jitter speed; the path
    must be stable to << lambda over the coherence time 1/df."""
    return wavelength_m * abs(freq_hz)
