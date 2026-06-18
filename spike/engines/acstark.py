"""
AC-Stark (light-shift) engine — the far-detuned single-beam light shift.

A laser beam of saturation parameter s = I/Isat, detuned by `delta` from a
transition of natural linewidth Gamma, shifts the levels by the AC-Stark (light)
shift

    delta_AC = Omega^2 / (4 delta) = s * Gamma^2 / (8 delta)      [Hz]

(Clos 2017 Eq. 2.2.24, far-detuned limit of E = +-(delta/2 + Omega^2/4delta);
Omega^2 = s Gamma^2 / 2). This is a COHERENT shift and is only meaningful FAR
from resonance (|delta| >> Gamma). A near-resonant beam SCATTERS instead (that is
the cooling engine's domain), so the formula must not be applied there. Hasse
2025 confirms the split experimentally: the far-detuned BDD beam (delta = -10
Gamma) shifts the cycling transition by ~2pi x 10 MHz, while the resonant
repumpers RD/RP "do not induce a significant ac Stark shift".
"""
from __future__ import annotations


def ac_stark_shift(s: float, gamma_hz: float, detuning_hz: float) -> float:
    """Far-detuned light shift of the driven ground state, delta_AC = s*Gamma^2/(8 delta)
    [Hz]. Sign follows `detuning_hz`: a red-detuned beam (delta < 0) shifts the GROUND
    LEVEL down (negative), which RAISES the transition frequency -- i.e. the resonance
    moves up, consistent with Hasse's "shifts the cycling transition to higher
    frequencies" for BDD. Validations compare magnitudes."""
    if detuning_hz == 0.0:
        raise ValueError("detuning must be non-zero (and far-detuned) for the AC-Stark shift")
    return s * gamma_hz * gamma_hz / (8.0 * detuning_hz)


def is_far_detuned(gamma_hz: float, detuning_hz: float, ratio: float = 3.0) -> bool:
    """Whether the coherent AC-Stark formula applies (|delta| >= ratio*Gamma). A
    detuning-only heuristic (ratio=3): it cleanly separates the five beams here (BDD
    at 10 Gamma vs the rest at <= 0.5 Gamma); near resonance the beam scatters rather
    than coherently shifting (-> the cooling engine)."""
    return abs(detuning_hz) >= ratio * abs(gamma_hz)
