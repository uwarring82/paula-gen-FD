"""
Raman off-resonant-scattering + differential-AC-Stark engine for the TPSR (two-
photon stimulated-Raman) carrier drive.

The PAULA OC/CC/AC/ROC "carriers" are driven by two beams (one blue, one red),
both detuned by Delta_R from the 3P_3/2 manifold (raman_detuning_from_p32 = 20
GHz). Two things follow from that single detuning, and BOTH must be in a twin of
a Raman flop:

1. COHERENT two-photon carrier Rabi rate (single dominant fine-structure level,
   valid here because Delta_R = 20 GHz << the P_1/2-P_3/2 fine-structure splitting
   ~2.75 THz):

       Omega = Omega_B Omega_R / (2 Delta_R)                              (rad/s)

   with Omega_B, Omega_R the single-photon Rabi frequencies of the two beams.

2. OFF-RESONANT (spontaneous) photon SCATTERING. Each beam parks a small virtual
   population rho_ee = Omega_single^2/(4 Delta_R^2) in 3P_3/2, which decays at the
   natural linewidth Gamma. Summed over both beams the total scattering rate is

       Gamma_sc = Gamma (Omega_B^2 + Omega_R^2) / (4 Delta_R^2).          (1/s)

   For BALANCED beams (Omega_B = Omega_R) this collapses onto the coherent rate:

       Gamma_sc = (Gamma / Delta_R) Omega          (balanced),

   so the dimensionless scattering-per-pi-pulse is detuning-only:

       P_SE(pi) = Gamma_sc * t_pi = pi * Gamma / Delta_R     (balanced),

   i.e. ~pi * (41.8 MHz / 20 GHz) = 6.6e-3 here. THIS IS THE WHOLE POINT OF A
   LARGE RAMAN DETUNING (Ozeri 2007, Wineland): Gamma_sc/Omega ~ Gamma/Delta_R
   falls as 1/Delta_R, so spontaneous emission is a *small* error at 20 GHz. A
   carrier flop that decays over only a few pi-rotations is therefore NOT
   scattering-limited -- its envelope is dominated by motional/thermal dephasing
   (the carrier Debye-Waller spread) or technical noise, which this engine does
   NOT claim to model. The engine quantifies the scattering floor so the twin can
   separate it from the rest.

DECOHERENCE convention. We model each scattering event as fully destroying the
qubit coherence and depolarising the qubit toward the maximally-mixed state
(P_flip -> 1/2). With equal longitudinal/transverse relaxation T1 = T2 = 1/Gamma_sc
the resonant damped-Rabi solution of the optical Bloch equations has the contrast
envelope e^{-(3/4) Gamma_sc t}; the 3/4 is recorded as CONTRAST_DECAY_FACTOR. This
is the leading-order picture: it ignores the Rayleigh/Raman branching (elastic
Rayleigh dephases without changing state; inelastic Raman can pump the ion OUT of
the qubit manifold -> leakage, not depolarisation) and any beam imbalance beyond
the `balance` handle. There is no measured Raman-scattering decoherence rate in
the theses to anchor it; the BDD benchmark (engines/acstark.py) validates only the
single-beam light-shift SCALE.

The DIFFERENTIAL AC-Stark shift between |down>=|3,3> and |up>=|2,2> (they see Raman
detunings differing by the qubit splitting omega_HF) is delegated to
sideband.raman_differential_stark_factor (leading order omega_HF/Delta_R) and
re-exported here as differential_stark_shift(...) so a Raman-flop twin has one
import.

Pure Python; all public functions take/return ordinary frequencies in Hz (rates in
1/s). Frequencies may be passed as Omega/2pi (Hz) consistently -- the ratios are
scale-free, so Hz vs rad/s cancels except in t_pi (handled explicitly).
"""
from __future__ import annotations

import math

from .sideband import raman_differential_stark_factor

# Contrast-envelope exponent prefactor for a resonant damped Rabi flop with
# T1 = T2 = 1/Gamma_sc (full depolarisation): P_flip ~ (1 - e^{-k Gamma_sc t}cos),
# k = 3/4. Documented in the module docstring.
CONTRAST_DECAY_FACTOR = 0.75


def two_photon_rabi(omega_b_hz: float, omega_r_hz: float, raman_detuning_hz: float) -> float:
    """Coherent TPSR carrier Rabi frequency Omega/2pi [Hz] = Omega_B Omega_R/(2 Delta_R)
    (single dominant fine-structure level). Inputs are the single-beam Rabi
    frequencies Omega_B/2pi, Omega_R/2pi [Hz] and the Raman detuning [Hz]."""
    if raman_detuning_hz == 0.0:
        raise ValueError("Raman detuning must be non-zero (and far-detuned)")
    return omega_b_hz * omega_r_hz / (2.0 * abs(raman_detuning_hz))


def scatter_rate(omega_b_hz: float, omega_r_hz: float, raman_detuning_hz: float,
                 gamma_hz: float) -> float:
    """Total off-resonant photon-scattering RATE [1/s] from both Raman beams:
    R = Gamma_decay * rho_ee = (2 pi gamma) (Omega_B^2 + Omega_R^2)/(4 Delta_R^2),
    rho_ee = far-detuned virtual excited-state population per beam, summed. Inputs
    are Omega/2pi values [Hz] (gamma_hz = Gamma/2pi); the 2pi turns the linewidth
    into the decay rate so the result is a genuine 1/s (consistent with
    scatter_rate_from_rabi). The Delta_R^2-vs-Omega^2 ratio is scale-free, so only
    that one 2pi is needed."""
    if raman_detuning_hz == 0.0:
        raise ValueError("Raman detuning must be non-zero (and far-detuned)")
    return 2.0 * math.pi * gamma_hz * (omega_b_hz ** 2 + omega_r_hz ** 2) / (4.0 * raman_detuning_hz ** 2)


def _imbalance_factor(balance: float) -> float:
    """(1 + r^2)/(2 r) with r = Omega_R/Omega_B (intensity/field ratio of the two
    beams). 1 for balanced beams; >1 otherwise -> MORE scattering per unit coherent
    Rabi. Minimum at r = 1."""
    r = balance
    if r <= 0.0:
        raise ValueError("beam balance ratio must be positive")
    return (1.0 + r * r) / (2.0 * r)


def scatter_rate_from_rabi(rabi_hz: float, raman_detuning_hz: float, gamma_hz: float,
                           balance: float = 1.0) -> float:
    """The PRACTICAL handle: the scattering rate [1/s] expressed through the MEASURED
    two-photon Rabi frequency rabi_hz (= Omega/2pi), since the single-beam Rabi
    frequencies are not directly recorded. For balanced beams

        Gamma_sc = (Gamma / Delta_R) Omega = 2 pi (Gamma/Delta_R) rabi_hz;

    `balance` = Omega_R/Omega_B carries any beam imbalance via (1+r^2)/(2r) >= 1.
    Note rabi_hz is Omega/2pi but Gamma_sc is a true RATE, so the 2pi is restored."""
    if raman_detuning_hz == 0.0:
        raise ValueError("Raman detuning must be non-zero (and far-detuned)")
    omega = 2.0 * math.pi * rabi_hz
    return (gamma_hz / abs(raman_detuning_hz)) * omega * _imbalance_factor(balance)


def se_probability_per_pi(raman_detuning_hz: float, gamma_hz: float,
                          balance: float = 1.0) -> float:
    """Mean number of scattered photons during one pi-pulse, Gamma_sc * t_pi. For
    balanced beams this is detuning-ONLY: P_SE(pi) = pi * Gamma / Delta_R (the
    coherent Rabi rate cancels). The dimensionless figure of merit for a Raman
    gate's spontaneous-emission error floor (Ozeri 2007)."""
    if raman_detuning_hz == 0.0:
        raise ValueError("Raman detuning must be non-zero (and far-detuned)")
    return math.pi * (gamma_hz / abs(raman_detuning_hz)) * _imbalance_factor(balance)


def contrast_decay_rate(rabi_hz: float, raman_detuning_hz: float, gamma_hz: float,
                        balance: float = 1.0) -> float:
    """Flop CONTRAST decay rate [1/s] from scattering alone: CONTRAST_DECAY_FACTOR *
    Gamma_sc (the e-folding rate of the Rabi envelope under full depolarisation,
    T1 = T2 = 1/Gamma_sc). 1/this is the scattering-limited coherence time."""
    return CONTRAST_DECAY_FACTOR * scatter_rate_from_rabi(rabi_hz, raman_detuning_hz, gamma_hz, balance)


def differential_stark_shift(rabi_hz: float, omega_hf_hz: float, raman_detuning_hz: float) -> float:
    """Differential AC-Stark shift [Hz] of the |down>-|up> qubit splitting under the
    Raman beams, as an effective carrier detuning: delta_AC = (omega_HF/Delta_R) *
    Omega (leading order; see sideband.raman_differential_stark_factor for the
    caveats -- factor-of-2 convention, CG / P_1/2 / beam-balance corrections)."""
    return raman_differential_stark_factor(omega_hf_hz, raman_detuning_hz) * rabi_hz


def flip_probability(t_s: float, rabi_hz: float, gamma_sc_hz: float = 0.0,
                     stark_detuning_hz: float = 0.0) -> float:
    """Forward twin of a Raman carrier flop: P(flipped) after a square pulse of
    duration t_s, including the differential AC-Stark detuning and the scattering
    decoherence.

        P_flip(t) = (A/2) [1 - e^{-k Gamma_sc t} cos(2 pi Omega_eff t)],

    with Omega_eff = hypot(Omega, delta_AC) the generalised Rabi frequency (the
    AC-Stark shift speeds the oscillation and CAPS its amplitude at A = Omega^2/
    Omega_eff^2 < 1), k = CONTRAST_DECAY_FACTOR, and Gamma_sc the scattering rate.
    The decohered population settles to A/2 (mixed within the addressed two-level
    subspace). Leading order: the detuning x damping cross-terms are dropped (valid
    for delta_AC, Gamma_sc << Omega, as here). All frequencies in Hz, t in s."""
    eff = math.hypot(rabi_hz, stark_detuning_hz)
    if eff == 0.0:
        return 0.0
    amp = (rabi_hz * rabi_hz) / (eff * eff)
    env = math.exp(-CONTRAST_DECAY_FACTOR * gamma_sc_hz * t_s) if gamma_sc_hz > 0.0 else 1.0
    return 0.5 * amp * (1.0 - env * math.cos(2.0 * math.pi * eff * t_s))


class RamanScatter:
    """Bundle the Raman-flop loss channels around one detuning + linewidth, the way
    a twin consumes them. Construct directly or via from_ledger."""

    def __init__(self, raman_detuning_hz: float, gamma_hz: float, omega_hf_hz: float):
        self.delta_r = float(raman_detuning_hz)
        self.gamma = float(gamma_hz)
        self.omega_hf = float(omega_hf_hz)

    @classmethod
    def from_ledger(cls, ledger,
                    detuning_name: str = "raman_detuning_from_p32",
                    gamma_name: str = "mg_p32_natural_linewidth",
                    hf_name: str = "hyperfine_splitting_25mg_f2_f3"):
        """Consume the Raman detuning, the 3P_3/2 linewidth, and the qubit hyperfine
        splitting -- all `input` (wall-enforced via input_quantity)."""
        return cls(
            raman_detuning_hz=ledger.input_quantity(detuning_name).value,
            gamma_hz=ledger.input_quantity(gamma_name).value,
            omega_hf_hz=ledger.input_quantity(hf_name).value,
        )

    def scatter_rate(self, rabi_hz: float, balance: float = 1.0) -> float:
        return scatter_rate_from_rabi(rabi_hz, self.delta_r, self.gamma, balance)

    def se_per_pi(self, balance: float = 1.0) -> float:
        return se_probability_per_pi(self.delta_r, self.gamma, balance)

    def contrast_decay_rate(self, rabi_hz: float, balance: float = 1.0) -> float:
        return contrast_decay_rate(rabi_hz, self.delta_r, self.gamma, balance)

    def stark_detuning(self, rabi_hz: float) -> float:
        return differential_stark_shift(rabi_hz, self.omega_hf, self.delta_r)

    def flip_probability(self, t_s: float, rabi_hz: float, balance: float = 1.0,
                         with_stark: bool = True) -> float:
        """Forward flop including this object's scattering + (optionally) AC-Stark."""
        g = self.scatter_rate(rabi_hz, balance)
        d = self.stark_detuning(rabi_hz) if with_stark else 0.0
        return flip_probability(t_s, rabi_hz, gamma_sc_hz=g, stark_detuning_hz=d)
