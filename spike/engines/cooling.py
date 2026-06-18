"""
Cooling engine — two-level photon scattering and the Doppler cooling limit for
the 25Mg+ Blue Doppler transition (S_1/2 -> P_3/2).

Unlike the microwave drive, this is clean textbook physics with no apparatus
model: given the natural linewidth Gamma, the detuning Delta and the saturation
parameter s = I/Isat (all in the ledger from the laser table), the steady-state
scattering rate is

    R_sc = (Gamma/2) * s / (1 + s + (2 Delta/Gamma)^2)

and the Doppler cooling limit (reached at Delta = -Gamma/2) is

    T_D = hbar Gamma / (2 k_B) = h * (Gamma/2pi) / (2 k_B) ~ 1 mK for 25Mg+.

Gamma here is the decay rate (= 2pi * gamma_hz); detuning and linewidth are
passed as ordinary frequencies (the dimensionless ratio 2*Delta/Gamma equals
2*detuning_hz/gamma_hz, so the 2pi cancels).
"""
from __future__ import annotations

import math

from .. import constants as C


def scatter_rate(detuning_hz: float, s: float, gamma_hz: float) -> float:
    """Steady-state photon scattering rate [photons/s].
    R = (Gamma/2) s / (1 + s + (2 detuning/Gamma)^2), Gamma = 2pi*gamma_hz."""
    gamma = 2.0 * math.pi * gamma_hz
    x = 2.0 * detuning_hz / gamma_hz
    return (gamma / 2.0) * s / (1.0 + s + x * x)


def doppler_limit_temperature(gamma_hz: float) -> float:
    """Doppler cooling limit T_D = hbar*Gamma/(2 k_B) = h*gamma_hz/(2 k_B) [K]."""
    return C.H_PLANCK * gamma_hz / (2.0 * C.K_BOLTZMANN)


def optimal_cooling_detuning(gamma_hz: float) -> float:
    """Detuning [Hz] that reaches the Doppler limit: Delta = -Gamma/2 = -gamma_hz/2."""
    return -gamma_hz / 2.0


def mean_occupation(omega_hz: float, temperature: float) -> float:
    """Thermal (Bose-Einstein) mean phonon number of a mode of frequency omega_hz
    at temperature `temperature` [K]: n_bar = 1 / (exp(h*omega/(k_B*T)) - 1)."""
    x = C.H_PLANCK * omega_hz / (C.K_BOLTZMANN * temperature)
    return 1.0 / math.expm1(x)


class DopplerCooling:
    """Doppler cooling / scattering on the S_1/2 -> P_3/2 transition, parameterised
    by the natural linewidth Gamma (= gamma_hz)."""

    def __init__(self, gamma_hz: float):
        self.gamma_hz = float(gamma_hz)

    @classmethod
    def from_ledger(cls, ledger, gamma_name: str = "mg_p32_natural_linewidth"):
        """Build from the natural-linewidth `input` record (wall-enforced)."""
        return cls(gamma_hz=ledger.input_quantity(gamma_name).value)

    def doppler_limit(self) -> float:
        return doppler_limit_temperature(self.gamma_hz)

    def optimal_detuning(self) -> float:
        return optimal_cooling_detuning(self.gamma_hz)

    def doppler_limit_occupation(self, omega_hz: float) -> float:
        """Mean phonon number of a mode of frequency omega_hz when Doppler-cooled
        to the limit T_D (a thermal state at T_D). Equals 1/(exp(2*omega/Gamma)-1),
        so it depends only on the ratio omega/Gamma."""
        return mean_occupation(omega_hz, self.doppler_limit())

    def scatter_rate(self, detuning_hz: float, s: float) -> float:
        return scatter_rate(detuning_hz, s, self.gamma_hz)

    def max_scatter_rate(self) -> float:
        """Saturated, on-resonance scattering rate Gamma/2 [photons/s]."""
        return math.pi * self.gamma_hz
