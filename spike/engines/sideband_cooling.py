"""
Resolved-sideband cooling engine — the ground-state cooling limit and its
consistency with the achieved mean occupations.

Pulsed RSB cooling (Thomm 2021: red-sideband pulses + RD/RP repump, after Doppler
cooling) is limited by OFF-RESONANT excitation of the carrier (and the blue
sideband) while driving the red sideband: during a red-sideband pi-pulse the
carrier, detuned by the mode frequency omega, is weakly excited with amplitude
~ Omega/(2 omega), and the subsequent spontaneous emission heats. Balancing this
against the per-pulse cooling gives a steady-state floor that scales as

    n_bar_min  ~  alpha * (kappa / (2 omega))^2 ,

where kappa is the effective cooling linewidth/rate (the carrier Rabi rate Omega
for pulsed cooling, or the repump-broadened linewidth Gamma_eff for continuous
cooling) and alpha is an O(1) geometry/recoil factor (resolved-sideband
leading order, LBMW03). Higher mode frequency -> deeper cooling.

This engine is a DIAGNOSTIC, not a sigma test: the achieved n_bar (Thomm 0.07-0.11)
are protocol-limited (finite pulse number, per-mode heating), not a single common
floor, so the cleanest move is to INVERT each measured n_bar for the implied
kappa/(2 omega) and check the picture is consistent (kappa ~ a fraction of omega,
i.e. genuinely in the resolved-sideband regime). Pure Python (math only).
"""
from __future__ import annotations

import math
from dataclasses import dataclass


def offres_carrier_excitation(rabi_carrier_hz: float, omega_hz: float) -> float:
    """Off-resonant carrier excitation probability while driving a sideband:
    ~ (Omega / (2 omega))^2 (the leading generalised-Rabi term for detuning omega).
    This is the per-pulse heating channel that limits pulsed RSB cooling."""
    if omega_hz <= 0.0:
        return float("nan")
    return (rabi_carrier_hz / (2.0 * omega_hz)) ** 2


def rsb_cooling_limit(omega_hz: float, kappa_hz: float, alpha: float = 1.0) -> float:
    """Resolved-sideband cooling floor n_bar_min ~ alpha*(kappa/(2 omega))^2, with
    kappa the effective cooling linewidth/rate and alpha an O(1) factor (default the
    leading-order 1). Valid in the resolved regime kappa < omega."""
    if omega_hz <= 0.0:
        return float("nan")
    return alpha * (kappa_hz / (2.0 * omega_hz)) ** 2


def kappa_from_nbar(nbar: float, omega_hz: float, alpha: float = 1.0) -> float:
    """Invert the floor: the effective cooling linewidth/rate kappa implied by an
    achieved n_bar at mode frequency omega — kappa = 2 omega sqrt(n_bar/alpha). If the
    same kappa explains all modes the cooling is at a common RSB floor; a spread means
    the achieved n_bar are protocol/heating-limited per mode."""
    if nbar < 0.0 or alpha <= 0.0:
        return float("nan")
    return 2.0 * omega_hz * math.sqrt(nbar / alpha)


def mean_from_ground_state_prob(p0: float) -> float:
    """For a (near-)thermal residual, n_bar = 1/P(n=0) - 1 (Bose-Einstein), a quick
    cross-check of a reported ground-state probability against a reported n_bar."""
    if not 0.0 < p0 <= 1.0:
        return float("nan")
    return 1.0 / p0 - 1.0


@dataclass
class SidebandCooling:
    """Per-mode RSB cooling diagnostic from the ledger's RSB-cooled occupations."""
    modes: tuple   # ((label, omega_hz, nbar), ...)

    @classmethod
    def from_ledger(cls, ledger):
        spec = [("axial lf", "omega_z_axial_com_25mg", "mg_rsb_cooled_nbar_axial_lf_25mg"),
                ("radial mf", "omega_radial_mf_25mg", "mg_rsb_cooled_nbar_radial_mf_25mg"),
                ("radial hf", "omega_radial_hf_25mg", "mg_rsb_cooled_nbar_radial_hf_25mg")]
        modes = tuple((lab, ledger.value(o), ledger.value(n))
                      for lab, o, n in spec if o in ledger and n in ledger)
        return cls(modes=modes)

    def inferred_kappa(self):
        """[(label, omega_hz, nbar, kappa_hz, kappa_over_omega)]."""
        out = []
        for lab, omega, nbar in self.modes:
            kap = kappa_from_nbar(nbar, omega)
            out.append((lab, omega, nbar, kap, kap / omega if omega else float("nan")))
        return out
