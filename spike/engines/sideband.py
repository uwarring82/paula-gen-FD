"""
Sideband engine — absolute Lamb-Dicke parameters and first-order sideband Rabi
rates for the Raman (TPSR) beam combinations, plus the Raman differential
AC-Stark shift that moves the sideband resonance.

Absolute Lamb-Dicke parameter for a combination onto a mode:

    eta(comb, mode, omega) = |Delta_k/k . e_mode| * k*z_bar(omega),

where Delta_k/k is the magnitude-bearing combination vector (raman_*_combination_25mg),
e_mode the mode axis (mode_axes from the projection engine), and
z_bar = sqrt(hbar/(2 m omega)) the mode ground-state extent (k*z_bar ~ omega^-1/2).
Anchored to the measured axial value eta = 0.32 for OC (|Delta_k/k| = sqrt2) onto
the lf mode at omega_ref/(2pi) = 1.92 MHz (raman_axial_lamb_dicke_25mg):

    k*z_bar(omega) = (eta_ref / |Delta_k_OC/k|) * sqrt(omega_ref / omega).

First-order sideband Rabi rate (Lamb-Dicke regime):

    Omega_{n,n+1} = eta sqrt(n+1) Omega_0   (blue);   Omega_{n,n-1} = eta sqrt(n) Omega_0 (red).

Raman differential AC-Stark shift (ORDER-OF-MAGNITUDE). The far-detuned Raman
beams (Delta_R = 20 GHz) light-shift both qubit states; because |down>=|3,3> and
|up>=|2,2> see Raman detunings differing by the qubit splitting (~1.775 GHz, ~0.8%
below the zero-field 3|A| = 1.789 GHz used here as omega_HF), the large common shift
largely cancels and only a fraction ~ omega_HF/Delta_R survives:

    delta_AC_diff ~ (omega_HF / Delta_R) * Omega_0 ~ 0.05-0.09 * Omega_0.

This is a LEADING-ORDER estimate good to ~a factor of 2: the prefactor (0.5-1.0) is
convention-dependent (a plain balanced two-state model gives omega_HF/(2 Delta_R) ~
0.045), and it also ignores Clebsch-Gordan differences, the opposite-sign P_1/2
fine-structure contribution (~Delta_R/Delta_FS), and beam imbalance. There is NO
measured Raman differential shift to anchor it; the BDD benchmark validates only the
single-beam light-shift SCALE (Omega^2/4delta), not this differential ratio.
"""
from __future__ import annotations

import math

from .projection import _dot, mode_axes

_COMBO_RECORDS = {
    "CC": "raman_cc_combination_25mg",
    "OC": "raman_oc_combination_25mg",
    "AC": "raman_ac_combination_25mg",
    "ROC": "raman_roc_combination_25mg",
}
_REF_OMEGA_HZ = 1.92e6     # reference axial frequency where the eta_ref anchor is quoted


class Sideband:
    """Absolute Lamb-Dicke + sideband Rabi for the four TPSR combinations."""

    def __init__(self, dk_vectors: dict, tilt_deg: float, ref_eta: float,
                 ref_omega_hz: float = _REF_OMEGA_HZ,
                 ref_combo: str = "OC", ref_mode: str = "lf"):
        self.axes = mode_axes(tilt_deg)
        self.dk = {k: tuple(float(c) for c in v) for k, v in dk_vectors.items()}
        # anchor projection derived from the actual ledger vectors (not a hardcoded
        # sqrt2): eta_ref = |Delta_k_ref . e_ref| * k*z_bar(omega_ref).
        ref_proj = abs(_dot(self.dk[ref_combo], self.axes[ref_mode]))
        self.kzbar_ref = ref_eta / ref_proj
        self.ref_omega_hz = float(ref_omega_hz)

    @classmethod
    def from_ledger(cls, ledger, tilt_name: str = "radial_mode_tilt_25mg",
                    eta_name: str = "raman_axial_lamb_dicke_25mg", combos=None):
        """Consume the Delta_k/k vectors + the radial tilt + the axial Lamb-Dicke
        anchor, all `input` (wall-enforced). NOTE: the eta = 0.32 anchor is from the
        Wittemer/Clos era (axial ~1.92 MHz); transferring it to the freddy lf mode
        (~1.3 MHz) via omega^-1/2 assumes the OC (orthogonal) Raman beam geometry,
        hence |Delta_k|, is unchanged across generations."""
        combos = combos or _COMBO_RECORDS
        tilt = ledger.input_quantity(tilt_name).value
        dk = {label: ledger.input_quantity(name).value for label, name in combos.items()}
        ref_eta = ledger.input_quantity(eta_name).value
        return cls(dk, tilt, ref_eta=ref_eta)

    def lamb_dicke(self, comb: str, mode: str, omega_hz: float) -> float:
        """Absolute Lamb-Dicke parameter eta = |Delta_k/k . e_mode| * k*z_bar(omega)."""
        if omega_hz <= 0.0:
            raise ValueError(f"mode frequency must be positive, got {omega_hz}")
        proj = abs(_dot(self.dk[comb], self.axes[mode]))
        return proj * self.kzbar_ref * math.sqrt(self.ref_omega_hz / omega_hz)

    def sideband_rabi(self, comb: str, mode: str, omega_hz: float,
                      carrier_rabi: float, n: int = 0, order: str = "blue") -> float:
        """First-order sideband Rabi rate Omega_{n,n+-1} = eta sqrt(n+1 | n) Omega_0."""
        eta = self.lamb_dicke(comb, mode, omega_hz)
        factor = math.sqrt(n + 1) if order == "blue" else math.sqrt(n)
        return eta * factor * carrier_rabi


def raman_differential_stark_factor(omega_hf_hz: float, raman_detuning_hz: float) -> float:
    """Leading-order Raman differential AC-Stark shift in units of the carrier Rabi
    Omega_0: omega_HF / Delta_R. This is an ORDER-OF-MAGNITUDE (upper) estimate good
    to ~a factor of 2 -- a plain balanced two-state model gives omega_HF/(2 Delta_R);
    the prefactor and CG / P_1/2 / beam-balance corrections are unanchored (no
    measured Raman differential shift exists)."""
    return omega_hf_hz / raman_detuning_hz
