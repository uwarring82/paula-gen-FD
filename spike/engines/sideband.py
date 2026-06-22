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

    def carrier_thermal_flip(self, comb: str, mode: str, omega_mode_hz: float,
                             carrier_rabi_hz: float, nbar: float, t_s: float,
                             detuning_hz: float = 0.0) -> float:
        """Thermal carrier-flop P(flipped) at duration t_s, using this combination's
        Lamb-Dicke parameter on `mode` (eta = lamb_dicke(comb, mode, omega_mode_hz)).
        See thermal_carrier_flip for the physics."""
        eta = self.lamb_dicke(comb, mode, omega_mode_hz)
        return thermal_carrier_flip(t_s, carrier_rabi_hz, eta, nbar, detuning_hz)


# --- carrier Debye-Waller thermal dephasing --------------------------------
# A CARRIER (Delta n = 0) Raman flop is NOT motion-free at finite Lamb-Dicke eta:
# the carrier Rabi rate depends on the phonon number n through the Debye-Waller
# factor, Omega_{n,n} = Omega_0 <n|e^{i eta (a+a^dag)}|n> = Omega_0 e^{-eta^2/2}
# L_n(eta^2), with L_n the Laguerre polynomial. A THERMAL motional state (mean
# nbar, P_n = nbar^n/(nbar+1)^{n+1}) is then a spread of Rabi frequencies, so the
# flop DEPHASES: F(t) = sum_n P_n cos(2 pi Omega_n t) decays even with no
# scattering and no technical noise. This is the standard thermal-Rabi-flop decay
# (Meekhof 1996, Wineland 1998); it is the LEADING motional contribution to a
# carrier-flop envelope and is ledger-anchorable from eta (sideband) + nbar
# (cooling benchmark). NOTE: for nbar < 1 the distribution is ~all in n=0, so the
# dephasing is SMALL (a few-% beat from the n=1 tail), NOT a clean decay -- the
# Gaussian-spread rate thermal_dephasing_rate OVERSTATES it; use the exact sum
# thermal_carrier_flip for the curve.


def _laguerre(n: int, x: float) -> float:
    """Laguerre polynomial L_n(x), pure-Python three-term recurrence
    (k+1)L_{k+1} = (2k+1-x)L_k - k L_{k-1}, L_0=1, L_1=1-x."""
    if n <= 0:
        return 1.0
    lkm1, lk = 1.0, 1.0 - x
    for k in range(1, n):
        lkm1, lk = lk, ((2 * k + 1 - x) * lk - k * lkm1) / (k + 1)
    return lk


def thermal_pn(nbar: float, n: int) -> float:
    """Thermal (Bose-Einstein) occupation probability P_n = nbar^n/(nbar+1)^{n+1}."""
    if nbar < 0.0:
        raise ValueError("nbar must be non-negative")
    if nbar == 0.0:
        return 1.0 if n == 0 else 0.0
    return nbar ** n / (nbar + 1.0) ** (n + 1)


def carrier_rabi_factor(eta: float, n: int) -> float:
    """Carrier (Delta n = 0) Rabi reduction for phonon number n:
    Omega_{n,n}/Omega_0 = e^{-eta^2/2} L_n(eta^2) (Debye-Waller x Laguerre). At
    n=0 this is the Debye-Waller factor e^{-eta^2/2}; eta->0 gives 1."""
    return math.exp(-eta * eta / 2.0) * _laguerre(n, eta * eta)


def _thermal_nmax(nbar: float, tol: float = 1e-5, cap: int = 5000) -> int:
    """Smallest n such that the thermal tail P(>n) = (nbar/(nbar+1))^{n+1} < tol."""
    if nbar <= 0.0:
        return 0
    r = nbar / (nbar + 1.0)
    n = int(math.ceil(math.log(tol) / math.log(r))) if r > 0 else 0
    return max(1, min(cap, n))


def thermal_carrier_flip(t_s: float, rabi0_hz: float, eta: float, nbar: float,
                         detuning_hz: float = 0.0, n_max: int | None = None) -> float:
    """Thermal carrier-flop P(flipped) after a square pulse of duration t_s:

        P_flip(t) = (1/2) sum_n P_n [1 - (Omega_n^2/Omega_{n,eff}^2)
                                         cos(2 pi Omega_{n,eff} t)],

    P_n thermal, Omega_n = Omega_0 e^{-eta^2/2} L_n(eta^2) the n-dependent carrier
    Rabi rate, Omega_{n,eff} = hypot(Omega_n, detuning) the generalised Rabi (the
    detuning carries the differential AC-Stark shift; each component is capped at
    Omega_n^2/Omega_{n,eff}^2). The sum is renormalised by the included weight.
    nbar=0 reduces to a single (Debye-Waller-reduced) coherent flop."""
    nmax = _thermal_nmax(nbar) if n_max is None else n_max
    num = wsum = 0.0
    for n in range(nmax + 1):
        pn = thermal_pn(nbar, n)
        if pn <= 0.0:
            continue
        om = rabi0_hz * carrier_rabi_factor(eta, n)
        eff = math.hypot(om, detuning_hz)
        amp = (om * om) / (eff * eff) if eff > 0.0 else 0.0
        # generalised-Rabi flip of this component: amp*sin^2(pi*eff*t) (=0 at t=0)
        num += pn * 0.5 * amp * (1.0 - math.cos(2.0 * math.pi * eff * t_s))
        wsum += pn
    return num / wsum if wsum > 0.0 else 0.0


def sideband_rabi_factor(eta: float, n: int, order: str = "red") -> float:
    """First-order (Lamb-Dicke) sideband Rabi reduction Omega_{n,n-+1}/Omega_0: red
    |n>->|n-1> = eta*sqrt(n) (ZERO at n=0 -- the ground state has no phonon to remove),
    blue |n>->|n+1> = eta*sqrt(n+1). The red sideband's vanishing at n=0 is the
    thermometer: a cold ion barely flops on the red."""
    return eta * math.sqrt(n if order == "red" else n + 1)


def thermal_sideband_flip(t_s: float, rabi0_hz: float, eta: float, nbar: float,
                          order: str = "red", gamma_common_hz: float = 0.0,
                          n_max: int | None = None) -> float:
    """Thermal RED/BLUE sideband-flop P(flipped) at duration t_s: sum_n P_n *
    sin^2(pi Omega_{n,n-+1} t), Omega_{n,n-+1} = sideband_rabi_factor * Omega_0. On the
    RED sideband the n=0 population CANNOT flop (no phonon to subtract), so the flop
    saturates at (1 - P_0)/2 -- a direct motional thermometer. `gamma_common_hz` damps
    the oscillation (a dephasing SHARED with the carrier: scattering + Raman-beam
    phase noise); the sideband flop's own decay therefore BOUNDS that common channel."""
    nmax = _thermal_nmax(nbar) if n_max is None else n_max
    num = wsum = 0.0
    for n in range(nmax + 1):
        pn = thermal_pn(nbar, n)
        wsum += pn
        f = rabi0_hz * sideband_rabi_factor(eta, n, order)
        if f <= 0.0:
            continue                                   # red n=0: no flop (contributes 0)
        env = math.exp(-gamma_common_hz * t_s) if gamma_common_hz > 0.0 else 1.0
        num += pn * 0.5 * (1.0 - env * math.cos(2.0 * math.pi * f * t_s))
    return num / wsum if wsum > 0.0 else 0.0


def thermal_coherence(t_s: float, rabi0_hz: float, eta: float, nbar: float,
                      n_max: int | None = None) -> float:
    """Magnitude of the carrier dephasing envelope, |sum_n P_n e^{i 2pi Omega_n t}|
    (= 1 at t=0, decaying as the n-components dephase). This is the factor that
    multiplies the coherent oscillation; -ln(|.|)/t is the effective dephasing rate
    over the window [0, t]. O(n_max), so cheap enough for an inversion sweep (unlike
    fitting the full flop curve)."""
    nmax = _thermal_nmax(nbar) if n_max is None else n_max
    re = im = wsum = 0.0
    for n in range(nmax + 1):
        pn = thermal_pn(nbar, n)
        if pn <= 0.0:
            continue
        ph = 2.0 * math.pi * rabi0_hz * carrier_rabi_factor(eta, n) * t_s
        re += pn * math.cos(ph)
        im += pn * math.sin(ph)
        wsum += pn
    return math.hypot(re, im) / wsum if wsum > 0.0 else 0.0


def thermal_dephasing_rate(rabi0_hz: float, eta: float, nbar: float) -> float:
    """LEADING-ORDER Gaussian spread of the carrier Rabi frequency from the thermal
    phonon distribution [1/s]: sigma_Omega = Omega_0 eta^2 sqrt(nbar(nbar+1))
    (small-eta: Omega_n ~ Omega_0(1 - eta^2(n+1/2)), Var(n) = nbar(nbar+1)). This is
    the characteristic dephasing rate; it OVERSTATES the loss for nbar < 1 (the
    geometric distribution is far from Gaussian there) -- thermal_carrier_flip is the
    faithful curve. Omega_0 = 2 pi rabi0_hz, so the result is a true rate."""
    return 2.0 * math.pi * rabi0_hz * eta * eta * math.sqrt(nbar * (nbar + 1.0))


def raman_differential_stark_factor(omega_hf_hz: float, raman_detuning_hz: float) -> float:
    """Leading-order Raman differential AC-Stark shift in units of the carrier Rabi
    Omega_0: omega_HF / Delta_R. This is an ORDER-OF-MAGNITUDE (upper) estimate good
    to ~a factor of 2 -- a plain balanced two-state model gives omega_HF/(2 Delta_R);
    the prefactor and CG / P_1/2 / beam-balance corrections are unanchored (no
    measured Raman differential shift exists)."""
    return omega_hf_hz / raman_detuning_hz
