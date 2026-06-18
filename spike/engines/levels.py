"""
Levels engine — hyperfine + Zeeman structure of an alkali-like ``2S_1/2`` ground
state (J = 1/2) with nuclear spin I and magnetic-dipole hyperfine constant A, in
a static magnetic field B.

For J = 1/2 the Breit-Rabi problem is exactly solvable in closed form (each
m_F = m_J + m_I sector is at most a 2x2 block), so no linear-algebra dependency
is needed. The implementation is sign-safe: A may be negative (the 25Mg+ ground
state is inverted, F = I+1/2 lying *below* F = I-1/2), because A enters the
transition frequencies only through A^2 and the diagonal m_J m_I products.

Hamiltonian (in frequency units, H/h):

    H/h = A (I . J) + (g_J mu_B/h) B J_z + (g_I mu_N/h) B I_z

with g_I the *signed* nuclear g-factor. I.J = I_z J_z + (1/2)(I_+ J_- + I_- J_+).

Convention note: the nuclear Zeeman term is written ``+ g_I mu_N B I_z`` with
g_I carrying its sign (g_I < 0 for 25Mg). The m_F=0 clock transition is only
*nearly* independent of g_I — it enters via the combination
(g_J mu_B - g_I mu_N) at the ~12 Hz level at 5.5 G, which the engine retains
exactly; g_I dominates the m_F != 0 Zeeman slopes.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .. import constants as C


@dataclass(frozen=True)
class Level:
    F: float
    mF: float
    energy_hz: float  # E/h relative to the unperturbed fine-structure centroid


class GroundStateZeeman:
    def __init__(self, A_hz: float, I: float, g_J: float = C.G_J_2S12, g_I: float = C.G_I_25MG):
        A_hz, I = float(A_hz), float(I)
        if not math.isfinite(A_hz):
            raise ValueError(f"hyperfine constant A must be finite, got {A_hz}")
        if I < 0.5:
            raise ValueError(
                f"levels engine requires nuclear spin I >= 1/2; got I={I}. "
                "I=0 isotopes (24Mg, 26Mg) have no hyperfine manifold — they are "
                "coolant/co-trapped, not qubits."
            )
        self.A = A_hz                 # magnetic-dipole hyperfine constant [Hz], signed
        self.I = I                    # nuclear spin
        self.J = 0.5
        self.g_J = float(g_J)
        self.g_I = float(g_I)

    # ------------------------------------------------------------------ #
    @classmethod
    def from_ledger(cls, ledger, a_name="hyperfine_a_constant_25mg",
                    i_name="nuclear_spin_25mg", gj_name="g_factor_electron_2s12",
                    gi_name="g_factor_nuclear_25mg", g_J=None, g_I=None):
        """Build the engine from `input` ledger records (the wall: it must not
        read benchmarks). Every consumed quantity — A, I, and the g-factors —
        goes through ledger.input_quantity(), which enforces kind:input.
        g_J/g_I are read from the ledger when present, else fall back to the
        CODATA/tabulated values in constants.py (backward compatible)."""
        A = ledger.input_quantity(a_name).value
        I = ledger.input_quantity(i_name).value
        if g_J is None:
            g_J = ledger.input_quantity(gj_name).value if gj_name in ledger else C.G_J_2S12
        if g_I is None:
            g_I = ledger.input_quantity(gi_name).value if gi_name in ledger else C.G_I_25MG
        return cls(A_hz=A, I=I, g_J=g_J, g_I=g_I)

    # ------------------------------------------------------------------ #
    def zero_field_splitting(self) -> float:
        """|3A| = |A|(I+1/2): the zero-field F=I+1/2 <-> F=I-1/2 splitting [Hz]."""
        return abs(self.A) * (self.I + 0.5)

    def levels(self, B: float) -> dict:
        """All ground-state sublevels at field B [Tesla]. Returns {(F, mF): E/h [Hz]}."""
        A, I = self.A, self.I
        gjB = self.g_J * C.MU_B_OVER_H * B          # Hz
        giB = self.g_I * C.MU_N_OVER_H * B          # Hz  (g_I signed)
        Fhi, Flo = I + 0.5, I - 0.5
        mF_max = I + 0.5
        out: dict = {}
        mF = -mF_max
        while mF <= mF_max + 1e-9:
            if abs(abs(mF) - mF_max) < 1e-9:
                # stretched state: a single |m_J=±1/2, m_I=±I> eigenstate
                s = 1.0 if mF > 0 else -1.0
                mJ, mI = s * 0.5, s * I
                out[(Fhi, mF)] = A * mI * mJ + gjB * mJ + giB * mI
            else:
                # 2x2 block in {|+1/2, mF-1/2>, |-1/2, mF+1/2>}
                mI1, mI2 = mF - 0.5, mF + 0.5
                a = A * mI1 * 0.5 + gjB * 0.5 + giB * mI1
                d = A * mI2 * (-0.5) + gjB * (-0.5) + giB * mI2
                b = 0.5 * A * math.sqrt(I * (I + 1) - (mF * mF - 0.25))
                mean = 0.5 * (a + d)
                half = math.hypot(0.5 * (a - d), b)
                e_minus, e_plus = mean - half, mean + half
                # For A < 0 (inverted) the F=I+1/2 manifold is the lower one.
                if A < 0:
                    out[(Fhi, mF)], out[(Flo, mF)] = e_minus, e_plus
                else:
                    out[(Flo, mF)], out[(Fhi, mF)] = e_minus, e_plus
            mF += 1.0
        return out

    def energy(self, F: float, mF: float, B: float) -> float:
        return self.levels(B)[(F, mF)]

    def transition(self, upper, lower, B: float) -> float:
        """Frequency [Hz] of |upper> <- |lower>, each given as (F, mF)."""
        lv = self.levels(B)
        return lv[tuple(upper)] - lv[tuple(lower)]

    def clock_transition(self, B: float) -> float:
        """The field-insensitive (first order) |I+1/2, 0> <-> |I-1/2, 0> clock
        transition [Hz], reported as a positive frequency. Requires half-integer
        I, so that F is integer and an m_F=0 sublevel exists."""
        if abs((self.I + 0.5) - round(self.I + 0.5)) > 1e-9:
            raise ValueError(
                f"clock |I+1/2,0> <-> |I-1/2,0> requires half-integer I (integer F, "
                f"so m_F=0 exists); got I={self.I}"
            )
        lv = self.levels(B)
        return abs(lv[(self.I + 0.5, 0.0)] - lv[(self.I - 0.5, 0.0)])

    def _require_half_integer_I(self):
        if abs((self.I + 0.5) - round(self.I + 0.5)) > 1e-9:
            raise ValueError(f"integer F (m_F=0 manifold) requires half-integer I; got I={self.I}")

    def hyperfine_transitions(self, B: float) -> dict:
        """All |F=I+1/2, m_F> <-> |F=I-1/2, m_F'> microwave hyperfine transitions
        [Hz], keyed by (m_F_upper, m_F_lower), with |delta m_F| <= 1 (magnetic
        dipole). The field-insensitive clock is the (0, 0) entry; the others
        carry first-order Zeeman shifts."""
        self._require_half_integer_I()
        lv = self.levels(B)
        hi = {mF: E for (F, mF), E in lv.items() if F == self.I + 0.5}
        lo = {mF: E for (F, mF), E in lv.items() if F == self.I - 0.5}
        return {
            (mFu, mFl): abs(Eu - El)
            for mFu, Eu in hi.items()
            for mFl, El in lo.items()
            if abs(mFu - mFl) <= 1.0
        }

    def clock_transition_from_spectrum(self, B: float) -> float:
        """The clock as the (0,0) entry of the hyperfine spectrum (cross-check of
        clock_transition())."""
        return self.hyperfine_transitions(B)[(0.0, 0.0)]

    def zeeman_splitting(self, F: float, mF_low: float, B: float) -> float:
        """Adjacent Zeeman splitting E(F, mF_low+1) - E(F, mF_low) [Hz] within a
        hyperfine manifold (~ g_F mu_B B at low field; sign follows g_F)."""
        lv = self.levels(B)
        return lv[(F, mF_low + 1.0)] - lv[(F, mF_low)]

    def quadratic_zeeman_coeff(self) -> float:
        """Leading quadratic Zeeman coefficient K [Hz/T^2] of the clock
        transition: nu(B) ~ nu0 + K B^2, K = (g_J mu_B/h - g_I mu_N/h)^2 / (2 nu0)."""
        nu0 = self.zero_field_splitting()
        slope = self.g_J * C.MU_B_OVER_H - self.g_I * C.MU_N_OVER_H   # Hz/T
        return slope * slope / (2.0 * nu0)
