"""
Polarization- and power-resolved optical (Raman-beam) light shifts and off-resonant
scattering on the 25Mg+ 3s 2S_1/2 -> 3p 2P_{1/2,3/2} D lines.

WHY this engine. The far-detuned Raman beams (Delta_R ~ 20 GHz red of 3P_3/2) both
AC-Stark-shift the qubit |down>=|3,3>, |up>=|2,2> AND scatter photons; BOTH depend
on each beam's POLARIZATION (which ground sublevels it couples, and to which fine-
structure line) and POWER (linearly). The scalar engine (engines/scatter.py) and the
leading-order differential shift (engines/sideband.raman_differential_stark_factor,
~omega_HF/Delta_R) ignore this. Here we do the full angular-momentum sum.

BASIS (the physically honest choice). At Delta_R ~ 20 GHz the excited 3P hyperfine +
Zeeman structure is UNRESOLVED, so F', mF' are NOT good quantum numbers for the P
levels. We therefore work in the EXCITED FINE-STRUCTURE basis |P_{J'} mJ'> (no
hyperfine), and decompose the GROUND |S_{1/2} F mF> -- where F IS good (hyperfine
1.79 GHz >> Zeeman at 5.5 G) -- into electronic-spin components |mJ, mI> with the
optical dipole acting only on mJ (the nuclear spin mI is a SPECTATOR, conserved):

    |F mF> = sum_{mJ} <S_{1/2} mJ, I (mF-mJ) | F mF> |mJ, mI=mF-mJ>.

The electronic dipole is <P_{J'} mJ'| d_q | S_{1/2} mJ> = CG(1/2 mJ; 1 q | J' mJ')*d0
with a COMMON reduced d0 -- the line strengths S(D2):S(D1) = 2:1 come from the
(2J'+1) multiplicity, not d0. Single-beam shift and scattering of |F mF>:

    delta_AC = sum_mJ |<mJ,mI|F mF>|^2 * (1/4) sum_{q,J'} (E xi_q)^2 |d_{mJ,q,J'}|^2/Delta(J',F)
    Gamma_sc = sum_mJ |<mJ,mI|F mF>|^2 *       sum_{q,J'} (E xi_q)^2 |d|^2 Gamma_{J'}/(4 Delta^2)

with Delta(J',F) = Delta_R + (J'=1/2 ? Delta_FS : 0) + ground-hyperfine(F) (|down> vs
|up> differ by omega_HF -> the SCALAR differential). The excited hyperfine is
degenerate (the honest far-detuned limit; an optional p_hyperfine hook is kept off).

The F'-basis sum (coupling/line_coupling_sq, via a Wigner 6j) is retained ONLY as a
degenerate-limit CROSS-CHECK: it equals the |mJ> sum to machine precision (basis
independence), which the tests assert. The transparent summary is the scalar+vector
polarizability (scalar_vector_shift): tensor = 0 for 2S_1/2, the vector part is the
mJ-ODD (circular-polarization) shift. Since Delta_R << Delta_FS the laser sits
essentially on P_3/2 alone, so the vector shift is NOT fine-structure-suppressed --
~0.5x the scalar for circular light: polarization has a LARGE effect on the
differential shift and the spin-flip scattering.

ANCHORING. The absolute field scale E is unknown, but it CANCELS in the ratios
delta_AC/Omega_2gamma and Gamma_sc/Omega_2gamma, where Omega_2gamma is the two-photon
carrier Rabi computed in the SAME relative dipole units (two_photon_rabi_relative).
So a twin multiplies these dimensionless ratios by the MEASURED flop Rabi -- the
engine refines the polarization/power factors while staying pinned to the data.

PUBLIC API (the |mJ> basis is the PRODUCTION path; the F'-basis is a cross-check):
  * RamanOptics.light_shift / scatter_rate / two_photon_rabi / scalar_vector_shift  -- production
  * RamanOptics.differential_stark_per_rabi -> DIMENSIONLESS (Hz/Hz): multiply by the
    measured rabi_hz to get delta_AC [Hz].
  * RamanOptics.scatter_per_rabi           -> DIMENSIONLESS (Gamma/Delta-like):
    multiply by 2*pi*rabi_hz to get Gamma_sc [1/s] (the 2pi turns the /2pi Rabi into a
    rate, matching engines/scatter.py's convention).
  * coupling / line_coupling_sq (F'-basis, via wigner_6j) -- CROSS-CHECK ONLY; equals
    the |mJ> sum in the degenerate limit (tested), not used by the production path.

Pure Python; reuses drive.clebsch_gordan, adds a Racah wigner_6j (cross-checked
against sympy in the tests). Frequencies in Hz.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .drive import _f, clebsch_gordan

# Reduced dipole matrix elements d_J' in relative units (line strengths 1 : 2).
_RED = {0.5: 1.0, 1.5: math.sqrt(2.0)}
# Excited fine-structure manifolds: J' -> the F' it spans for I = 5/2.
_J_PRIME = (0.5, 1.5)


def wigner_6j(j1, j2, j3, j4, j5, j6) -> float:
    """Wigner 6j symbol {j1 j2 j3; j4 j5 j6} via the Racah formula. Returns 0 if any
    triad fails the triangle rule (half-integer arguments allowed)."""
    def tri(a, b, c):
        # triangle inequality AND integer perimeter (a 6j triad must sum to an integer;
        # a half-integer perimeter, e.g. (4,3,2.5), is forbidden -> 6j = 0)
        return (a + b - c) >= 0 and (a - b + c) >= 0 and (-a + b + c) >= 0 and \
               abs(round(a + b + c) - (a + b + c)) < 1e-9

    triads = ((j1, j2, j3), (j1, j5, j6), (j4, j2, j6), (j4, j5, j3))
    if not all(tri(*t) for t in triads):
        return 0.0

    def delta(a, b, c):
        return math.sqrt(_f(a + b - c) * _f(a - b + c) * _f(-a + b + c) / _f(a + b + c + 1))

    pre = delta(j1, j2, j3) * delta(j1, j5, j6) * delta(j4, j2, j6) * delta(j4, j5, j3)
    # sum over t; the four "t - sum" args and three "sum - t" args must be >= 0
    s1, s2, s3, s4 = j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3
    p1, p2, p3 = j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j1 + j3 + j4 + j6
    tmin = int(round(max(s1, s2, s3, s4)))
    tmax = int(round(min(p1, p2, p3)))
    total = 0.0
    for t in range(tmin, tmax + 1):
        total += (-1) ** t * _f(t + 1) / (
            _f(t - s1) * _f(t - s2) * _f(t - s3) * _f(t - s4)
            * _f(p1 - t) * _f(p2 - t) * _f(p3 - t)
        )
    return pre * total


def f_levels(J_prime: float, I: float = 2.5):
    """The F' = |J'-I| .. J'+I spanned by an excited fine-structure manifold."""
    lo, hi = abs(J_prime - I), J_prime + I
    n = int(round(hi - lo))
    return [lo + k for k in range(n + 1)]


def coupling(F, mF, q, J_prime, F_prime, mF_prime, J=0.5, I=2.5) -> float:
    """Relative dipole amplitude C(F mF -> J' F' mF', q), the PHYSICAL per-state
    matrix element |<e|d_q|g>| (up to the overall reduced scale). By Wigner-Eckart
    (3j convention) and the hyperfine reduction,

        |C|^2 = CG(F mF; 1 q | F' mF')^2 * (2F+1) * {J' F' I; F J 1}^2 * d_J'^2,

    so C = CG * sqrt(2F+1) * {6j} * d_J' (the reduced element's (2F'+1) cancels
    against the 3j<->CG conversion -- WITHOUT this the total |C|^2 would spuriously
    depend on F, violating the F-independence of <g|d^2|g>). reduced ratio
    d_{1/2}:d_{3/2} = 1:sqrt2 (line strengths 1:2). Zero unless mF' = mF + q."""
    if abs((mF + q) - mF_prime) > 1e-9:
        return 0.0
    cg = clebsch_gordan(F, mF, 1, q, F_prime, mF_prime)
    if cg == 0.0:
        return 0.0
    phase = (-1) ** int(round(J_prime + I + F + 1))
    rme_f = math.sqrt(2 * F + 1) * wigner_6j(J_prime, F_prime, I, F, J, 1)
    return cg * phase * rme_f * _RED[J_prime]


def line_coupling_sq(F, mF, q, J_prime, J=0.5, I=2.5) -> float:
    """sum over F' of |C(F mF -> J' F' mF=mF+q, q)|^2 for one fine-structure line --
    the total |coupling|^2 from |F,mF> via polarization q into manifold J'."""
    return sum(coupling(F, mF, q, J_prime, Fp, mF + q, J, I) ** 2 for Fp in f_levels(J_prime, I))


# --- fine-structure (|J', mJ'>) basis: the physically honest far-detuned picture --
# At Delta_R ~ 20 GHz the excited 3P hyperfine + Zeeman structure is UNRESOLVED, so
# F', mF' are NOT good quantum numbers for the P levels. We therefore work in the
# fine-structure |P_{J'} mJ'> basis (no excited hyperfine), and decompose the GROUND
# |S_{1/2} F mF> -- where F IS good (hyperfine 1.79 GHz >> Zeeman at 5.5 G) -- into
# |mJ, mI> with the optical dipole acting ONLY on mJ (the nuclear spin mI is a
# SPECTATOR, conserved). The F'-basis sum above is kept only as a degenerate-limit
# cross-check; the two bases agree (basis independence), which the tests assert.

def hf_amplitude(F, mF, mJ, I=2.5) -> float:
    """<S_{1/2} mJ, I mI | F mF> with mI = mF - mJ: the weight of electronic-spin
    projection mJ in the ground hyperfine state |F, mF>."""
    return clebsch_gordan(0.5, mJ, I, mF - mJ, F, mF)


# Common electronic reduced dipole <S_1/2 || d || P_J'>: J'-INDEPENDENT (same radial
# integral). The D2:D1 = 2:1 line strength then emerges from the (2J'+1) multiplicity
# of the CG sum, NOT from the reduced element -- the physically correct statement in
# this basis (and where the F'-basis instead carries an explicit sqrt2 in _RED). The
# 1/sqrt2 fixes the overall scale to match the F'-basis cross-check exactly.
_D0 = 1.0 / math.sqrt(2.0)


def electronic_dipole(mJ, q, J_prime, mJ_prime) -> float:
    """Relative electronic dipole <P_{J'} mJ'| d_q | S_{1/2} mJ> = CG(1/2 mJ; 1 q | J'
    mJ') * d0, with d0 COMMON to both lines (the 2:1 line strength comes from the
    (2J'+1) CG multiplicity). Nuclear spin is a spectator. Zero unless mJ' = mJ + q."""
    if abs((mJ + q) - mJ_prime) > 1e-9 or abs(mJ_prime) > J_prime + 1e-9:
        return 0.0
    return clebsch_gordan(0.5, mJ, 1, q, J_prime, mJ_prime) * _D0


_QMAP = {1: 0, 0: 1, -1: 2}      # spherical q (sigma+, pi, sigma-) -> pol-tuple index
_MJ = (0.5, -0.5)


@dataclass
class RamanBeam:
    """One Raman beam: a relative INTENSITY (e.g. the .dat per-beam power, arbitrary
    units -- only ratios matter) and a spherical-polarization decomposition pol =
    (sigma+, pi, sigma-) as INTENSITY fractions (summing to 1) along the quantization
    axis B. The field amplitude in component q is sqrt(intensity * pol_q)."""
    intensity: float
    pol: tuple                                   # (f_sigma+, f_pi, f_sigma-), sum = 1
    label: str = ""

    def amp_sq(self, q: int) -> float:
        """|E_q|^2 = intensity * polarization fraction in component q."""
        return self.intensity * self.pol[_QMAP[q]]


# canonical polarizations (intensity fractions along B); see ledger records.
POL_PI = (0.0, 1.0, 0.0)                          # pure pi (B1)
POL_LINEAR_PERP = (0.5, 0.0, 0.5)                 # sigma+ + sigma- or sigma+ - sigma- (R1, R2):
#   both are LINEAR (perp to B) -> equal |sigma+|^2, |sigma-|^2, no pi. They differ only by a
#   relative PHASE, which does not enter single-beam |C|^2 sums (intensity-level model).
POL_SIGMA_PLUS = (1.0, 0.0, 0.0)
POL_SIGMA_MINUS = (0.0, 0.0, 1.0)

# combination -> (blue beam, red beam) (Doerr 2024 Sec. 2.1.4)
_COMBO_BEAMS = {"CC": ("b1", "r1"), "OC": ("b1", "r2"),
                "AC": ("b3", "r1"), "ROC": ("b3", "r2")}


def beams_from_ledger(ledger, combo: str = "OC", powers=None):
    """(blue, red) RamanBeam for a combination: polarizations from the ledger
    (raman_<beam>_polarization_25mg, input), relative INTENSITIES from `powers`
    (a dict beam-label -> relative power, e.g. the .dat pwr_b1/pwr_r2; default 1)."""
    powers = powers or {}
    bl, rd = _COMBO_BEAMS[combo]

    def mk(name):
        pol = tuple(ledger.input_quantity(f"raman_{name}_polarization_25mg").value)
        return RamanBeam(float(powers.get(name, 1.0)), pol, name)

    return mk(bl), mk(rd)


class RamanOptics:
    """Polarization+power-resolved optical light shifts and scattering for the Raman
    beams on the 25Mg+ D lines, in the fine-structure |J', mJ'> excited basis.
    Detunings are SIGNED (red-detuned < 0): the Raman laser sits ~20 GHz RED of
    3P_3/2 and ~2.7 THz BLUE of 3P_1/2 (opposite sign, ~Delta_R/Delta_FS smaller) --
    that interplay IS the fine-structure dependence of the vector light shift /
    spin-flip scattering. The excited hyperfine is unresolved (degenerate, the honest
    far-detuned limit); an optional p_hyperfine hook is kept for a future F'-resolved
    refinement but defaults off."""

    def __init__(self, raman_detuning_hz: float, fine_structure_hz: float,
                 omega_hf_hz: float, gamma_p32_hz: float, gamma_p12_hz: float,
                 p_hyperfine=None, I: float = 2.5):
        self.delta_p32 = -abs(raman_detuning_hz)     # SIGNED (red, negative)
        self.delta_fs = abs(fine_structure_hz)
        self.omega_hf = abs(omega_hf_hz)
        self.gamma = {0.5: abs(gamma_p12_hz), 1.5: abs(gamma_p32_hz)}
        self.I = I
        # ground hyperfine energies (F=3 lower for inverted 25Mg+): E(F=2)-E(F=3)=+omega_hf
        self.ground_energy = {I + 0.5: 0.0, I - 0.5: self.omega_hf}
        self.p_hyperfine = p_hyperfine or {}         # optional excited-hf offsets (default off)

    @classmethod
    def from_ledger(cls, ledger, p_hyperfine=None):
        """Consume the Raman detuning, the (derived) fine-structure splitting, the
        qubit hyperfine splitting, and the two P linewidths -- all `input`."""
        return cls(
            raman_detuning_hz=ledger.input_quantity("raman_detuning_from_p32").value,
            fine_structure_hz=ledger.input_quantity("mg_fine_structure_splitting_3p_25mg").value,
            omega_hf_hz=ledger.input_quantity("hyperfine_splitting_25mg_f2_f3").value,
            gamma_p32_hz=ledger.input_quantity("mg_p32_natural_linewidth").value,
            gamma_p12_hz=ledger.input_quantity("mg_p12_natural_linewidth").value,
            p_hyperfine=p_hyperfine,
        )

    def detuning(self, J_prime: float, F: float) -> float:
        """Signed laser detuning [Hz] from |S_{1/2} F> -> |P_{J'}>: delta_p32 + (ground
        hyperfine of F) + (Delta_FS if P_1/2). Red of P3/2 -> negative; far blue of
        P1/2 -> large positive. (Excited hyperfine unresolved -> no F' dependence.)"""
        d = self.delta_p32 + self.ground_energy.get(F, 0.0)
        if J_prime == 0.5:
            d += self.delta_fs
        return d

    def _electronic_shift(self, beam, mJ, F) -> float:
        """AC-Stark shift of the electronic state |S_{1/2} mJ> in ground manifold F:
        (1/4) sum_{q,J'} |E_q|^2 |<P_{J'} mJ+q|d_q|mJ>|^2 / Delta(J',F). Signed."""
        s = 0.0
        for q in (1, 0, -1):
            eq2 = beam.amp_sq(q)
            if eq2 == 0.0:
                continue
            for Jp in _J_PRIME:
                d = electronic_dipole(mJ, q, Jp, mJ + q)
                if d == 0.0:
                    continue
                s += eq2 * d * d / (4.0 * self.detuning(Jp, F))
        return s

    def light_shift(self, beam: RamanBeam, F: float, mF: float) -> float:
        """Single-beam AC-Stark shift of |F, mF> = sum_{mJ} |<mJ,mI|F mF>|^2 *
        electronic_shift(mJ). Signed (red beam -> negative)."""
        return sum(hf_amplitude(F, mF, mJ, self.I) ** 2 * self._electronic_shift(beam, mJ, F)
                   for mJ in _MJ)

    def _electronic_scatter(self, beam, mJ, F) -> float:
        s = 0.0
        for q in (1, 0, -1):
            eq2 = beam.amp_sq(q)
            if eq2 == 0.0:
                continue
            for Jp in _J_PRIME:
                d = electronic_dipole(mJ, q, Jp, mJ + q)
                if d == 0.0:
                    continue
                dl = self.detuning(Jp, F)
                s += eq2 * d * d * self.gamma[Jp] / (4.0 * dl * dl)
        return s

    def scatter_rate(self, beam: RamanBeam, F: float, mF: float) -> float:
        """Single-beam total photon-scattering rate from |F, mF> [intensity units]:
        sum_{mJ} |<mJ,mI|F mF>|^2 sum_{q,J'} |E_q|^2 |d|^2 Gamma_{J'}/(4 Delta^2). Same
        relative dipole units as light_shift / two_photon_rabi (absolute field scale
        cancels in the twin's ratios)."""
        return sum(hf_amplitude(F, mF, mJ, self.I) ** 2 * self._electronic_scatter(beam, mJ, F)
                   for mJ in _MJ)

    def two_photon_rabi(self, beam_b: RamanBeam, beam_r: RamanBeam,
                        lower=(3.0, 3.0), upper=(2.0, 2.0)) -> float:
        """Relative two-photon (carrier) Rabi |Omega_2gamma| for lower<->upper, summed
        over intermediate |P_{J'} mJ'>: beam_b absorbs (q_b) from the |mJ_lo,mI> part of
        `lower`, beam_r stimulates (q_r) to the |mJ_hi,mI> part of `upper` with the SAME
        nuclear spin mI (the optical dipole conserves mI). Same dipole units as the
        shifts -- so delta_AC/Omega_2gamma and Gamma_sc/Omega_2gamma are independent of
        the absolute field scale (only the intensity RATIO and polarizations enter)."""
        Fl, mFl = lower
        Fu, mFu = upper
        amp = 0.0
        for mJl in _MJ:
            al = hf_amplitude(Fl, mFl, mJl, self.I)
            if al == 0.0:
                continue
            for mJu in _MJ:
                au = hf_amplitude(Fu, mFu, mJu, self.I)
                if au == 0.0 or abs((mFl - mJl) - (mFu - mJu)) > 1e-9:   # mI conserved
                    continue
                for q_b in (1, 0, -1):
                    eb = math.sqrt(beam_b.amp_sq(q_b))
                    if eb == 0.0:
                        continue
                    mJp = mJl + q_b
                    q_r = int(round(mJp - mJu))
                    if q_r not in (-1, 0, 1):
                        continue
                    er = math.sqrt(beam_r.amp_sq(q_r))
                    if er == 0.0:
                        continue
                    for Jp in _J_PRIME:
                        dl = electronic_dipole(mJl, q_b, Jp, mJp)
                        du = electronic_dipole(mJu, q_r, Jp, mJp)
                        if dl == 0.0 or du == 0.0:
                            continue
                        amp += al * au * eb * dl * er * du / (2.0 * self.detuning(Jp, Fl))
        return abs(amp)

    # --- scalar / vector light-shift summary (basis-independent, versatile) ----
    def scalar_vector_shift(self, beam: RamanBeam, F: float):
        """(scalar, vector) decomposition of the electronic light shift in manifold F:
        scalar = [shift(mJ=+1/2) + shift(-1/2)]/2  (mF-independent, set by total
        intensity), vector = [shift(+1/2) - shift(-1/2)]/2  (the mJ-ODD part, present
        ONLY for circular polarization and proportional to (1/Delta_{P3/2} -
        1/Delta_{P1/2}) -- the fine-structure-sensitive piece). For 2S_1/2 the tensor
        polarizability vanishes, so these two terms are complete."""
        up = self._electronic_shift(beam, 0.5, F)
        dn = self._electronic_shift(beam, -0.5, F)
        return 0.5 * (up + dn), 0.5 * (up - dn)

    # --- ABSOLUTE anchoring: relative engine units -> Hz via the saturation -----
    def cycling_scale_hz(self) -> float:
        """Engine(relative)->Hz scale, anchored so that a sigma+ beam of saturation
        s=1 on the cycling |3,3>->|P3/2 4,4> gives the textbook scalar light shift
        s*Gamma^2/(8|Delta|). Multiply any light_shift (with beam.intensity set to the
        saturation parameter s) by this to get the absolute shift in Hz."""
        ref = self.light_shift(RamanBeam(1.0, POL_SIGMA_PLUS), 3.0, 3.0)
        phys = self.gamma[1.5] ** 2 * 0.5 / (4.0 * self.delta_p32)   # SIGNED: red -> negative
        return phys / ref if ref else float("nan")                  # +ve scale (preserves sign)

    def differential_stark_hz(self, beam_at_saturation: RamanBeam,
                              up=(2.0, 2.0), dn=(3.0, 3.0)) -> float:
        """ABSOLUTE single-beam differential AC-Stark shift delta_AC = LS(up)-LS(dn)
        [Hz] of the qubit (|up>=|2,2>, |dn>=|3,3>). The beam's .intensity must be its
        saturation parameter s (= I/I_sat = 2P/(pi w^2 I_sat)). Signed (red beam)."""
        sc = self.cycling_scale_hz()
        return sc * (self.light_shift(beam_at_saturation, *up)
                     - self.light_shift(beam_at_saturation, *dn))

    @staticmethod
    def saturation(power_w: float, waist_m: float, i_sat_w_m2: float) -> float:
        """On-axis saturation parameter s = I/I_sat for a Gaussian beam, I = 2P/(pi w^2)."""
        return (2.0 * power_w / (math.pi * waist_m * waist_m)) / i_sat_w_m2

    def two_photon_rabi_hz(self, beam_b: RamanBeam, beam_r: RamanBeam,
                           lower=(3.0, 3.0), upper=(2.0, 2.0)) -> float:
        """ABSOLUTE two-photon (carrier) Rabi |Omega_2gamma|/2pi [Hz] for beams whose
        .intensity are saturation parameters s -- anchored by the SAME cycling scale as
        the light shifts (both are dipole^2/Delta in the engine's units). Lets the twin
        check the NOMINAL beam intensities against the OBSERVED flop rate (the real
        measure of the intensity at the ion) and re-anchor the AC-Stark shifts to it."""
        return abs(self.two_photon_rabi(beam_b, beam_r, lower, upper)) * self.cycling_scale_hz()

    def circular_fraction(self, beam: RamanBeam) -> float:
        """Net circular polarization C = (|sigma+|^2 - |sigma-|^2)/(total) along B --
        the quantity that controls the vector light shift (0 for pi or balanced
        linear, +-1 for pure sigma+/-)."""
        sp, _, sm = beam.pol
        tot = sp + sm
        return (sp - sm) / tot if tot > 0 else 0.0

    # --- twin-facing dimensionless ratios (E-scale cancels) -----------------
    def differential_stark_per_rabi(self, beams, beam_b, beam_r,
                                   lower=(3.0, 3.0), upper=(2.0, 2.0)) -> float:
        """[delta_AC(upper) - delta_AC(lower)] / Omega_2gamma -- the differential
        carrier light shift in units of the two-photon Rabi. Polarization+power
        resolved; replaces the leading-order omega_HF/Delta_R factor."""
        d_lo = sum(self.light_shift(b, *lower) for b in beams)
        d_hi = sum(self.light_shift(b, *upper) for b in beams)
        rabi = self.two_photon_rabi(beam_b, beam_r, lower, upper)
        return (d_hi - d_lo) / rabi if rabi else float("nan")

    def scatter_per_rabi(self, beams, beam_b, beam_r, state=(3.0, 3.0),
                        lower=(3.0, 3.0), upper=(2.0, 2.0)) -> float:
        """Total scattering rate of `state` (summed over beams) / Omega_2gamma --
        scattering per unit coherent Rabi, polarization+power resolved. Generalises the
        scalar Gamma_sc/Omega = Gamma/Delta_R x balance factor."""
        g = sum(self.scatter_rate(b, *state) for b in beams)
        rabi = self.two_photon_rabi(beam_b, beam_r, lower, upper)
        return g / rabi if rabi else float("nan")
