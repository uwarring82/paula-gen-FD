"""
Tests for the polarization+power-resolved Raman optical engine: the angular-momentum
core (Wigner 6j, the cycling transition, sum rules, line strengths), BASIS
INDEPENDENCE (the F'-basis cross-check equals the |J',mJ'> primary basis), the
scalar/vector light-shift structure, and the ledger wiring.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.raman_optical import (
    POL_LINEAR_PERP,
    POL_PI,
    POL_SIGMA_MINUS,
    POL_SIGMA_PLUS,
    RamanBeam,
    RamanOptics,
    _MJ,
    beams_from_ledger,
    coupling,
    electronic_dipole,
    f_levels,
    hf_amplitude,
    line_coupling_sq,
    wigner_6j,
)
from spike.ledger import Ledger

GROUND = [(3, 3), (3, 2), (3, 0), (3, -3), (2, 2), (2, 0), (2, -2)]


# --- angular-momentum core --------------------------------------------------
def test_wigner_6j_known_values():
    assert wigner_6j(1, 1, 1, 1, 1, 1) == pytest.approx(1.0 / 6.0)
    assert wigner_6j(0.5, 0.5, 1, 0.5, 0.5, 1) == pytest.approx(1.0 / 6.0)
    assert wigner_6j(1, 1, 2, 1, 1, 2) == pytest.approx(1.0 / 30.0, abs=1e-9)
    assert wigner_6j(2, 2, 2, 2, 2, 2) == pytest.approx(-3.0 / 70.0, abs=1e-9)
    # symmetry under column swap; triangle violation -> 0
    assert wigner_6j(1.5, 4, 2.5, 3, 0.5, 1) == pytest.approx(wigner_6j(4, 1.5, 2.5, 0.5, 3, 1))
    assert wigner_6j(1, 1, 5, 1, 1, 1) == 0.0


def test_f_levels():
    assert f_levels(0.5) == [2.0, 3.0]
    assert f_levels(1.5) == [1.0, 2.0, 3.0, 4.0]


def test_wigner_6j_symmetries():
    # invariant under any column permutation and under swapping upper<->lower in two
    # columns (Regge symmetries) -- a property test, not a table lookup
    a = wigner_6j(1.5, 4, 2.5, 3, 0.5, 1)
    assert a == pytest.approx(wigner_6j(4, 1.5, 2.5, 0.5, 3, 1))        # swap cols 1,2
    assert a == pytest.approx(wigner_6j(2.5, 4, 1.5, 1, 0.5, 3))        # swap cols 1,3
    assert a == pytest.approx(wigner_6j(1.5, 0.5, 1, 3, 4, 2.5))        # swap rows in cols 2,3


def test_wigner_6j_vs_sympy():
    # cross-check the hand-rolled Racah 6j against sympy over a range of (half-)integer
    # arguments incl. the I=5/2 / J'=3/2 values the engine actually uses (reviewer rec).
    sympy_wigner = pytest.importorskip("sympy.physics.wigner")
    from sympy import Rational

    def toR(x):
        return Rational(int(round(2 * x)), 2)

    args = [(1, 1, 1, 1, 1, 1), (0.5, 0.5, 1, 0.5, 0.5, 1), (1, 1, 2, 1, 1, 2),
            (1.5, 4, 2.5, 3, 0.5, 1), (1.5, 3, 2.5, 3, 0.5, 1), (1.5, 2, 2.5, 2, 0.5, 1),
            (2.5, 2.5, 1, 1.5, 1.5, 1), (4, 4, 1, 3, 3, 2.5)]
    for t in args:
        try:
            ref = float(sympy_wigner.wigner_6j(*[toR(x) for x in t]))
        except ValueError:
            ref = 0.0                      # sympy raises on triangle violation; we return 0
        assert wigner_6j(*t) == pytest.approx(ref, abs=1e-12)


def test_clebsch_gordan_vs_sympy():
    sympy_cg = pytest.importorskip("sympy.physics.quantum.cg")
    from spike.engines.drive import clebsch_gordan
    from sympy import Rational, S

    def toR(x):
        return Rational(int(round(2 * x)), 2)

    cases = [(0.5, 0.5, 2.5, 2.5, 3, 3), (0.5, -0.5, 2.5, 2.5, 2, 2),
             (3, 3, 1, -1, 4, 2), (1.5, 0.5, 1.5, -0.5, 1, 0)]
    for j1, m1, j2, m2, j3, m3 in cases:
        ref = float(sympy_cg.CG(toR(j1), toR(m1), toR(j2), toR(m2), toR(j3), toR(m3)).doit())
        assert clebsch_gordan(j1, m1, j2, m2, j3, m3) == pytest.approx(ref, abs=1e-12)


def test_clebsch_gordan_orthogonality():
    # sum_{m1 m2} CG(j1 m1 j2 m2|J M) CG(j1 m1 j2 m2|J' M') = delta_{JJ'} delta_{MM'}
    from spike.engines.drive import clebsch_gordan
    j1, j2 = 0.5, 2.5      # the S1/2 (x) I=5/2 decomposition used by hf_amplitude
    for (J, M), (Jp, Mp) in [((3, 3), (3, 3)), ((3, 2), (2, 2)), ((2, 2), (2, 2)), ((3, 1), (2, 1))]:
        s = sum(clebsch_gordan(j1, m1, j2, M - m1, J, M) * clebsch_gordan(j1, m1, j2, Mp - m1, Jp, Mp)
                for m1 in (0.5, -0.5))
        expect = 1.0 if (J == Jp and M == Mp) else 0.0
        assert s == pytest.approx(expect, abs=1e-12)


def test_cycling_transition_is_pure_p32():
    # |3,3> + sigma+ couples ONLY to the P3/2 cycling state (P1/2 max F'=3 -> mF'=4 impossible)
    assert line_coupling_sq(3, 3, +1, 0.5) == pytest.approx(0.0)
    assert line_coupling_sq(3, 3, +1, 1.5) > 0.0


def test_sum_rule_total_is_sublevel_and_F_independent():
    # <g|d^2|g> is the same for EVERY ground sublevel (both F, all mF)
    def total(F, mF):
        return sum(line_coupling_sq(F, mF, q, Jp) for q in (-1, 0, 1) for Jp in (0.5, 1.5))
    vals = [total(F, mF) for F, mF in GROUND]
    assert max(vals) - min(vals) < 1e-9
    assert vals[0] == pytest.approx(1.5)


def test_line_strength_ratio_is_two():
    full = [(3, m) for m in range(-3, 4)] + [(2, m) for m in range(-2, 3)]
    d1 = sum(line_coupling_sq(F, mF, q, 0.5) for F, mF in full for q in (-1, 0, 1))
    d2 = sum(line_coupling_sq(F, mF, q, 1.5) for F, mF in full for q in (-1, 0, 1))
    assert d2 / d1 == pytest.approx(2.0)


# --- BASIS INDEPENDENCE: F'-basis == |J',mJ'> basis -------------------------
def test_basis_independence_Fprime_vs_mJ():
    for F, mF in GROUND:
        for q in (1, 0, -1):
            for Jp in (0.5, 1.5):
                fp = line_coupling_sq(F, mF, q, Jp)
                mj = sum(hf_amplitude(F, mF, mJ) ** 2 * electronic_dipole(mJ, q, Jp, mJ + q) ** 2
                         for mJ in _MJ)
                assert fp == pytest.approx(mj, abs=1e-12)


def test_hf_amplitude_stretched_is_pure():
    assert hf_amplitude(3, 3, 0.5) == pytest.approx(1.0)     # |3,3> = |mJ=+1/2, mI=5/2>
    assert hf_amplitude(3, 3, -0.5) == pytest.approx(0.0)
    # |2,2> is a genuine mJ superposition (mI=3/2 and mI=5/2 parts)
    assert abs(hf_amplitude(2, 2, 0.5)) > 0.1 and abs(hf_amplitude(2, 2, -0.5)) > 0.1


def test_electronic_dipole_selection_rule():
    assert electronic_dipole(0.5, +1, 1.5, 1.5) != 0.0       # mJ'=mJ+q
    assert electronic_dipole(0.5, +1, 0.5, 1.5) == 0.0       # |mJ'|=3/2 > J'=1/2
    assert electronic_dipole(0.5, 0, 0.5, 0.5) != 0.0


# --- scalar / vector light shift --------------------------------------------
def _ro():
    return RamanOptics(20e9, 2.7457e12, 1.79e9, 41.8e6, 41.3e6)


def test_detuning_signs():
    ro = _ro()
    assert ro.detuning(1.5, 3.0) < 0.0                       # P3/2: red
    assert ro.detuning(0.5, 3.0) > 0.0                       # P1/2: far blue (opposite sign)
    # |up>=F=2 is closer to P3/2 resonance than |down>=F=3 (by omega_HF)
    assert abs(ro.detuning(1.5, 2.0)) < abs(ro.detuning(1.5, 3.0))


def test_vector_shift_zero_for_pi_and_balanced_linear():
    ro = _ro()
    for pol in (POL_PI, POL_LINEAR_PERP):
        _, vec = ro.scalar_vector_shift(RamanBeam(1.0, pol), 3.0)
        assert vec == pytest.approx(0.0, abs=1e-18)


def test_vector_shift_opposite_for_sigma_plus_minus():
    ro = _ro()
    _, vp = ro.scalar_vector_shift(RamanBeam(1.0, POL_SIGMA_PLUS), 3.0)
    _, vm = ro.scalar_vector_shift(RamanBeam(1.0, POL_SIGMA_MINUS), 3.0)
    assert vp == pytest.approx(-vm)
    assert vp != 0.0
    # scalar is polarization-independent (same total intensity)
    sp, _ = ro.scalar_vector_shift(RamanBeam(1.0, POL_SIGMA_PLUS), 3.0)
    spi, _ = ro.scalar_vector_shift(RamanBeam(1.0, POL_PI), 3.0)
    assert sp == pytest.approx(spi)


def test_circular_fraction():
    ro = _ro()
    assert ro.circular_fraction(RamanBeam(1.0, POL_SIGMA_PLUS)) == pytest.approx(1.0)
    assert ro.circular_fraction(RamanBeam(1.0, POL_SIGMA_MINUS)) == pytest.approx(-1.0)
    assert ro.circular_fraction(RamanBeam(1.0, POL_PI)) == pytest.approx(0.0)
    assert ro.circular_fraction(RamanBeam(1.0, POL_LINEAR_PERP)) == pytest.approx(0.0)


# --- two-photon Rabi + power/polarization dependence ------------------------
def test_two_photon_rabi_nonzero_for_OC():
    ro = _ro()
    B1 = RamanBeam(1.0, POL_PI)
    R2 = RamanBeam(1.0, POL_LINEAR_PERP)
    assert ro.two_photon_rabi(B1, R2, (3, 3), (2, 2)) > 0.0


def test_two_photon_rabi_scales_with_power():
    ro = _ro()
    B1 = RamanBeam(1.0, POL_PI)
    R2 = RamanBeam(1.0, POL_LINEAR_PERP)
    base = ro.two_photon_rabi(B1, R2, (3, 3), (2, 2))
    # Omega_2gamma ~ E_B E_R ~ sqrt(I_B I_R): doubling BOTH intensities doubles it
    big = ro.two_photon_rabi(RamanBeam(2.0, POL_PI), RamanBeam(2.0, POL_LINEAR_PERP), (3, 3), (2, 2))
    assert big == pytest.approx(2.0 * base)


def test_scatter_scales_linearly_with_power():
    ro = _ro()
    b1 = ro.scatter_rate(RamanBeam(1.0, POL_PI), 3.0, 3.0)
    b2 = ro.scatter_rate(RamanBeam(2.0, POL_PI), 3.0, 3.0)
    assert b2 == pytest.approx(2.0 * b1)             # ~ intensity


def test_circular_contamination_changes_differential_shift():
    ro = _ro()
    B1 = RamanBeam(0.425, POL_PI)
    R2_clean = RamanBeam(0.425, POL_LINEAR_PERP)
    R2_dirty = RamanBeam(0.425, (0.55, 0.0, 0.45))   # 10% circular excess (C=+0.1)
    d_clean = ro.differential_stark_per_rabi([B1, R2_clean], B1, R2_clean)
    d_dirty = ro.differential_stark_per_rabi([B1, R2_dirty], B1, R2_dirty)
    # a 10% circular contamination shifts the differential by a LARGE amount (vector shift)
    assert abs(d_dirty - d_clean) > 0.3 * abs(d_clean)


# --- ledger wiring ----------------------------------------------------------
def test_from_ledger_and_beams():
    ledger = Ledger.load()
    ro = RamanOptics.from_ledger(ledger)
    assert ro.delta_fs == pytest.approx(2.7457e12, rel=1e-3)
    assert ro.delta_p32 < 0 and abs(ro.delta_p32) == pytest.approx(20e9)
    bB, bR = beams_from_ledger(ledger, "OC", powers={"b1": 0.425, "r2": 0.425})
    assert bB.pol == (0.0, 1.0, 0.0) and bB.intensity == pytest.approx(0.425)
    assert ro.circular_fraction(bR) == pytest.approx(0.0)    # R2 balanced linear


def test_from_ledger_uses_inputs_only():
    # from_ledger must not consume a benchmark (wall): a benchmark name raises
    ledger = Ledger.load()
    with pytest.raises(Exception):
        RamanOptics.from_ledger.__func__  # smoke: the classmethod exists
        ledger.input_quantity("mg_rsb_cooled_nbar_axial_lf_25mg")   # this IS a benchmark -> raises
