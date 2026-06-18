"""
Tests for the drive engine — the relative magnetic-dipole (Clebsch-Gordan) Rabi
couplings of the 25Mg+ ground-state hyperfine transitions. The reference |CG|
values were cross-checked against sympy.physics.quantum.cg.CG.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.drive import HyperfineDrive, clebsch_gordan

# |CG(3,mF; 1,q | 2,mF')| for the 8 Doerr/Kaufmann transitions (verified vs sympy)
KNOWN = {
    (3, 2): 0.845154, (1, 2): 0.218218, (1, 0): 0.534522, (0, 0): 0.654654,
    (0, -1): 0.377964, (-1, -1): 0.617213, (-1, -2): 0.218218, (-3, -2): 0.845154,
}


def test_couplings_match_reference_cg():
    d = HyperfineDrive(3, 2)
    for (a, b), v in KNOWN.items():
        assert d.coupling(a, b) == pytest.approx(v, abs=1e-5)


def test_clebsch_gordan_known_identities():
    assert abs(clebsch_gordan(1, 0, 1, 0, 0, 0)) == pytest.approx(1 / math.sqrt(3))
    assert abs(clebsch_gordan(1, 1, 1, -1, 2, 0)) == pytest.approx(1 / math.sqrt(6))
    assert abs(clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0)) == pytest.approx(1 / math.sqrt(2))


def test_mirror_symmetry():
    d = HyperfineDrive(3, 2)
    for (a, b) in KNOWN:
        assert d.coupling(a, b) == pytest.approx(d.coupling(-a, -b))


def test_polarization_labels():
    d = HyperfineDrive(3, 2)
    assert d.polarization(0, 0) == "pi"        # the clock
    assert d.polarization(3, 2) == "sigma-"
    assert d.polarization(1, 2) == "sigma+"


def test_forbidden_transitions_are_zero():
    d = HyperfineDrive(3, 2)
    assert d.coupling(3, 0) == 0.0    # |dmF| = 3, beyond a rank-1 (dipole) photon
    assert d.coupling(3, -2) == 0.0


def test_relative_couplings_normalised():
    d = HyperfineDrive(3, 2)
    rc = d.relative_couplings([(3, 2), (1, 2), (0, 0)])
    assert max(rc.values()) == pytest.approx(1.0)
    assert rc[(1, 2)] < rc[(0, 0)] < rc[(3, 2)]


def test_absolute_rabi_factorisation():
    # Omega = |CG| * apparatus_factor; the empirical factor round-trips
    d = HyperfineDrive(3, 2)
    measured = 59453.0
    af = d.apparatus_factor(3, 2, measured)
    assert af == pytest.approx(measured / d.coupling(3, 2))
    assert d.absolute_rabi(3, 2, af) == pytest.approx(measured)
