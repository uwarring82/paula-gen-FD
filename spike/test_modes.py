"""
Tests for the modes engine + the pure-Python linear algebra. Pins the axial
normal-mode eigenvalues against known values (James 1998: N=2 -> 1,3;
N=3 -> 1,3,5.8), the sqrt(3) stretch relation, the wall, and reproduction of the
Wittemer 2-ion stretch benchmark.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.modes import (
    AxialModes,
    RadialModes,
    axial_mode_eigenvalues,
    equilibrium_positions,
)
from spike.ledger import Ledger
from spike.linalg import eigvalsh, solve


# --- linear algebra ---------------------------------------------------------
def test_solve_2x2():
    x = solve([[2.0, 1.0], [1.0, 3.0]], [3.0, 5.0])
    assert x[0] == pytest.approx(0.8) and x[1] == pytest.approx(1.4)


def test_eigvalsh_known():
    assert eigvalsh([[2.0, -1.0], [-1.0, 2.0]]) == pytest.approx([1.0, 3.0])


# --- equilibrium + eigenvalues ----------------------------------------------
def test_equilibrium_two_ion():
    u = equilibrium_positions(2)
    half = (1.0 / 4.0) ** (1.0 / 3.0)   # = 0.62996
    assert u[0] == pytest.approx(-half, abs=1e-6)
    assert u[1] == pytest.approx(+half, abs=1e-6)


def test_axial_eigenvalues_universal_first_two():
    for N in (2, 3, 4, 5):
        lam = axial_mode_eigenvalues(N)
        assert lam[0] == pytest.approx(1.0, abs=1e-9)   # COM
        assert lam[1] == pytest.approx(3.0, abs=1e-6)   # stretch


def test_axial_eigenvalues_three_ion():
    assert axial_mode_eigenvalues(3) == pytest.approx([1.0, 3.0, 5.8], abs=1e-3)


# --- engine ------------------------------------------------------------------
def test_stretch_is_sqrt3_com():
    e = AxialModes(1.30e6)
    assert e.stretch_frequency(2) == pytest.approx(math.sqrt(3) * 1.30e6, rel=1e-9)


def test_com_mode_equals_single_ion_frequency():
    e = AxialModes(2.0e6)
    for N in (2, 3, 4):
        assert e.mode_frequencies(N)[0] == pytest.approx(2.0e6, rel=1e-9)


def test_stretch_requires_two_ions():
    with pytest.raises(ValueError):
        AxialModes(1.30e6).stretch_frequency(1)


def test_from_ledger_reproduces_stretch_benchmark():
    ledger = Ledger.load()
    e = AxialModes.from_ledger(ledger)
    predicted = e.stretch_frequency(2)
    bench = ledger.benchmark_quantity("omega_z_axial_stretch_2ion_25mg")
    # sqrt(3) relation holds to ~1% in the 3-sig-fig data (~22 kHz)
    assert abs(predicted - bench.value) < 30e3


def test_from_ledger_refuses_benchmark_as_com():
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        AxialModes.from_ledger(ledger, com_name="omega_z_axial_stretch_2ion_25mg")


# --- robustness (review hardening) ------------------------------------------
def test_solve_rejects_nonsquare():
    with pytest.raises(ValueError):
        solve([[1.0, 2.0]], [1.0])


def test_eigvalsh_raises_when_not_converged():
    A = [[2.0, 1.0, 0.5], [1.0, 3.0, 1.0], [0.5, 1.0, 4.0]]
    with pytest.raises(RuntimeError):
        eigvalsh(A, max_sweeps=1)


def test_equilibrium_raises_when_not_converged():
    with pytest.raises(RuntimeError):
        equilibrium_positions(10, max_iter=2)


def test_large_chain_converges():
    # N >= 40 used to stall at the 2-norm noise floor; the inf-norm criterion fixes it
    lam = axial_mode_eigenvalues(40)
    assert lam[0] == pytest.approx(1.0, abs=1e-9)   # COM
    assert lam[1] == pytest.approx(3.0, abs=1e-7)   # stretch
    assert all(x >= 0.0 for x in lam)


# --- radial (transverse) modes ----------------------------------------------
def test_radial_rocking_two_ion():
    r = RadialModes(1.30e6, 2.88e6)
    expected = math.sqrt((2.88e6) ** 2 - (1.30e6) ** 2)   # sqrt(omega_r^2 - omega_z^2)
    assert r.rocking_frequency(2) == pytest.approx(expected, rel=1e-9)


def test_radial_com_is_omega_r():
    r = RadialModes(1.30e6, 2.88e6)
    assert r.com_frequency() == 2.88e6
    assert r.mode_frequencies(2)[-1] == pytest.approx(2.88e6, rel=1e-9)   # COM is the highest
    assert r.mode_frequencies(2)[0] < r.mode_frequencies(2)[1]            # rocking below COM


def test_radial_instability_raises_for_soft_radial():
    with pytest.raises(ValueError):
        RadialModes(2.0e6, 1.0e6).rocking_frequency(2)   # omega_r < omega_z -> zigzag


def test_radial_from_ledger_reproduces_calc():
    ledger = Ledger.load()
    r = RadialModes.from_ledger(ledger)
    bench = ledger.benchmark_quantity("omega_radial_rocking_2ion_25mg")
    assert abs(r.rocking_frequency(2) - bench.value) < 1e3   # matches the 2.57 MHz calc to ~0.1 kHz


def test_radial_from_ledger_refuses_benchmark():
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        RadialModes.from_ledger(ledger, radial_com_name="omega_radial_rocking_2ion_25mg")
