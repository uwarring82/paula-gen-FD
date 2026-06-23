"""
Tests for the stroboscopic spin-motion simulator (the active-phase-grating detuning
scan): the displacement matrix D(i*eta), and the detuning spectrum -- a comb of
resonances at delta = k*f_lf (carrier + first sidebands), narrow (set by N*DELTA_t),
with the off-resonant regions dark.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.strobo_sim import displacement_matrix, strobo_detuning_scan

ETA = 0.389
F_LF = 1.3001e6


def test_displacement_matrix_elements():
    M = displacement_matrix(ETA, 8)
    dw = math.exp(-ETA * ETA / 2.0)
    assert abs(M[0][0]) == pytest.approx(dw, rel=1e-6)            # carrier = Debye-Waller
    assert abs(M[1][0]) == pytest.approx(ETA * dw, rel=1e-3)      # blue ~ eta * DW
    assert abs(M[0][1]) == pytest.approx(abs(M[1][0]))            # |red| = |blue| magnitude
    # the ground-state column is normalised (D unitary; |0> couples only to low n)
    col0 = sum(abs(M[m][0]) ** 2 for m in range(8))
    assert col0 == pytest.approx(1.0, abs=1e-6)


def _scan(span=1.6, npts=321, omega=4.99e5, delta_t=0.02):
    dets = [(-span + 2 * span * k / (npts - 1)) * F_LF for k in range(npts)]
    P = strobo_detuning_scan(ETA, omega, delta_t, 0.769172, 50, F_LF, dets, F=10)
    return dets, P


def test_carrier_and_first_sidebands_present():
    dets, P = _scan()
    def at(target_flf):
        k = min(range(len(dets)), key=lambda i: abs(dets[i] - target_flf * F_LF))
        return P[k]
    assert at(0.0) > 0.9                                          # carrier full flop (delta_t~pi/2)
    assert at(+1.0) > 0.9                                         # first blue-side comb tooth
    assert at(-1.0) > 0.9                                         # first red-side comb tooth (aliased)


def test_off_resonance_is_dark():
    dets, P = _scan()
    def at(target_flf):
        k = min(range(len(dets)), key=lambda i: abs(dets[i] - target_flf * F_LF))
        return P[k]
    assert at(0.5) < 0.2                                          # between teeth -> dark
    assert at(1.5) < 0.2


def test_probability_bounded():
    _dets, P = _scan()
    assert all(-1e-9 <= p <= 1.0 + 1e-6 for p in P)              # valid populations (no leak)


def test_no_pulse_no_flip():
    # delta_t -> 0 (no pulse area) leaves the spin in |up> -> P_flip ~ 0 everywhere
    _dets, P = _scan(delta_t=0.0)
    assert max(P) < 1e-6
