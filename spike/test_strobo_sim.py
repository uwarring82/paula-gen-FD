"""
Tests for the stroboscopic spin-motion simulator (the active-phase-grating detuning
scan): the displacement matrix D(i*eta), and the detuning spectrum -- a comb of
resonances at delta = k*f_lf (carrier + first sidebands), narrow (set by N*DELTA_t),
with the off-resonant regions dark.

Run:  pytest spike/
"""
import math

import pytest

from spike.engines.strobo_sim import (
    displacement_matrix, strobo_detuning_scan, strobo_population_vs_cycles)

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


def test_carrier_and_first_floquet_sidebands_present():
    dets, P = _scan()
    def at(target_flf):
        k = min(range(len(dets)), key=lambda i: abs(dets[i] - target_flf * F_LF))
        return P[k]
    assert at(0.0) > 0.9                                          # carrier full flop (delta_t~pi/2)
    assert at(+1.0) > 0.9                                         # +1st FLOQUET tooth (NOT a blue motional SB)
    assert at(-1.0) > 0.9                                         # -1st FLOQUET tooth (symmetric, NOT red SB)


def test_comb_is_symmetric_floquet_independent_of_eta():
    # At the EXACT strobe (DELTA_t = 1/f_lf) the teeth are FLOQUET sidebands of the
    # pulsed drive, not motional ones: they are full-contrast even at eta=0 (no motion)
    # and do not change with eta. (This is the corrected claim; the old docstring wrongly
    # called the comb red/blue asymmetric.)
    dets = [0.0, F_LF, -F_LF]
    vals = []
    for eta in (0.0, 0.2, 0.389):
        P = strobo_detuning_scan(eta, 4.99e5, 0.02, 0.769172, 50, F_LF, dets, F=12)
        assert P[1] == pytest.approx(P[2], abs=1e-3)             # +f_lf tooth == -f_lf tooth (symmetric)
        vals.append((P[0], P[1]))
    # eta-independence: carrier and +1 tooth barely move from eta=0 to eta=0.389
    assert vals[0][0] == pytest.approx(vals[-1][0], abs=2e-3)
    assert vals[0][1] == pytest.approx(vals[-1][1], abs=2e-3)


def test_motional_coupling_appears_when_strobe_is_detuned():
    # The eta (motional) coupling is hidden at the exact strobe (the motion wraps to
    # identity each cycle); detune the strobe off the motional period and it shows up.
    f_mis = F_LF * 1.03                                          # motion 3% off the strobe period
    P0 = strobo_detuning_scan(0.0, 4.99e5, 0.02, 0.769172, 50, f_mis, [0.0], F=12)[0]
    Pe = strobo_detuning_scan(0.389, 4.99e5, 0.02, 0.769172, 50, f_mis, [0.0], F=12)[0]
    assert P0 > 0.99                                             # eta=0: carrier still ~full
    assert Pe < P0 - 5e-3                                        # eta>0: motional coupling now bites


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


def test_heterodyne_cycle_domain_beat():
    # The train is a sampling mixer: ON a comb tooth (f_IF=0) the cycle-domain population
    # builds up to the pi flop; detuned by f_IF it NUTATES, and the first turning point
    # tracks the half-beat ~ 1/(2 f_IF DELTA_t) -- the down-converted intermediate freq.
    DT = 0.769172
    on = strobo_population_vs_cycles(0.389, 4.99e5, 0.02, DT, 60, F_LF, 0.0, F=10)
    assert on[49] > 0.9 and on[0] < 0.05                        # monotone build to ~pi at N~50
    for f_IF, n_pred in ((50e3, 13), (100e3, 6)):
        P = strobo_population_vs_cycles(0.389, 4.99e5, 0.02, DT, 60, F_LF, f_IF, F=10)
        turn = next(n for n in range(2, 60) if P[n] < P[n - 1])  # first turning point
        assert abs(turn - n_pred) <= 2                          # ~ 1/(2 f_IF DELTA_t)
