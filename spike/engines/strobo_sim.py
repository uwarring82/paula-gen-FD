"""
Stroboscopic spin-motion simulator for the active-phase-grating DETUNING scan.

The "active phase grating" applies N identical OC carrier pulses, one per period
DELTA_t (locked to the lf motional period), each of width delta_t. Scanning the DRIVE
detuning delta (the offset of the strobo frequency from the qubit carrier) probes the
Floquet/stroboscopic resonance comb: the spin responds whenever delta matches a comb
tooth, delta = k / DELTA_t = k * f_lf -- the carrier (k=0) and the +-1, +-2, ... teeth
at multiples of the strobe (= motional) frequency.

    U_pulse = exp(-i (Omega_strobo/2)(s+ (x) D(i eta) + s- (x) D(i eta)^dag) delta_t)
    U_free(delta) = diag exp(-i (E_s + n omega_lf) DELTA_t),  E_up=-delta/2, E_down=+delta/2
    P_flip(delta) = sum_n |<down,n| (U_free(delta) U_pulse)^N |up,0>|^2.

NATURE OF THE COMB (corrected -- see the 2026-06-23 review). The teeth are FLOQUET
sidebands of the periodic pulsed drive, NOT ordinary motional red/blue sidebands. The
comb is SYMMETRIC: at the exact strobe (DELTA_t = 1/f_lf) the +-k teeth are full-
contrast and INDEPENDENT of eta -- they survive even at eta = 0 (no motion at all). The
reason: with DELTA_t = the motional period the free motional phase wraps to 2*pi each
cycle, so D(i eta) is a coherent passenger that cancels in the spin-flip sum sum_n. The
ordinary sideband-thermometry asymmetry (red Delta_n=-1 forbidden from |0>) does NOT
apply, because the train coherently drives multi-phonon Floquet paths that alias onto
every tooth. The MOTIONAL coupling (eta) becomes visible only when the strobe is
DETUNED off the motional period (the phase stops wrapping) or the motion is DISPLACED
(displ != 0) -- see strobo_population_vs_cycles and the eta tests.

HETERODYNE / MIXER VIEW. The train is a sampling mixer: the strobe is a local
oscillator (comb at 1/DELTA_t = f_lf), the spin-motion coupling is the signal, and the
coherent spin accumulation is the mixer + integrator. The detuning delta is the
intermediate frequency: ON a tooth (delta = k*f_lf) the signal down-converts to DC and
builds up (homodyne, the pi flop); off a tooth by f_IF the cycle-domain population
NUTATES at f_IF (heterodyne) -- strobo_population_vs_cycles. The narrow tooth width
1/(N*DELTA_t) is the finite-integration IF bandwidth.

APPROXIMATIONS. (1) INSTANTANEOUS pulse: U_pulse is applied at a point and U_free runs
the full DELTA_t. The per-period DETUNING phase is correctly delta*DELTA_t (so the comb
POSITIONS are exact -- the delta*delta_t accrued during the pulse is lumped into
U_free); what is neglected is the motion's evolution DURING the pulse (the pulse samples
the motional phase over a window omega_lf*delta_t ~ 0.16 rad at delta_t = 0.02 us,
growing to ~0.8 rad at delta_t = 0.1 us). (2) No decoherence/noise model -- this engine
is the ideal coherent propagator only. Pure Python (built-in complex), no numpy: a few
small (2F x 2F) matrices.
"""
from __future__ import annotations

import cmath
import math


def _laguerre_gen(n: int, k: int, x: float) -> float:
    """Generalized Laguerre L_n^{(k)}(x) (three-term recurrence)."""
    if n <= 0:
        return 1.0
    lkm1, lk = 1.0, 1.0 + k - x
    for j in range(1, n):
        lkm1, lk = lk, ((2 * j + 1 + k - x) * lk - (j + k) * lkm1) / (j + 1)
    return lk


def displacement_matrix(eta: float, F: int):
    """<m|D(i eta)|n> for m,n in 0..F-1 (F x F complex). D(i eta) is unitary, so the
    returned matrix satisfies M M^dag = I (to the Fock truncation)."""
    a = 1j * eta
    a2 = eta * eta
    e = math.exp(-a2 / 2.0)
    fact = [1.0] * (F + 1)
    for i in range(1, F + 1):
        fact[i] = fact[i - 1] * i
    M = [[0j] * F for _ in range(F)]
    for m in range(F):
        for n in range(F):
            if m >= n:
                k = m - n
                M[m][n] = math.sqrt(fact[n] / fact[m]) * (a ** k) * e * _laguerre_gen(n, k, a2)
            else:
                k = n - m
                M[m][n] = math.sqrt(fact[m] / fact[n]) * ((-a.conjugate()) ** k) * e * _laguerre_gen(m, k, a2)
    return M


def _matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    C = [[0j] * p for _ in range(n)]
    for i in range(n):
        Ai, Ci = A[i], C[i]
        for k in range(m):
            a = Ai[k]
            if a == 0:
                continue
            Bk = B[k]
            for j in range(p):
                Ci[j] += a * Bk[j]
    return C


def _expm(A, terms: int = 18):
    """Matrix exponential via Taylor series (A here has small norm: Omega*delta_t<<1)."""
    n = len(A)
    R = [[1j * 0 + (1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]  # I
    term = [row[:] for row in R]
    for k in range(1, terms):
        term = _matmul(term, A)
        inv = 1.0 / k
        for i in range(n):
            ti, ri = term[i], R[i]
            for j in range(n):
                ti[j] *= inv
                ri[j] += ti[j]
    return R


def _matpow(A, p: int):
    n = len(A)
    R = [[(1.0 if i == j else 0.0) + 0j for j in range(n)] for i in range(n)]  # I
    base = [row[:] for row in A]
    while p > 0:
        if p & 1:
            R = _matmul(R, base)
        p >>= 1
        if p:
            base = _matmul(base, base)
    return R


def _matvec(A, v):
    n = len(A)
    return [sum(A[i][k] * v[k] for k in range(len(v))) for i in range(n)]


def strobo_population_vs_cycles(eta: float, omega_strobo_hz: float, delta_t_us: float,
                                DELTA_t_us: float, n_cycles: int, f_lf_hz: float,
                                detuning_hz: float, F: int = 10):
    """P_flip after each of n = 1..n_cycles stroboscopic cycles at a FIXED drive detuning
    -- the CYCLE-DOMAIN (slow) signal. ON a comb tooth (detuning = k*f_lf) it builds up
    monotonically (the pi flop); detuned by f_IF = detuning - k*f_lf from a tooth it
    NUTATES, and for f_IF above the per-cycle Rabi the nutation frequency IS f_IF: the
    train acts as a sampling MIXER that down-converts the strobe<->transition beat into
    the slow cycle domain (homodyne at f_IF=0, heterodyne otherwise). Returns the list
    [P_flip(after 1 cycle), ..., P_flip(after n_cycles)]."""
    M = displacement_matrix(eta, F)
    D = 2 * F
    om = 2.0 * math.pi * omega_strobo_hz
    H = [[0j] * D for _ in range(D)]
    for m in range(F):
        for n in range(F):
            H[F + m][n] += 0.5 * om * M[m][n]
            H[n][F + m] += 0.5 * om * M[m][n].conjugate()
    U_pulse = _expm([[-1j * H[i][j] * (delta_t_us * 1e-6) for j in range(D)] for i in range(D)])
    wlf = 2.0 * math.pi * f_lf_hz
    dt = DELTA_t_us * 1e-6
    dang = 2.0 * math.pi * detuning_hz
    Ufree = [[0j] * D for _ in range(D)]
    for s in range(2):
        E = (-dang / 2.0) if s == 0 else (+dang / 2.0)
        for n in range(F):
            Ufree[s * F + n][s * F + n] = cmath.exp(-1j * (E + n * wlf) * dt)
    Ucyc = _matmul(Ufree, U_pulse)
    psi = [0j] * D
    psi[0] = 1.0 + 0j                            # |up, 0>
    out = []
    for _ in range(n_cycles):
        psi = _matvec(Ucyc, psi)
        out.append(sum((psi[F + n] * psi[F + n].conjugate()).real for n in range(F)))
    return out


def strobo_detuning_scan(eta: float, omega_strobo_hz: float, delta_t_us: float,
                         DELTA_t_us: float, n_cycles: int, f_lf_hz: float,
                         detunings_hz, F: int = 8):
    """P_flip(delta) for the stroboscopic pulse train vs the DRIVE detuning delta [Hz]
    (relative to the carrier). Start |up, 0>. Returns a list aligned with detunings_hz.
    delta_t fixed; eta, Omega_strobo, DELTA_t, N, f_lf the apparatus parameters."""
    M = displacement_matrix(eta, F)
    D = 2 * F                                   # spin (0=up,1=down) x Fock; idx = s*F + n
    om = 2.0 * math.pi * omega_strobo_hz
    # H_pulse coupling (Hermitian): (Omega/2)(|down,m><up,n| M[m][n] + h.c.)
    H = [[0j] * D for _ in range(D)]
    for m in range(F):
        for n in range(F):
            H[F + m][n] += 0.5 * om * M[m][n]
            H[n][F + m] += 0.5 * om * M[m][n].conjugate()
    U_pulse = _expm([[-1j * H[i][j] * (delta_t_us * 1e-6) for j in range(D)] for i in range(D)])
    wlf = 2.0 * math.pi * f_lf_hz
    dt = DELTA_t_us * 1e-6
    out = []
    for det in detunings_hz:
        dang = 2.0 * math.pi * det
        Ufree = [[0j] * D for _ in range(D)]
        for s in range(2):
            E = (-dang / 2.0) if s == 0 else (+dang / 2.0)
            for n in range(F):
                Ufree[s * F + n][s * F + n] = cmath.exp(-1j * (E + n * wlf) * dt)
        Utot = _matpow(_matmul(Ufree, U_pulse), n_cycles)
        p = 0.0
        for n in range(F):
            amp = Utot[F + n][0]                # <down,n| U |up,0>
            p += (amp * amp.conjugate()).real
        out.append(p)
    return out
