"""
Tests for engines.grating_tomography -- the characteristic-function / Wigner transfer
functions of the stroboscopic phase grating (docs/notes/strobo_grating_transfer_function.md).

Checks: chi(beta) identities; the exact eta=0 kicked-two-level formula vs the full Floquet
propagator (machine precision, all orders, arbitrary eta on the exact strobe); the
weak-pulse double-sum kernel vs Floquet; chi/Wigner Fourier consistency; and that the
chi->W reconstruction recovers the analytic Wigner of vacuum / coherent / cat states.

Run:  pytest spike/
"""
import cmath
import math

import pytest

from spike.engines import grating_tomography as gt
from spike.engines.strobo_sim import strobo_detuning_scan

DT = 0.769172
OM = 4.99e5
N = 50


# ---- characteristic functions ----------------------------------------------
def test_chi_values_at_origin_and_vacuum():
    assert gt.chi_vacuum(0j) == pytest.approx(1.0)               # Tr rho = 1
    assert gt.chi_coherent(0j, 1 + 0.5j) == pytest.approx(1.0)
    assert gt.chi_fock(0j, 3) == pytest.approx(1.0)
    assert gt.chi_thermal(0j, 2.0) == pytest.approx(1.0)
    assert gt.chi_cat(0j, 1.2 + 0j, +1) == pytest.approx(1.0)
    # vacuum is the nbar=0 thermal and the n=0 Fock
    b = 0.7 - 0.3j
    assert gt.chi_thermal(b, 0.0) == pytest.approx(gt.chi_vacuum(b).real)
    assert gt.chi_fock(b, 0) == pytest.approx(gt.chi_vacuum(b))


def test_chi_coherent_is_displaced_vacuum_modulus():
    # |chi_coherent(beta,gamma)| = |chi_vacuum(beta)| (the displacement is a pure phase)
    b, g = 0.8 + 0.2j, 1.0 - 0.6j
    assert abs(gt.chi_coherent(b, g)) == pytest.approx(abs(gt.chi_vacuum(b)))


# ---- exact eta=0 transfer function -----------------------------------------
def test_exact_eta0_matches_floquet_machine_precision():
    f_ex = 1.0 / (DT * 1e-6)
    for dt in (0.02, 0.05):
        th = gt.theta_of(OM, dt)
        for dfac in (0.0, 0.3, 0.7, 1.0):
            det = dfac * f_ex
            p_ana = gt.exact_eta0_probability(det, th, N, DT)
            p_num = strobo_detuning_scan(0.0, OM, dt, DT, N, f_ex, [det], F=4)[0]
            assert abs(p_ana - p_num) < 1e-9


def test_exact_strobe_motion_blind_all_orders():
    # All-orders result: on the exact strobe P_flip == eta=0 formula for ANY eta, even at
    # the pi pulse (theta not small).
    f_ex = 1.0 / (DT * 1e-6)
    th = gt.theta_of(OM, 0.02)                                   # ~pi pulse over N=50
    for det in (0.0, 0.3 * f_ex, 1.0 * f_ex):
        p0 = gt.exact_eta0_probability(det, th, N, DT)
        for eta in (0.2, 0.389):
            p_num = strobo_detuning_scan(eta, OM, 0.02, DT, N, f_ex, [det], F=16)[0]
            assert p_num == pytest.approx(p0, abs=1e-6)


def test_exact_eta0_on_tooth_is_full_flop():
    # delta=0, N*theta = pi  ->  P = sin^2(N theta/2) = 1
    th = math.pi / N
    assert gt.exact_eta0_probability(0.0, th, N, DT) == pytest.approx(1.0, abs=1e-12)


# ---- weak-pulse double-sum kernel ------------------------------------------
def test_kernel_probability_matches_floquet_weak_pulse():
    f_lf = 1.03 / (DT * 1e-6)                                    # off-strobe 3%
    dt = 0.001                                                  # weak pulse
    th = gt.theta_of(OM, dt)
    p_ana = gt.kernel_probability(gt.chi_vacuum, 0.389, 0.0, 0.0, N, DT, f_lf, th)
    p_num = strobo_detuning_scan(0.389, OM, dt, DT, N, f_lf, [0.0], F=16)[0]
    assert abs(p_ana - p_num) / p_num < 5e-3                     # ~0.18% at this theta


def test_kernel_probability_motion_blind_on_exact_strobe():
    # On the exact strobe the kernel is eta-independent (chi(0)=1 only).
    f_ex = 1.0 / (DT * 1e-6)
    th = gt.theta_of(OM, 0.001)
    p_eta = gt.kernel_probability(gt.chi_vacuum, 0.389, 0.0, 0.0, N, DT, f_ex, th)
    p_thermal = gt.kernel_probability(lambda b: gt.chi_thermal(b, 3.0), 0.389, 0.0, 0.0, N, DT, f_ex, th)
    assert p_eta == pytest.approx(p_thermal, rel=1e-9)          # state-independent on strobe


# ---- coherence channel = DFT of chi ----------------------------------------
def test_kernel_coherence_on_strobe_reads_chi_at_one_point():
    # On the exact strobe Tr[rho A] = -i theta/2 chi(i eta e^{i phi}) S_N(delta).
    f_ex = 1.0 / (DT * 1e-6)
    eta, phi = 0.389, 0.4
    th = gt.theta_of(OM, 0.001)
    coh = gt.kernel_coherence(lambda b: gt.chi_coherent(b, 0.5 + 0.2j), eta, phi, 0.0, N, DT, f_ex, th)
    assert isinstance(coh, complex)                             # API: coherence is complex (Re,Im)
    S_N = N                                                     # delta=0 -> S_N = N
    expect = -0.5j * th * gt.chi_coherent(1j * eta * cmath.exp(1j * phi), 0.5 + 0.2j) * S_N
    assert coh == pytest.approx(expect, rel=1e-9)


# ---- chi -> W reconstruction -----------------------------------------------
def _recon_err(chi, wig, half=5.0, nb=41, ahalf=1.5, na=13):
    beta_pts, db = gt.square_grid(half, nb)
    chi_vals = [chi(b) for b in beta_pts]
    step = 2 * ahalf / (na - 1)
    alpha_pts = [complex(-ahalf + i * step, -ahalf + j * step) for j in range(na) for i in range(na)]
    W_rec = gt.wigner_from_samples(beta_pts, chi_vals, alpha_pts, db * db)
    W_ana = [wig(a) for a in alpha_pts]
    return max(abs(r - a) for r, a in zip(W_rec, W_ana)), max(abs(a) for a in W_ana)


def test_reconstruct_wigner_vacuum():
    err, peak = _recon_err(gt.chi_vacuum, gt.wigner_vacuum)
    assert err < 1e-3 * peak + 1e-4


def test_reconstruct_wigner_coherent():
    g = 0.8 + 0.5j
    err, _ = _recon_err(lambda b: gt.chi_coherent(b, g), lambda a: gt.wigner_coherent(a, g))
    assert err < 5e-4


def test_reconstruct_wigner_cat_has_negativity():
    # The ODD cat (|g>-|-g>) has W(0) = -2/pi < 0 -- a Wigner-negativity witness the
    # chi->W reconstruction must reproduce. (The EVEN cat has W(0)=+2/pi; both show fringes.)
    g = 1.6 + 0j
    beta_pts, db = gt.square_grid(7.0, 81)
    chi_vals = [gt.chi_cat(b, g, -1) for b in beta_pts]
    W0_rec = gt.wigner_from_samples(beta_pts, chi_vals, [0j], db * db)[0]
    assert gt.wigner_cat(0j, g, -1) == pytest.approx(-2.0 / math.pi, abs=1e-6)
    assert W0_rec < 0.0                                          # reconstruction reproduces it
    assert W0_rec == pytest.approx(gt.wigner_cat(0j, g, -1), abs=5e-3)


def test_chi_wigner_fourier_consistency_fock():
    # chi_fock and wigner_fock must be a 2-D Fourier pair (reconstruct W_1 from chi_1).
    err, _ = _recon_err(lambda b: gt.chi_fock(b, 1), lambda a: gt.wigner_fock(a, 1))
    assert err < 5e-3


# ---- Ramsey characteristic-function interferometer -------------------------
def test_ramsey_identity_exact_vs_analytic():
    # Two recoil-dressed pi/2 pulses: the exact unitary population == the boxed identity
    # P_down(phi) = 1/2[1 + Re(e^{i[phi+Im(bg br*)]} chi(bg-br))], for pure states.
    br, bg = 0.15 + 0.05j, 0.45 + 0.25j
    cases = [(gt.chi_vacuum, gt.coherent_state_vec(0j, 30)),
             (lambda b: gt.chi_coherent(b, 0.9 - 0.4j), gt.coherent_state_vec(0.9 - 0.4j, 30)),
             (lambda b: gt.chi_fock(b, 1), gt.fock_state_vec(1, 30))]
    for chi, psi in cases:
        for phi in (0.0, math.pi / 2, 1.3):
            assert gt.ramsey_population_exact(br, bg, phi, psi) == pytest.approx(
                gt.ramsey_population(chi, br, bg, phi), abs=1e-10)


def test_ramsey_recovers_chi_from_two_populations():
    g = 0.9 - 0.4j
    chi = lambda b: gt.chi_coherent(b, g)               # noqa: E731
    br, bg = 0.15 + 0.05j, 0.45 + 0.25j
    p0 = gt.ramsey_population(chi, br, bg, 0.0)
    p90 = gt.ramsey_population(chi, br, bg, math.pi / 2)
    assert gt.ramsey_chi_from_populations(p0, p90, br, bg) == pytest.approx(chi(bg - br), abs=1e-9)


def test_ramsey_population_valid_and_fringe_contrast_is_chi():
    chi = lambda b: gt.chi_coherent(b, 0.7j)            # noqa: E731
    br, bg = 0.1 + 0j, 0.3 + 0.2j
    db = bg - br
    Ps = [gt.ramsey_population(chi, br, bg, i * 2 * math.pi / 64) for i in range(64)]
    assert all(-1e-12 <= p <= 1 + 1e-12 for p in Ps)   # valid probabilities
    # peak-to-trough fringe (at the extremal phases) = |chi(Delta beta)|
    phi_pk = -(bg * br.conjugate()).imag - cmath.phase(chi(db))
    p_pk = gt.ramsey_population(chi, br, bg, phi_pk)
    p_tr = gt.ramsey_population(chi, br, bg, phi_pk + math.pi)
    assert p_pk - p_tr == pytest.approx(abs(chi(db)), rel=1e-9)


def test_ramsey_disk_reach_is_two_eta():
    eta = 0.389
    vals = [abs(gt.delta_beta(eta, pr, pg))
            for pr in [i * math.pi / 12 for i in range(24)]
            for pg in [i * math.pi / 12 for i in range(24)]]
    assert all(v <= 2 * eta + 1e-9 for v in vals)      # confined to the disk
    assert max(vals) == pytest.approx(2 * eta, rel=1e-6)
    assert abs(gt.delta_beta(eta, 0.0, math.pi)) == pytest.approx(2 * eta)  # diameter
