"""
Characteristic-function / Wigner transfer-function kernels for the stroboscopic phase
grating -- the analytics of docs/notes/strobo_grating_transfer_function.md, packaged so
they can be evaluated and NUMERICALLY SELF-CHECKED against the full Floquet propagator
(engines.strobo_sim).

Conventions (cyclic frequencies in Hz, matching the note and strobo_sim):

    beta_k = i*eta*exp(i(phi_g - k*Phi)),   Phi = 2*pi*f_lf*Delta_t,
    theta  = 2*pi*Omega_strobo[Hz]*delta_t[s]   (the per-cycle pulse area).

Two transfer functions on the SAME amplitude operator A = -i(theta/2) sum_k
exp(-i2pi k delta Delta_t) D(beta_k):

  * PROBABILITY  P_down = <A^dag A>  -- a quadratic DOUBLE sum sampling chi at the kick
    differences; motion-blind on the exact strobe (the Dirichlet comb).
  * COHERENCE    Tr[rho A]           -- a single-sum discrete Fourier transform of chi.

Plus the exact eta=0 kicked-two-level result (valid at the N*theta~pi operating point),
analytic chi(beta)=Tr[rho D(beta)] and Wigner W(alpha) for standard states, and a direct
(pure-Python) chi->W transform for reconstruction. No numpy.
"""
from __future__ import annotations

import cmath
import math


# ---------- generalized Laguerre L_n(x) (Fock chi / Wigner) ------------------
def laguerre(n: int, x: float) -> float:
    if n <= 0:
        return 1.0
    lkm1, lk = 1.0, 1.0 - x
    for k in range(1, n):
        lkm1, lk = lk, ((2 * k + 1 - x) * lk - k * lkm1) / (k + 1)
    return lk


# ---------- characteristic functions  chi(beta) = Tr[rho D(beta)] ------------
def chi_vacuum(beta: complex) -> complex:
    return cmath.exp(-abs(beta) ** 2 / 2.0)


def chi_coherent(beta: complex, gamma: complex) -> complex:
    return cmath.exp(-abs(beta) ** 2 / 2.0 + beta * gamma.conjugate() - beta.conjugate() * gamma)


def chi_fock(beta: complex, n: int) -> complex:
    return cmath.exp(-abs(beta) ** 2 / 2.0) * laguerre(n, abs(beta) ** 2)


def chi_thermal(beta: complex, nbar: float) -> complex:
    return cmath.exp(-(nbar + 0.5) * abs(beta) ** 2)


def _braket_D(a: complex, beta: complex, b: complex) -> complex:
    """<a| D(beta) |b> for coherent states a, b (used to build the cat)."""
    return cmath.exp(1j * (beta * b.conjugate()).imag) * cmath.exp(
        -abs(a) ** 2 / 2.0 - abs(beta + b) ** 2 / 2.0 + a.conjugate() * (beta + b))


def chi_cat(beta: complex, gamma: complex, parity: int = +1) -> complex:
    """Even (parity=+1) / odd (-1) cat |psi> ~ |gamma> + parity|-gamma>."""
    g, p = gamma, parity
    norm = 2.0 + 2.0 * p * math.exp(-2.0 * abs(g) ** 2)        # <psi|psi> (unnormalised: 1/norm)
    val = (_braket_D(g, beta, g) + _braket_D(-g, beta, -g)
           + p * _braket_D(-g, beta, g) + p * _braket_D(g, beta, -g))
    return val / norm


# ---------- analytic Wigner functions  W(alpha)  (for validation) ------------
def wigner_vacuum(alpha: complex) -> float:
    return (2.0 / math.pi) * math.exp(-2.0 * abs(alpha) ** 2)


def wigner_coherent(alpha: complex, gamma: complex) -> float:
    return (2.0 / math.pi) * math.exp(-2.0 * abs(alpha - gamma) ** 2)


def wigner_fock(alpha: complex, n: int) -> float:
    return (2.0 / math.pi) * (-1) ** n * laguerre(n, 4.0 * abs(alpha) ** 2) * math.exp(-2.0 * abs(alpha) ** 2)


def wigner_thermal(alpha: complex, nbar: float) -> float:
    s = 2.0 * nbar + 1.0
    return (2.0 / (math.pi * s)) * math.exp(-2.0 * abs(alpha) ** 2 / s)


def wigner_cat(alpha: complex, gamma: complex, parity: int = +1) -> float:
    """Even/odd cat Wigner: two Gaussians + the interference fringe."""
    p = parity
    norm = math.pi * (1.0 + p * math.exp(-2.0 * abs(gamma) ** 2))
    return (math.exp(-2.0 * abs(alpha - gamma) ** 2) + math.exp(-2.0 * abs(alpha + gamma) ** 2)
            + 2.0 * p * math.exp(-2.0 * abs(alpha) ** 2) * math.cos(4.0 * (alpha * gamma.conjugate()).imag)
            ) / norm


# ---------- grating geometry -------------------------------------------------
def theta_of(omega_strobo_hz: float, delta_t_us: float) -> float:
    """Per-cycle pulse area theta = 2*pi*Omega_strobo[Hz]*delta_t[s]."""
    return 2.0 * math.pi * omega_strobo_hz * delta_t_us * 1e-6


def beta_samples(eta: float, phi_g: float, n_cycles: int, DELTA_t_us: float, f_lf_hz: float):
    """The displacements sampled by the grating, beta_k = i*eta*exp(i(phi_g - k*Phi))."""
    Phi = 2.0 * math.pi * f_lf_hz * DELTA_t_us * 1e-6
    return [1j * eta * cmath.exp(1j * (phi_g - k * Phi)) for k in range(n_cycles)]


# ---------- the two transfer-function kernels --------------------------------
def kernel_probability(chi, eta, phi_g, det_hz, n_cycles, DELTA_t_us, f_lf_hz, theta):
    """Weak-pulse spin-flip PROBABILITY P_down(phi_g, delta) = (theta/2)^2 sum_{k,k'}
    e^{i2pi(k-k')delta Dt} e^{i Im(b_k' b_k*)} chi(b_k' - b_k). `chi` is a callable."""
    b = beta_samples(eta, phi_g, n_cycles, DELTA_t_us, f_lf_hz)
    dt = DELTA_t_us * 1e-6
    s = 0j
    for k in range(n_cycles):
        bk = b[k]
        for kp in range(n_cycles):
            ph = 2.0 * math.pi * (k - kp) * det_hz * dt + (b[kp] * bk.conjugate()).imag
            s += cmath.exp(1j * ph) * chi(b[kp] - bk)
    return (theta / 2.0) ** 2 * s.real


def kernel_coherence(chi, eta, phi_g, det_hz, n_cycles, DELTA_t_us, f_lf_hz, theta):
    """Weak-pulse COHERENCE Tr[rho A] = -i(theta/2) sum_k e^{-i2pi k delta Dt} chi(b_k)
    -- a discrete Fourier transform of chi over the sampled ring. Returns a complex number
    (Re,Im map to <sigma_x>,<sigma_y> up to the reference coefficient)."""
    b = beta_samples(eta, phi_g, n_cycles, DELTA_t_us, f_lf_hz)
    dt = DELTA_t_us * 1e-6
    s = sum(cmath.exp(-1j * 2.0 * math.pi * k * det_hz * dt) * chi(b[k]) for k in range(n_cycles))
    return -0.5j * theta * s


def exact_eta0_probability(det_hz, theta, n_cycles, DELTA_t_us):
    """EXACT (all orders in theta) eta=0 kicked-two-level transfer function:
    P = sin^2(theta/2)/sin^2(lam) * sin^2(N lam), cos(lam)=cos(theta/2)cos(pi delta Dt).
    cos(lam) = (1/2) Tr U_cycle with U_cycle = U_free R_x(theta)."""
    cl = math.cos(theta / 2.0) * math.cos(math.pi * det_hz * DELTA_t_us * 1e-6)
    cl = max(-1.0, min(1.0, cl))
    lam = math.acos(cl)
    sl = math.sin(lam)
    if abs(sl) < 1e-15:
        return (n_cycles * math.sin(theta / 2.0)) ** 2          # on-tooth limit
    return math.sin(theta / 2.0) ** 2 / sl ** 2 * math.sin(n_cycles * lam) ** 2


# ---------- chi -> W reconstruction (direct 2-D transform) -------------------
def wigner_from_samples(beta_pts, chi_vals, alpha_pts, dbeta_area):
    """Reconstruct W on `alpha_pts` from chi sampled at `beta_pts` (a list of complex),
    W(alpha) = (1/pi^2) sum_beta chi(beta) e^{alpha beta* - alpha* beta} * dbeta_area.
    This is the 2-D Fourier inversion that turns measured chi(beta) into the Wigner
    function (real part returned)."""
    out = []
    pref = dbeta_area / math.pi ** 2
    for a in alpha_pts:
        acc = 0j
        for b, c in zip(beta_pts, chi_vals):
            acc += c * cmath.exp(a * b.conjugate() - a.conjugate() * b)
        out.append((pref * acc).real)
    return out


def square_grid(half_width, n):
    """n x n square grid of complex points on [-half_width, half_width]^2 (flattened),
    returns (points, spacing)."""
    step = 2.0 * half_width / (n - 1)
    pts = [complex(-half_width + i * step, -half_width + j * step)
           for j in range(n) for i in range(n)]
    return pts, step


# ---------- numerical self-checks vs the full Floquet propagator -------------
def check_probability_vs_floquet(eta=0.389, omega_strobo_hz=4.99e5, delta_t_us=0.001,
                                 N=50, DELTA_t_us=0.769172, f_lf_factor=1.03, det_hz=0.0):
    """Relative error of kernel_probability (weak pulse, state |0>) vs strobo_sim."""
    from .strobo_sim import strobo_detuning_scan
    f_lf = f_lf_factor / (DELTA_t_us * 1e-6)
    theta = theta_of(omega_strobo_hz, delta_t_us)
    p_ana = kernel_probability(chi_vacuum, eta, 0.0, det_hz, N, DELTA_t_us, f_lf, theta)
    p_num = strobo_detuning_scan(eta, omega_strobo_hz, delta_t_us, DELTA_t_us, N, f_lf, [det_hz], F=16)[0]
    return p_ana, p_num, abs(p_ana - p_num) / p_num if p_num else float("nan")


def check_exact_eta0_vs_floquet(omega_strobo_hz=4.99e5, delta_t_us=0.02, N=50,
                                DELTA_t_us=0.769172, det_factor=0.3):
    """Max abs error of exact_eta0_probability vs strobo_sim (eta=0)."""
    from .strobo_sim import strobo_detuning_scan
    f_lf = 1.0 / (DELTA_t_us * 1e-6)
    theta = theta_of(omega_strobo_hz, delta_t_us)
    det = det_factor * f_lf
    p_ana = exact_eta0_probability(det, theta, N, DELTA_t_us)
    p_num = strobo_detuning_scan(0.0, omega_strobo_hz, delta_t_us, DELTA_t_us, N, f_lf, [det], F=4)[0]
    return p_ana, p_num, abs(p_ana - p_num)
