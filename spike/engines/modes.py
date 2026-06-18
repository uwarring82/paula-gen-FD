"""
Modes engine — axial normal modes of a linear chain of N equal-mass ions in a
harmonic trap (motion subsystem).

Given the single-ion axial secular frequency ``omega_z`` (the centre-of-mass
mode), the N axial normal modes are ``omega_p = sqrt(lambda_p) * omega_z`` where
``lambda_p`` are the eigenvalues of the dimensionless axial Hessian evaluated at
the ion equilibrium positions (James, Appl. Phys. B 66, 181 (1998)). The first
two eigenvalues are universal: ``lambda_1 = 1`` (COM) and ``lambda_2 = 3``
(stretch/breathing), so for N=2 the stretch mode is ``sqrt(3)*omega_z``
(analytically exact; computed numerically here to ~1e-15 relative).

Lengths are in units of the characteristic scale
``l = (e^2 / (4 pi eps0 M omega_z^2))^(1/3)``; only the dimensionless ratios
``omega_p / omega_z`` are computed, so the scale never needs evaluating.

Pure Python (uses spike.linalg), to stay numpy-free like the levels engine.
"""
from __future__ import annotations

import math

from ..linalg import eigvalsh, solve


def equilibrium_positions(N: int, tol: float = 1e-12, max_iter: int = 200):
    """Dimensionless axial equilibrium positions u_1 < ... < u_N (centred at 0).
    Solves the force balance u_m = sum_{n<m} 1/(u_m-u_n)^2 - sum_{n>m} 1/(u_m-u_n)^2
    by Newton's method (Jacobian = the axial Hessian). Convergence is judged on
    the max-component |g| (N-independent) and RAISES if not reached."""
    if N < 1:
        raise ValueError(f"need N >= 1 ions, got {N}")
    if N == 1:
        return [0.0]
    # initial guess: evenly spaced, centred
    u = [float(i) - (N - 1) / 2.0 for i in range(N)]
    converged = False
    for _ in range(max_iter):
        g = _grad(u)
        if max(abs(gi) for gi in g) < tol:
            converged = True
            break
        du = solve(axial_hessian(u), g)
        u = [u[i] - du[i] for i in range(N)]
    if not converged:
        raise RuntimeError(f"equilibrium_positions(N={N}) did not converge (|g|_inf > {tol})")
    return u


def _grad(u):
    """g_m = dV/du_m = u_m - sum_{n!=m} sign(u_m-u_n)/(u_m-u_n)^2."""
    N = len(u)
    g = []
    for m in range(N):
        s = u[m]
        for n in range(N):
            if n == m:
                continue
            d = u[m] - u[n]
            s -= math.copysign(1.0, d) / (d * d)
        g.append(s)
    return g


def axial_hessian(u):
    """Dimensionless axial Hessian A: A_mm = 1 + sum_{n!=m} 2/|u_m-u_n|^3,
    A_mn = -2/|u_m-u_n|^3 (m!=n)."""
    N = len(u)
    A = [[0.0] * N for _ in range(N)]
    for m in range(N):
        diag = 1.0
        for n in range(N):
            if n == m:
                continue
            inv3 = 2.0 / abs(u[m] - u[n]) ** 3
            diag += inv3
            A[m][n] = -inv3
        A[m][m] = diag
    return A


def axial_mode_eigenvalues(N: int):
    """The dimensionless eigenvalues lambda_p (ascending); omega_p = sqrt(lambda_p)*omega_z.
    The harmonic+Coulomb model gives lambda >= 1 for all modes; a roundoff-negative
    value is clamped to 0, but a significantly negative one is raised as a sign of
    an unstable/ill-conditioned configuration (before it reaches sqrt())."""
    out = []
    for lam in eigvalsh(axial_hessian(equilibrium_positions(N))):
        if lam < -1e-9:
            raise ValueError(
                f"negative axial eigenvalue {lam:.3e}: unstable/ill-conditioned configuration"
            )
        out.append(max(lam, 0.0))
    return out


def axial_mode_frequencies(omega_z: float, N: int):
    """The N axial normal-mode frequencies [Hz], ascending (omega_z is the COM)."""
    return [omega_z * math.sqrt(lam) for lam in axial_mode_eigenvalues(N)]


class AxialModes:
    """Axial normal modes of an N-ion chain, parameterised by the single-ion
    (COM) axial secular frequency."""

    def __init__(self, omega_z_com: float):
        self.omega_z = float(omega_z_com)

    @classmethod
    def from_ledger(cls, ledger, com_name: str = "omega_z_axial_com_25mg"):
        """Build from the COM axial frequency `input` (the wall enforces kind:input)."""
        return cls(omega_z_com=ledger.input_quantity(com_name).value)

    def mode_frequencies(self, N: int):
        return axial_mode_frequencies(self.omega_z, N)

    def stretch_frequency(self, N: int = 2) -> float:
        """The first non-COM (stretch/breathing) axial mode [Hz]. For N=2 this is
        sqrt(3)*omega_z (analytically exact; numerically to ~1e-15 relative)."""
        if N < 2:
            raise ValueError("the stretch mode requires N >= 2 ions")
        return self.mode_frequencies(N)[1]


# --------------------------------------------------------------------------- #
# Radial (transverse) modes                                                   #
# --------------------------------------------------------------------------- #
# Transverse displacements couple to Coulomb with HALF the magnitude and the
# OPPOSITE sign of the axial direction, so the radial modes soften below the
# single-ion radial frequency omega_r. The COM (all ions in phase) sits at
# omega_r; lower modes are reduced, and if omega_r is too small the lowest mode
# goes imaginary (the linear chain buckles to a zigzag).
def transverse_hessian(u, alpha: float):
    """Dimensionless transverse Hessian B (alpha = (omega_r/omega_z)^2):
    B_mm = alpha - sum_{n!=m} 1/|u_m-u_n|^3,  B_mn = +1/|u_m-u_n|^3 (m!=n)."""
    N = len(u)
    B = [[0.0] * N for _ in range(N)]
    for m in range(N):
        diag = alpha
        for n in range(N):
            if n == m:
                continue
            inv3 = 1.0 / abs(u[m] - u[n]) ** 3
            diag -= inv3
            B[m][n] = inv3
        B[m][m] = diag
    return B


def radial_mode_eigenvalues(N: int, alpha: float):
    """Transverse eigenvalues mu_p (ascending); omega_radial,p = sqrt(mu_p)*omega_z.
    A negative mu means the linear chain is unstable (zigzag) — raised loudly."""
    out = []
    for mu in eigvalsh(transverse_hessian(equilibrium_positions(N), alpha)):
        if mu < -1e-9:
            raise ValueError(
                f"radial mode unstable (mu={mu:.3e} < 0): the linear chain buckles "
                "to a zigzag; omega_r is too small relative to omega_z"
            )
        out.append(max(mu, 0.0))
    return out


def radial_mode_frequencies(omega_z: float, omega_r: float, N: int):
    """The N radial normal-mode frequencies [Hz], ascending (the COM, = omega_r,
    is the highest; lower modes are Coulomb-softened)."""
    alpha = (omega_r / omega_z) ** 2
    return [omega_z * math.sqrt(mu) for mu in radial_mode_eigenvalues(N, alpha)]


class RadialModes:
    """Radial (transverse) normal modes of an N-ion chain, parameterised by the
    single-ion axial (omega_z) and radial (omega_r) secular frequencies."""

    def __init__(self, omega_z: float, omega_r: float):
        self.omega_z = float(omega_z)
        self.omega_r = float(omega_r)

    @classmethod
    def from_ledger(cls, ledger, axial_com_name: str = "omega_z_axial_com_25mg",
                    radial_com_name: str = "omega_radial_com_25mg"):
        """Build from the axial + radial COM `input` records (wall-enforced)."""
        return cls(omega_z=ledger.input_quantity(axial_com_name).value,
                   omega_r=ledger.input_quantity(radial_com_name).value)

    def mode_frequencies(self, N: int):
        return radial_mode_frequencies(self.omega_z, self.omega_r, N)

    def com_frequency(self) -> float:
        """The radial COM mode [Hz] = omega_r (the highest radial mode)."""
        return self.omega_r

    def rocking_frequency(self, N: int = 2) -> float:
        """The lowest (out-of-phase / rocking) radial mode [Hz]. For N=2 this is
        sqrt(omega_r^2 - omega_z^2)."""
        if N < 2:
            raise ValueError("the rocking mode requires N >= 2 ions")
        return self.mode_frequencies(N)[0]
