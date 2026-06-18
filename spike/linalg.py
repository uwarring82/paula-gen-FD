"""
Tiny pure-Python linear algebra for the modes engine (keeps the spike numpy-free,
like the levels engine). Matrices are lists of lists of floats.

- solve(A, b): dense linear solve via Gaussian elimination with partial pivoting.
- eigvalsh(A): eigenvalues of a real SYMMETRIC matrix via the cyclic Jacobi
  algorithm, returned ascending.

These target the small (N up to ~30 ions) symmetric systems the normal-mode
engine produces; they are not meant as a general LA library.
"""
from __future__ import annotations

import math


def solve(A, b):
    """Solve A x = b for a square A (Gaussian elimination, partial pivoting)."""
    n = len(A)
    if any(len(row) != n for row in A) or len(b) != n:
        raise ValueError("solve() requires a square matrix A and a matching b")
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            raise ValueError("singular matrix in solve()")
        M[col], M[piv] = M[piv], M[col]
        pivot = M[col][col]
        for r in range(n):
            if r != col and M[r][col] != 0.0:
                f = M[r][col] / pivot
                for c in range(col, n + 1):
                    M[r][c] -= f * M[col][c]
    return [M[i][n] / M[i][i] for i in range(n)]


def eigvalsh(A, tol: float = 1e-14, max_sweeps: int = 200):
    """Ascending eigenvalues of a real symmetric matrix (cyclic Jacobi).

    Convergence is judged on the off-diagonal norm RELATIVE to the matrix scale,
    and the routine RAISES if it does not converge within ``max_sweeps`` (rather
    than silently returning an unconverged diagonal). Intended for the small,
    well-conditioned Hessians the modes engine produces; it is not robust to
    extreme dynamic range (near-coincident ions), which the engine never makes.
    """
    n = len(A)
    a = [list(row) for row in A]
    scale = math.sqrt(sum(a[i][i] ** 2 for i in range(n))) or 1.0
    converged = False
    for _ in range(max_sweeps):
        off = math.sqrt(sum(a[p][q] ** 2 for p in range(n) for q in range(p + 1, n)))
        if off <= tol * scale:
            converged = True
            break
        for p in range(n):
            for q in range(p + 1, n):
                apq = a[p][q]
                if apq == 0.0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * apq)
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                app, aqq = a[p][p], a[q][q]
                a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
                a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
                a[p][q] = a[q][p] = 0.0
                for i in range(n):
                    if i != p and i != q:
                        aip, aiq = a[i][p], a[i][q]
                        a[i][p] = a[p][i] = c * aip - s * aiq
                        a[i][q] = a[q][i] = s * aip + c * aiq
    if not converged:
        raise RuntimeError(f"eigvalsh did not converge in {max_sweeps} sweeps")
    return sorted(a[i][i] for i in range(n))
