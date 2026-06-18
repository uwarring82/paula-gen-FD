"""
Drive engine — relative magnetic-dipole (microwave) Rabi couplings of the
25Mg+ ground-state hyperfine transitions |F, mF> <-> |F', mF'>.

By the Wigner-Eckart theorem the Rabi rate of a magnetic-dipole transition is
``Omega ∝ |<F' mF'| T^1_q |F mF>| = |<F mF; 1 q | F' mF'>| * <F'||T^1||F>`` with
``q = mF' - mF`` the photon polarization (sigma-: -1, pi: 0, sigma+: +1). Within
one F<->F' manifold the reduced matrix element is common, so the engine predicts
the *relative* couplings — the ATOMIC (Clebsch-Gordan) part only.

IMPORTANT: the measured absolute Rabi rates also fold in the microwave antenna's
polarization geometry and its frequency response, which DOMINATE the variation
across the ~26 MHz manifold (see docs/LOGBOOK.md). So this engine reproduces the
atomic structure, not the apparatus; ``measured / coupling`` is the apparatus
factor (a diagnostic, not a clean benchmark).

Pure Python (a Racah-formula Clebsch-Gordan), to stay dependency-light.
"""
from __future__ import annotations

import math


def _f(n):
    return math.factorial(int(round(n)))


def clebsch_gordan(j1, m1, j2, m2, j3, m3) -> float:
    """<j1 m1; j2 m2 | j3 m3> via the Racah formula (integer-step arguments)."""
    if abs((m1 + m2) - m3) > 1e-9:
        return 0.0
    if j3 < abs(j1 - j2) - 1e-9 or j3 > j1 + j2 + 1e-9:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    tri = _f(j1 + j2 - j3) * _f(j1 - j2 + j3) * _f(-j1 + j2 + j3) / _f(j1 + j2 + j3 + 1)
    pre = math.sqrt((2 * j3 + 1) * tri)
    pre *= math.sqrt(
        _f(j3 + m3) * _f(j3 - m3) * _f(j1 - m1) * _f(j1 + m1) * _f(j2 - m2) * _f(j2 + m2)
    )
    kmin = int(round(max(0.0, j2 - j3 - m1, j1 + m2 - j3)))
    kmax = int(round(min(j1 + j2 - j3, j1 - m1, j2 + m2)))
    s = 0.0
    for k in range(kmin, kmax + 1):
        s += (-1) ** k / (
            _f(k) * _f(j1 + j2 - j3 - k) * _f(j1 - m1 - k) * _f(j2 + m2 - k)
            * _f(j3 - j2 + m1 + k) * _f(j3 - j1 - m2 + k)
        )
    return pre * s


_POL = {1: "sigma+", 0: "pi", -1: "sigma-"}


class HyperfineDrive:
    """Relative microwave Rabi couplings for |F_a, mF> <-> |F_b, mF'> magnetic-
    dipole transitions (F_a, F_b are the two ground-state hyperfine manifolds)."""

    def __init__(self, F_a: float = 3.0, F_b: float = 2.0):
        self.F_a = float(F_a)
        self.F_b = float(F_b)

    def coupling(self, mF_a: float, mF_b: float) -> float:
        """|<F_a mF_a; 1 q | F_b mF_b>|, q = mF_b - mF_a — the relative Rabi
        amplitude (up to the common reduced matrix element)."""
        q = mF_b - mF_a
        return abs(clebsch_gordan(self.F_a, mF_a, 1, q, self.F_b, mF_b))

    def polarization(self, mF_a: float, mF_b: float) -> str:
        return _POL.get(int(round(mF_b - mF_a)), f"q={mF_b - mF_a:+g}")

    def relative_couplings(self, transitions) -> dict:
        """Couplings for a list of (mF_a, mF_b) transitions, normalised to the
        strongest in the set."""
        c = {t: self.coupling(*t) for t in transitions}
        peak = max(c.values()) or 1.0
        return {t: v / peak for t, v in c.items()}

    # --- absolute rate = atomic |CG| x apparatus calibration ---------------- #
    # The absolute Rabi rate factorises as Omega = |CG| * apparatus_factor, where
    # apparatus_factor = (MW drive amplitude) x (antenna polarization gain) x
    # (frequency response). The atomic |CG| is predicted above; the apparatus
    # factor must come from CALIBRATION.
    #
    # NOTE (verified against Doerr 2024 + Kaufmann 2022): a smooth, predictive
    # apparatus model is NOT extractable from the theses. The antenna response is
    # asymmetric, peaked (~1786 MHz), and polarization-dependent; a 5-parameter
    # physical model (3 polarization gains + a quadratic frequency response) fits
    # the 8 Doerr rates to only ~25% RMS, and the detailed characterization lives
    # in Doerr's figures (3.2-3.7), not tables. The apparatus also changed between
    # generations (Doerr is ~2-3x faster than Kaufmann, non-uniformly). So
    # apparatus_factor is an empirical PER-TRANSITION calibration, not a model.
    def apparatus_factor(self, mF_a: float, mF_b: float, measured_rabi_hz: float) -> float:
        """The empirical apparatus factor = measured_rabi / |CG| (calibration)."""
        cg = self.coupling(mF_a, mF_b)
        return measured_rabi_hz / cg if cg else math.nan

    def absolute_rabi(self, mF_a: float, mF_b: float, apparatus_factor: float) -> float:
        """Absolute Rabi rate [Hz] = |CG| * apparatus_factor, with the apparatus
        factor supplied from calibration (see the note above)."""
        return self.coupling(mF_a, mF_b) * apparatus_factor
