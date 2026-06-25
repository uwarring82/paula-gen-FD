"""
Finite acoustic-transit model of an AOM-switched (Raman) drive.

Canonical, testable implementation behind
[docs/notes/aom_finite_sound_velocity_rabi.md] and the figures in
docs/figures/make_aom_rise_figs.py. Pure Python (``math`` only, no numpy) to match the
other spike engines and to keep the dependency surface small.

Physics. A single rising acoustic edge sweeps the Gaussian beam, so the diffracted FIELD
rises as an error function. For a nominal rectangular RF gate [0, dt] the delivered
single-pass field envelope (normalised to 1 at full diffraction) is a rounded box

    a(t) = 1/2 [ erf(t/tau_f) - erf((t-dt)/tau_f) ],     tau_f = w/V = D/(2 V)

with w the beam 1/e^2 intensity radius (= D/2). The two-photon Rabi rate carries one
field factor per SWITCHED diffraction event:  Omega(t) ~ a(t)^n, with
    n = 1  current scheme (only the single-pass R2 is chopped; B1 left on),
    n = 2  e.g. both beams single-pass and switched,
    n = 3  alternative (R2 single-pass + B1 double-pass, both switched).

Device default: IntraAction ASM-2202B3 (registry key ``intraaction_asm2202b3``),
fused silica, acoustic velocity V = 5.95 mm/us.

UNITS. The model is unit-agnostic in time: t, dt and tau_f must share ONE time unit, and
the comb detuning ``delta`` must be the reciprocal of that unit (MHz with us, GHz with ns).
``tau_f(D, V)`` returns us for D in mm and V in mm/us.
"""
from __future__ import annotations

from math import erf, exp, pi, sin, sqrt

V_FUSED_SILICA = 5.95  # mm/us (ASM-2202B3 datasheet)


def tau_f(D_mm: float, V: float = V_FUSED_SILICA) -> float:
    """Field 1/e rise time constant [us] for beam 1/e^2 diameter ``D_mm`` [mm]; w = D/2."""
    return (D_mm / 2.0) / V


def field_envelope(t: float, dt: float, tf: float) -> float:
    """Single-pass diffracted field a(t) for a gate [0, dt] (t, dt, tf in the same unit)."""
    return 0.5 * (erf(t / tf) - erf((t - dt) / tf))


def peak(dt: float, tf: float, n: int = 1) -> float:
    """Peak of a(t)^n = erf(dt/2 tau_f)^n  (reached at t = dt/2)."""
    return erf(dt / (2.0 * tf)) ** n


def _integral_pow(dt: float, tf: float, n: int, span: float = 16.0, npts: int = 8001) -> float:
    """Composite-trapezoid integral of a(t)^n over a wide window [us·(field units)]."""
    lo, hi = -span * tf, dt + span * tf
    h = (hi - lo) / (npts - 1)
    s = 0.0
    for i in range(npts):
        a = field_envelope(lo + i * h, dt, tf)
        s += (0.5 if (i == 0 or i == npts - 1) else 1.0) * a ** n
    return s * h


def area_ratio(dt: float, tf: float, n: int = 1) -> float:
    """Rotation-area efficiency  R_area = (1/dt) integral a^n dt, relative to an ideal
    instantaneous-switch rectangle. For n=1 this is **exactly 1** (the area theorem:
    the slow turn-on is cancelled by the equally slow turn-off tail)."""
    if n == 1:
        return 1.0
    return _integral_pow(dt, tf, n) / dt


def w_equiv(dt: float, tf: float, n: int = 1) -> float:
    """Equivalent-rectangle width  (integral a^n)/max(a^n) [same unit as dt].
    For n=1: dt/erf(dt/2 tau_f), which SATURATES at the floor sqrt(pi) tau_f as dt->0."""
    return area_ratio(dt, tf, n) * dt / peak(dt, tf, n)


def w_equiv_floor(tf: float, n: int = 1) -> float:
    """Small-dt floor of the equivalent width: sqrt(pi/n) * tau_f."""
    return sqrt(pi / n) * tf


def fwhm(dt: float, tf: float, n: int = 1, span: float = 16.0, npts: int = 24001) -> float:
    """Full width at half maximum of a(t)^n [same unit as dt]."""
    lo, hi = -span * tf, dt + span * tf
    h = (hi - lo) / (npts - 1)
    pk = peak(dt, tf, n)
    half = pk / 2.0
    first = last = None
    for i in range(npts):
        if field_envelope(lo + i * h, dt, tf) ** n >= half:
            t = lo + i * h
            if first is None:
                first = t
            last = t
    return (last - first) if first is not None else 0.0


def comb_envelope(delta: float, dt: float, tf: float, aom: bool = True) -> float:
    """Single-pulse spectral envelope |a~(delta)|^2 that weights the detuning-comb teeth.

    Weak-pulse limit: P_down(delta) ~ |array factor(delta)|^2 * comb_envelope(delta).
    The array factor (teeth at delta = k f_lf, width 1/(N Dt)) is unchanged by the AOM;
    the AOM enters ONLY here, multiplying the ideal sinc^2 by a Gaussian roll-off:

        |a~|^2 = sinc^2(pi delta dt) * [ exp(-2 pi^2 tau_f^2 delta^2)  if aom ]

    ``delta`` reciprocal to the time unit of dt, tf (MHz with us)."""
    x = pi * delta * dt
    sinc2 = 1.0 if x == 0.0 else (sin(x) / x) ** 2
    g = exp(-2.0 * (pi * tf * delta) ** 2) if aom else 1.0
    return sinc2 * g


def comb_halfwidth(dt: float, tf: float, aom: bool = True, fmax_units: float = 200.0,
                   npts: int = 40001) -> float:
    """1/e half-width of the comb envelope [reciprocal of the dt/tf time unit].
    Without the AOM this scales ~1/dt; with it, it SATURATES near 1/(sqrt(2) pi tau_f)."""
    h = fmax_units / (npts - 1)
    thr = exp(-1.0)
    for i in range(npts):
        d = i * h
        if comb_envelope(d, dt, tf, aom) < thr:
            return d
    return fmax_units
