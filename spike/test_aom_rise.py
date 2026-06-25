"""
Tests for the AOM finite-acoustic-transit model (spike/engines/aom_rise.py): the
error-function envelope, the area theorem, the peak, the effective-width floor, and the
detuning-comb spectral envelope. Backs docs/notes/aom_finite_sound_velocity_rabi.md.

Run:  pytest spike/
"""
import math

from spike.engines.aom_rise import (
    V_FUSED_SILICA,
    area_ratio,
    comb_envelope,
    comb_halfwidth,
    field_envelope,
    fwhm,
    peak,
    tau_f,
    w_equiv,
    w_equiv_floor,
)

_TF = tau_f(1.0)               # us, D = 1 mm -> 0.084 us


def test_tau_f_matches_datasheet():
    # tau_f = w/V = D/(2V); for D=1 mm, V=5.95 mm/us -> 0.0840 us
    assert math.isclose(tau_f(1.0), 0.5 / V_FUSED_SILICA, rel_tol=1e-12)
    assert math.isclose(tau_f(1.0) * 1e3, 84.03, abs_tol=0.1)
    # spec 10-90% INTENSITY rise: a(t)^2; cross-check ~0.75 D/V (~125 ns at D=1 mm)
    # (datasheet quotes 0.64 D/V = 110 ns; the ~15% is the coherent-vs-geometric convention)


def test_area_theorem_single_pass():
    # n=1: the rotation area is EXACTLY preserved (= dt) for any gate width.
    for dt in (0.005, 0.06, 0.5, 5.0):           # us
        assert area_ratio(dt, _TF, 1) == 1.0
    # and the numerical integral agrees with the analytic theorem
    from spike.engines.aom_rise import _integral_pow
    for dt in (0.01, 0.06, 0.3):
        assert math.isclose(_integral_pow(dt, _TF, 1), dt, rel_tol=2e-4)


def test_peak_is_erf():
    for dt in (0.01, 0.06, 0.2, 1.0):
        assert math.isclose(peak(dt, _TF, 1), math.erf(dt / (2 * _TF)), rel_tol=1e-12)
    # multi-pass peak is the single-pass peak to the n-th power
    assert math.isclose(peak(0.06, _TF, 3), peak(0.06, _TF, 1) ** 3, rel_tol=1e-12)


def test_field_envelope_peaks_at_centre_and_is_bounded():
    dt = 0.06
    centre = field_envelope(dt / 2, dt, _TF)
    assert math.isclose(centre, peak(dt, _TF, 1), rel_tol=1e-9)
    # always within [0, 1]
    for k in range(-50, 51):
        a = field_envelope(dt / 2 + k * 0.01, dt, _TF)
        assert -1e-12 <= a <= 1.0 + 1e-12


def test_effective_width_floor():
    # w_eff = dt/erf(dt/2 tau_f) for n=1, saturating at sqrt(pi) tau_f as dt -> 0
    assert math.isclose(w_equiv(0.06, _TF, 1), 0.06 / math.erf(0.06 / (2 * _TF)), rel_tol=1e-9)
    floor1 = w_equiv_floor(_TF, 1)
    assert math.isclose(floor1, math.sqrt(math.pi) * _TF, rel_tol=1e-12)
    assert math.isclose(w_equiv(1e-4, _TF, 1), floor1, rel_tol=1e-3)   # tiny gate -> floor
    # n=3 floor is sqrt(pi/3) tau_f (sharper) and below the n=1 floor
    floor3 = w_equiv_floor(_TF, 3)
    assert math.isclose(floor3, math.sqrt(math.pi / 3) * _TF, rel_tol=1e-12)
    assert floor3 < floor1
    assert math.isclose(w_equiv(1e-4, _TF, 3), floor3, rel_tol=3e-3)


def test_fwhm_floor_and_growth():
    # short gate -> FWHM near the kernel floor; long gate -> FWHM ~ dt
    assert math.isclose(fwhm(2.0, _TF, 1), 2.0, rel_tol=2e-3)
    assert 0.12 < fwhm(0.005, _TF, 1) < 0.16        # ~140 ns floor (in us)


def test_area_collapse_double_pass_regression():
    # regression against the note table (D=1 mm)
    assert math.isclose(area_ratio(0.06, _TF, 3), 0.086, abs_tol=0.002)
    assert math.isclose(area_ratio(0.06, _TF, 2), 0.273, abs_tol=0.003)
    # each extra switched factor costs ~one power of (dt/tau_f) at short dt: R2 < R3-ratio etc.
    assert area_ratio(0.03, _TF, 3) < area_ratio(0.03, _TF, 2) < area_ratio(0.03, _TF, 1)


def test_comb_envelope():
    # carrier (delta=0) is full; AOM only ever suppresses; matches sinc^2 * gaussian
    assert comb_envelope(0.0, 0.06, _TF, aom=True) == 1.0
    for d in (0.5, 1.3, 2.6, 5.2):                  # MHz
        e_ideal = comb_envelope(d, 0.06, _TF, aom=False)
        e_aom = comb_envelope(d, 0.06, _TF, aom=True)
        assert e_aom < e_ideal
        x = math.pi * d * 0.06
        ref = (math.sin(x) / x) ** 2 * math.exp(-2 * (math.pi * _TF * d) ** 2)
        assert math.isclose(e_aom, ref, rel_tol=1e-9)


def test_comb_width_saturates_with_aom():
    f_lf = 1.3001
    # WITHOUT the effect the comb keeps broadening as 1/dt: 10 ns wider than 100 ns
    assert comb_halfwidth(0.010, _TF, aom=False) > 3 * comb_halfwidth(0.100, _TF, aom=False)
    # WITH the effect the width saturates: 10 ns ~ 30 ns, both near the ~2.1 f_lf floor
    w10 = comb_halfwidth(0.010, _TF, aom=True) / f_lf
    w30 = comb_halfwidth(0.030, _TF, aom=True) / f_lf
    assert math.isclose(w10, w30, rel_tol=0.05)
    assert 1.8 < w10 < 2.4
