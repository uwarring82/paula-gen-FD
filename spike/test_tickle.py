"""
Tests for the tickle engine (Kalis 2016 sinc + Bessel motional-mode spectroscopy).

Run:  pytest spike/
"""
import glob
import math

import pytest

from spike.datfile import DatFile
from spike.engines.tickle import (
    _amp_shape,
    besselj,
    excitation_amplitude,
    fit_tickle,
    fluorescence,
)


def test_besselj_known_values():
    assert besselj(0, 0.0) == 1.0
    assert besselj(1, 0.0) == 0.0
    assert besselj(0, 1.0) == pytest.approx(0.7651976866, abs=1e-6)
    assert besselj(1, 1.0) == pytest.approx(0.4400505857, abs=1e-6)
    assert besselj(0, 15.0) == pytest.approx(-0.0142244728, abs=1e-5)     # large-x stability


def test_excitation_amplitude_sinc_shape():
    f0, texc = 1.3e6, 200e-6
    assert _amp_shape(f0, f0, texc) == pytest.approx(1.0)                 # max at resonance
    # the finite-pulse sinc has its first zero at Delta f = 1/texc (FWHM-ish), NOT a Lorentzian
    assert abs(_amp_shape(f0 + 1.0 / texc, f0, texc)) < 1e-6
    assert abs(_amp_shape(f0 + 0.5 / texc, f0, texc)) > 0.5               # still high half a lobe in


def test_fluorescence_dips_with_modulation():
    # On-resonance detection (Delta_BD=0): any sideband spread leaves the line peak -> DIP.
    f = lambda b: fluorescence(b, detuning_bd_hz=0.0, omega_hz=1.3e6, gamma_w_hz=41.8e6)
    assert f(0.0) == pytest.approx(1.0)                                   # no motion -> full fluorescence
    assert f(3.0) < 1.0 and f(6.0) < f(3.0)                               # more modulation -> deeper


def test_fit_recovers_injected_f0():
    f0, texc = 1.30e6, 200e-6
    freqs = [1.29e6 + 1e3 * k for k in range(21)]                         # +-10 kHz, 0.5 kHz steps
    counts = [3.0 - 2.0 * _amp_shape(f, f0, texc) ** 2 for f in freqs]    # a sinc^2 dip
    fit = fit_tickle(freqs, counts, texc_s=texc)
    assert fit["f0_hz"] == pytest.approx(f0, abs=300.0)                   # within 0.3 kHz
    assert fit["fwhm_hz"] == pytest.approx(1.0 / texc)                    # Fourier-limited width


def test_fit_real_lf_tickle_data():
    files = sorted(glob.glob("sources/data/Tickle/PDQ_LF_FScan/*.dat"))
    if not files:
        pytest.skip("tickle data not present")
    pts = DatFile(files[0]).point_shots()
    xs = [x * 1e6 for x, _ in pts]
    ys = [sum(c) / len(c) for _, c in pts]
    fit = fit_tickle(xs, ys, texc_s=200e-6)
    assert 1.297e6 < fit["f0_hz"] < 1.300e6              # axial mode ~1.299 MHz (ledger nominal 1.30)
    assert fit["resolved"]                               # a genuine, bracketed, significant dip


def test_significance_gate_rejects_miss_scan():
    # The MF directory mixes a NARROW calibration scan that misses the mode (not
    # resolved) with scans that bracket the real dip at ~3.224 MHz (resolved).
    files = sorted(glob.glob("sources/data/Tickle/PDQ_MF_FScan/*.dat"))
    if len(files) < 2:
        pytest.skip("MF tickle data not present")
    res = []
    for fn in files:
        pts = DatFile(fn).point_shots()
        order = sorted(range(len(pts)), key=lambda i: pts[i][0])
        xs = [pts[i][0] * 1e6 for i in order]
        ys = [sum(pts[i][1]) / len(pts[i][1]) for i in order]
        res.append(fit_tickle(xs, ys, texc_s=200e-6))
    resolved = [f["f0_hz"] / 1e6 for f in res if f["resolved"]]
    assert resolved                                      # at least one file resolves
    assert len(resolved) < len(files)                    # the narrow miss-scan is rejected
    assert all(3.20 < f < 3.25 for f in resolved)        # real mode ~3.224, not nominal 3.0/edge 3.124
