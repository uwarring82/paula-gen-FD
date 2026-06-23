"""
Tests for the end-to-end Wigner-tomography twin (spike/twin_wigner_tomography.py): the
raw-data round-trip, the noise-free reconstruction (exact pipeline), and recovery of the
displaced state's amplitude and phase from shot-noisy simulated data.

Run:  pytest spike/
"""
import cmath
import math

import pytest

from spike import twin_wigner_tomography as tw
from spike.engines import grating_tomography as gt

GAMMA = 1.3 * cmath.exp(1j * math.pi / 4.0)


def test_raw_file_roundtrip(tmp_path):
    raw, h = tw.simulate(GAMMA, beta_n=11, seed=1)
    p = tw.write_raw(raw, GAMMA, h, path=tmp_path / "raw.dat")
    rows, h_read = tw.read_raw(p)
    assert h_read == pytest.approx(h, abs=1e-4)
    assert len(rows) == len(raw)
    for (b0, nx0, ny0), (b1, M, nx1, ny1) in zip(raw, rows):
        assert b1 == pytest.approx(b0, abs=1e-4)
        assert (nx1, ny1, M) == (nx0, ny0, tw.SHOTS)


def test_noise_free_reconstruction_recovers_gamma():
    # Bypass shot noise: build chi_hat directly from the analytic chi, reconstruct W, and
    # check the centroid <alpha> = gamma to high precision (grid-robust).
    pts, h = gt.square_grid(4.0, 25)
    beta = [b for b in pts if abs(b) <= 4.0]
    chi = [gt.chi_coherent(b, GAMMA) for b in beta]
    step = 2 * 2.5 / 40
    ax = [-2.5 + i * step for i in range(41)]
    alpha = [complex(x, y) for y in ax for x in ax]
    W = gt.wigner_from_samples(beta, chi, alpha, h * h)
    g = tw.recovered_gamma(W, alpha)
    assert abs(g - GAMMA) < 0.05                       # centroid ~ gamma
    W_ana = [gt.wigner_coherent(a, GAMMA) for a in alpha]
    assert max(abs(r - a) for r, a in zip(W, W_ana)) < 0.02   # exact pipeline, fine grid


def test_full_pipeline_recovers_amplitude_and_phase():
    # Seeded shot-noisy run -> recover |gamma| and arg(gamma) within the shot-noise band.
    raw, h = tw.simulate(GAMMA, beta_n=21, seed=tw.SEED)
    rows, _ = (raw_rows(raw), h)                        # use raw directly (skip file I/O)
    _chi, _b, W, alpha, _axl, _s = tw.reconstruct(rows, tw.CONTRAST, alpha_n=41, h=h)
    g = tw.recovered_gamma(W, alpha)
    assert abs(g) == pytest.approx(1.3, abs=0.1)        # amplitude
    assert cmath.phase(g) == pytest.approx(math.pi / 4, abs=0.1)   # phase
    assert abs(g - GAMMA) < 0.1


def raw_rows(raw):
    """(beta, n_x, n_y) -> (beta, M, n_x, n_y) rows, as read_raw would return."""
    return [(b, tw.SHOTS, nx, ny) for (b, nx, ny) in raw]


def test_chi_recovery_is_unbiased_on_average():
    # The reconstructed chi_hat is an UNBIASED estimator of chi: the SIGNED mean residual
    # averages down over the grid (~ sigma/sqrt(N)), confirming the (2n/M-1)/C calibration.
    raw, h = tw.simulate(GAMMA, beta_n=21, seed=3)
    rows = raw_rows(raw)
    chi_hat, beta, _W, _a, _x, _s = tw.reconstruct(rows, tw.CONTRAST, alpha_n=11, h=h)
    res = [ch - gt.chi_coherent(b, GAMMA) for ch, b in zip(chi_hat, beta)]
    assert abs(sum(res) / len(res)) < 0.01             # unbiased: signed mean ~ 0
    assert sum(abs(r) for r in res) / len(res) < 0.09  # per-point spread ~ the shot level
