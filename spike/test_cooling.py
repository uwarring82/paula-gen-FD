"""
Tests for the cooling engine — two-level scattering and the Doppler limit.

Run:  pytest spike/
"""
import math

import pytest

from spike import constants as C
from spike.engines.cooling import (
    DopplerCooling,
    doppler_limit_temperature,
    mean_occupation,
    optimal_cooling_detuning,
    scatter_rate,
)
from spike.ledger import Ledger

GAMMA = 41.8e6   # 25Mg+ P3/2 linewidth /2pi (Hz)


def test_doppler_limit_is_about_1mK():
    assert doppler_limit_temperature(GAMMA) == pytest.approx(1.003e-3, abs=5e-6)
    # = h*gamma/(2 kB)
    assert doppler_limit_temperature(GAMMA) == pytest.approx(
        C.H_PLANCK * GAMMA / (2 * C.K_BOLTZMANN))


def test_optimal_detuning_is_minus_half_gamma():
    assert optimal_cooling_detuning(GAMMA) == pytest.approx(-GAMMA / 2)


def test_max_scatter_rate_is_gamma_over_two():
    # on resonance, fully saturated -> Gamma/2 (angular) = pi*gamma_hz
    assert DopplerCooling(GAMMA).max_scatter_rate() == pytest.approx(math.pi * GAMMA)


def test_scatter_rate_lorentzian_limits():
    # at resonance, s=1: R = (Gamma/2)/2
    r0 = scatter_rate(0.0, 1.0, GAMMA)
    assert r0 == pytest.approx(math.pi * GAMMA / 2)
    # far off resonance -> strongly suppressed
    far = scatter_rate(100 * GAMMA, 1.0, GAMMA)
    assert far < r0 / 1000
    # symmetric in detuning sign
    assert scatter_rate(GAMMA, 1.0, GAMMA) == pytest.approx(scatter_rate(-GAMMA, 1.0, GAMMA))


def test_scatter_rate_at_bd_setting():
    # BD: detuning -0.5*Gamma, s=0.5 -> ~26 Mphotons/s
    assert scatter_rate(-0.5 * GAMMA, 0.5, GAMMA) == pytest.approx(2.626e7, rel=1e-3)


def test_from_ledger_reproduces_doppler_limit_benchmark():
    ledger = Ledger.load()
    pred = DopplerCooling.from_ledger(ledger).doppler_limit()
    bench = ledger.benchmark_quantity("doppler_cooling_limit_25mg")
    assert abs(pred - bench.value) < 5e-5   # within 0.05 mK of the ~1 mK benchmark


def test_from_ledger_refuses_benchmark():
    ledger = Ledger.load()
    with pytest.raises(ValueError):
        DopplerCooling.from_ledger(ledger, gamma_name="doppler_cooling_limit_25mg")


# --- Doppler-cooled mean occupation (Clos measured benchmark) ---------------
def test_mean_occupation_bose_einstein_and_high_T_limit():
    # high-temperature limit: n ~ kB*T/(h*omega) - 1/2
    n = mean_occupation(1.0e6, 1.0e-3)
    approx = C.K_BOLTZMANN * 1.0e-3 / (C.H_PLANCK * 1.0e6) - 0.5
    assert n == pytest.approx(approx, rel=1e-3)


def test_doppler_limit_occupation_matches_clos():
    c = DopplerCooling(41.8e6)
    # Clos: n_bar ~ 10 at 2 MHz (his stated estimate); 10.4 at the 1.915 MHz measurement
    assert c.doppler_limit_occupation(2.0e6) == pytest.approx(9.96, abs=0.1)
    assert c.doppler_limit_occupation(1.915e6) == pytest.approx(10.42, abs=0.05)
    # closed form: depends only on omega/Gamma -> 1/(exp(2 omega/Gamma) - 1)
    assert c.doppler_limit_occupation(1.915e6) == pytest.approx(
        1.0 / math.expm1(2 * 1.915e6 / 41.8e6))


def test_from_ledger_reproduces_occupation_benchmark():
    ledger = Ledger.load()
    eng = DopplerCooling.from_ledger(ledger)
    omega = ledger.input_quantity("omega_z_axial_clos_25mg").value
    bench = ledger.benchmark_quantity("doppler_cooled_occupation_25mg")
    pred = eng.doppler_limit_occupation(omega)
    assert abs(pred - bench.value) < bench.sigma   # within the n_bar = 10(1) error bar
