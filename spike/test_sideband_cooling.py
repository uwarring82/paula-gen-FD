"""
Tests for the resolved-sideband-cooling diagnostic engine.

Run:  pytest spike/
"""
import pytest

from spike.engines.sideband_cooling import (
    SidebandCooling,
    kappa_from_nbar,
    mean_from_ground_state_prob,
    offres_carrier_excitation,
    rsb_cooling_limit,
)
from spike.ledger import Ledger


def test_cooling_limit_scaling_and_inversion():
    assert rsb_cooling_limit(2e6, 0.7e6) < rsb_cooling_limit(1e6, 0.7e6)   # higher omega -> colder
    nb = rsb_cooling_limit(1.3e6, 0.7e6)
    assert kappa_from_nbar(nb, 1.3e6) == pytest.approx(0.7e6, rel=1e-9)    # inverts the floor


def test_offres_carrier_excitation():
    assert offres_carrier_excitation(0.5e6, 1.3e6) == pytest.approx((0.5 / (2 * 1.3)) ** 2, rel=1e-9)


def test_ground_state_prob_consistency():
    assert mean_from_ground_state_prob(0.94) == pytest.approx(1 / 0.94 - 1, abs=1e-9)   # ~0.064


def test_from_ledger_resolved_regime():
    sc = SidebandCooling.from_ledger(Ledger.load())
    assert len(sc.modes) == 3
    for _lab, _omega, _nbar, _kap, kappa_over_omega in sc.inferred_kappa():
        assert 0.0 < kappa_over_omega < 1.0                   # all in the resolved-sideband regime
