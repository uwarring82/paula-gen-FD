"""
Tests for the readout engine — single-shot fidelity + ML P_down precision.

Run:  pytest spike/
"""
import pytest

from spike.engines.readout import (
    ReadoutModel,
    fisher_information_p_down,
    p_down_uncertainty,
    qpn_uncertainty,
    single_shot_fidelity,
)
from spike.ledger import Ledger


def test_single_shot_fidelity_separation_and_depump():
    assert single_shot_fidelity(20.0, 0.05) > 0.999             # well separated -> ~1
    f0 = single_shot_fidelity(2.682, 0.036)                     # Thomm levels -> ~0.95
    assert 0.93 < f0 < 0.96
    assert single_shot_fidelity(2.682, 0.036, depump_bright=0.3) < f0   # depumping degrades it


def test_fisher_information_limits():
    p = 0.5
    i_sep = fisher_information_p_down(p, 30.0, 0.01)            # well separated -> 1/(p(1-p))=4
    assert i_sep == pytest.approx(1.0 / (p * (1.0 - p)), rel=0.05)
    assert fisher_information_p_down(p, 2.682, 0.036) < i_sep   # overlap reduces information


def test_crb_at_least_qpn():
    for p in (0.3, 0.5, 0.9):                                   # readout overhead >= 1
        assert p_down_uncertainty(p, 1000, 2.682, 0.036) >= qpn_uncertainty(p, 1000) - 1e-12


def test_from_ledger():
    rm = ReadoutModel.from_ledger(Ledger.load())
    assert rm.lam_bright == pytest.approx(2.682)
    assert rm.lam_dark == pytest.approx(0.036)
    assert 0.93 < rm.single_shot_fidelity() < 0.96
