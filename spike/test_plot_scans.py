"""
Tests for the (transition-generic) scan plotter: parsing + discovery + a smoke
render. The physics (generalized_rabi, qpn, fit_rabi) is tested elsewhere.

Run:  pytest spike/
"""
from spike.plot_scans import (
    _base_transition,
    _discover,
    _label,
    _levels_key,
    _secs,
    plot_pair,
)


def test_base_transition_strips_prefix_and_variant():
    assert _base_transition("fr_mw_3p1_2p2_s") == "3p1_2p2"
    assert _base_transition("t_mw_3p3_2p2") == "3p3_2p2"
    assert _base_transition("fr_mw_3m1_2m2") == "3m1_2m2"


def test_levels_key_maps_to_F3_F2_mf():
    assert _levels_key("3p3_2p2") == (3.0, 2.0)
    assert _levels_key("3p1_2p2") == (1.0, 2.0)
    assert _levels_key("3m1_2m1") == (-1.0, -1.0)
    assert _levels_key("3m3_2m2") == (-3.0, -2.0)


def test_label_has_both_kets():
    lab = _label("3p1_2p2")
    assert "3,+1" in lab and "2,+2" in lab


def test_secs_from_filename():
    assert _secs("09_40_09_12_06_2026.dat") == 9 * 3600 + 40 * 60 + 9


def test_discover_finds_both_mw_control_transitions():
    pairs = {p[0] for p in _discover("sources/data/MW_Control")}
    assert {"3p3_2p2", "3p1_2p2"} <= pairs


def test_plot_pair_smoke(tmp_path):
    base, fpath, dpath = next(p for p in _discover("sources/data/MW_Control") if p[0] == "3p1_2p2")
    out = tmp_path / "t.png"
    plot_pair(base, fpath, dpath, out, ledger=None)
    assert out.exists() and out.stat().st_size > 2000
