"""
Tests for the .dat reader.

Run:  pytest spike/
"""
import pytest

from spike.datfile import DatFile

_DUR = "sources/data/microwave/13_28_39_15_06_2026.dat"
_FREQ = "sources/data/microwave/13_28_34_15_06_2026.dat"


def test_scan_definition_and_settings():
    d = DatFile(_DUR)
    assert d.scan["name"] == "t_mw_3p3_2p2"
    assert d.scan["points"] == 17 and d.scan["shots"] == 75
    assert d.scan["upper"] == pytest.approx(34.875)
    assert d.settings["EU_beam_fr"] == pytest.approx(1774.83)   # MW base frequency


def test_signal_is_the_flop_not_the_reference():
    d = DatFile(_DUR)
    x, y, s = d.signal()
    assert len(x) == 17 and x[0] == 0.0
    assert y[0] == pytest.approx(5.76, abs=0.05)                # bright at t=0
    assert min(y) < 0.6 and max(y) > 5.0                        # full flop, not the 0.013 reference


def test_histograms_mean_matches_signal():
    d = DatFile(_DUR)
    hists = d.histograms()
    # every histogram has 75 shots; the bright one reproduces the t=0 signal
    assert all(sum(h.values()) == 75 for h in hists)
    means = sorted(DatFile.hist_mean(h) for h in hists)
    assert means[0] == pytest.approx(0.013, abs=0.01)           # dark reference
    assert means[-1] == pytest.approx(5.8, abs=0.1)             # bright


def test_timestamp_from_filename():
    assert DatFile(_DUR).timestamp == "2026-06-15 13:28:39"
    assert DatFile(_FREQ).timestamp == "2026-06-15 13:28:34"


def test_freq_scan_resonance_dip():
    d = DatFile(_FREQ)
    x, y, _ = d.signal()
    f_res = min(zip(x, y), key=lambda p: p[1])[0]
    assert f_res == pytest.approx(1775.49, abs=0.02)
