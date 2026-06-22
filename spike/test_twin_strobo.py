"""
Tests for the stroboscopic OC carrier-flop twin (phase-grating n=0 baseline): the
per-cycle Rabi extraction (freq/N), the N=50 / motional-period strobo structure, and
the stroboscopic dephasing-decoupling (the flop keeps far more contrast over the fixed
38.5-us train than the continuous carrier flop's T_phi would allow).

Run:  pytest spike/
"""
import pytest

from spike import twin_strobo


_HAS = twin_strobo._DATAFILE.exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="Strobo2.0 .dat not present")


@pytest.fixture(scope="module")
def info():
    return twin_strobo.analyze(n_boot=60)


def test_strobo_structure(info):
    assert info["N"] == pytest.approx(50)
    assert info["DELTA_t"] == pytest.approx(0.769, abs=0.01)     # = 1/(1.30 MHz) motional period
    assert info["det_khz"] == pytest.approx(-40, abs=5)          # fr_oc_strobo near the carrier
    assert info["total_us"] == pytest.approx(38.5, rel=0.05)     # fixed N*DELTA_t


def test_per_cycle_rabi(info):
    # Omega_strobo = fit freq / N; the N=50 amplification gives a fast (~25 cyc/us) flop
    assert 300e3 < info["omega_strobo"] < 700e3
    assert info["omega_strobo_err"] >= 0
    assert info["turns"] > 1.0                                   # several flop turns over the scan
    assert info["omega0"] > info["omega_strobo"]                 # bare rate > Debye-Waller-reduced


def test_stroboscopic_dephasing_decoupling(info):
    assert 0.0 < info["contrast"] < 1.0
    # the strobo keeps MUCH more contrast than a continuous flop at T_phi~15us over 38.5us
    assert info["contrast"] > 3.0 * info["cont_contrast"]
    assert info["tphi_eff_us"] > 15.0                            # effective T_phi longer -> decoupled


def test_report_mentions_grating_and_decoupling(info):
    txt = twin_strobo.report(info)
    assert "phase-grating" in txt or "phase grating" in txt
    assert "DECOUPLES" in txt or "decoupl" in txt.lower()
    assert "n=0" in txt and "displ" in txt
