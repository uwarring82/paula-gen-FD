"""
Tests for the sideband-thermometry discriminator: DatFile.counter_blocks returns both
sideband flops, the RSB/BSB ratio gives a COLD nbar, and the carrier decay decomposed
at that nbar is dominated by Raman-beam dephasing (resolving the carrier degeneracy).

Run:  pytest spike/
"""
import pytest

from spike import twin_sideband
from spike.datfile import DatFile

_HAS = twin_sideband._DATAFILE.exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="RSB+BSB sideband .dat not present")


@pytest.fixture(scope="module")
def info():
    return twin_sideband.analyze(n_boot=60)             # one analysis shared by all tests


def test_counter_blocks_returns_both_sidebands():
    dat = DatFile(twin_sideband._DATAFILE)
    blocks = dat.counter_blocks()
    assert len(blocks) == 2
    for x, y, s in blocks:
        assert len(x) == len(y) == len(s) == dat.scan["points"]
    # block 0 (BSB) has much larger swing than block 1 (RSB, near-constant)
    swing = lambda b: max(b[1]) - min(b[1])
    assert swing(blocks[0]) > 2.0 * swing(blocks[1])


def test_sideband_thermometry_is_cold(info):
    assert info["a_bsb"] > 3.0 * info["a_rsb"]          # BSB full, RSB suppressed
    assert info["ratio"] < 0.4
    assert 0.0 < info["nbar"] < 0.6                      # COLD (not the carrier's ~1)
    assert info["nbar_err"] > 0
    assert info["nbar"] < 0.5 * info["nbar_eff_carrier"]   # below the all-motional reading


def test_carrier_decay_dominated_by_raman_dephasing(info):
    # at the COLD measured nbar, the carrier motional floor is a minority of the decay
    assert info["g_floor"] < info["g_obs"]
    assert info["g_raman"] > info["g_floor"]            # Raman dephasing dominates
    assert info["dnu"] > 0 and info["tphi_us"] > 0
    txt = twin_sideband.report(info)
    assert "RAMAN-BEAM DEPHASING" in txt and "COLD" in txt
