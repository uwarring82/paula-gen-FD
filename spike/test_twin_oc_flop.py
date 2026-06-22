"""
Tests for the OC axial carrier-flop twin: it builds from the real .dat + ledger,
composes the four loss channels (coherent / AC-Stark / scattering / motional), and
the decomposition finds the ledger floor (scattering + thermal at the RSB-cooled
nbar) to be a SMALL part of the observed decay — with the effective-nbar inversion
recovering a near-unity nbar, far above the cooled benchmark.

Run:  pytest spike/
"""
import math
import random

import pytest

from spike.engines.scatter import RamanScatter
from spike.engines.sideband import Sideband
from spike.ledger import Ledger
from spike.twin_oc_flop import OCFlopTwin, build, report

_HAS_DATA = __import__("spike.twin_oc_flop", fromlist=["_DATAFILE"])._DATAFILE.exists()
pytestmark = pytest.mark.skipif(not _HAS_DATA, reason="OC carrier-flop .dat not present")


def _twin(rabi=2.0e5, nbar=0.07, **kw):
    ledger = Ledger.load()
    return OCFlopTwin(rabi_hz=rabi, scatter=RamanScatter.from_ledger(ledger),
                      sideband=Sideband.from_ledger(ledger),
                      omega_lf_hz=1.30e6, nbar=nbar, mu_bright=2.9, mu_dark=0.03, **kw)


def test_flip_prob_limits_coherent():
    tw = _twin(nbar=0.0)
    t_pi_us = 1e6 / (2.0 * tw.omega_eff)        # AC-Stark speed-up
    assert tw.flip_prob(0.0, with_scatter=False) == pytest.approx(0.0)
    # cold + no scatter + no AC-Stark detuning beyond the cap -> near-full pi flop.
    # eta reduces the rate (Debye-Waller), so use the engine's own peak time.
    assert tw.flip_prob(t_pi_us, with_scatter=False, with_motional=False) == pytest.approx(
        tw.amp_cap, rel=1e-6)


def test_warmer_nbar_decays_faster():
    tw = _twin()
    late = 8.0
    centre = 0.5
    p_cold = tw.flip_prob(late, nbar=0.0)
    p_warm = tw.flip_prob(late, nbar=1.0)
    # the warm flop has collapsed toward 1/2 more than the cold one
    assert abs(p_warm - centre) < abs(p_cold - centre)


def test_effective_decay_monotone_in_nbar():
    tw = _twin()
    g = [tw.effective_decay(nb) for nb in (0.0, 0.07, 0.3, 1.0)]
    assert g[0] == pytest.approx(0.0, abs=tw.gamma_sc_contrast + 1.0) or g[0] < g[1]
    assert g[1] < g[2] < g[3]                   # hotter -> faster effective decay


def test_invert_nbar_roundtrips():
    tw = _twin()
    target = tw.effective_decay(0.5)
    nbar_eff, saturated = tw.invert_nbar(target)
    assert not saturated
    assert nbar_eff == pytest.approx(0.5, abs=0.1)


def test_signal_curve_maps_to_count_levels():
    tw = _twin()
    ts = [0.2 * k for k in range(50)]
    sig = tw.signal_curve(ts, nbar=1.0)
    assert min(sig) >= tw.mu_dark - 1e-9
    assert max(sig) <= tw.mu_bright + 1e-9


def test_simulate_counts_nonneg_ints():
    tw = _twin()
    rng = random.Random(0)
    counts = tw.simulate_counts(2.4, rng, nbar=1.0)
    assert len(counts) == tw.n_shots
    assert all(isinstance(c, int) and c >= 0 for c in counts)
    assert sum(counts) / len(counts) > tw.mu_dark


def test_build_motional_dominates_scattering_but_floor_is_small():
    twin, fit, info = build()
    # at the RSB-cooled nbar, motional dephasing > scattering, but both are small
    assert info["g_thermal"] > info["g_sc"]
    assert info["g_floor"] < 0.25 * info["g_obs"]      # ledger floor < ~25% of observed
    assert info["eta"] == pytest.approx(0.389, rel=0.03)
    assert info["nbar_rsb"] == pytest.approx(0.07, rel=1e-6)


def test_build_nbar_eff_is_near_unity_and_hotter():
    twin, fit, info = build()
    assert info["nbar_eff"] > 10.0 * info["nbar_rsb"]  # much hotter than the cooled benchmark
    assert 0.5 < info["nbar_eff"] < 1.5                # ~ unity for this flop
    # the inversion reproduces the observed decay at nbar_eff
    assert twin.effective_decay(info["nbar_eff"]) == pytest.approx(info["g_obs"], rel=0.12)


def test_build_rabi_is_physical():
    twin, fit, info = build()
    assert 80e3 < info["rabi_hz"] < 260e3
    assert not fit["grid_edge"]


def test_report_mentions_inversion_and_channels():
    twin, fit, info = build()
    text = report(info)
    assert "motional" in text and "scattering" in text
    assert "nbar_eff" in text and "hotter" in text


# --- Raman-beam dephasing as the alternative residual explanation -----------
def test_gamma_raman_channel_adds_decay():
    tw0 = _twin()
    tw1 = _twin(gamma_raman_hz=5e4)
    late = 8.0
    centre = 0.5
    # the Raman-phase channel collapses the oscillation toward 1/2 faster
    assert abs(tw1.flip_prob(late, nbar=0.0) - centre) < abs(tw0.flip_prob(late, nbar=0.0) - centre)
    # and it raises the effective decay at fixed nbar
    assert tw1.effective_decay(0.07) > tw0.effective_decay(0.07)


def test_raman_explanation_reproduces_observed_decay():
    from spike.engines.raman_dephasing import mutual_linewidth_from_rate
    twin, fit, info = build()
    # at the COOLED nbar, adding gamma_raman = residual reproduces the observed decay
    tw_raman = _twin(rabi=info["rabi_hz"], gamma_raman_hz=info["g_resid"])
    assert tw_raman.effective_decay(info["nbar_rsb"]) == pytest.approx(info["g_obs"], rel=0.05)
    # the reported mutual linewidth is the inversion of the residual
    assert info["dnu_req"] == pytest.approx(mutual_linewidth_from_rate(info["g_resid"]))
    assert info["dnu_req"] > 0 and info["tphi_req"] > 0


def test_dual_explanation_in_report():
    twin, fit, info = build()
    text = report(info)
    assert "MOTIONAL" in text and "RAMAN-BEAM DEPHASING" in text
    assert "DEGENERATE" in text and "dnu" in text
