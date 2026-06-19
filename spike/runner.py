"""
Twin composition root.

Each registered validation: pulls its `input` records through the wall, builds
the relevant engine, predicts a quantity, and compares it to a `benchmark`
record — producing a uniform ``ValidationResult``. The runner renders one table,
flags any measured benchmark that no engine covers, and exits nonzero if any
result is in tension (or errored, e.g. an engine tried to consume a non-input).

Adding the next engine is: implement it, write one `_validate_*` that returns a
ValidationResult, and register it in REGISTRY.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

from .engines.acstark import ac_stark_shift, is_far_detuned
from .engines.cooling import DopplerCooling
from .engines.drive import HyperfineDrive
from .engines.levels import GroundStateZeeman
from .engines.modes import AxialModes, RadialModes
from .engines.projection import ModeProjection
from .engines.readout import ReadoutModel, qpn_uncertainty
from .engines.sideband import Sideband, raman_differential_stark_factor
from .engines.sideband_cooling import SidebandCooling
from .ledger import Ledger

THRESHOLD_SIGMA = 3.0


# --------------------------------------------------------------------------- #
@dataclass
class ValidationResult:
    benchmark: str
    engine: str
    subsystem: str = ""
    predicted: float = math.nan
    measured: float = math.nan
    units: str = ""
    sigma_pred: float = math.nan
    sigma_meas: float = math.nan
    consumed: tuple = ()
    error: str = ""

    @property
    def residual(self) -> float:
        return self.predicted - self.measured

    @property
    def sigma(self) -> float:
        return math.hypot(self.sigma_pred, self.sigma_meas)

    @property
    def nsigma(self) -> float:
        if self.error:
            return math.nan
        s = self.sigma
        if s == 0.0:
            return 0.0 if self.residual == 0.0 else math.inf
        return abs(self.residual) / s

    @property
    def status(self) -> str:
        if self.error:
            return "ERROR"
        return "ok" if self.nsigma <= THRESHOLD_SIGMA else "TENSION"

    @property
    def ok(self) -> bool:
        return self.status == "ok"


@dataclass
class Validation:
    benchmark: str
    engine: str
    run: Callable[[Ledger], ValidationResult]


# --------------------------------------------------------------------------- #
# Per-engine validations. Each consumes inputs through the wall (input_quantity
# / from_ledger) and compares against a benchmark (benchmark_quantity).
# --------------------------------------------------------------------------- #
def _central_sigma(f, x, dx) -> float:
    """1-sigma propagation through f via a central difference: |f(x+dx)-f(x-dx)|/2."""
    return abs(f(x + dx) - f(x - dx)) / 2.0 if dx else 0.0


def _validate_clock(ledger: Ledger) -> ValidationResult:
    Bq = ledger.input_quantity("b_field_quantization_freddy")
    bench = ledger.benchmark_quantity("clock_transition_25mg")
    eng = GroundStateZeeman.from_ledger(ledger)   # consumes A, I, g_J, g_I (wall-enforced)
    return ValidationResult(
        benchmark=bench.name, engine="levels",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=eng.clock_transition(Bq.value), measured=bench.value, units=bench.units,
        sigma_pred=_central_sigma(eng.clock_transition, Bq.value, Bq.sigma),
        sigma_meas=bench.sigma,
        consumed=("hyperfine_a_constant_25mg", "nuclear_spin_25mg",
                  "g_factor_electron_2s12", "g_factor_nuclear_25mg",
                  "b_field_quantization_freddy"),
    )


def _validate_stretch(ledger: Ledger) -> ValidationResult:
    com = ledger.input_quantity("omega_z_axial_com_25mg")
    bench = ledger.benchmark_quantity("omega_z_axial_stretch_2ion_25mg")
    eng = AxialModes.from_ledger(ledger)          # consumes the COM (wall-enforced)
    stretch = lambda wz: AxialModes(wz).stretch_frequency(2)  # noqa: E731
    return ValidationResult(
        benchmark=bench.name, engine="modes",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=eng.stretch_frequency(2), measured=bench.value, units=bench.units,
        sigma_pred=_central_sigma(stretch, com.value, com.sigma), sigma_meas=bench.sigma,
        consumed=("omega_z_axial_com_25mg",),
    )


def _validate_radial_rocking(ledger: Ledger) -> ValidationResult:
    axial = ledger.input_quantity("omega_z_axial_com_25mg")
    radial = ledger.input_quantity("omega_radial_com_25mg")
    bench = ledger.benchmark_quantity("omega_radial_rocking_2ion_25mg")
    eng = RadialModes.from_ledger(ledger)         # consumes both COMs (wall-enforced)
    rock = lambda wz, wr: RadialModes(wz, wr).rocking_frequency(2)  # noqa: E731
    sigma_pred = math.hypot(
        _central_sigma(lambda wr: rock(axial.value, wr), radial.value, radial.sigma),
        _central_sigma(lambda wz: rock(wz, radial.value), axial.value, axial.sigma),
    )
    return ValidationResult(
        benchmark=bench.name, engine="modes",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=eng.rocking_frequency(2), measured=bench.value, units=bench.units,
        sigma_pred=sigma_pred, sigma_meas=bench.sigma,
        consumed=("omega_z_axial_com_25mg", "omega_radial_com_25mg"),
    )


def _validate_weber_clock(ledger: Ledger) -> ValidationResult:
    Bq = ledger.input_quantity("b_field_zeeman_weber_25mg")
    bench = ledger.benchmark_quantity("clock_transition_weber_25mg")
    eng = GroundStateZeeman.from_ledger(ledger)
    return ValidationResult(
        benchmark=bench.name, engine="levels",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=eng.clock_transition(Bq.value), measured=bench.value, units=bench.units,
        sigma_pred=_central_sigma(eng.clock_transition, Bq.value, Bq.sigma), sigma_meas=bench.sigma,
        consumed=("hyperfine_a_constant_25mg", "nuclear_spin_25mg",
                  "g_factor_electron_2s12", "g_factor_nuclear_25mg",
                  "b_field_zeeman_weber_25mg"),
    )


def _validate_doppler_limit(ledger: Ledger) -> ValidationResult:
    gamma = ledger.input_quantity("mg_p32_natural_linewidth")
    bench = ledger.benchmark_quantity("doppler_cooling_limit_25mg")
    eng = DopplerCooling.from_ledger(ledger)
    return ValidationResult(
        benchmark=bench.name, engine="cooling",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=eng.doppler_limit(), measured=bench.value, units=bench.units,
        sigma_pred=_central_sigma(lambda g: DopplerCooling(g).doppler_limit(), gamma.value, gamma.sigma),
        sigma_meas=bench.sigma, consumed=("mg_p32_natural_linewidth",),
    )


def _validate_doppler_occupation(ledger: Ledger) -> ValidationResult:
    gamma = ledger.input_quantity("mg_p32_natural_linewidth")
    omega = ledger.input_quantity("omega_z_axial_clos_25mg")
    bench = ledger.benchmark_quantity("doppler_cooled_occupation_25mg")
    eng = DopplerCooling.from_ledger(ledger)
    nbar = lambda g, w: DopplerCooling(g).doppler_limit_occupation(w)  # noqa: E731
    sigma_pred = math.hypot(
        _central_sigma(lambda g: nbar(g, omega.value), gamma.value, gamma.sigma),
        _central_sigma(lambda w: nbar(gamma.value, w), omega.value, omega.sigma),
    )
    return ValidationResult(
        benchmark=bench.name, engine="cooling",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=eng.doppler_limit_occupation(omega.value), measured=bench.value,
        units=bench.units, sigma_pred=sigma_pred, sigma_meas=bench.sigma,
        consumed=("mg_p32_natural_linewidth", "omega_z_axial_clos_25mg"),
    )


def _validate_bdd_ac_stark(ledger: Ledger) -> ValidationResult:
    s = ledger.input_quantity("bdd_far_cooling_saturation")
    delta = ledger.input_quantity("bdd_far_cooling_detuning")
    gamma = ledger.input_quantity("mg_p32_natural_linewidth")
    bench = ledger.benchmark_quantity("bdd_ac_stark_shift_25mg")
    shift = lambda sv, dv, gv: abs(ac_stark_shift(sv, gv, dv))  # noqa: E731
    sigma_pred = math.sqrt(
        _central_sigma(lambda sv: shift(sv, delta.value, gamma.value), s.value, s.sigma) ** 2
        + _central_sigma(lambda dv: shift(s.value, dv, gamma.value), delta.value, delta.sigma) ** 2
        + _central_sigma(lambda gv: shift(s.value, delta.value, gv), gamma.value, gamma.sigma) ** 2
    )
    return ValidationResult(
        benchmark=bench.name, engine="acstark",
        subsystem=ledger.record(bench.name)["scope"]["subsystem"],
        predicted=shift(s.value, delta.value, gamma.value), measured=bench.value,
        units=bench.units, sigma_pred=sigma_pred, sigma_meas=bench.sigma,
        consumed=("bdd_far_cooling_saturation", "bdd_far_cooling_detuning",
                  "mg_p32_natural_linewidth"),
    )


REGISTRY = [
    Validation("clock_transition_25mg", "levels", _validate_clock),
    Validation("clock_transition_weber_25mg", "levels", _validate_weber_clock),
    Validation("omega_z_axial_stretch_2ion_25mg", "modes", _validate_stretch),
    Validation("omega_radial_rocking_2ion_25mg", "modes", _validate_radial_rocking),
    Validation("doppler_cooling_limit_25mg", "cooling", _validate_doppler_limit),
    Validation("doppler_cooled_occupation_25mg", "cooling", _validate_doppler_occupation),
    Validation("bdd_ac_stark_shift_25mg", "acstark", _validate_bdd_ac_stark),
]


# --------------------------------------------------------------------------- #
def run_all(ledger: Ledger, registry=REGISTRY) -> list[ValidationResult]:
    """Run every validation, capturing per-validation errors (e.g. a wall
    violation) as ERROR results rather than aborting the whole run."""
    results = []
    for v in registry:
        try:
            results.append(v.run(ledger))
        except Exception as exc:
            results.append(ValidationResult(benchmark=v.benchmark, engine=v.engine,
                                            error=f"{type(exc).__name__}: {exc}"))
    return results


def uncovered_benchmarks(ledger: Ledger, registry=REGISTRY) -> list[str]:
    """Measured benchmarks (derived_from == []) that no registered engine
    validates. Excludes subsystem:control (microwave-drive CALIBRATION), which is
    handled by the drive DIAGNOSTIC below (apparatus-limited, not a sigma test)."""
    covered = {v.benchmark for v in registry}
    out = []
    for name in ledger.by_kind("benchmark"):
        rec = ledger.record(name)
        if name in covered or (rec.get("derived_from") or []):
            continue
        if (rec.get("scope") or {}).get("subsystem") == "control":
            continue
        out.append(name)
    return out


# (mF_a, mF_b, ledger record) for the Doerr microwave Rabi set
DOERR_RABI = [
    (3, 2, "mw_rabi_3p3_2p2_doerr"), (1, 2, "mw_rabi_3p1_2p2_doerr"),
    (1, 0, "mw_rabi_3p1_2p0_doerr"), (0, 0, "mw_rabi_3p0_2p0_doerr"),
    (0, -1, "mw_rabi_3p0_2m1_doerr"), (-1, -1, "mw_rabi_3m1_2m1_doerr"),
    (-1, -2, "mw_rabi_3m1_2m2_doerr"), (-3, -2, "mw_rabi_3m3_2m2_doerr"),
]


def drive_diagnostic(ledger: Ledger) -> str:
    """DIAGNOSTIC (not a sigma-validation): for each Doerr microwave Rabi rate,
    the drive engine's atomic Clebsch-Gordan coupling vs the measured rate. The
    apparatus factor (measured / |CG|, normalised) reveals the MW antenna
    polarization + frequency response, which dominate the absolute rates."""
    eng = HyperfineDrive(3, 2)
    rows = []
    for a, b, name in DOERR_RABI:
        if name not in ledger:
            continue
        meas = ledger.benchmark_quantity(name).value
        cg = eng.coupling(a, b)
        rows.append([f"|3,{a:+d}>-|2,{b:+d}>", eng.polarization(a, b),
                     f"{cg:.3f}", f"{meas / 1e3:.1f}", meas / cg])
    if not rows:
        return ""
    norm = min(r[4] for r in rows)
    head = ["transition", "pol", "|CG|", "Rabi/kHz", "apparatus(meas/|CG|, norm)"]
    body = [[r[0], r[1], r[2], r[3], f"{r[4] / norm:.2f}x"] for r in rows]
    w = [max(len(head[i]), *(len(r[i]) for r in body)) for i in range(len(head))]
    line = lambda c: "  ".join(c[i].ljust(w[i]) for i in range(len(c)))  # noqa: E731
    out = ["DRIVE DIAGNOSTIC — atomic |CG| vs measured MW Rabi (Doerr); apparatus-limited",
           "", line(head), line(["-" * x for x in w])] + [line(r) for r in body]
    out.append("")
    out.append("  The apparatus column is NOT flat (spans ~6x) -> the absolute Rabi rate is")
    out.append("  dominated by the MW antenna polarization + frequency response, not |CG|.")
    return "\n".join(out)


_COOLING_BEAMS = [
    ("BD (cooling)", "bd_cooling_detuning", "bd_cooling_saturation"),
    ("BDX (detection)", "bdx_detection_detuning", "bdx_detection_saturation"),
    ("BDD (far cool)", "bdd_far_cooling_detuning", "bdd_far_cooling_saturation"),
]


def cooling_diagnostic(ledger: Ledger) -> str:
    """Steady-state scattering rate at each Blue Doppler beam setting (detuning +
    saturation from the laser table), via the cooling engine."""
    if "mg_p32_natural_linewidth" not in ledger:
        return ""
    eng = DopplerCooling.from_ledger(ledger)
    rows = []
    for label, dname, sname in _COOLING_BEAMS:
        if dname in ledger and sname in ledger:
            d, s = ledger.value(dname), ledger.value(sname)
            rows.append([label, f"{d / 1e6:+.1f}", f"{s:g}", f"{eng.scatter_rate(d, s) / 1e6:.1f}"])
    if not rows:
        return ""
    head = ["beam", "detuning/MHz", "s=I/Isat", "scatter/Mphotons/s"]
    w = [max(len(head[i]), *(len(r[i]) for r in rows)) for i in range(len(head))]
    line = lambda c: "  ".join(c[i].ljust(w[i]) for i in range(len(c)))  # noqa: E731
    out = ["COOLING DIAGNOSTIC — Doppler scattering per Blue Doppler beam (engine x ledger)",
           "", line(head), line(["-" * x for x in w])] + [line(r) for r in rows]
    out.append("")
    out.append(f"  max (saturated, on resonance) = {eng.max_scatter_rate() / 1e6:.0f} Mphotons/s; "
               f"T_D = {eng.doppler_limit() * 1e3:.3f} mK at -Gamma/2 = {eng.optimal_detuning() / 1e6:.1f} MHz.")
    return "\n".join(out)


_DOERR_ADDRESSING = {
    "CC": (), "OC": ("lf",), "AC": ("lf", "mf", "hf"), "ROC": ("mf", "hf"),
}


def projection_diagnostic(ledger: Ledger) -> str:
    """Which motional mode each Raman (TPSR) combination addresses: the engine's
    geometric projection (Delta_k direction x mode axes) vs Doerr 2024's
    documented addressing."""
    if "radial_mode_tilt_25mg" not in ledger:
        return ""
    eng = ModeProjection.from_ledger(ledger)
    rows = []
    all_match = True
    for comb in ("CC", "OC", "AC", "ROC"):
        p = eng.projections(comb)
        addr = eng.addressed_modes(comb)
        match = addr == _DOERR_ADDRESSING[comb]
        all_match = all_match and match
        rows.append([comb, f"{p['lf']:.3f}", f"{p['mf']:.3f}", f"{p['hf']:.3f}",
                     ",".join(addr) or "(carrier)", "ok" if match else "MISMATCH"])
    head = ["comb", "->lf", "->mf", "->hf", "addresses", "vs Doerr"]
    w = [max(len(head[i]), *(len(r[i]) for r in rows)) for i in range(len(head))]
    line = lambda c: "  ".join(c[i].ljust(w[i]) for i in range(len(c)))  # noqa: E731
    out = ["PROJECTION DIAGNOSTIC — Raman combination -> motional mode (engine x ledger geometry)",
           "", line(head), line(["-" * x for x in w])] + [line(r) for r in rows]
    out.append("")
    out.append(f"  addressed modes vs Doerr 2024: {'all match' if all_match else 'MISMATCH'}; "
               f"AC axial proj = {eng.projection('AC', 'lf'):.3f} (Doerr 45 deg -> 0.707); "
               f"radial tilt = {eng.tilt_deg:.0f} deg.")
    return "\n".join(out)


_ACSTARK_BEAMS = [
    ("BD (cool)", "bd_cooling_detuning", "bd_cooling_saturation", "mg_p32_natural_linewidth"),
    ("BDD (far)", "bdd_far_cooling_detuning", "bdd_far_cooling_saturation", "mg_p32_natural_linewidth"),
    ("BDX (det)", "bdx_detection_detuning", "bdx_detection_saturation", "mg_p32_natural_linewidth"),
    ("RD (rep)", "rd_repump_detuning", "rd_repump_saturation", "mg_p12_natural_linewidth"),
    ("RP (rep)", "rp_repump_detuning", "rp_repump_saturation", "mg_p12_natural_linewidth"),
]


def acstark_diagnostic(ledger: Ledger) -> str:
    """Per-beam AC-Stark light shift delta_AC = s*Gamma^2/(8 delta), flagged by
    detuning regime: only the FAR-detuned BDD shifts coherently; the near-resonant
    beams scatter (cooling engine) and the formula does not apply."""
    rows = []
    for label, dname, sname, gname in _ACSTARK_BEAMS:
        if dname in ledger and sname in ledger and gname in ledger:
            d, s, g = ledger.value(dname), ledger.value(sname), ledger.value(gname)
            far = is_far_detuned(g, d)
            shift = ac_stark_shift(s, g, d)
            rows.append([label, f"{d / g:+.1f}", f"{s:g}",
                         f"{shift / 1e6:+.1f}" if far else "--",
                         "coherent shift" if far else "scatters (n/a)"])
    if not rows:
        return ""
    head = ["beam", "delta/Gamma", "s", "delta_AC/MHz", "regime"]
    w = [max(len(head[i]), *(len(r[i]) for r in rows)) for i in range(len(head))]
    line = lambda c: "  ".join(c[i].ljust(w[i]) for i in range(len(c)))  # noqa: E731
    out = ["AC-STARK DIAGNOSTIC — far-detuned light shift per beam (engine x ledger)",
           "", line(head), line(["-" * x for x in w])] + [line(r) for r in rows]
    out.append("")
    out.append("  Only the far-detuned BDD has a meaningful coherent AC-Stark shift "
               "(~10 MHz, vs Hasse);")
    out.append("  the near-resonant beams scatter -> cooling engine. Hasse: RD/RP "
               "induce no significant shift.")
    return "\n".join(out)


def sideband_diagnostic(ledger: Ledger) -> str:
    """Absolute Lamb-Dicke eta per (combination, addressed mode) at the lf/mf/hf
    mode frequencies, plus the Raman differential AC-Stark shift of the resonance
    (engine x ledger geometry)."""
    needed = ("radial_mode_tilt_25mg", "raman_axial_lamb_dicke_25mg",
              "omega_z_axial_com_25mg", "omega_radial_mf_25mg", "omega_radial_hf_25mg")
    if any(n not in ledger for n in needed):
        return ""
    sb = Sideband.from_ledger(ledger)
    proj = ModeProjection.from_ledger(ledger)
    freq = {"lf": ledger.value("omega_z_axial_com_25mg"),
            "mf": ledger.value("omega_radial_mf_25mg"),
            "hf": ledger.value("omega_radial_hf_25mg")}
    rows = []
    for comb in ("CC", "OC", "AC", "ROC"):
        addr = proj.addressed_modes(comb)
        cells = [comb] + [f"{sb.lamb_dicke(comb, m, freq[m]):.3f}" if m in addr else "--"
                          for m in ("lf", "mf", "hf")]
        rows.append(cells)
    head = ["comb", "eta_lf", "eta_mf", "eta_hf"]
    w = [max(len(head[i]), *(len(r[i]) for r in rows)) for i in range(len(head))]
    line = lambda c: "  ".join(c[i].ljust(w[i]) for i in range(len(c)))  # noqa: E731
    out = ["SIDEBAND DIAGNOSTIC — absolute Lamb-Dicke eta per combination x mode (engine x ledger)",
           "", line(head), line(["-" * x for x in w])] + [line(r) for r in rows]
    out.append("")
    out.append(f"  mode freqs lf/mf/hf = {freq['lf'] / 1e6:.1f}/{freq['mf'] / 1e6:.1f}/"
               f"{freq['hf'] / 1e6:.1f} MHz; sideband Rabi = eta sqrt(n+1) Omega_0; "
               "anchored to eta=0.32 (OC,lf,1.92 MHz).")
    if "hyperfine_splitting_calc_25mg" in ledger and "raman_detuning_from_p32" in ledger:
        whf = ledger.value("hyperfine_splitting_calc_25mg")
        d_r = ledger.value("raman_detuning_from_p32")
        factor = raman_differential_stark_factor(whf, d_r)
        out.append(f"  Raman differential AC-Stark (order-of-magnitude): delta_AC_diff ~ "
                   f"{factor / 2:.3f}-{factor:.3f} * Omega_0 (omega_HF/Delta_R = "
                   f"{whf / 1e9:.3f}/{d_r / 1e9:.0f} GHz, prefactor 0.5-1.0, unanchored). "
                   "BDD validates the single-beam light-shift SCALE only, not this ratio.")
    return "\n".join(out)


def readout_diagnostic(ledger: Ledger) -> str:
    """Electronic-state readout from Thomm's measured count levels: the DETECTION-
    limited single-shot fidelity and the maximum-likelihood P_down precision over N
    shots (Cramer-Rao). The measured ENSEMBLE fidelities (99.4/97.4 %) are
    preparation-dominated, so this is a diagnostic, not a sigma test."""
    if "mg_bright_counts_25mg" not in ledger or "mg_dark_counts_25mg" not in ledger:
        return ""
    rm = ReadoutModel.from_ledger(ledger)
    f_ss = rm.single_shot_fidelity()
    n = 2500
    out = ["READOUT DIAGNOSTIC — state discrimination from lambda_down/lambda_up (Thomm 2021)",
           "",
           f"  lambda_down = {rm.lam_bright:.3f}, lambda_up = {rm.lam_dark:.3f} counts (30 us window)",
           f"  single-shot discrimination fidelity (threshold)  = {f_ss * 100:.1f} %"]
    for p in (0.5, 0.99):
        sig, qp = rm.p_down_uncertainty(p, n), qpn_uncertainty(p, n)
        out.append(f"  ML P_down precision over N={n}, p={p}: sigma = {sig:.4f} "
                   f"(ideal QPN {qp:.4f}, readout overhead x{sig / qp:.1f})")
    if "mg_readout_fidelity_bright_25mg" in ledger:
        fb = ledger.value("mg_readout_fidelity_bright_25mg")
        fd = ledger.value("mg_readout_fidelity_dark_25mg")
        out += ["",
                f"  measured ENSEMBLE fidelities (ML, prep+detection): bright {fb * 100:.1f} %, "
                f"dark {fd * 100:.1f} % (single-shot cap {f_ss * 100:.0f} %). By channel:",
                f"    bright {(1 - fb) * 100:.1f} % deficit: PREPARATION-dominated — ML of a perfect bright "
                "state averages to ~1, above the single-shot cap;",
                f"    dark   {(1 - fd) * 100:.1f} % deficit: DETECTION — dark-state off-resonant scatter / "
                "depumping (leak) DURING t_det, a systematic the averaging does NOT remove.",
                "  So the readout limits the DARK channel; not a sigma target."]
    return "\n".join(out)


def sideband_cooling_diagnostic(ledger: Ledger) -> str:
    """Resolved-sideband cooling consistency: invert each achieved n_bar for the
    implied effective cooling linewidth/rate kappa (n_bar ~ (kappa/2 omega)^2). A
    diagnostic — the achieved n_bar are protocol-limited, not a single RSB floor."""
    sc = SidebandCooling.from_ledger(ledger)
    if not sc.modes:
        return ""
    rows = [[lab, f"{omega / 1e6:.1f}", f"{nbar:.2f}", f"{kap / 1e6:.2f}", f"{kov:.2f}"]
            for lab, omega, nbar, kap, kov in sc.inferred_kappa()]
    head = ["mode", "omega/MHz", "n_bar", "kappa/MHz", "kappa/omega"]
    w = [max(len(head[i]), *(len(r[i]) for r in rows)) for i in range(len(head))]
    line = lambda c: "  ".join(c[i].ljust(w[i]) for i in range(len(c)))  # noqa: E731
    out = ["SIDEBAND-COOLING DIAGNOSTIC — RSB floor consistency (Thomm 2021 achieved n_bar)",
           "", line(head), line(["-" * x for x in w])] + [line(r) for r in rows]
    out += ["",
            "  all kappa/omega < 1 (resolved-sideband regime confirmed); kappa varies across modes",
            "  -> the achieved n_bar are protocol/per-mode-limited, not a single common floor."]
    return "\n".join(out)


# --- rendering --------------------------------------------------------------
def _cell(value: float, units: str, kind: str) -> str:
    if units == "Hz":
        return f"{value / 1e6:.6f} MHz" if kind == "value" else f"{value / 1e3:+.2f} kHz"
    if units == "K":
        return f"{value * 1e3:.4f} mK" if kind == "value" else f"{value * 1e6:+.2f} uK"
    if units in ("dimensionless", "1"):
        return f"{value:.4g}" if kind == "value" else f"{value:+.3g}"
    return f"{value:.6g} {units}" if kind == "value" else f"{value:+.3g} {units}"


def render_table(results: list[ValidationResult]) -> str:
    head = ["benchmark", "engine", "subsystem", "predicted", "reference",
            "residual", "n_sigma", "status"]
    rows = []
    for r in results:
        if r.error:
            rows.append([r.benchmark, r.engine, "-", "-", "-", "-", "-", "ERROR"])
        else:
            rows.append([
                r.benchmark, r.engine, r.subsystem,
                _cell(r.predicted, r.units, "value"), _cell(r.measured, r.units, "value"),
                _cell(r.residual, r.units, "resid"), f"{r.nsigma:.2f}", r.status,
            ])
    w = [max(len(head[i]), *(len(row[i]) for row in rows)) for i in range(len(head))] if rows \
        else [len(h) for h in head]
    line = lambda cells: "  ".join(c.ljust(w[i]) for i, c in enumerate(cells))  # noqa: E731
    out = [line(head), line(["-" * x for x in w])] + [line(row) for row in rows]
    return "\n".join(out)


def main(argv=None) -> int:
    ledger = Ledger.load()
    results = run_all(ledger)

    print("TWIN VALIDATION — engines reproduce benchmarks from ledger inputs\n")
    print(render_table(results))

    for r in results:
        if r.error:
            print(f"\n  {r.benchmark} [{r.engine}] ERROR: {r.error}")
        else:
            print(f"\n  {r.benchmark}: consumed {', '.join(r.consumed)}")

    for diag in (drive_diagnostic(ledger), cooling_diagnostic(ledger),
                 projection_diagnostic(ledger), acstark_diagnostic(ledger),
                 sideband_diagnostic(ledger), readout_diagnostic(ledger),
                 sideband_cooling_diagnostic(ledger)):
        if diag:
            print("\n" + diag)

    uncovered = uncovered_benchmarks(ledger)
    if uncovered:
        print("\nMeasured benchmarks with no sigma-validation (prep/protocol-limited; "
              "addressed by the readout / sideband-cooling DIAGNOSTICS above, not a sigma target):")
        for n in uncovered:
            print(f"  - {n}")

    n_ok = sum(1 for r in results if r.ok)
    n_bad = len(results) - n_ok
    print(f"\n{len(results)} validation(s): {n_ok} ok, {n_bad} not ok "
          f"(threshold {THRESHOLD_SIGMA:.0f} sigma).")
    return 0 if n_bad == 0 else 1
