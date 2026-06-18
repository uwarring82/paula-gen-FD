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

from .engines.levels import GroundStateZeeman
from .engines.modes import AxialModes, RadialModes
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


REGISTRY = [
    Validation("clock_transition_25mg", "levels", _validate_clock),
    Validation("clock_transition_weber_25mg", "levels", _validate_weber_clock),
    Validation("omega_z_axial_stretch_2ion_25mg", "modes", _validate_stretch),
    Validation("omega_radial_rocking_2ion_25mg", "modes", _validate_radial_rocking),
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
    """Measured benchmarks (derived_from == []) that no registered engine validates."""
    covered = {v.benchmark for v in registry}
    out = []
    for name in ledger.by_kind("benchmark"):
        rec = ledger.record(name)
        if name not in covered and not (rec.get("derived_from") or []):
            out.append(name)
    return out


# --- rendering --------------------------------------------------------------
def _cell(value: float, units: str, kind: str) -> str:
    if units == "Hz":
        return f"{value / 1e6:.6f}" if kind == "value" else f"{value / 1e3:+.2f}"
    return f"{value:.6g}" if kind == "value" else f"{value:+.3g}"


def render_table(results: list[ValidationResult]) -> str:
    head = ["benchmark", "engine", "subsystem", "predicted/MHz", "benchmark/MHz",
            "residual/kHz", "n_sigma", "status"]
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

    uncovered = uncovered_benchmarks(ledger)
    if uncovered:
        print("\nUNCOVERED measured benchmarks (no engine validates these):")
        for n in uncovered:
            print(f"  - {n}")

    n_ok = sum(1 for r in results if r.ok)
    n_bad = len(results) - n_ok
    print(f"\n{len(results)} validation(s): {n_ok} ok, {n_bad} not ok "
          f"(threshold {THRESHOLD_SIGMA:.0f} sigma).")
    return 0 if n_bad == 0 else 1
