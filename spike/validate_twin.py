"""
Twin validation runner (spike).

Builds the levels engine from the ledger's `input` records, predicts the 25Mg+
clock transition at the freddy operating field, and compares it to the measured
`benchmark` record clock_transition_25mg — the input/benchmark wall, executed.

Run from the repo root (a source checkout):  python -m spike.validate_twin
"""
from __future__ import annotations

import math
import sys

from .constants import GAUSS, MHZ
from .engines.levels import GroundStateZeeman
from .engines.modes import AxialModes
from .ledger import Ledger


def _fmt_mhz(hz: float) -> str:
    return f"{hz / MHZ:.6f} MHz"


def check_clock(ledger: Ledger) -> bool:
    # every consumed quantity goes through the wall (kind:input); the comparison
    # target must be kind:benchmark
    A = ledger.input_quantity("hyperfine_a_constant_25mg")
    I = ledger.input_quantity("nuclear_spin_25mg")
    Bq = ledger.input_quantity("b_field_quantization_freddy")
    bench = ledger.benchmark_quantity("clock_transition_25mg")

    eng = GroundStateZeeman.from_ledger(ledger)
    predicted = eng.clock_transition(Bq.value)

    # propagate the dominant input uncertainty (the field) into the prediction
    dpred_dB = (
        eng.clock_transition(Bq.value + Bq.sigma) - eng.clock_transition(Bq.value - Bq.sigma)
    ) / 2.0
    sigma_pred = abs(dpred_dB)
    residual = predicted - bench.value
    sigma_res = (sigma_pred ** 2 + bench.sigma ** 2) ** 0.5
    if not math.isfinite(sigma_res):
        # a required uncertainty is missing/unparseable — don't silently call it tension
        raise ValueError(
            "cannot judge consistency: a required uncertainty is missing or non-finite "
            f"(sigma_pred={sigma_pred}, bench.sigma={bench.sigma})"
        )
    if sigma_res == 0.0:
        nsig = 0.0 if residual == 0.0 else float("inf")
    else:
        nsig = abs(residual) / sigma_res
    ok = nsig <= 3.0

    print("LEVELS ENGINE — 25Mg+ ground-state clock transition |F=3,mF=0> <-> |F=2,mF=0>")
    print("  inputs (consumed):")
    print(f"    A   = {A.value / MHZ:.6f} MHz        (hyperfine_a_constant_25mg)")
    print(f"    I   = {I.value}                       (nuclear_spin_25mg)")
    print(f"    B   = {Bq.value / GAUSS:.3f} +/- {Bq.sigma / GAUSS:.3f} G   (b_field_quantization_freddy)")
    print(f"  zero-field splitting 3|A|      = {_fmt_mhz(eng.zero_field_splitting())}")
    # K: Hz/T^2 -> Hz/G^2  (1 T^2 = 1e8 G^2)
    print(f"  quadratic Zeeman coefficient K = {eng.quadratic_zeeman_coeff() / 1e8:.1f} Hz/G^2")
    print(f"  predicted clock @ {Bq.value / GAUSS:.2f} G        = {_fmt_mhz(predicted)}  (+/- {sigma_pred:.0f} Hz)")
    print(f"  benchmark (measured, Doerr)     = {_fmt_mhz(bench.value)}  (+/- {bench.sigma:.0f} Hz)")
    print(f"  residual (pred - meas)          = {residual:+.0f} Hz  = {nsig:.2f} sigma")
    print(f"  --> {'CONSISTENT' if ok else 'TENSION'} within combined uncertainty\n")
    return ok


def check_modes(ledger: Ledger) -> bool:
    com = ledger.input_quantity("omega_z_axial_com_25mg")
    bench = ledger.benchmark_quantity("omega_z_axial_stretch_2ion_25mg")

    eng = AxialModes.from_ledger(ledger)
    predicted = eng.stretch_frequency(2)            # sqrt(3) * COM
    sigma_pred = abs(
        AxialModes(com.value + com.sigma).stretch_frequency(2)
        - AxialModes(com.value - com.sigma).stretch_frequency(2)
    ) / 2.0
    residual = predicted - bench.value
    sigma_res = (sigma_pred ** 2 + bench.sigma ** 2) ** 0.5
    nsig = abs(residual) / sigma_res if sigma_res else float("inf")
    ok = nsig <= 3.0

    print("MODES ENGINE — 2-ion 25Mg+ axial stretch mode (out-of-phase)")
    print("  inputs (consumed):")
    print(f"    omega_z (COM) = {com.value / MHZ:.4f} +/- {com.sigma / 1e3:.1f} kHz   (omega_z_axial_com_25mg)")
    print(f"  axial eigenvalue ratios (N=2)   = 1, 3  -> stretch = sqrt(3)*COM")
    print(f"  predicted stretch               = {predicted / MHZ:.4f} MHz  (+/- {sigma_pred / 1e3:.1f} kHz)")
    print(f"  benchmark (measured, Wittemer)  = {bench.value / MHZ:.4f} MHz  (+/- {bench.sigma / 1e3:.1f} kHz)")
    print(f"  residual (pred - meas)          = {residual / 1e3:+.1f} kHz  = {nsig:.2f} sigma")
    print(f"  --> {'CONSISTENT' if ok else 'TENSION'} within combined uncertainty\n")
    return ok


def main(argv=None) -> int:
    ledger = Ledger.load()
    results = [check_clock(ledger), check_modes(ledger)]
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
