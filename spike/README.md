# Twin spike

A throwaway-ish, end-to-end proof that the structured parameter layer can drive
**physics engines**. Each engine consumes `input` records from the
[ledger](../records/) and is checked against `benchmark` records — the
input/benchmark wall, executed.

> **Boundary.** Per the task card, solver/physics code is *out of scope* for the
> substrate (`records/` + `validator/`). This package is kept separate and may
> later graduate to its own repo(s) — `iontrap-levels` / `-modes` / `-optics`
> plus a twin composition root. Pure Python, no numpy: the `2S_1/2` Breit-Rabi
> problem is closed-form, and the small normal-mode linear algebra is in
> [`linalg.py`](linalg.py).

## Composition root

[`runner.py`](runner.py) is the twin's composition root: it runs every
registered engine validation against its benchmark and renders one table.

```
$ python -m spike.validate_twin
TWIN VALIDATION — engines reproduce benchmarks from ledger inputs

benchmark                        engine  subsystem       predicted/MHz  measured/MHz  residual/kHz  n_sigma  status
clock_transition_25mg            levels  internal_state  1788.829549    1788.832200   -2.65         1.09     ok
omega_z_axial_stretch_2ion_25mg  modes   motion          2.251666       2.230000      +21.67        2.17     ok

2 validation(s): 2 ok, 0 not ok (threshold 3 sigma).
```

It exits nonzero if any result is in tension (> 3σ) or errors (e.g. an engine
tries to consume a non-input), and flags any measured benchmark that no engine
covers. **Adding an engine** is: implement it, write one `_validate_*` returning
a `ValidationResult`, and add it to `REGISTRY`.

## Engine: `levels` (Breit-Rabi)

[`engines/levels.py`](engines/levels.py) models the ²⁵Mg⁺ ground-state
(3s ²S₁/₂) hyperfine + Zeeman structure from the magnetic-dipole constant `A`,
nuclear spin `I`, electronic/nuclear g-factors, and the field `B`. Exact for
`J = 1/2` and **sign-safe**: `A < 0` (the inverted ²⁵Mg⁺ manifold, F=3 below F=2)
enters only through `A²` and diagonal `m_J m_I` terms. Reproduces the measured
`clock_transition_25mg` benchmark (residual −2.65 kHz = 1.09σ).

## Engine: `modes` (axial normal modes)

[`engines/modes.py`](engines/modes.py) computes the axial normal modes of an
N-ion chain from the COM secular frequency: `omega_p = sqrt(lambda_p) * omega_z`
with `lambda_p` the eigenvalues of the dimensionless axial Hessian at the ion
equilibrium positions (James 1998; `lambda_1 = 1` COM, `lambda_2 = 3` stretch).
Reproduces the 2-ion axial stretch (`omega_z_axial_stretch_2ion_25mg`):
predicted √3·COM = 2.2517 MHz vs measured 2.23 MHz, +21.7 kHz = 2.17σ (the √3
relation holds to ~1% in the 3-sig-fig data).

## Layout

```
spike/
  constants.py      CODATA constants + 25Mg atomic g-factors (sourced)
  ledger.py         loads records/*.yaml into a queryable Ledger (pyyaml only)
  linalg.py         tiny pure-Python solve() + eigvalsh() (Jacobi)
  engines/
    levels.py       2S_1/2 hyperfine+Zeeman engine (closed-form Breit-Rabi)
    modes.py        axial normal modes (equilibrium + Hessian eigenvalues)
  runner.py         composition root: registry, ValidationResult, table
  validate_twin.py  CLI shim -> runner.main
  test_levels.py    levels physics limits + benchmark + wall
  test_modes.py     modes eigenvalues + benchmark + wall + linalg robustness
  test_runner.py    runner: result math, tension detection, wall refusal, coverage
```

## Run

```bash
pytest spike/                  # all engine + runner tests
python -m spike.validate_twin  # the twin validation table
```

## Wall discipline

The wall is enforced at the ledger boundary: `Ledger.input_quantity()` refuses a
`benchmark`, `benchmark_quantity()` refuses an `input`. Each engine's
`from_ledger(...)` consumes only inputs (`levels` reads A/I/g_J/g_I; `modes`
reads the COM), so it cannot "predict" a number it was handed. The runner reads
each benchmark once, at comparison time, via `benchmark_quantity()`. If a
validation violates this, the runner surfaces it as an `ERROR` row and exits
nonzero.

## Not yet (follow-ups)

- `beams` / `optics` engine (AC Stark shifts, scattering) — blocked on a measured
  differential-AC-Stark benchmark, which is not yet in the ledger.
- Extend `modes` to the radial spectrum and the full N-ion mode table.
- Graduate the spike to its own repo(s) once the schema/engines stabilise.
