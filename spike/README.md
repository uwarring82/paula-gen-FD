# Twin spike

A throwaway-ish, end-to-end proof that the structured parameter layer can drive
**physics engines**. Each engine consumes `input` records from the
[ledger](../records/) and is checked against `benchmark` records — the
input/benchmark wall, executed.

> **Boundary.** Per the task card, solver/physics code is *out of scope* for the
> substrate (`records/` + `validator/`). This package is kept separate and may
> later graduate to its own repo(s) — `iontrap-levels` / `-fields` / `-optics`
> plus a twin composition root. It depends only on `pyyaml` (no numpy): for a
> `2S_1/2` ground state the Breit-Rabi problem is solvable in closed form.

## First engine: `levels` (Breit-Rabi)

[`engines/levels.py`](engines/levels.py) models the ²⁵Mg⁺ ground-state
(3s ²S₁/₂) hyperfine + Zeeman structure from the magnetic-dipole constant `A`,
nuclear spin `I`, electronic/nuclear g-factors, and the field `B`. It is exact
for `J = 1/2` and **sign-safe**: `A < 0` (the ²⁵Mg⁺ manifold is inverted, F=3
below F=2) enters only through `A²` and diagonal `m_J m_I` terms.

It reproduces the `clock_transition_25mg` benchmark we validated by hand:

```
$ python -m spike.validate_twin
LEVELS ENGINE — 25Mg+ ground-state clock transition |F=3,mF=0> <-> |F=2,mF=0>
  inputs (consumed):
    A   = -596.254376 MHz   (hyperfine_a_constant_25mg)
    I   = 2.5               (nuclear_spin_25mg)
    B   = 5.500 +/- 0.100 G (b_field_quantization_freddy)
  zero-field splitting 3|A|      = 1788.763128 MHz
  quadratic Zeeman coefficient K = 2195.x Hz/G^2
  predicted clock @ 5.50 G        = 1788.8296 MHz  (+/- ~2400 Hz)
  benchmark (measured, Doerr)     = 1788.8322 MHz  (+/- 200 Hz)
  residual (pred - meas)          = -2.6e+03 Hz  = 1.1 sigma
  --> CONSISTENT within combined uncertainty
```

## Layout

```
spike/
  constants.py      CODATA constants + 25Mg atomic g-factors (sourced)
  ledger.py         loads records/*.yaml into a queryable Ledger (pyyaml only)
  engines/
    levels.py       2S_1/2 hyperfine+Zeeman engine (closed-form Breit-Rabi)
  validate_twin.py  predict from inputs, compare to the benchmark, report residual
  test_levels.py    physics limits + benchmark reproduction + wall enforcement
```

## Run

```bash
pytest spike/                  # engine tests (physics limits + benchmark)
python -m spike.validate_twin  # the clock-transition validation report
```

## Wall discipline

`GroundStateZeeman.from_ledger(...)` refuses to read a `benchmark` record — the
engine may consume only `input`s, so it cannot "predict" a number it was handed.
The benchmark is touched once, at comparison time, in `validate_twin`.

## Not yet (follow-ups)

- The atomic g-factors `g_J` / `g_I` live in `constants.py`; graduate them to
  `input` ledger records so the wall covers them too.
- More engines: `fields` (trap pseudopotential → secular frequencies),
  `beams`/`optics` (AC Stark shifts, scattering), and a composition root.
