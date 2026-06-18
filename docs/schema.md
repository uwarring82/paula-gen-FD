# The record contract

This is the **per-quantity contract** for `iontrap-reference`. It is the actual
deliverable; the numbers are seeded on top of it. The machine-checkable form
lives in [`schema/record.schema.json`](../schema/record.schema.json) (field
level) and [`validator/validate.py`](../validator/validate.py) (graph level).

## Requirement levels

| Level | Meaning |
|-------|---------|
| **[req]**   | Mandatory for **any** record. |
| **[req-b]** | Required for `kind: benchmark`. |
| **[conf]**  | Required to reach `status: confirmed`. |
| **[opt]**   | Optional / recommended. |

## Fields

| Field | Level | Notes |
|-------|-------|-------|
| `name` | req | Unique identifier; referenced by other records' `derived_from`. |
| `value` | req | Scalar (symmetric default), vector of numbers (e.g. a k-vector), or a categorical label. |
| `units` | req | SI contract. `dimensionless` / `1` for pure numbers; a category name for categorical values. |
| `kind` | req | `input` (twin consumes) \| `benchmark` (twin must reproduce). **The train/test wall.** |
| `uncertainty.value` | req | Scalar (symmetric) **or** `{lower, upper}` (asymmetric / one-sided). Consumers must not assume symmetry. |
| `uncertainty.type` | req | `statistical` \| `systematic` \| `exact` (fundamental constant). |
| `source.ref` | req | Key into [`registries/sources.yaml`](../registries/sources.yaml); must resolve. |
| `source.loc` | req | Precise location (eq / table / fig / §) — precise enough to recover the number without interpretation. |
| `source.extraction_note` | opt | What the quoted number refers to / how it was read. |
| `source.extracted_by` | req | Who extracted it (may auto-fill from git author). |
| `source.extracted_on` | req | `YYYY-MM-DD` the number was typed in. **Distinct from `measured_on`.** |
| `scope.isotope` | req | `24Mg` \| `25Mg` \| `26Mg` \| `apparatus`. |
| `scope.subsystem` | req | `fields` \| `motion` \| `internal_state` \| `optics` \| `detection` \| `control`. **Calibration is not a subsystem.** |
| `observation_type` | req | `direct` \| `fitted` \| `inferred` \| `derived` \| `simulated`. `derived` ⇒ non-empty `derived_from`. |
| `derived_from` | req | List of `name`s this quantity is computed from; `[]` for directly obtained values. |
| `generation` | req | Key into [`registries/generations.yaml`](../registries/generations.yaml); lineage stage. |
| `configuration` | conf | Compositional slot keys (`trap` / `beams` / `b_field`); `null` only if genuinely apparatus-independent. |
| `measured_on` | req-b | `YYYY-MM-DD` the **physical quantity** was valid (not the extraction date). |
| `conditions` | opt | Physical state of validity not captured by the configuration key. |
| `caveats` | opt | Applicability / breakdown warnings. Distinct from `extraction_note` (referent) and `conditions` (state). |
| `status` | req | `provisional` \| `confirmed`. |

### `observation_type` glosses

- `direct` — read from electronics / known exactly.
- `fitted` — free parameters optimised to match data.
- `inferred` — deduced from a model with fixed parameters.
- `derived` — computed from other ledger quantities (requires `derived_from`).
- `simulated` — output of a simulation.

## Schema invariants (non-negotiable)

1. `kind` is mandatory and must be `input` or `benchmark`.
2. `benchmark` records must never be loaded into parameterisation / fitting.
3. Any quantity whose **transitive** `derived_from` closure contains a
   `benchmark` **inherits `benchmark` status**, unless re-derived exclusively
   from `input` records.
4. **Benchmark by default:** differential AC Stark shifts, measured
   decoherence / dephasing / heating rates, final transition-frequency
   residuals.
5. **Input by default:** beam powers, waists, detunings, k-vectors,
   polarisations, trap-drive parameters, externally calibrated constants —
   unless deliberately held out as validation data.
6. `generation` (lineage) and `configuration` (physical state) must not be
   conflated; both resolve to registry keys.
7. `source.loc` must be precise enough to recover the number without
   interpretation; `source.ref` must resolve in the sources registry.
8. `observation_type` of `fitted` / `inferred` / `derived` / `simulated` must
   be declared explicitly; `derived` requires non-empty `derived_from`.

### What enforces what

| Invariant | Enforced by |
|-----------|-------------|
| 1, 8 (field level) | `schema/record.schema.json` (enum, conditional `if/then`) |
| `measured_on` for benchmarks, `configuration` for confirmed | `schema/record.schema.json` (`if/then`) |
| 2, 4, 5 (the wall + default-kind lints) | `validator/validate.py` |
| 3 (transitive benchmark inheritance) | `validator/validate.py` |
| 6, 7 (reference resolution) | `validator/validate.py` |
| no `derived_from` cycles | `validator/validate.py` |
| dates are real `YYYY-MM-DD` | `validator/validate.py` |

## Fixed vocabularies (expensive to re-key — change deliberately)

- **`scope.subsystem`**: `fields`, `motion`, `internal_state`, `optics`,
  `detection`, `control`.
- **`configuration` top-level slots**: `trap`, `beams`, `b_field`. Adding a
  slot is a deliberate schema change — see
  [`docs/decisions/0005-subsystem-and-slot-vocabulary.md`](decisions/0005-subsystem-and-slot-vocabulary.md).

See the [README](../README.md) for the worked `input`/`benchmark` wall example.
