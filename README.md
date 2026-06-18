# iontrap-reference

A **citable source-of-truth** for the lineage Magnesium-ion (²⁵Mg⁺) trap
apparatus at the University of Freiburg. The lineage theses and papers are the
human-readable source; on top sits a **structured, schema-validated parameter
layer** that cites *into* them, so every number is traceable to where it was
measured or derived.

> **Status:** schema frozen at task-card revision 3; seeded with a first
> automated extraction pass (provisional, awaiting review). See
> [`docs/LOGBOOK.md`](docs/LOGBOOK.md).

## The one idea: the `input` / `benchmark` wall

Every record carries a mandatory `kind`:

- **`input`** — the digital twin *consumes* it (beam powers, detunings,
  trap-drive parameters, known atomic constants …).
- **`benchmark`** — the twin must *reproduce* it and it is **never** fed into
  parameterisation/fitting (differential AC Stark shifts, measured
  decoherence/heating rates, transition-frequency residuals …).

Without this wall, calibration data leaks into the parameterisation set and the
validation proves nothing. It is mandatory and **machine-enforced** (ADR-0004).

### Worked example (real, from the seed)

The ²⁵Mg⁺ ground-state hyperfine structure, split across the wall:

```yaml
# INPUT — the twin consumes the known atomic structure
- name: hyperfine_splitting_25mg_f2_f3
  value: 1.79e9            # Hz; literature value via Clos 2017 §3.1 p.33
  kind: input             #     (orig. Itano & Wineland 1981)

# BENCHMARK — the twin must REPRODUCE the in-house precision measurement
- name: clock_transition_25mg
  value: 1.7888322e9      # Hz; measured in-house, Doerr 2024 Fig. 2.13 p.33
  kind: benchmark         #     ± 200 Hz; never fed in
  measured_on: 2024-01-01 # the physical-validity date (≠ extraction date)
```

Given the atomic structure (`input`) and the magnetic field, the twin must
predict the measured clock transition (`benchmark`). The canonical case in the
physics is the **differential AC Stark shift**: predictable from beam parameters
(`input`) yet routinely calibrated — so it is `benchmark`, *never* `input`.

## Repository layout

```
schema/         JSON Schemas — field-level contract (record + 3 registries)
registries/     sources.yaml · generations.yaml · configuration_slots.yaml
records/        the parameter layer: fields.yaml · ions.yaml · beams.yaml
validator/      validate.py (schema + graph invariants) + self-tests
docs/           schema.md (the contract) · LOGBOOK.md · decisions/ (ADRs)
sources/pdf/    local-only source PDFs (untracked; see ADR-0002)
```

The full contract — every field, requirement level, and invariant — is
[`docs/schema.md`](docs/schema.md).

### Three registries (kept apart so records stay lean)

- **`sources`** — citation key → resolvable link/DOI. Every `source.ref` resolves
  here; the PDFs themselves are not committed (ADR-0002).
- **`generations`** — apparatus **lineage** stages (e.g. `hasse`, `freddy`).
- **`configuration_slots`** — physical **state**, compositional but shallow:
  fixed slots `trap` / `beams` / `b_field`, each a key into its own namespace.

`generation` (lineage) and `configuration` (state) are orthogonal and must not
be conflated (invariant 6).

## Validate

Requires **Python ≥ 3.10** and **jsonschema ≥ 4.18** (Draft 2020-12). A system
Python with an older `jsonschema` (e.g. Anaconda's 3.x) will fail — use a
virtualenv:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"     # pyyaml, jsonschema, pytest
pytest validator/           # self-tests: the wall, inheritance, cycle, and
                            #   malformed-input checks must all fire
python validator/validate.py
```

Two layers run in CI ([`.github/workflows/validate.yml`](.github/workflows/validate.yml))
and pre-commit:

1. **Field level** — JSON Schema (required fields, enums, types, conditionals).
2. **Graph level** — reference resolution; the **wall test**; the
   **inheritance test** (no `kind: input` whose `derived_from` closure touches a
   `benchmark`); `derived_from` cycle detection; real-date checks.

See ADR-0006. The validator ships self-tests that prove the graph invariants
(the wall, transitive inheritance, cycle detection) reject real violations and
that malformed records are reported rather than crashing the run.

## FAIR & scientific practice

- **Findable / Accessible** — resolvable links/DOIs in `sources.yaml`;
  [`CITATION.cff`](CITATION.cff) for machine-readable citation.
- **Interoperable** — open YAML + JSON Schema; SI units on every value.
- **Reusable** — explicit licensing: data/docs **CC-BY-4.0**
  ([`LICENSE-DATA`](LICENSE-DATA)), code **MIT** ([`LICENSE`](LICENSE)); ADR-0003.
- **Provenance** — every number carries `source.ref`+`loc`, `extracted_by/on`,
  and (for benchmarks) `measured_on`; decisions are logged in
  [`docs/LOGBOOK.md`](docs/LOGBOOK.md) and [`docs/decisions/`](docs/decisions/README.md).

## Scope

Seed coverage across fields / ions / beams — **coverage, not completeness**. No
physics computation, no twin orchestration; this is the substrate they cite.

## License

Data & documentation: [CC-BY-4.0](LICENSE-DATA). Code: [MIT](LICENSE). The cited
theses/papers remain under their own copyright (referenced, not redistributed).
