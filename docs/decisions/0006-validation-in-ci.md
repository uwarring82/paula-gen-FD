# 0006 — Two-layer validation (JSON Schema + graph invariants) enforced in CI

- **Status:** Accepted
- **Date:** 2026-06-18
- **Deciders:** UW (task card author), Claude (implementation)

## Context

The contract is only as good as its enforcement. Manual review does not scale
and does not protect the `input`/`benchmark` wall (ADR-0004). Some invariants
are field-level (required fields, enums, types, simple conditionals) and fit
JSON Schema; others are **graph-level** (reference resolution, transitive
benchmark inheritance, cycle-freedom) and JSON Schema *cannot* express them.

The `uncertainty` field is also made **forward-compatible** now: it accepts a
scalar (symmetric, the seed default) *or* `{lower, upper}` (asymmetric /
one-sided limits), so a later migration to asymmetric errors is unnecessary.
Consumers must not assume symmetry.

## Decision

Validation is **two layers**, both run in CI and pre-commit:

1. **Field level** — `schema/*.json` (JSON Schema 2020-12): required fields,
   enums, types, and conditional rules (`derived` ⇒ non-empty `derived_from`;
   `benchmark` ⇒ `measured_on`; `confirmed` ⇒ `configuration`).
2. **Graph level** — `validator/validate.py`: `source.ref`/`generation`/slot/
   `derived_from`/`supersedes` resolution; the **wall test** (default-kind lints
   + the hard differential-AC-Stark-shift rule); the **inheritance test**
   (no `kind: input` whose `derived_from` closure touches a benchmark);
   `derived_from` cycle detection; real-date checks.

The validator ships **self-tests** (`validator/test_validate.py`) that prove the
graph invariants *fire* on constructed violations — the wall (AC-Stark hard rule
+ default-kind lints), transitive benchmark inheritance, and cycle detection —
and that malformed records are *reported, not crashed on* (field-level errors
short-circuit the graph pass). So the enforcement itself is regression-protected,
not just the data. CI runs `pytest` then the
validator over the whole substrate ([`.github/workflows/validate.yml`](../../.github/workflows/validate.yml));
[`.pre-commit-config.yaml`](../../.pre-commit-config.yaml) runs the same locally
plus a large-file guard (ADR-0002).

## Consequences

- **+** The wall and the contract are machine-enforced on every commit/PR.
- **+** `uncertainty` won't need a migration when asymmetric errors arrive.
- **−** Contributors need `pyyaml` + `jsonschema` (and `pytest` for dev); pinned
  in `pyproject.toml`.
- **−** Graph checks are custom code that must itself be trusted — hence the
  self-tests.
