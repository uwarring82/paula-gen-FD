# 0003 — CC-BY-4.0 for data/docs, MIT for code

- **Status:** Accepted
- **Date:** 2026-06-18
- **Deciders:** UW (decision), Claude (implementation)

## Context

FAIR's *Reusable* principle (R1.1) requires a clear, machine-findable usage
license. A repository that mixes **data** (parameter records, registries, prose)
with **code** (the validator, schemas, CI) is best served by licensing each
under the instrument appropriate to it; a single software license over data is a
poor fit, and a single data license over code is unconventional. Options
weighed: CC-BY-4.0 + MIT; defer/private; CC0.

## Decision

- **Data and prose documentation** (`records/`, `registries/`, `docs/`) →
  **CC-BY-4.0** ([`LICENSE-DATA`](../../LICENSE-DATA)). Attribution is required;
  reuse (incl. commercial) is permitted. CC0 was rejected because thesis-derived
  scientific data should retain an attribution requirement.
- **Software/code** (`validator/`, `schema/`, CI) → **MIT** ([`LICENSE`](../../LICENSE)).
- Machine-readable citation metadata in [`CITATION.cff`](../../CITATION.cff).

Applied now even though the repository stays private under `uwarring82` while the
schema is unstable (task card, "Home / migration") — licensing is cheap to set
now and awkward to retrofit across many files later.

## Consequences

- **+** Satisfies FAIR R1.1; unambiguous reuse terms for each artifact class.
- **+** The CC-BY grant explicitly does **not** relicense the underlying cited
  theses (see `LICENSE-DATA`), avoiding an overclaim.
- **−** Two licenses means contributors must know which covers what; the README
  and each LICENSE file state the split.
