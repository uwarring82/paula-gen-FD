# 0001 — Records and registries in YAML, validated by JSON Schema

- **Status:** Accepted
- **Date:** 2026-06-18
- **Deciders:** UW (decision), Claude (implementation)

## Context

The structured parameter layer is hand-extracted from theses by physicists. The
substrate must be (a) easy to author and review by humans, (b) diff-friendly in
git, (c) comment-able (extraction notes, provenance), and (d) machine-validated.
The task card's worked examples are already written in YAML. Candidate formats:
YAML, JSON, TOML.

## Decision

Store records (`records/*.yaml`) and registries (`registries/*.yaml`) as
**YAML**, and validate them against **JSON Schema** (`schema/*.json`). JSON
Schema validates the parsed data structure regardless of the surface syntax, so
we keep YAML's ergonomics and JSON Schema's tooling.

One subtlety is pinned in the loader (`validator/validate.py`): physicists write
scientific notation as `2.213e6` (unsigned exponent), which PyYAML's default
resolver treats as a **string**, and ISO dates as `2026-06-18`, which it treats
as a `date` object. The custom `RecordLoader` adds a YAML-1.2-style float
resolver and keeps dates as strings, so `value: 2.213e6` is a number and
`extracted_on: 2026-06-18` is a string the schema's `format: date` can check.

## Consequences

- **+** Comments carry provenance inline; reviews read naturally; numbers in
  natural scientific notation.
- **+** Field-level validity is portable JSON Schema.
- **−** The custom loader is a small piece of non-obvious machinery; it is unit-
  tested (`test_unsigned_exponent_parses_as_float`, `test_iso_date_stays_string`)
  so the behaviour can't silently regress.
- **−** YAML has footguns (e.g. the Norway problem); mitigated because the
  schema constrains every field to known types/enums.
