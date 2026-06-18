# Decision records (ADRs)

Architecture/▸ design decisions for `iontrap-reference`, in the lightweight
[ADR](https://adr.github.io/) style. Each file is one decision: the **context**
(forces at play), the **decision**, and its **consequences**. They are
append-only — a superseded decision is not deleted but marked `Superseded by
NNNN`, so the reasoning trail stays intact (good scientific practice: the
record of *why* survives even when the *what* changes).

The running chronological narrative lives separately in
[`../LOGBOOK.md`](../LOGBOOK.md); ADRs capture the load-bearing decisions that
the logbook references.

| # | Title | Status |
|---|-------|--------|
| [0001](0001-record-format-yaml.md) | Records and registries in YAML, validated by JSON Schema | Accepted |
| [0002](0002-source-pdfs-out-of-vcs.md) | Source PDFs kept out of version control; cite resolvable links | Accepted |
| [0003](0003-licensing-cc-by-and-mit.md) | CC-BY-4.0 for data/docs, MIT for code | Accepted |
| [0004](0004-input-benchmark-wall.md) | The mandatory `input`/`benchmark` wall, designed in from the start | Accepted |
| [0005](0005-subsystem-and-slot-vocabulary.md) | Fixed `subsystem` enum and shallow-compositional `configuration` slots | Accepted |
| [0006](0006-validation-in-ci.md) | Two-layer validation (JSON Schema + graph invariants) enforced in CI | Accepted |
