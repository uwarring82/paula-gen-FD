# 0005 — Fixed `subsystem` enum and shallow-compositional `configuration` slots

- **Status:** Accepted
- **Date:** 2026-06-18
- **Deciders:** UW (task card author), Claude (implementation)

## Context

Two vocabularies are **query axes** and are expensive to re-key once entries
accumulate (task card, "Open questions"): `scope.subsystem` and the top-level
`configuration` slot set. They must be fixed *early*, before the seed grows. A
competing pressure is dense parameter-space exploration: if `configuration` were
a free-form nested object, every swept axis would explode the registry.

## Decision

**`scope.subsystem`** is a fixed enum: `fields`, `motion`, `internal_state`,
`optics`, `detection`, `control`. Calibration is deliberately **not** a
subsystem — a calibrated quantity takes the subsystem of the physics it measures
(calibrated trap frequency → `motion`; calibrated qubit transition →
`internal_state`).

**`configuration`** is *compositional but shallow*: a fixed set of named
top-level slots — `trap`, `beams`, `b_field` — each holding a key into its own
slot namespace in [`registries/configuration_slots.yaml`](../../registries/configuration_slots.yaml).
Changing one parameter axis means swapping one slot key. The JSON Schema sets
`additionalProperties: false` on both the slot set and each record's
`configuration`, so an unknown slot is a **validation error** — adding a slot is
a deliberate, reviewed schema change, never an accident. Per-slot lineage is
carried by `supersedes` (resolved within the slot by the validator).

`generation` (lineage) and `configuration` (physical state) are kept strictly
orthogonal (invariant 6).

## Consequences

- **+** Stable query axes; dense sweeps cost one slot-key swap, not a registry
  explosion.
- **+** Vocabulary drift is impossible silently — it trips CI.
- **−** Genuinely new physics axes require a schema edit + this ADR's successor.
  That friction is the point.
- **Open:** the slot *values* seeded so far (`mg_linear_trap_v3`, `raman_R3`,
  `nominal`) are placeholder labels from the task card, not yet confirmed
  against the apparatus; records avoid asserting them until confirmed (see
  `LOGBOOK.md`).
