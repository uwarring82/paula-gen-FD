# 0004 — The mandatory `input`/`benchmark` wall, designed in from the start

- **Status:** Accepted
- **Date:** 2026-06-18
- **Deciders:** UW (task card author), Claude (implementation)

## Context

The eventual digital twin needs parameters to *consume* and, separately,
quantities it must *reproduce* to be validated. If both are fed in
indiscriminately, calibration data leaks into the parameterisation set and the
validation proves nothing — a train/test leak. The task card's central thesis:
the `kind` flag is *nearly free to add while extracting numbers, and very
expensive to reconstruct after the fact*. The same logic motivates making
`uncertainty` forward-compatible and `configuration` compositional from day one
(see ADR-0005, ADR-0006): harden the contract before the first number, not after
200 records exist.

## Decision

`kind ∈ {input, benchmark}` is **mandatory on every record** (no default). The
wall is enforced mechanically, not by convention:

- **Benchmark by default:** differential AC Stark shifts, measured
  decoherence/dephasing/heating rates, transition-frequency residuals.
- **Input by default:** beam powers, waists, detunings, k-vectors,
  polarisations, trap-drive parameters, externally calibrated constants.
- The **differential AC Stark shift is `benchmark`, never `input`** — it is both
  routinely calibrated *and* predictable from beam parameters, so feeding it in
  closes the loop. The validator hard-errors on a `*stark*shift*` named `input`.
- **Transitive inheritance:** any quantity whose `derived_from` closure touches
  a `benchmark` inherits `benchmark` status. The validator rejects a `kind:
  input` record whose closure reaches a benchmark.

## Schema-freeze note

Per the task card, the contract is **frozen at revision 3**: review passes had
moved from surfacing gaps to surfacing generically-good ideas, so the dominant
remaining risk is *extraction discipline*, not the contract. We built the
validator and did one real extraction pass (Wittemer / Clos / Doerr theses; see
[`../LOGBOOK.md`](../LOGBOOK.md)) to stress-test the contract empirically before
seeding broadly.

## Consequences

- **+** The validation set is structurally protected from leakage.
- **+** The seed already carries a real worked pair: literature hyperfine
  splitting (`input`) vs. the in-house measured clock transition (`benchmark`).
- **−** Every extraction must make a kind decision up front; the default-kind
  lints reduce that cost by flagging likely mistakes.
