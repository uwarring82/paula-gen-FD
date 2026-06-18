# Task Card — Source-of-truth repository + structured parameter layer

| | |
|---|---|
| **Repo** | `iontrap-reference` |
| **Status** | proposed |
| **Revision** | 3 — config granularity resolved (compositional / shallow); `measured_on`; sources registry; `observation_type` glosses; CI-enforced invariants; uncertainty made forward-compatible |
| **Depends on** | — |
| **Feeds** | end-to-end spike → `iontrap-fields` / `-optics` / `-levels` → twin composition root |

> **Schema freeze recommendation.** Three review passes have moved from
> surfacing gaps to surfacing generically-good ideas. The dominant remaining
> risk is no longer the contract but *extraction discipline* — which no further
> review will surface. Freeze at this revision, build the validator, and do one
> real extraction pass against an actual thesis page to stress-test the contract
> empirically before seeding broadly.

---

## Why

The twin needs one citable substrate. The lineage theses and papers are the
human-readable source; on top sits a **structured parameter layer** that cites
*into* them, so every number is traceable to where it was measured or derived.

The **train/test flag (`kind`)** must be designed in now. It is nearly free to
add while extracting numbers, and very expensive to reconstruct once the twin
exists and everything has been fed in indiscriminately. Without it, calibration
data leaks into the parameterisation set and the validation proves nothing.

Schema richness is *not* population completeness: hardening the contract before
the first number is cheap; doing it after 200 records exist is not. The same
logic is why the uncertainty field is forward-compatible and the configuration
model is compositional from the start — to make a later migration unnecessary.

---

## Scope

1. **Document corpus** — collect the lineage theses (Friedenauer, Clos,
   Wittemer, Hasse, + others) and the key journal papers; register each in the
   **sources registry** with a stable citation key + resolvable link, not just a
   dumped PDF.
2. **Record schema** — define the per-quantity contract *before* populating.
   This is the actual deliverable; the numbers come later.
3. **Registries** — `sources`, `generations`, and the per-slot
   `configuration_slots`, all separate from per-record entries.
4. **Seed population** — enough entries to exercise the schema across all three
   domains (fields / ions / beams). Coverage, not completeness.

---

## Record schema (the contract)

Requirement levels: **[req]** mandatory for any record · **[req-b]** required for
`kind: benchmark` · **[conf]** required to reach `status: confirmed` ·
**[opt]** optional/recommended.

- `name` **[req]** — identifier
- `value` **[req]**
- `units` **[req]** — SI contract
- `kind` **[req]** — `input` (twin consumes) | `benchmark` (twin must
  reproduce) ← the train/test wall
- `uncertainty` **[req]**
  - `value` — scalar (symmetric) **or** `{lower, upper}` (asymmetric, e.g.
    one-sided limits); scalar is the seed default. Consumers must not assume
    symmetry.
  - `type` — `statistical` / `systematic` / `exact` (fundamental constant)
- `source` **[req]**
  - `ref` **[req]** — key into the sources registry
  - `loc` **[req]** — precise location (eq / table / fig / §), precise enough
    that another person recovers the number without interpretation
  - `extraction_note` **[opt]** — what the quoted number refers to / how it was
    read (intake interpretation)
  - `extracted_by` **[req]** — may be auto-filled from git author
  - `extracted_on` **[req]** — date typed in; may be auto-filled from commit
    date. *Distinct from `measured_on`.*
- `scope` **[req]**
  - `isotope` — `24Mg` / `25Mg` / `26Mg` / `apparatus`
  - `subsystem` — `fields` / `motion` / `internal_state` / `optics` /
    `detection` / `control` *(enum to be fixed early — it is a query axis).*
    **Calibration is not a subsystem**: a calibrated quantity takes the
    subsystem of the physics it measures (calibrated trap frequency → `motion`;
    calibrated qubit transition → `internal_state`).
- `observation_type` **[req]** — how the value was obtained; `derived` requires
  non-empty `derived_from`:
  - `direct` — read from electronics / known exactly
  - `fitted` — free parameters optimised to match data
  - `inferred` — deduced from a model with fixed parameters
  - `derived` — computed from other ledger quantities
  - `simulated` — output of a simulation
- `derived_from` **[req]** — list of `name`s this quantity is computed from;
  `[]` for directly obtained values
- `generation` **[req]** — key into the `generations` registry (lineage stage at
  which the quantity entered the apparatus record)
- `configuration` **[conf]** — compositional record of slot keys (see below);
  `null` only if genuinely apparatus-independent
- `measured_on` **[req-b]** — date the *physical quantity* was valid (not the
  extraction date). Mandatory for benchmarks, where temporal drift matters.
- `conditions` **[opt]** — physical state of validity not captured by the
  configuration key
- `caveats` **[opt]** — applicability / breakdown warnings (e.g. "assumes linear
  trap; breaks down near the RF node"). Distinct remit from `extraction_note`
  (referent) and `conditions` (state).
- `status` **[req]** — `provisional` / `confirmed`

### Worked examples

```yaml
- name: omega_z_com
  value: 2.213e6
  units: Hz
  kind: input                         # measured; fed in
  uncertainty: {value: 1.0e3, type: statistical}
  source:
    ref: hasse2021
    loc: "Tab. 4.2"
    extraction_note: "Axial COM mode frequency for stated trap setting"
    extracted_by: UW
    extracted_on: 2026-06-18
  scope: {isotope: apparatus, subsystem: motion}
  observation_type: fitted
  derived_from: []
  generation: hasse
  configuration: {trap: mg_linear_trap_v3, beams: null, b_field: nominal}
  status: confirmed

- name: differential_ac_stark_shift_qubit
  value: 1.2e3
  units: Hz
  kind: benchmark                     # predicted from beam params; never fed in
  uncertainty: {value: 2.0e2, type: statistical}
  source:
    ref: freddy_logbook
    loc: "run 2026-05-12, fig. 3"
    extraction_note: "Shift of 25Mg+ qubit transition under Raman setting R3"
    extracted_by: UW
    extracted_on: 2026-06-18
  scope: {isotope: 25Mg, subsystem: internal_state}
  observation_type: fitted
  derived_from: []
  generation: freddy
  configuration: {trap: mg_linear_trap_v3, beams: raman_R3, b_field: nominal}
  measured_on: 2026-05-12             # physical validity, not extraction date
  status: provisional
```

The twin may consume the beam parameters, detunings, polarisation geometry, and
magnetic-field direction as `input`; it must then reproduce the Stark shift,
held out as `benchmark`. This worked pair should anchor the README — it
communicates the wall faster than any abstract rule.

---

## Registries

`generation` records lineage; `configuration` records physical state; they are
orthogonal and kept apart so per-record entries stay lean.

**Configuration is compositional but shallow** — a fixed set of named slots,
each a key into its own slot registry. This keeps dense parameter-space
exploration from exploding the registry (change one axis → swap one slot) while
deliberately *not* becoming a free-form nesting DSL.

```yaml
sources:
  hasse2021:
    title: "..."
    author: "Hasse"
    year: 2021
    link: "https://freidok.uni-freiburg.de/..."   # DOI / repository permalink
  clos2016:
    doi: "10.1103/PhysRevLett.117.170401"

generations:
  hasse:
    description: "Hasse-thesis generation of the Mg+ apparatus"
    changes:
      - "Raman beam geometry updated"
      - "Trap-drive electronics revised"
      - "New magnetic-field calibration"
  freddy:
    description: "Current in-lab generation (Freddy PhD)"
    changes: [...]

configuration_slots:
  trap:
    mg_linear_trap_v3:
      description: "Linear trap, axial confinement setting A"
      supersedes: null
  beams:
    raman_R3:
      description: "Raman beam setting R3"
      supersedes: raman_R2          # configuration lineage, per slot
  b_field:
    nominal:
      description: "Nominal B-field calibration"
```

---

## Physics inventory to seed (contents, not modules)

- **Fields** — analytical RF E-field, DC E-field, B-field.
  *B defines the quantisation axis; every polarisation and k-vector elsewhere is
  expressed relative to it.*
- **Ions** — Mg ²⁴ / ²⁵ / ²⁶ ground-state structure. Anchor the schema on
  **²⁵Mg⁺ hyperfine** (F = 2 / F = 3, ~1.8 GHz, the qubit) — the seed must
  express the individual F states and the splitting as distinct records. ²⁴ and
  ²⁶ are I = 0 and express a *role* (coolant / co-trapped), never a faked
  hyperfine manifold.
- **Beams** — Raman pair (k₁, k₂, Δk, co- / counter-propagating), cooling,
  preparation, detection. Each carries power, detuning, k-vector, polarisation,
  waist at the ion.
- **Calibration** — frequencies, timings, phase offsets, electronic and
  motional dephasing / decoherence rates. *Predominantly `benchmark`;
  cross-cutting across subsystems, not a subsystem of its own.*

---

## Out of scope (boundaries)

- No physics computation / solver code — that is the spike and the modules.
  (The validator below is schema tooling, not physics.)
- No twin orchestration.
- Not a complete parameter set — seed only.
- The **differential AC Stark shift** is `benchmark`, never `input` — it is both
  routinely calibrated in the lab *and* predictable from beam parameters, so
  feeding it in closes the loop and the validation proves nothing.
- The spike build itself — separate sibling card.

---

## Schema invariants (non-negotiable)

1. `kind` is mandatory and must be `input` or `benchmark`.
2. `benchmark` records must never be loaded into parameterisation or fitting
   routines.
3. Any quantity whose **transitive** `derived_from` closure contains a
   `benchmark` record **inherits `benchmark` status**, unless re-derived
   exclusively from `input` records.
4. **Benchmark by default:** differential AC Stark shifts, measured
   decoherence / dephasing rates, heating rates, final transition-frequency
   residuals.
5. **Input by default:** beam powers, waists, detunings, k-vectors,
   polarisations, trap-drive parameters, externally calibrated constants —
   unless deliberately held out as validation data.
6. `generation` (lineage) and `configuration` (physical state) must not be
   conflated; both resolve to registry keys.
7. `source.loc` must be precise enough to recover the number without
   interpretation; `source.ref` must resolve in the sources registry.
8. `observation_type` of `fitted` / `inferred` / `derived` / `simulated` must be
   declared explicitly; `derived` requires non-empty `derived_from`.

---

## Definition of done

- [ ] Schema documented as a contract, with `kind` **mandatory** and the
      requirement levels stated.
- [ ] `sources`, `generations`, and `configuration_slots` registries scaffolded
      with at least the Hasse and Freddy generations, the lineage theses, and
      one key per slot.
- [ ] Seed entries across all three domains, each fully populated (incl. `kind`,
      `scope`, `observation_type`, `derived_from`, precise `source.loc`, and
      `measured_on` for benchmarks).
- [ ] At least one `input` and one `benchmark` entry demonstrating the wall
      end-to-end; the ²⁵Mg⁺ hyperfine manifold expressed without I = 0 faking.
- [ ] **Validation in CI**, not manual:
  - JSON Schema for field-level validity (required fields, enums, types).
  - A small Python validator for graph-level invariants JSON Schema cannot
    express: the **wall test** (no `input` leakage into the benchmark set and
    converse) and the **inheritance test** (no `kind: input` record whose
    transitive `derived_from` closure touches a `benchmark`).
  - Both run as a pre-commit / CI check.

---

## Open questions

- **Slot vocabulary** — fix the `configuration_slots` top-level set
  (`trap` / `beams` / `b_field` / …) and the `scope.subsystem` enum early; both
  are query axes and expensive to re-key once entries accumulate.
- **Home / migration** — stays under `uwarring82` while the schema is unstable;
  move to an org (and consider a public-facing alias) only after the schema has
  survived the end-to-end spike and one real extraction pass.
