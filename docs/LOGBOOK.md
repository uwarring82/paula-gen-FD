# Logbook

A chronological lab notebook for `iontrap-reference`. Newest entries on top.
Load-bearing decisions are captured as ADRs under
[`decisions/`](decisions/README.md) and linked from here. The goal is that the
*reasoning* behind the repository state is recoverable, not just the state.

---

## 2026-06-18 — Repository bootstrap and first extraction pass

**Goal.** Stand up the source-of-truth repository and structured parameter layer
defined in `task card/iontrap-reference-task-card.md` (revision 3, schema-freeze
recommended), following FAIR principles and good-scientific-practice logging.

**Decisions taken** (each with an ADR):

- Records + registries in **YAML**, validated by **JSON Schema** — ADR-0001.
  Confirmed by UW.
- Source PDFs (~148 MB, 12 docs) kept **out of version control**; the sources
  registry cites resolvable permalinks/DOIs instead — ADR-0002. Confirmed by UW.
  Rationale also: copyright, lean history, FAIR Findable/Accessible.
- Licensing: **CC-BY-4.0** for data/docs, **MIT** for code — ADR-0003. Confirmed
  by UW.
- The mandatory `input`/`benchmark` **wall**, with transitive benchmark
  inheritance, designed in now — ADR-0004 (from the task card).
- Fixed `subsystem` enum and shallow-compositional `configuration` slots —
  ADR-0005.
- Two-layer validation (JSON Schema + graph invariants) in CI + pre-commit, with
  validator self-tests — ADR-0006.

**Built.** Directory scaffold; `schema/*.json`; registries
(`generations.yaml`, `configuration_slots.yaml`); `validator/validate.py` with
self-tests (`pytest`: 8 passed); CI workflow; pre-commit config; `docs/schema.md`
(the contract); licensing + `CITATION.cff`.

**First real extraction pass** (task card's pre-seed-broadly recommendation —
stress-test the contract against actual thesis pages). An automated pass over
five theses (`pdftotext -layout` + grep + context reading) produced the seed
records. Sources and precise locations:

| Record | Value | Source · loc | kind |
|--------|-------|--------------|------|
| `omega_z_axial_com_25mg` | 1.30 MHz | Wittemer, Table 3.2, p. 53 | input |
| `omega_radial_com_25mg` | 2.88 MHz | Wittemer, Table 3.2, p. 53 | input |
| `b_field_quantization` | 0.585 mT | Clos 2017, §3.1, p. 34 | input |
| `hyperfine_splitting_25mg_f2_f3` | 2π·1.79 GHz | Clos 2017, §3.1, p. 33 (orig. Itano & Wineland 1981) | **input** |
| `hyperfine_a_constant_25mg` | \|A\| = Δ/3 | derived | input |
| `clock_transition_25mg` | 1788.8322(2) MHz | Doerr 2024, Fig. 2.13, p. 33 | **benchmark** |
| `raman_detuning_from_p32` | 2π·20 GHz | Doerr 2024, §2.1.4, p. 11 | input |
| `nuclear_spin_25mg / 24mg / 26mg` | 5/2, 0, 0 | atomic constants | input |

The worked **wall pair** is real: the literature hyperfine splitting is `input`
(the twin consumes the atomic structure); the in-house precision-measured clock
transition is `benchmark` (the twin must reproduce it, never consume it).

**Data-quality caveats — KNOWN GAPS, do not treat seed as confirmed:**

1. **All seed records are auto-extracted and NOT yet human-reviewed.**
   `extracted_by: claude-agent` records this honestly; nearly all are
   `status: provisional`. UW should review each and re-attribute on confirmation.
2. **`clock_transition_25mg.measured_on` is a PLACEHOLDER** (`2024-01-01`, the
   Doerr-2024 period). The benchmark contract requires a `measured_on`; the
   exact logbook date must replace the placeholder before `status: confirmed`.
3. **`observation_type` for the Wittemer secular frequencies** is recorded as
   `inferred` provisionally — the table does not state whether the values are
   measured or modelled. Confirm.
4. **Uncertainties tagged "quotation precision"** are derived from the rounding
   of the quoted value (±½ last digit), not from a stated experimental error.
5. **Configuration slot labels** (`mg_linear_trap_v3`, `raman_R3`, `nominal`) are
   placeholder vocabulary from the task card; records avoid asserting them on
   real measurements until the actual setting labels are confirmed (ADR-0005).
6. **24Mg/26Mg I=0 records:** the `source.loc` points to the section introducing
   the isotopes; the exact line was not verified in this automated pass. The
   coolant/co-trapped *role* is noted in caveats but not source-verified.

**Source-link verification** (FAIR Findable/Accessible; resolvability checked by
HTTP 200 / DOI / URN resolution, 2026-06-18). Of 12 corpus documents, **5 have
resolvable archival copies**; **7 do not** (the Diplom/MSc theses + Schmitz) and
are registered honestly with `archived: false`/`verified: false` rather than
hidden. Corrections to the initial guesses, now in `sources.yaml`:

- **Friedenauer 2010 is an LMU München PhD, not Freiburg** (DOI 10.5282/edoc.11595).
- **Wittemer is 2019** (not ~2020) and has a DOI (10.6094/UNIFR/151582) — key
  re-set to `wittemer2019` and the two field records updated.
- **Enderlein is 2013, not 2012**, FreiDok 8886, URN only (no DOI) — key
  `enderlein2013`. (FreiDok 9632 is a *different* thesis — do not confuse.)
- **Clos 2017** (FreiDok 12400) and **Hasse 2025** (FreiDok 274764) verified.
- **There is no `hasse2021` paper.** The task card's `hasse2021` key has no
  matching publication; Hasse's earliest paper is Palani, Hasse, …, Warring,
  Schätz, *PRA* 107, L050601 (2023) — registered as `palani2023`. The task-card
  worked example using `hasse2021` is illustrative only; no seed record uses it.
- **Doerr 2024 thesis is not publicly archived**; existence confirmed via the
  author's RTG-DynCAM profile (secondary source). The clock-transition benchmark
  and Raman detuning cite the local PDF; the validator WARNS that this referenced
  source lacks a resolvable identifier — an accepted, visible accessibility gap.

This drove a schema change: `sources.schema.json` no longer hard-requires
`link`/`doi`; instead the validator *warns* on referenced sources that are
unresolvable or unverified (ADR-0006 rationale: register-and-flag beats hide).

**PAULA Mathematica notebooks.** Five own primary-analysis notebooks (~4.2 MB
ASCII, 2014–2018: trap simulations, BEM shim optimisation, mode orientations, Mg
details/scatter rate) were found under `sources/Mathematica/`. Decision (UW):
**keep local, register as sources** — gitignored like the PDFs, but each given a
`paula_*` key in `sources.yaml` (degree `misc`, `archived: false`) so records can
cite them as the origin of `simulated`/`derived` quantities. This also exercises
the internal-source path the schema relaxation opened up. Not yet referenced by
any record.

**Validator status.** `pytest` green (8/8); full-substrate validation green
(`python validator/validate.py` → exit 0, 2 intentional warnings: `doerr2024`
unresolvable, `itano_wineland_1981` unverified).

**Next steps.**
- [x] Finalise `sources.yaml` with verified FreiDok/DOI permalinks (done; 5/12
  resolvable, 7 flagged unarchived).
- [ ] Click-through-verify `itano_wineland_1981` DOI and find/confirm the exact
  Doerr-2024 thesis record if it is later archived.
- [ ] UW review pass: confirm/expand the auto-extracted records; replace the
  placeholder `measured_on`; set confirmed records' `configuration`.
- [ ] Break the `legacy` generation into the real lineage stages; confirm the
  `hasse` generation change-list (currently the task-card placeholder).
- [ ] Expand seed coverage (beam waists/powers were not locatable with a precise
  source this pass — `raman_beam_waist` reported NOT_FOUND).
