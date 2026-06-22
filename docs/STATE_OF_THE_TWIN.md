# State of the Twin

A maintained, structured summary of what is **validated**, what is **diagnostic**,
and what is **open** — so a reader does not have to grep the >1300-line chronological
[LOGBOOK](LOGBOOK.md). Regenerate the tables from the code:
`python -m spike.validate_twin` (validations) and `python -m validator.validate`
(source warnings). Last reconciled: 2026-06-22.

> **Glossary.** UW = Ulrich Warring (group lead, decisions). PAULA = the Freiburg
> ²⁵Mg⁺ apparatus. *freddy* / *legacy* = apparatus generations (freddy = F. Dörr's
> current setup). TPSR = two-photon stimulated Raman. CC/OC/AC/ROC = the four Raman
> beam combinations (carrier / orthogonal / antilinear / radial-orthogonal). *The
> wall* = the input/benchmark separation ([ADR-0004](decisions/0004-input-benchmark-wall.md)):
> an engine that makes a prediction to be σ-validated may consume only `input`
> records, never the `benchmark` it is checked against.

## 1. Validated benchmarks (σ-validations, wall-enforced)

These run through [`spike/runner.py`](../spike/runner.py); each engine consumes only
`input` records and is compared against a `benchmark` it did not see. CI fails on any
tension > 3σ.

| benchmark | engine | subsystem | residual | n_sigma |
|---|---|---|---|---|
| clock_transition_25mg | levels | internal_state | −2.65 kHz | 1.09 |
| clock_transition_weber_25mg | levels | internal_state | +0.31 kHz | 0.51 |
| omega_z_axial_stretch_2ion_25mg | modes | motion | +21.67 kHz | 2.17 |
| omega_radial_rocking_2ion_25mg | modes | motion | −0.10 kHz | 0.01 |
| doppler_cooling_limit_25mg | cooling | motion | +3.04 µK | 0.03 |
| doppler_cooled_occupation_25mg | cooling | motion | +0.421 | 0.42 |
| bdd_ac_stark_shift_25mg | acstark | optics | +450 kHz | 0.15 |

**7 σ-validations, all ok.** Only `levels`, `modes`, `cooling`, `acstark` are formal
validations; everything below is a diagnostic or a raw-data tool.

## 2. Engine registry (validation / diagnostic / raw-data; benchmark consumption)

| engine | role | consumes benchmarks? | note |
|---|---|---|---|
| levels | **validation** | no (wall) | Breit-Rabi, 2 clock benchmarks |
| modes | **validation** | no (wall) | axial+radial normal modes |
| cooling | **validation** | no (wall) | Doppler limit + occupation |
| acstark | **validation** | no (wall) | far-detuned BDD light shift |
| drive | diagnostic | no | relative MW couplings (CG); apparatus-dominated |
| projection | diagnostic | no | Raman comb → mode geometry |
| sideband | diagnostic | no (input_quantity) | η + sideband Rabi + carrier Debye-Waller |
| scatter | capability | no | scalar Raman scatter + differential AC-Stark |
| raman_optical | capability | no (input_quantity) | polarization+power-resolved shifts/scatter |
| readout | diagnostic | **no** (count levels are `input`) | single-shot fidelity + Fisher info |
| sideband_cooling | diagnostic | **YES, by design** | inverts the *measured* n̄ (benchmark) for κ — uses `benchmark_quantity`, documented |
| tickle | raw-data | n/a | secular-freq spectroscopy on `.dat` |
| rabi | raw-data | n/a | damped-flop fit on `.dat` |
| detection | raw-data | n/a | bright/dark discrimination on `.dat` |
| spin / twin* | composition | n/a | Bloch state; the integrated cycle |

\* `twin_oc_flop` deliberately consumes the RSB-cooled n̄ benchmark via
`benchmark_quantity` (the motional channel is an *inversion*, not a zero-parameter
prediction) — see [ADR-0007](decisions/0007-raman-scatter-vs-dephasing-in-flop-twin.md).

There is **no measured benchmark** for: the Raman *differential* AC-Stark shift, the
absolute MW Rabi rate, the sideband-cooling floor, or the Raman-scattering decoherence
rate. Those engines are capability/diagnostic by necessity, not choice.

## 3. Open physics flags / risks (ranked)

1. **Raman beam polarizations are provisional** ([ADR-0008](decisions/0008-polarization-power-resolved-raman-optics.md)).
   `raman_{b1,b3,r1,r2}_polarization_25mg` are seeded from Clos Tab. 3.2; the
   lab-frame → quantization-axis (B) projection geometry is **not fully tabulated**,
   and **B3's polarization is not stated anywhere** (placeholder, large σ). Because
   Δ_R ≪ Δ_FS the vector shift is unsuppressed, so real circular contamination would
   move δ_AC a lot. **Needs UW: the actual Stokes vectors along B.**
2. **OC flop decay is unexplained at the cooled n̄.** Scattering (~4%) + carrier
   Debye-Waller at n̄=0.07 (~9%) explain only ~13% of the observed ~60% contrast loss;
   the n̄_eff inversion gives n̄_eff ≈ 1 (≈15× the cooled benchmark) → cooling
   underperformance/heating this run and/or technical (intensity/B-field) dephasing.
   Not closed; flagged, not fabricated.
3. **Radial mode nominals need refining.** Tickle measures 3.22/4.71 MHz vs ledger
   3.0/4.5 (+7.5%/+4.7%); also differs from Thomm 2.7/4.4 (three epochs). Not recorded.
4. **Excited-state hyperfine treated as degenerate** in `raman_optical` (~1–2%
   detuning spread at 20 GHz; configurable hook off by default).
5. **Fits report no uncertainties.** `rabi.fit_rabi` / `tickle.fit_tickle` are
   grid searches returning point estimates (χ²_red only); the n̄_eff ≈ 1 result is
   therefore order-of-magnitude. A covariance/bootstrap wrapper is a follow-up.

## 4. Source-traceability warnings (validator, non-fatal)

6 `WARN`s, all FAIR-Accessible traceability (not correctness): `doerr2024`,
`kaufmann2022`, `thomm2021`, `weber2025` have no resolvable DOI/permalink;
`codata2018`, `stone2005` identifiers are unverified. Records count: **77 provisional,
6 confirmed** — provisional is the honest default until a value is independently
confirmed.

## 5. Known code-entropy follow-ups (from the 2026-06-22 review)

Tracked, not yet done (deferred as spike-stage gold-plating until graduation):
shared square-pulse propagator (duplicated damped-flop math across rabi/scatter/
sideband/spin); a shared `angular` module (CG in drive, +6j in raman_optical);
fit-parameter uncertainties; centralizing fit-heuristic constants. The angular and
linalg numerics are hand-rolled but **validated** (CG/6j vs known values + basis
independence to 1e-12; eigensolver vs numpy to 1e-14) and operate only on small
quantum numbers / N — not a near-term risk.
