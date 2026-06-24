# State of the Twin

A maintained, structured summary of what is **validated**, what is **diagnostic**,
and what is **open** — so a reader does not have to grep the >1300-line chronological
[LOGBOOK](LOGBOOK.md). Regenerate the tables from the code:
`python -m spike.validate_twin` (validations) and `python -m validator.validate`
(source warnings). Last reconciled: 2026-06-23.

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
| raman_optical | capability | no (input_quantity) | polarization+power-resolved shifts/scatter; absolute differential_stark_hz (saturation-anchored) |
| strobo_sim | capability | n/a | spin⊗Fock stroboscopic propagator → symmetric Floquet comb + heterodyne IF beat; transfer function in [notes/strobo_grating_transfer_function.md](notes/strobo_grating_transfer_function.md) |
| grating_tomography | capability | n/a | char.-function / Wigner transfer-function kernels (chi, double/single-sum kernels, exact η=0, χ→W reconstruction) + Ramsey 2-pulse χ-interferometer (population→χ(Δβ), disk \|Δβ\|≤2η) self-checked vs exact sims; tutorial [notebooks/strobo_grating_tomography](notebooks/strobo_grating_tomography.ipynb) |
| twin_wigner_tomography | raw-data demo | n/a | end-to-end Wigner tomography of a displaced state — reconstruction PIPELINE + measurement MODEL only (analytic χ + idealized SDF readout + shot noise → raw data → reconstructed W); **does NOT run the grating propagator** (strobo_sim). Worked example [examples/wigner_tomography](examples/wigner_tomography/README.md), walkthrough [notes/wigner_tomography_walkthrough](notes/wigner_tomography_walkthrough.md) |
| raman_dephasing | capability | no | relative-phase noise of the 2 beams → Δν/T_φ readout of the residual |
| sideband (thermal) | capability | no | RSB/BSB thermal flops → nbar thermometry (twin_sideband discriminator) |
| twin_strobo | raw-data | n/a | stroboscopic OC carrier flop (phase-grating n=0 baseline; strobo dephasing-decoupling; detuning-scan + AC-Stark-vs-N) |
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
2. **OC flop decay — RESOLVED by sideband thermometry: mostly Raman-beam dephasing.**
   The carrier flop alone was degenerate: its ~87% residual decay read as either
   (A) motional n̄_eff = 1.06 ± 0.27 or (B) Raman-beam dephasing Δν ≈ 28 kHz. The
   **RSB+BSB sideband scan** (`OC_Axial/1_1R_LF_MA`, `twin_sideband`) breaks it: the
   red sideband is near-constant while the blue flops fully, giving **n̄ = 0.27 ± 0.13
   (COLD)** from the RSB/BSB ratio — far below the carrier's apparent 1.06. So the ion
   is cold and the carrier loss is **~64% Raman-beam dephasing** (Δν ≈ 21 kHz, T_φ ≈
   15 µs) + ~36% motional. The apparent "hot ion" was Raman dephasing posing as motion.
   That dephasing ⇔ a **differential beam-path instability of ~λ/2π ≈ 45 nm** (path
   stable to ≪ λ on the µs timescale, ~6 mm/s jitter; `raman_dephasing.path_*`) — a
   STATIC imbalance is excluded, so it points to dynamic (acoustic/vibration) path
   jitter and/or RF-phase / pointing noise.
   Remaining: n̄ = 0.27 still exceeds the RSB-cooled benchmark 0.07 (this run/sequence
   adds motion); and the Raman mutual linewidth has no direct measurement yet.
3. **Radial mode nominals need refining.** Tickle measures 3.22/4.71 MHz vs ledger
   3.0/4.5 (+7.5%/+4.7%); also differs from Thomm 2.7/4.4 (three epochs). Not recorded.
4. **Excited-state hyperfine treated as degenerate** in `raman_optical` (~1–2%
   detuning spread at 20 GHz; configurable hook off by default).

## 4. Source-traceability warnings (validator, non-fatal)

6 `WARN`s, all FAIR-Accessible traceability (not correctness): `doerr2024`,
`kaufmann2022`, `thomm2021`, `weber2025` have no resolvable DOI/permalink;
`codata2018`, `stone2005` identifiers are unverified. Records count: **77 provisional,
6 confirmed** — provisional is the honest default until a value is independently
confirmed.

## 5. Known code-entropy follow-ups (from the 2026-06-22 review)

DONE since the review: **fit-parameter uncertainties** — `spike/bootstrap.py`
(parametric + shot resamplers); `fit_rabi`/`fit_tickle` take `n_boot` and emit
`<key>_err`; the twin propagates a 300-replica shot bootstrap to **n̄_eff = 1.06 ±
0.27**. Also: explicit wall in the diagnostics; sympy cross-checks for the
angular-momentum code (which caught a latent 6j integer-perimeter bug).

Still tracked, deferred as spike-stage gold-plating until graduation: a shared
square-pulse propagator (duplicated damped-flop math across rabi/scatter/sideband/
spin); a shared `angular` module (CG in drive, +6j in raman_optical); centralizing
fit-heuristic constants. The angular and linalg numerics are hand-rolled but
**validated** (CG/6j vs sympy + basis independence to 1e-12; eigensolver vs numpy to
1e-14) and operate only on small quantum numbers / N — not a near-term risk.
