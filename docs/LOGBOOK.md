# Logbook

A chronological lab notebook for `iontrap-reference`. Newest entries on top.
Load-bearing decisions are captured as ADRs under
[`decisions/`](decisions/README.md) and linked from here. The goal is that the
*reasoning* behind the repository state is recoverable, not just the state.

> **Looking for the current state, not the history?** See
> [`STATE_OF_THE_TWIN.md`](STATE_OF_THE_TWIN.md): the validated-benchmark table, the
> engine registry (validation / diagnostic / raw-data, and which consume benchmarks),
> the ranked open flags, and the source warnings — maintained, with a glossary. This
> log is the chronological appendix; that file is the index.

---

## 2026-06-22 (later 2) — Polarization+power-resolved Raman optics (fine-structure basis)

UW: carefully distinguish the polarization and power of the individual Raman beams —
both significantly affect the shifts and the spontaneous emission. Mid-build UW added
the key physics steer: **F, mF are not good quantum numbers for the P levels** — set
it up more versatile.

Built `spike/engines/raman_optical.py`. The differential AC-Stark shift was a leading-
order scalar (omega_HF/Delta_R) and the scattering had a single scalar `balance`
handle — both polarization- and power-blind. The new engine resolves each beam's
(sigma+, pi, sigma-) polarization and relative power.

**Basis (the versatile/honest choice, per UW's steer).** At Delta_R = 20 GHz the
excited 3P hyperfine + Zeeman is UNRESOLVED, so F', mF' are not good for P. The engine
works in the EXCITED FINE-STRUCTURE |P_{J'} mJ'> basis (no excited hyperfine), with
the GROUND |S_{1/2} F mF> (F good: hf 1.79 GHz >> Zeeman) decomposed into |mJ, mI>
with the nuclear spin mI a SPECTATOR of the optical dipole. The dipole is CG(1/2 mJ;
1 q | J' mJ')*d0 with a COMMON reduced d0 — the D2:D1 = 2:1 line strength comes from
the (2J'+1) multiplicity. The old |F',mF'>-basis sum (Wigner 6j) is kept ONLY as a
degenerate-limit cross-check: it equals the |mJ> sum to **machine precision (1e-12)**
— basis independence, asserted in the tests.

**A normalization bug the validation caught.** The first F'-basis coupling carried an
extra sqrt(2F'+1) (CG vs 3j convention), so the sum rule <g|d^2|g> came out F-
dependent (11.07 for F=3 vs 8.30 for F=2) — unphysical. Fixed to C = CG*sqrt(2F+1)*
{6j}*d_J': the sum rule is now 1.5 for EVERY sublevel and D2:D1 = 2.000 exactly.

**Anchoring.** The absolute field scale cancels in the dimensionless ratios
differential_stark_per_rabi / scatter_per_rabi (delta_AC and Gamma_sc per two-photon
Rabi); the twin multiplies by the measured flop Rabi (delta_AC = ratio*rabi;
Gamma_sc = ratio*2pi*rabi, matching the scalar engine).

**Records (provisional, seeded+flagged).** mg_fine_structure_splitting_3p_25mg =
2.7457 THz DERIVED from the two D-line wavelengths (solid). Per-beam polarizations
raman_{b1,b3,r1,r2}_polarization_25mg from Clos Tab. 3.2 (B1 pi; R1/R2 balanced
linear, C=0), FLAGGED: the lab-frame->B geometry isn't fully tabulated and **B3's
polarization is not stated anywhere** (large-uncertainty placeholder). POWER is a
per-run .dat setting (pwr_b1/pwr_r2), read by the twin.

**KEY FINDING (why polarization matters so much here).** Since Delta_R << Delta_FS
(2.7 THz) the laser sits essentially on P_3/2 alone, so the VECTOR light shift is NOT
fine-structure-suppressed (~0.5x the scalar for circular light). A **10% circular
contamination of R2 changes the differential shift by ~47%**. For the OC run (B1 pi +
R2 balanced linear, equal 0.425 powers) the resolved differential is delta_AC =
**-44.5 kHz** (vs the scalar +18.6 kHz — now signed, ~4.4% amplitude cap), and
Gamma_sc ~2x the scalar. The crude factors were explicitly order-of-magnitude; this
is the documented refinement. The flop's dominant decay is still motional/technical
(ADR-0007), so the ledger floor rises only ~11% -> ~13%; n_bar_eff ~ 1 unchanged.

Wired into twin_oc_flop (report shows resolved-vs-scalar). ADR-0008. Tests 204 -> 222
(+18 raman_optical: 6j, basis independence, cycling, sum rule, line strength,
scalar/vector, two-photon, power/polarization scaling, ledger). Validator green.

---

## 2026-06-22 (later) — Motional channel for the OC flop twin: carrier Debye-Waller + n̄_eff inversion

UW: continue with the suggested follow-up — turn the (fitted) motional/thermal
dephasing residual of the OC carrier-flop twin into a third **predicted**,
ledger-anchored channel.

Done, and the result is more informative than just "explaining" the decay. Added the
**carrier Debye-Waller thermal dephasing** to `spike/engines/sideband.py`: a Δn=0
carrier is NOT motion-free at finite η — Ω_{n,n} = Ω₀·e^(−η²/2)·L_n(η²)
(`carrier_rabi_factor`, pure-Python Laguerre `_laguerre`), so a thermal state (mean
n̄, `thermal_pn`) is a spread of Rabi frequencies and the flop dephases. New:
`thermal_carrier_flip` (exact flip-probability sum, with the AC-Stark detuning folded
in per component), `thermal_coherence` (the envelope |Σ P_n e^{iΩ_n t}| — the cheap
inversion handle), and `thermal_dephasing_rate` (leading-order σ_Ω = Ω₀ η²√(n̄(n̄+1)),
documented as OVERSTATING for n̄<1).

`spike/twin_oc_flop.py` now composes **four ledger-anchored channels** (coherent,
AC-Stark, scattering, motional) and **drops the empirical exponential residual**. The
motional channel is η = 0.389 (sideband, OC→lf at `omega_z_axial_com_25mg` = 1.30
MHz) × the RSB-cooled n̄ = 0.07 (`mg_rsb_cooled_nbar_axial_lf_25mg`, a measured
benchmark — consumed via `benchmark_quantity`, honestly "anchored to a measured n̄").

**KEY FINDING (the honest payoff).** At the cooled n̄ = 0.07 the ledger floor
(scattering ~2 % + carrier Debye-Waller ~9 %) explains only ~**11 %** of the observed
~60 % contrast loss. So the follow-up does NOT make the residual disappear — it shows
the cooled motional state is *also* too cold to cause the decay. The twin then adds an
**effective-n̄ inversion** (solve the same thermal model for the n̄ reproducing the
observed decay): **n̄_eff ≈ 1.1**, ~16× the cooled benchmark. So this OC flop is
consistent with a near-unity-n̄ state — sideband-cooling underperformance / heating in
this run and/or technical Raman intensity-phase / B-field dephasing — **not** the
cooled 0.07. Flagged for the experimentalist, not fabricated away. Figure updated
(`docs/figures/twin_oc_axial_carrier_flop.png`): the bare ledger floor (n̄ = 0.07)
barely decays and sits above the data; the n̄_eff curve tracks it.

**Bug caught by the tests (good).** The first `thermal_carrier_flip` used ½(1−amp·cos)
instead of the generalised-Rabi ½·amp·(1−cos), so a detuned component gave P(0) =
½(1−amp) ≠ 0 — fixed (now exactly 0 at t = 0). Also replaced the inversion's
per-n̄ `fit_rabi` (42 s/`build()`) with the analytic `thermal_coherence` envelope
(γ_floor = (3/4)Γ_sc + (−ln|C(t_span)|/t_span), additive) → `build()` ~3 s, same
n̄_eff. ADR-0007 updated (2026-06-22 section). Tests 193 → 204 (+9 sideband carrier
thermal, twin rewritten to the 4-channel API); validator green (6 pre-existing
warnings).

---

## 2026-06-22 — OC axial Raman carrier-flop twin (AC-Stark + off-resonant scattering)

UW: take a look at the `OC_Axial/0_Car_Flop` data file and set up a twin simulation
including AC-Stark shifts and the spontaneous-emission effects of off-resonant
scattering.

The file is an **OC orthogonal-carrier axial Raman carrier flop** (TPSR, B1+R2,
Δk‖z; PAULA 2026-06-12 09:46:23): the pulse duration `t_oc` is scanned 0→9.75 µs
(21 pts × 40 shots) and the addressed |3,3⟩↔|2,2⟩ population Rabi-oscillates, its
contrast decaying over ~3 periods. Both effects UW asked for follow from the single
Raman detuning Δ_R = 20 GHz (`raman_detuning_from_p32`) + the linewidth Γ.

**`spike/engines/scatter.py`** — the Raman loss budget, all ledger-derived:
- coherent rate Ω = Ω_BΩ_R/(2Δ_R);
- **off-resonant scattering** Γ_sc = Γ(Ω_B²+Ω_R²)/(4Δ_R²), which for balanced beams
  collapses to the clean **Γ_sc = (Γ/Δ_R)·Ω**, so the SE floor is detuning-only:
  **P_SE(π) = π·Γ/Δ_R ≈ 0.66 %**;
- **differential AC-Stark** δ_AC ≈ (ω_HF/Δ_R)·Ω (reuses `sideband.raman_differential_
  stark_factor`), speeding Ω_eff and capping the amplitude at Ω²/Ω_eff²;
- `flip_probability(t,Ω,Γ_sc,δ_AC)` forward flop with the scattering envelope
  e^{−(3/4)Γ_sc t} (full-depolarisation leading order, `CONTRAST_DECAY_FACTOR`);
- `RamanScatter.from_ledger` consumes Δ_R, Γ, ω_HF (all `input`, wall-safe).

**`spike/twin_oc_flop.py`** (`python -m spike.twin_oc_flop` → `docs/figures/twin_oc_
axial_carrier_flop.png`) fits the flop (`rabi.fit_rabi` → Ω/2π ≈ 208 kHz, t_π ≈ 2.4
µs) and **separates the ledger-anchored floor from an empirical residual**.

**KEY, HONEST FINDING.** At Δ_R = 20 GHz the scattering floor removes only ~2 % of
contrast over the scan (Γ_sc-contrast ≈ 2.0×10³ /s, τ ≈ 0.5 ms), while the data lose
~60 % (observed ≈ 1.0×10⁵ /s, τ ≈ 10 µs). **Spontaneous scattering explains only
~2 % of the observed decay** — the *designed* payoff of a large Raman detuning
(Γ_sc/Ω ∝ Γ/Δ_R; Ozeri 2007). The residual (≈ 9.9×10⁴ /s) is **motional/thermal
dephasing** (carrier Debye-Waller spread over the thermal phonon distribution),
carried as a *labelled, fitted* term — NOT dressed up as a prediction. The AC-Stark
differential is δ_AC ≈ 18.6 kHz (0.8 % amplitude cap). The figure overlays the
per-shot count cloud, the full twin, and the bare scattering-only "floor" (which
visibly barely decays) so the separation is legible. Per the FAIR/integrity memo,
the two physical channels stay zero-parameter and wall-enforced; only the residual
is fitted, and it is named.

This is **capability + diagnostic, not a σ-validation**: there is no measured
Raman-scattering decoherence rate in the theses to anchor against (mirrors how
`sideband` treats the also-unanchored differential AC-Stark ratio). Decision +
caveats in **[ADR-0007](decisions/0007-raman-scatter-vs-dephasing-in-flop-twin.md)**.

A 2π unit slip in `scatter_rate` (linewidth → decay rate) was caught by the
balanced-collapse test (Γ_sc=(Γ/Δ_R)Ω must equal the two-beam form) and fixed.
Caveat: the flop fit is bound-sensitive on the noisy 40-shot data (χ²_red ≈ 2.9), so
Ω lands at ~208 kHz with an 80–260 kHz window; the decomposition conclusion
(scattering ≪ observed) is robust to that. Follow-up: predict the dephasing residual
from n̄ (sideband-cooling) × the carrier η² spread (sideband), turning it into a
third ledger-anchored channel.

Tests 172 → 193 (+21: `test_scatter` 14, `test_twin_oc_flop` 7); validator green
(6 pre-existing source-identifier warnings).

---

## 2026-06-19 (later 30) — Tickle engine (Kalis 2016 secular-frequency spectroscopy)

UW: implement a tickle engine — and corrected an initial Lorentzian-dip attempt: the
tickle lineshape follows the group's own published method (kalis2016, then kalis2017).

Added the published reference `kalis2016` (Kalis, Hakelberg, Wittemer, Mielenz, **Warring**,
Schaetz, PRA 94, 023401 (2016), DOI-verified) to `sources.yaml`/`SOURCES.md`, and built
`spike/engines/tickle.py`. A finite resonant excitation pulse (t_exc) drives a mode as a
classical oscillator → amplitude A ∝ sin([ω_exc−ω_i]t_exc/2)/(ω_exc²−ω_i²) — a **sinc** of
FWHM ≈ 1/t_exc (Eq. 3, NOT a Lorentzian). The Doppler modulation (β=⟨u,k_BD⟩A) shifts
carrier population to motional sidebands ±vω_i with Bessel weights J_v(β)² (pure-Python
`besselj`, verified J₀(1)=0.7652/J₀(15)=−0.0142), so F = Σ_v J_v(β)²·g(Δ_BD+vω) (Eq. 2)
dips at resonance — but only for on-resonance detection Δ_BD≈0 (red/blue Δ_BD on a slope
gives a peak). `fit_tickle` fits the robust leading-order **sinc² dip** (depth≥0
constrained) for f₀; the absolute depth/⟨n⟩ depends on the detection sensitivity
(Δ_BD/Γ_w/beam-waist) and is left free (f₀ is robust to <1 kHz).

Applied to `Tickle/PDQ_*_FScan` (`python -m spike.plot_tickle`,
`docs/figures/tickle_modes.png`): t_exc=200/200/100 µs → FWHM ~1/t_exc; **axial
1.299 MHz** (ledger 1.30, −0.1%), **radial 3.224/4.712 MHz** (ledger 3.0/4.5 →
**+7.5%/+4.7%**). The radial nominals need refining (also differ from Thomm's 2.7/4.4 —
three epochs; flagged, not yet recorded). Tests 166→171.

**Verification caught a real (high) bug, fixed.** A first pass used only `files[0]` per
mode and quoted +4.1%/+2.9% — but that file for MF/HF is a NARROW calibration scan that
MISSES the mode, so the "dip" was a sub-noise, edge-pinned fit to flat data. Fixed:
`fit_tickle` now returns an F-test significance + edge-pinned flag + `resolved` bool;
`plot_tickle` fits ALL files, keeps only significantly-resolved bracketed dips, and
reports the robust median (LF 10/10, MF 3/4, HF 4/4 files). The genuine radial offsets
are +7.5%/+4.7%, not the earlier +4.1%/+2.9%. Also corrected the readout diagnostic's
bright/dark decomposition (dark deficit is detection-leak, not preparation; see later 29).

## 2026-06-19 (later 29) — Readout + sideband-cooling diagnostic engines (Thomm benchmarks)

UW: "go with both" — a readout engine and a sideband-cooling engine for the pending
Thomm-2021 benchmarks.

Built both, but as **diagnostics, not σ-validations**: the benchmark values are
preparation/protocol-limited, not zero-parameter predictions, so a forced σ-test would
be circular or show a false tension (FAIR-honest call).
- `spike/engines/readout.py`: single-shot discrimination fidelity (94.8 % at Thomm's
  λ↓=2.682/λ↑=0.036), Fisher-information / Cramér-Rao precision of an ML P↓ estimate
  (overhead ×1.1 at p=0.5, ×3 near bright). The runner `readout` diagnostic decomposes
  the measured 99.4/97.4 % BY CHANNEL: bright 0.6 % deficit = preparation (ML of a
  perfect bright state → ~1, above the 95 % single-shot cap); **dark 2.6 % = DETECTION**
  (dark-state off-resonant scatter/depumping during t_det, a systematic averaging does
  not remove). [Initially mis-stated as all-preparation; corrected after the adversarial
  review flagged it against the ledger's own dark caveat.]
- `spike/engines/sideband_cooling.py`: the RSB floor n̄_min ∼ α(κ/2ω)² (off-resonant
  carrier, LBMW03; α O(1)) and a per-mode κ inversion. The `sideband_cooling` diagnostic
  confirms all κ/ω < 1 (resolved regime) but κ varies → the achieved n̄ (0.07/0.11/0.07)
  are protocol/per-mode-limited, not a common floor; cross-checks Thomm's P(n=0)↔n̄.

The 5 pending benchmarks are now "addressed by diagnostics (not a σ target)". Tests
153→166; removed the dead `optimal_detection_time`.

## 2026-06-19 (later 28) — Realistic detection (depumping tail) + optional ML readout; Raman detuning resolved

UW resolved the Raman-detuning conflict and asked to (a) make the detection engine
more realistic — "we observe repumping effects during detection via BD; bright
histograms show significant amounts of zero and few photons (cf. Thomm 2021)" — and
(b) add a maximum-likelihood state readout, **kept optional** ("in some experiments we
require only photon counts").

**Raman detuning resolved.** History 50→200→**20 GHz**: the final reduction (doerr2024
= current) was deliberate — a new laser system + BBO doubling-stage trouble, traded a
smaller detuning for larger Rabi rates (Ω∝1/Δ_R). The 20 GHz record stands; Thomm's
50/200 GHz are earlier epochs. Updated `raman_detuning_from_p32` caveats + SOURCES.md
(conflict → resolved history).

**Realistic detection.** Added `transition_count_pmf(k, lam_start, lam_other, decay)` to
`engines/detection.py`: a state scatters at `lam_start` but switches to `lam_other` at a
Poisson time (rate `decay`=Γ·t_det); switching at fraction u → Poisson(lam_start·u +
lam_other·(1−u)). So a **bright ion depumps** (lam_start=λ↓>λ↑) → Poisson core + a
**low-count tail** (the Thomm zero/few-photon events), and a **dark ion leaks** into the
cycle → a high-count tail. The twin's `simulate_counts` now applies it via
`MWModel.depump_bright` / `leak_dark` (both default 0 → pure Poisson, backward
compatible); the matching sampler `_detect_count` was verified against the PMF.

**Optional ML readout.** `ml_estimate_p_down(counts, λ↓, λ↑, depump_bright, leak_dark)`
— EM on the two-component mixture P(k)=P↓·B(k)+(1−P↓)·D(k) with the realistic B/D. It's
opt-in; the twin's primary output stays raw counts. Demo `spike/twin_detection.py`
(Thomm λ↓=2.682/λ↑=0.036, representative depump Γt=0.3): the bright histogram shows the
tail, and a Rabi flop read out by ML tracks the true P↓ (RMSE 0.028) while a fixed
threshold is biased low by the tail (RMSE 0.090); ignoring depumping in the ML
underestimates P↓. Figure `docs/figures/twin_detection_depumping.png`.

Tests 153→158 (4 detection: PMF normalize/tail/leak/ML-recovery; 1 twin: depump tail +
mean-lowering). Validator green; adversarial verification run.

---

## 2026-06-19 (later 27) — Thomm 2021 extracted: Freiburg state-detection reference (readout + motional)

UW added `sources/pdf/thomm2021.pdf` (Robin Thomm, MSc Freiburg 2021, "State
Detection of Trapped Magnesium Ions") to get ready to build more on motional effects
and electronic-state readout.

Extracted figures of merit via a 4-lens parallel workflow (readout / cooling+Lamb-Dicke
/ motional-prep / apparatus), each pulling FOMs with verbatim quotes + line numbers
and flagging ledger conflicts. Registered the source in `sources.yaml` and catalogued
it in `SOURCES.md`. **This is the actual Freiburg/PAULA detection reference** — it
supersedes the LMU `friedenauer2010` placeholders.

Added **8 records**:
- *Electronic readout* (`records/beams.yaml`, `detection`): `mg_detection_time_25mg`
  30 µs; `mg_bright_counts_25mg` λ↓=2.682(9) / `mg_dark_counts_25mg` λ↑=0.036(3) per
  30 µs (single ion, ~89 kHz bright); benchmarks `mg_readout_fidelity_bright_25mg`
  99.4(6)% / `mg_readout_fidelity_dark_25mg` 97.4(6)%. The fidelities are **maximum-
  likelihood state-probabilities** (2/3-Poissonian), not single-shot thresholds — a
  hard threshold on these means gives only ~93% bright (recorded in caveats).
- *Motional ground state* (`records/fields.yaml`, `motion`): RSB-cooled n̄ per mode
  `mg_rsb_cooled_nbar_{axial_lf,radial_mf,radial_hf}_25mg` = 0.07(2)/0.11(3)/0.07(4)
  (P(n=0)=94/90/94%, 3D 79(4)%, T≈26 µK) — the resolved-sideband counterpart to the
  Doppler `doppler_cooled_occupation_25mg` (n̄≈10).

**Conflicts handled (FAIR: documented, not overwritten):**
- ⚠️ **Raman detuning** — Thomm: virtual-level detuning raised **50→200 GHz** (new BBO
  cavity); ledger `raman_detuning_from_p32` = **20 GHz** (doerr2024). 10× mismatch
  flagged in that record's caveats + SOURCES.md; **needs UW** (Doerr typo? different
  epoch? re-detuning?). Matters: sideband ∝1/Δ_R, AC-Stark ∝1/Δ_R².
- ✅ **Lamb-Dicke** — Thomm η=0.39 at 1.3 MHz is *consistent* with the ledger's 0.32 at
  1.92 MHz via η∝ω^(−1/2) (0.32·√(1.92/1.30)=0.389); recorded as a third cross-check,
  not a conflict.

Catalogued (not yet records) for the next build: mode freqs 1.3/2.7/4.4 MHz, RF 56 MHz,
B≈5.5 G, hyperfine ω₀≈1774 MHz; and the motional-engineering set (Fock |1⟩ 89(3)%/|2⟩
76(2)%, coherent |α|=0.67–0.84, squeezing |ξ|=0.74(3), Wigner negativity, Ωₙ fit with
n_max=18) for a future Wigner/motional-tomography engine. Validator green (6 warnings).

---

## 2026-06-19 (later 26) — Rabi vs Ramsey FREQUENCY scans + adversarial verification

UW: do frequency scans Rabi vs Ramsey at τ = 100/300/600 µs, with the detuning range
adapted to the expected resolution.

Added `make_seq_ramsey_freq(tau)` (frequency-domain Ramsey) and `spike/twin_freqscan.py`
(`python -m spike.twin_freqscan` → `docs/figures/twin_freqscan_rabi_ramsey.png`). The
Rabi spectroscopy line (broad, ±25 kHz) is **pulled to f₀+δ_ACZ** (MW on), while the
Ramsey **fringe comb sits on the bare f₀** (free precession, MW off) with spacing
1/τ_eff, τ_eff = τ + 1/πΩ (the finite-pulse correction — measured 8.65/3.17/1.62 kHz
match exactly). Per UW, the window is **resolution-matched** (±2.5 fringe spacings), so
longer τ zooms in (RPE-style ladder). AC-Zeeman = Rabi dip − Ramsey comb centre.

A 3-lens adversarial workflow (physics / estimators / figure) returned 16 findings, **2
confirmed medium**, both fixed:
- The Ramsey comb centre is **envelope-pulled** by the AC-Zeeman shift (a uniform per-τ
  shift: +0.43/+0.16/+0.08 kHz at τ=100/300/600), biasing AC-Zeeman ~3% low — larger
  than its error bar and absent from the budget. Fixed: subtract the known per-τ
  `_model_pull` (noiseless comb×envelope; for real data one fits it). The symmetric-mean
  fix was tested and *rejected* (the pull shifts all teeth equally, so it doesn't cancel).
- The "17× sharper" headline rested on a lucky-low 12-seed std of a **heavy-tailed**
  estimator (rare catastrophic fringe mis-ID). Fixed: N_SEED 12→100 + **robust stats**
  (median, 16–84%). Now AC-Zeeman = **3.02 ± 0.41 kHz** (injected 3.00, unbiased),
  Ramsey comb ~13× sharper than the Rabi dip. The physics + 1/τ_eff law were confirmed
  correct; `_curve_fringe_spacing` relabelled as the model (noiseless) spacing.

---

## 2026-06-19 (later 25) — Integrated twin: prepare→drive→detect, AC-Zeeman + dephasing, Rabi-vs-Ramsey inference

UW: the engines must *weave together* state-preparation (incl. cooling),
manipulation and detection — and the microwave engine should consider AC-Zeeman
shifts and spin dephasing (individual runs **and** ensemble averages), with a first
test case inferring AC-Zeeman + dephasing from a Rabi-vs-Ramsey comparison.
Detection efficiency from Friedenauer or Schmitz.

Built the **integrated experiment-cycle twin** (distinct from the ledger-validation
engines): a single quantum **state** flows through prepare → drive → detect.

- `spike/engines/spin.py` — the qubit **state** as a Bloch vector (z=+1→|↓⟩,
  P_up=(1−z)/2) with coherent operations `prepare(eps)` / `pulse(Ω,δ,t,φ)`
  (Rodrigues rotation about the generalised-Rabi axis) / `free(δ,t)`.
- `spike/twin.py` — the cycle. Per shot: prepare (with prep error) → MW drive
  carrying an **AC-Zeeman shift** + **quasi-static Gaussian dephasing**
  (σ_δ=√2/2πT₂\*, drawn per shot) → **projective measurement** (QPN) →
  **Poisson detection** (μ_bright for |↓⟩ / μ_dark for |↑⟩). Produces both the
  individual-run **count cloud** (`simulate_scan`) and the **ensemble average**
  (`ensemble_p_up`, Gaussian-weighted over the dephasing). `fit_ramsey` fits the
  fringe (cosine × Gaussian envelope). The detection levels come from
  `detection_levels` = R_scatter × η × t_det, delegating to the existing
  `engines.detection.expected_bright_counts`.
- **Key physics for the test case:** the AC-Zeeman shift is present only while the
  MW is *on* (the pulses), absent during a Ramsey free gap. So the Rabi (MW-on)
  resonance *includes* it, while the Ramsey fringe (free precession) runs at the
  bare detuning δ_set — which **equals** the AC-Zeeman shift when the drive is tuned
  to the Rabi resonance (δ_set=acz). The fringe frequency then *is* the AC-Zeeman
  shift; the fringe's Gaussian decay is T₂\*.
- `spike/twin_demo.py` (`python -m spike.twin_demo`) — inject (Ω=50 kHz,
  AC-Zeeman=3 kHz, T₂\*=800 µs), simulate Rabi + Ramsey per shot, infer over **N=16
  Monte-Carlo replicas**: Ω = 50.15 ± 0.74 kHz (Rabi), AC-Zeeman = 3.00 ± 0.04 kHz,
  T₂\* = 791 ± 58 µs (Ramsey) — all consistent with injected within the scatter.
  Figure `docs/figures/twin_rabi_ramsey_inference.png`.

**Detection provenance.** Added `mg_detection_efficiency_25mg` = 5.6×10⁻³ and
`mg_fluorescence_count_rate_25mg` = 250 kHz (Friedenauer 2010, subsystem
`detection`). The efficiency is the three-factor product 4.7% solid-angle ×
(1−14% Schott-filter/mirror loss) × 14% PMT QE at 280 nm = 5.66×10⁻³; Schmitz
cross-references the resulting >99.99% readout fidelity. Both `provisional`, with a
caveat that this is the LMU/MPQ predecessor apparatus (lineage; the Freiburg value
may differ).

**Verification.** A 5-lens adversarial workflow (Bloch physics, AC-Zeeman/dephasing
model, inference, detection provenance, architecture) returned 18 findings, **1
confirmed material**: the demo originally reported a single seed's Ω=49.3 kHz (a
1.4%-low Monte-Carlo fluctuation, not a bias — the estimator is unbiased, 50.08 kHz
over 200 seeds) **without an uncertainty**. Fixed per the FAIR/integrity memo by
inferring over N=16 replicas and quoting mean ± std. Also tightened: wired
`detection_levels` into the demo (so the "Friedenauer detection" footnote is
literally true), stated the δ_set=acz precondition in the docstrings (the fringe
frequency is δ_set, = acz only on the Rabi resonance), renamed
`collection_eff`→`detection_eff`, and bumped the ensemble grid 21→41 for large-τ
accuracy. One finding (a `pulse_us` AttributeError) was refuted — already fixed.

Tests 139 → 151 (spin 4, twin 7, +1 detection rename); validator green.

**Next.** Real Ramsey/RPE data exists (`Strobo2.0/…RPE_Fdet…`, phase scans) — apply
the twin to it to extract a *measured* AC-Zeeman shift + T₂\*; weave the cooling
engine's scatter rate into the detection levels for the Raman/sideband cycle.

---

## 2026-06-19 (later 24) — MW workflow generalized to any transition (MW_Control)

UW: apply the MW twin-vs-data workflow to sources/data/MW_Control; a "next
transition" is now calibrated. Also: keep PHOTON COUNTS on the y-axis (no
normalisation).

Generalised `spike/plot_scans.py` to ANY |F,mF> <-> |F',mF'> transition: it parses
the transition from the scan parameter (fr_mw_<tr> / t_mw_<tr>, dropping variant
tags like _s), maps it to the levels-engine key + the t_mw_<tr>/fr_mw_<tr> ion
properties, auto-discovers + pairs freq/duration scans per transition (same
session), and writes docs/figures/mw_<tr>_<date>_twin_vs_data.png. Made `fit_rabi`
auto-derive its frequency grid from the time span (so a 7 kHz flop and a 54 kHz
flop both fit). Switched the y-axis to RAW photon counts (per UW): the twin curve
and QPN band are mapped into count space via the per-scan bright/dark levels
(used only to place the twin, not to rescale the data).

MW_Control (2026-06-12) has two transitions: the stretched |3,+3>-|2,+2> (t_pi 9 us,
~54 kHz) and the **newly-calibrated |3,+1>-|2,+2>** (~7x weaker, t_pi 70 us,
7.1 kHz). The well-calibrated 3p3 matches the twin tightly within QPN; the 3p1 twin
EXPOSES two real calibration issues: the actual t_pi (~66 us from the dip spacing)
is ~6% below the configured 70 us, and the t=0 prep sits ~1 count below the full-
bright 2pi revival (an early-time contrast transient). The levels-engine f0 offsets
(3p3 ~ -100/-120 kHz, 3p1 -61 kHz at the Weber field) differ per transition by their
field-sensitivity -> a 2-transition joint field determination is now possible.

Committed the MW_Control .dat data (group's own; format per kalis2017). Tests
133 -> 139; substrate green.

---

## 2026-06-19 (later 23) — Provenance correction: data origin vs .dat format

UW: label the data with the date/time from the .dat filenames; Kalis only described
the DAQ details, he did NOT take the data.

Correction applied across the repo: **kalis2017 is the .dat FORMAT / DAQ-method
reference ONLY.** The example data under `sources/data/microwave/` are the PAULA
group's OWN primary measurements, taken **2026-06-15** (13:28:34 frequency scan,
13:28:39 duration scan — parsed from the filenames `HH_MM_SS_DD_MM_YYYY`, added as
`DatFile.timestamp`). Fixed: the figure title/subtitles now read "PAULA, 2026-06-15"
with each scan's time and a "format: kalis2017" footnote (was mis-titled
"kalis2017"); `registries/sources.yaml` and `docs/SOURCES.md` clarify kalis2017 is
format-only and the data are the group's; the analyze_data header carries the
measurement date. Going forward, any record from this data cites the **.dat file
(date/time)** for the measurement and **kalis2017** for the format/method — not
kalis2017 for the data.

---

## 2026-06-19 (later 22) — Twin-vs-experiment plot with QPN (PAULA, 2026-06-15)

UW: plot the freq + duration scans with the digital-twin result, including QPN, all
relevant parameters read from the data-file ion properties.

`spike/plot_scans.py` (-> docs/figures/mw_3p3_2p2_twin_vs_data.png) overlays the
generalized-Rabi twin P = Omega^2/Omega_eff^2 sin^2(Omega_eff t/2) (new
`rabi.generalized_rabi`) on both scans, with the quantum-projection-noise band
(new `detection.qpn` = sqrt(P(1-P)/N)). EVERY twin parameter is read from the
`<ionproperties>`: N = exp_point = 75; t_pi = t_mw_3p3_2p2 = 9.30 us ->
Omega/2pi = 53.76 kHz; resonance fr_mw_3p3_2p2 = 1775.470 MHz; freq-scan pulse =
9.30 us. Counts -> P(|up>) via the bright/dark asymptotes from the duration fit.

Result: the twin tracks the experiment within QPN on BOTH scans -- the freq-scan
sinc lineshape (central lobe + null + sidelobe) and the duration sin^2 flop. The
data resonance sits ~17 kHz above the operating setting, and the levels-engine
prediction (1775.599 MHz, Weber B) is +129 kHz (the field cross-check, annotated).
A small duration-scan phase drift by ~30 us reflects configured 53.76 vs fitted
53.3 kHz. Tests 129 -> 132; substrate green.

---

## 2026-06-19 (later 21) — Raw-data engines: rabi + detection (kalis2017 .dat)

UW: "build more engines: Rabi rate, and detection" — consuming the group's own raw
microwave-Rabi data (|3,+3> <-> |2,+2>, taken 2026-06-15; .dat FORMAT per kalis2017
— see the later-23 provenance correction).

**`spike/datfile.py`** reads the PAULA DAQ `.dat` files (settings + scan signal +
per-shot count histograms). **`spike/engines/rabi.py`** fits a damped Rabi flop
(pure-Python weighted linear LS for c,a,b at each (f,gamma) grid point via
linalg.solve) — the duration scan gives **Omega/2pi = 53.3 kHz** (t_pi 9.38 us),
-10% vs mw_rabi_3p3_2p2_doerr (59.45 kHz): the expected MW-power/day dependence of
the apparatus-limited rate (a MEASUREMENT, not a prediction). **`spike/engines/
detection.py`** does Poissonian bright/dark discrimination (optimal threshold,
readout fidelity, Mandel Q). Run via **`spike/analyze_data.py`** (kept OUT of the
ledger-based runner — these read raw data, not ledger records, and produce
measurements, not engine-vs-benchmark validations). The freq scan's resonance
(1775.49 MHz) cross-checks the levels engine's field-sensitive (3,+3)<->(2,+2)
prediction (1775.60 MHz at the Weber field; ~0.05 G residual).

**Verified by a 4-lens workflow (fit math + Poisson stats CONFIRMED correct):**
independent period 18.85 us -> 53.05 kHz matches the fit; linalg.solve agrees with
numpy; threshold/fidelity verified by hand. Fixes applied -- two important: (1) the
detection fidelity gap (0.992 Poisson vs 0.973 empirical) was mis-attributed to
"dark-tail events"; it is ONE-SIDED, the BRIGHT histogram's low-count tail =
bright-state loss/depumping during detection (the dark/reference channel matches
Poisson); (2) the ~0.013 counter is the REFERENCE/normalisation channel (|up>-ion
proxy), NOT a measured dark-ion state, and the two counters are written as separate
BLOCKS, not interleaved -- corrected datfile.py, DATA_FORMAT.md, the report.
Plus: dof off-by-one (5 params -> chi2_red 2.53, >1 = residual scatter exceeds shot
noise); structural (variance-based) signal-counter selection; grid-edge flag;
empty-histogram + NaN/Inf guards. Tests 114 -> 129; substrate green.

---

## 2026-06-18 (later 20) — Source catalog (docs/SOURCES.md)

UW: carefully catalog all sources in a dedicated md file with metadata, own-words
summaries, and the digital-twin figures of merit each highlights.

Wrote **`docs/SOURCES.md`**: a FOM-at-a-glance table (the 7 sigma-validations and
their sources), a master index of all 33 registered sources, Part A (the 10 ACTIVE
sources supplying the twin's 65 records -- full metadata + own-words summary + the
exact input/benchmark FOMs each contributes), Part B (registered lineage corpus,
14 resolvable + 4 unresolved), Part C (5 internal PAULA notebooks), and an
open-provenance section. The FOM map was built mechanically by grouping records by
source.ref (clos2017 23, doerr2024 20, itano 6, wittemer2019 5, kaufmann2022 5,
weber2025 2, friedenauer2010/hasse2025/codata2018/stone2005 1).

Verified by a 2-lens workflow (registry/FOM fidelity + summary grounding vs the
thesis texts): all metadata, counts, values, and summary claims confirmed faithful.
Fixes applied -- **provenance correction**: the doerr2024 registry title was a
flagged placeholder paraphrase; read the real title from the thesis title page and
corrected BOTH `registries/sources.yaml` and the catalog: "Advanced Interferometer
Techniques for Measuring Near-Resonant Light Shifts and Superresolving Trapped-Ion
Dynamics". Plus: codata2018/stone2005 index icons corrected to "unverified" (they
are archived:true but verified:false); noted BDD = Clos's "BDdet"; listed the 20th
doerr record (clock_transition_residual). Substrate green.

---

## 2026-06-18 (later 19) — Sideband + AC-Stark engines (Raman couplings & light shifts)

UW: "go for it" (sideband engine) + "Should we include AC stark shifts of all
beams?" -> clarified: "I meant adding it for the Raman beams."

**Scope answer (evidence-backed): NOT all beams.** The AC-Stark shift is a
COHERENT far-detuned effect (delta_AC = s Gamma^2/(8 delta) = Omega^2/(4 delta),
Clos Eq. 2.2.24). Only far-detuned beams shift; near-resonant beams scatter (->
cooling engine). Hasse 2025 confirms it experimentally: the far-detuned BDD beam
shifts the cycling transition by ~2pi x 10 MHz, while the resonant repumpers RD/RP
"induce no significant ac Stark shift".

**Engine `acstark`** + benchmark `bdd_ac_stark_shift_25mg` (Hasse ~10 MHz): the
7th sigma-validation, predicted |delta_AC| = 10.45 MHz vs 10(2), 0.15 sigma. A
diagnostic flags per beam which regime applies (only BDD = coherent shift).

**Engine `sideband`**: absolute Lamb-Dicke eta(comb,mode,omega) = |Delta_k/k .
e_mode| * k z_bar(omega) (z_bar ~ omega^-1/2), anchored to the measured eta = 0.32
(OC->lf, 1.92 MHz). Gives sideband Rabi Omega_{n,n+-1} = eta sqrt(n+1|n) Omega_0 and
the **Raman differential AC-Stark shift** delta_AC_diff ~ (omega_HF/Delta_R) Omega_0.
Re-encoded the four Delta_k records to carry MAGNITUDE (Delta_k/k: CC 0, OC/ROC
sqrt2, AC 2); the projection engine normalises internally (cosines unchanged). Added
single-ion radial mode freqs omega_radial_mf/hf_25mg (Doerr 3.0/4.5 MHz). eta table:
OC->lf 0.389; AC->lf 0.389, mf 0.222, hf 0.105; ROC->mf 0.222, hf 0.105.

**Verified by a 4-lens workflow (physics correct):** every eta re-derived from
scratch, formula = Clos Eq. 2.2.24, BDD matches Hasse independently, sign/magnitude
handling honest. Applied findings -- the Raman differential factor was the soft spot:
reframed it as an ORDER-OF-MAGNITUDE estimate (~0.045-0.089 Omega_0; the prefactor
0.5-1.0 is convention-dependent -- a plain two-state model gives omega_HF/2Delta_R;
it also ignores CG / P_1/2 / beam-imbalance; NO measured anchor -- BDD validates only
the single-beam SCALE, not this ratio; and omega_HF=3|A| vs the qubit splitting
differ ~0.8%). Plus: derived the anchor projection from the ledger (removed a
hardcoded sqrt2); clarified the AC-Stark level-vs-resonance sign; noted the Gamma
source discrepancy (Clos 41.8 vs Hasse 42.7 MHz, both within band); cross-generation
|Delta_k| invariance note; zero-input guards; pinned the OC-vector magnitude in tests.

Records 62 -> 65; tests 98 -> 114; validations 6 -> 7. Substrate green.

---

## 2026-06-18 (later 18) — Measured Doppler-cooled occupation benchmark (Clos n_bar)

UW: "ok. go" -> the Clos measured Doppler-temperature / mean-phonon benchmark.

**Found the measurement** in Clos 2017, Fig. 3.14 (p. 51): a single Doppler-cooled
25Mg+ ion has **n_bar = 10(1)** at omega_1/(2pi) = 1.915(2) MHz (from Rabi flops
sensitive to the thermal motional distribution). Clos's own consistency estimate:
the Doppler limit (Gamma = 41.8(4) MHz, delta = -Gamma/2, s = 0.5) is <~1 mK ->
n_bar ~ 10 at 2 MHz.

**Cooling engine** gained `mean_occupation(omega, T)` (Bose-Einstein) and
`doppler_limit_occupation(omega)` = the thermal occupation at T_D, which reduces to
the parameter-free **n_bar = 1/(exp(2 omega/Gamma) - 1)** (depends only on
omega/Gamma, since T_D = hbar*Gamma/2kB). Records: `omega_z_axial_clos_25mg`
(1.915(2) MHz input) + `doppler_cooled_occupation_25mg` (n_bar = 10(1) benchmark).

**6th validation -- the first MEASURED motional test.** Predicted n_bar(1.915 MHz)
= 10.42 vs measured 10(1) -> +0.42 sigma. Unlike the theory-consistency Doppler
limit (both sides from Gamma), this is independent: the engine predicts n_bar from
Gamma + omega, the 10(1) is the measurement. Cross-check: n_bar(2.0 MHz) = 9.96
matches Clos's stated ~10. (Made the runner table render dimensionless cells.)

**Adversarial verification (3-lens workflow, ultracode): all "correct".** Physics
-- thermal-state-at-T_D is the standard model AND verbatim what Clos does (p_n =
n_bar^n/(n_bar+1)^(n+1)); the full Bose-Einstein form (via expm1) matters at ~5%
(classical kT/hw = 10.91 vs BE 10.42, and 10.42 is what matches Clos). Thesis --
n_bar=10(1), omega_1=1.915(2), the conditions, Fig. 3.14, and **page 51 all
confirmed** (page verified from the running headers); n_bar=10(1) is a genuine
fit-extracted measurement (observation_type: fitted). Wiring -- wall-clean
(consumes Gamma + omega, reads the benchmark separately), non-circular,
sigma_pred = 0.053 (Gamma-dominated). Applied low-severity refinements: noted T_D
is the optimal-detuning FLOOR (true steady-state n_bar sits >= it; the agreement
validates the BD setting is near-ideal); flagged Clos's own Fig. 3.14
carrier-vs-stretch text/caption ambiguity; added the Clos corroboration (explicit
Gamma/delta/s conditions) to doppler_cooling_limit_25mg.

Records 60 -> 62; tests 95 -> 98; substrate green.

---

## 2026-06-18 (later 17) — Mode-projection engine (Raman -> motional mode)

UW: "mode projection first" (before the Clos Doppler-temperature benchmark).

**`spike.engines.projection`** — predicts which motional mode each Raman (TPSR)
combination addresses, from the effective k-vector direction and the radial-mode
tilt, both consumed from the ledger. Coupling onto a mode = the direction cosine
|Delta_k_hat . e_mode| (relative Lamb-Dicke). Mode axes: lf along z; mf/hf in the
x-y plane tilted 30 deg. The geometry reproduces Doerr 2024's addressing exactly:

    comb  ->lf    ->mf    ->hf    addresses    vs Doerr
    CC    0.000   0.000   0.000   (carrier)    ok
    OC    1.000   0.000   0.000   lf           ok
    AC    0.707   0.612   0.354   lf,mf,hf     ok    (0.707 = cos45, Doerr's "45 deg")
    ROC   0.000   0.866   0.500   mf,hf        ok    (cos30 / sin30 on the tilted axes)

Wired as a runner DIAGNOSTIC (not a sigma-row: the projections have no
independent numeric benchmark beyond the tautological 45 deg; the teeth are the
addressed-modes reproduction). Re-encoded the four raman_*_combination_25mg
records from direction STRINGS to normalised Delta_k VECTORS so the engine
consumes them directly; extended `Ledger._coerce_value` to resolve vector
(-> tuple) and categorical values alongside scalars.

**Adversarial verification (4-lens workflow, ultracode).** Independent
re-derivation, convention-skeptic, thesis-grounding, and code review all
confirmed the PHYSICS correct (every direction cosine re-derived from scratch;
the two-step radial projection proven equal to the direct dot product because
Delta_k_y = 0 and the radial modes have e_z = 0). Findings applied:
  * MEDIUM (provenance): the per-beam laser parameters were cited as Clos 2017
    "Tab. 3.1" but Clos's laser table is **Tab. 3.2** ("Specifications of all
    laser systems", p. 39); Tab. 3.1 is a different linewidth matrix. Corrected
    all 14 per-beam citations + the header (Wittemer's table genuinely is 3.1).
  * LOW: tightened the 4 Raman records to Doerr "Sec. 2.1.4"; noted that only the
    Delta_k AXIS (not sign) is fixed and the engine uses |.| ; noted Doerr does
    not pin which of mf/hf sits on which tilted axis (labelling only); documented
    the two-step==dot precondition and that abs() discards the displacement-phase
    sign; added a 3-vector dimensionality guard + threshold justification.

Records unchanged in count (60); tests 82 -> 95; substrate green. A future
sideband engine can take |Delta_k| (CC 0, OC/ROC sqrt2, AC 2) + the mode
frequencies to turn these direction cosines into absolute eta per mode.

---

## 2026-06-18 (later 16) — Raman (TPSR) beam combinations + radial-mode tilt

UW: add the Raman beam-combination properties (CC/OC/AC/ROC) and "note our radial
mode tilt"; check the PhD thesis work to fill in the data. (UW first gave
AC = B1+R3, then corrected: **AC = B3+R1** — matching Doerr.)

**Extracted the Raman geometry** from Doerr 2024 (Sec. 2.1 + Fig. 2.4),
confirmed by the Hasse 2025 glossary, with polarisations from Clos 2017 (Tab.
3.2) and the Lamb-Dicke parameter from Wittemer 2019 / Clos 2017. The setup is a
clean 2x2 of two blue (B1, B3) x two red (R1, R2) beams, all in the trap x-z
plane; the effective k-vector Delta_k = k_B - k_R picks the mode:

| comb | beams | geometry | Delta_k | addresses |
|------|-------|----------|---------|-----------|
| CC   | B1+R1 | parallel     | ~0          | carrier only (no motion) |
| OC   | B1+R2 | orthogonal   | \|\| z       | axial lf mode (~1.3 MHz) |
| AC   | B3+R1 | antiparallel | \|\| -(x+z)  | ALL three modes (45 deg to lf axis) |
| ROC  | B3+R2 | orthogonal   | \|\| x       | radial mf/hf (~3.0 / ~4.5 MHz) |

Polarisations (Clos Tab. 3.2): B1 pi, R1 sigma+ + sigma-, R2 sigma+ - sigma-.

**Radial-mode tilt (flagged by UW).** Doerr Fig. 2.4: the two radial modes are
NOT in the x-z plane but lie at ~30 deg to the x/y axis. Recorded as
`radial_mode_tilt_25mg` = 30 deg. This is exactly why the AC combination (k in
the x-z plane) needs a TWO-step projection to reach the radial modes — captured
in the AC record's caveat.

**Lamb-Dicke.** `raman_axial_lamb_dicke_25mg` = 0.32 for Delta_k along axial at
omega_ax ~ 1.9 MHz (Wittemer eta=0.32 @ 1.920 MHz; Clos eta~0.32 @ ~2 MHz),
scaling as omega_ax^(-1/2).

**Records (6 new, 60 total).** Four categorical combination records (value =
Delta_k DIRECTION string in the trap frame; the schema's value oneOf already
allows strings + a category `units`), the axial Lamb-Dicke, and the radial tilt.
All `input` (configuration/geometry the twin consumes); no engine validates them
yet -> a future mode-projection engine could consume the Delta_k directions +
the 30 deg tilt to predict each combination's per-mode Lamb-Dicke projection.
Substrate green; 82 tests unchanged.

---

## 2026-06-18 (later 15) — Cooling engine (Doppler scattering + limit)

UW: "go for it" on the cooling/scattering engine.

**`spike.engines.cooling`** — clean two-level physics (no apparatus model): from
the natural linewidth Gamma, detuning Delta and saturation parameter s = I/Isat
(all in the ledger from the laser table), R_sc = (Gamma/2) s/(1 + s + (2Delta/
Gamma)^2) and T_D = hbar*Gamma/(2 k_B). Added k_B to constants.py.

**Validation (5th in the composition root).** The engine reproduces the Doppler
limit from Gamma: predicted **1.0030 mK** vs the theses' stated **~1 mK**
(`doppler_cooling_limit_25mg` benchmark), +3 uK = 0.03 sigma. The optimal
detuning -Gamma/2 = -20.9 MHz EQUALS the BD cooling setting (bd_cooling_detuning),
i.e. the lab cools at the textbook-optimal detuning. Honest framing: the limit is
T_D = hbar*Gamma/(2 k_B), so this is a theory-consistency check (both derive from
Gamma); a measured Doppler temperature / mean phonon number would be a stronger
benchmark (noted on the record).

**Capability + diagnostic.** Scatter rate per Blue Doppler beam (from detuning +
s): BD 26.3, BDX 21.2, BDD 6.2 Mphotons/s; max (saturated, on resonance) =
Gamma/2 = 131 Mphotons/s. Added a `cooling_diagnostic` to the runner. Made the
results table unit-aware (Hz -> MHz/kHz, K -> mK/uK).

Tests 75 -> 82; substrate green. The composition root now spans levels (clock x2),
modes (axial + radial), and cooling (Doppler limit) — three subsystems.

---

## 2026-06-18 (later 14) — Beams completed from the canonical laser table

UW: check clos2017, hasse2025, wittemer2019, friedenauer2010 (to fill the beam
gaps from later-13).

**Found the canonical laser table** — Clos 2017 Tab. 3.1 (confirmed by Wittemer
2019 Tab. 3.1): per beam, wavelength | polarisation | transition |
Gamma_nat/(2pi) | detuning Delta/Gamma | intensity I/Isat.

**RESOLVED the P1/2-vs-P3/2 ambiguity:** RD/RP operate on **3P_1/2** (280.353 nm),
distinct from BD/BDX/BDD on 3P_3/2 (279.635 nm). Confirmed by ALL FOUR theses
(Clos, Wittemer, Hasse 'coupling 2S1/2 and 2P1/2'; Friedenauer 'S1/2 F=2 ->
2P1/2'). Doerr's beam-section P1/2 was right; the stray 'P3/2 F=3' was wrong.

**Filled the gaps** and rewrote `records/beams.yaml` (17 optics/detection records,
53 total):
- detunings (Delta/Gamma * Gamma): BD -0.5G = -20.9 MHz, BDD -10G = -418 MHz,
  BDX -0.1G = -4.18 MHz (near resonance), RD/RP -0.5G(P1/2) = -20.65 MHz;
- intensities I/Isat (saturation parameter s): BD 0.5, BDD 20, BDX 0.2, RD 0.5,
  RP 1 — recorded as separate `*_saturation` records;
- constants: Gamma_P3/2 = 41.8 MHz, Gamma_P1/2 = 41.3 MHz, Isat = 255 mW/cm^2
  (2550 W/m^2, Friedenauer), BD wavelength 279.635 nm, RD wavelength 280.353 nm,
  BD detection waist ~50 um (Clos; Wittemer ~40 um; Friedenauer/LMU ~30 um).

Beam POWER is now recoverable (P = s * Isat * beam area) though not yet stored as
a value. This is exactly the input set a future cooling/scattering engine needs
(scatter rate R = (Gamma/2) s / (1 + s + (2 Delta/Gamma)^2)).

---

## 2026-06-18 (later 13) — Established the cooling/preparation/detection beams

UW: establish the cooling (BD), preparation (RD/RP) and detection (BDX) beams
(lab labels).

Extracted from Doerr 2024 + Kaufmann 2022: the Blue Doppler (BD) laser
(~279.635 nm, S1/2 -> P3/2, sigma+, natural linewidth Gamma ~ 42 MHz) is
AOM-split into three beams driving the cycling transition
|down>=|F=3,mF=3> <-> P3/2|F=4,mF=4>: **BDX** (near resonance, detection),
**BD** (red by ~Gamma/2 ~ 20 MHz, Doppler-limit cooling to ~1 mK), and **BDD**
(far red ~10-12 Gamma, cooling hot ions). **RD/RP** (sigma+) optically pump the
ion into |down> for preparation (RD: S1/2 F=2 -> P; RP: S1/2 F=3 -> P).

Established 7 input records in `records/beams.yaml`: `mg_p32_natural_linewidth`,
`bd_laser_wavelength`, `bd_doppler_cooling_detuning`, `bdd_far_cooling_detuning`,
`bdx_detection_detuning`, `rd_repump_polarization`, `rp_repump_polarization`
(joining the existing `raman_detuning_from_p32`). 44 records total; substrate
green.

Gaps flagged (coverage, not completeness): beam POWERS and WAISTS at the ion are
not given in the sources (omitted); the RD/RP exact detunings/powers are not
stated (recorded the polarisation categorically); and the RD/RP excited-state
target is ambiguous in Doerr (P_1/2 in the beam section vs P_3/2 F=3 elsewhere) —
flagged for confirmation. BDD detuning spans the quoted ~10-12 Gamma.

---

## 2026-06-18 (later 12) — Absolute-rate model: attempted, NOT achievable (honest)

UW: "go for it" — extract Doerr Section 3 and build the absolute-rate drive
engine (|CG| x antenna model).

**Extracted Doerr Section 3.** Omega ∝ sqrt(P_MW) confirmed. But the antenna
polarization + frequency-response characterization is in FIGURES (3.2-3.7), not
tables; Table 3.1 is a Rabi-fit-method comparison for the clock (two fit models),
not antenna data.

**Finding (negative, but important).** A rigorous absolute-rate model is NOT
extractable from the theses:
- A 5-parameter physical model (3 polarization gains + a quadratic frequency
  response) x |CG| fits the 8 Doerr rates to only **~25% RMS** (worst residual
  +53%).
- Even within ONE polarization (sigma-, 4 points), the response (rate/|CG|) is
  non-monotonic and asymmetric: 70 -> 85 (peak ~1786 MHz) -> 64 -> 32 — not
  reducible to a few parameters, and only 8 sparse points total.
- The apparatus changed between generations (Doerr ~2-3x faster than Kaufmann,
  NON-uniformly), so no shared model.

**Decision:** did NOT ship an overfit/circular fit. Added only the correct
STRUCTURE to the drive engine — `absolute_rabi = |CG| * apparatus_factor`, with
`apparatus_factor()` returning the empirical per-transition calibration
(measured/|CG|). The validated atomic deliverable stays the relative |CG| (+ the
drive diagnostic). Tests +1 (75); substrate green.

**Natural next step (offered):** the near-resonant AC ZEEMAN shift (Doerr's
actual thesis result) is the canonical light-shift benchmark and IS predictable
from the Rabi rates (calibration, already in the ledger) + the detunings (levels
engine) via delta ~ Omega^2/(4 Delta). But its MEASURED values are also in
Doerr's figures, so a clean validation needs the raw figure data or a dedicated
measurement.

---

## 2026-06-18 (later 11) — Drive engine (microwave Rabi couplings) + calibration

UW: check Kaufmann 2022 and Doerr 2024 for the microwave Rabi rates driving the
hyperfine transitions; then (1) build the CG relative-rate engine and (2) add the
rates as calibration records.

**Extracted.** Doerr Table 2.1 (B~5.5 G): all 8 |3,mF⟩↔|2,m′F⟩ π-times → Rabi
rates 7.0–59.5 kHz. Kaufmann Table 4.1 (B=5.6860 G): 5 transitions, 2.1–23.8 kHz
— the EARLIER apparatus, ~2-3x slower than Doerr (a real apparatus-evolution
finding).

**Key physics (the headline).** The measured Rabi rate folds the atomic
Clebsch-Gordan coupling AND the MW antenna polarization + frequency response —
and the APPARATUS DOMINATES. Mirror pairs with identical |CG|=0.845
(|3,+3⟩↔|2,+2⟩ vs |3,−3⟩↔|2,−2⟩) differ 5x in rate; the apparatus factor
(measured/|CG|) spans ~6x across the manifold. So these are CALIBRATION, not a
clean atomic benchmark.

**(1) Drive engine** `spike/engines/drive.py`: a pure-Python Racah-formula
Clebsch-Gordan + `HyperfineDrive`, predicting the *relative* atomic magnetic-
dipole couplings |⟨F mF; 1 q|F′ m′F⟩|. Verified against sympy EXACTLY (all 8
transitions + CG identities); mirror symmetry and σ±/π labels tested.

**(2) Calibration records** `records/control.yaml`: 13 Rabi rates (Doerr 8 =
generation freddy; Kaufmann 5 = legacy), `kind: benchmark`, `subsystem: control`,
`observation_type: fitted`. No `microwave` configuration slot exists yet, so
configuration is omitted (flagged).

**Runner.** `drive_diagnostic` prints |CG| vs measured Rabi → the apparatus
factor, honestly labelled NOT a sigma-validation. `uncovered_benchmarks` excludes
`subsystem: control` (drive calibration, covered by the diagnostic). The 4 clean
sigma-validations (levels×2, modes×2) are unchanged.

Follow-up: an absolute-rate drive engine needs an MW antenna model (polarization
geometry + frequency response, cf. Doerr Sec. 3); this engine captures the atomic
part only. Tests 68 -> 74; substrate green (5 honest unresolved-source warnings).

---

## 2026-06-18 (later 10) — Levels engine vs Weber 2025 measured manifold

UW: compare the ground-state manifold splittings against weber2025.pdf.

Extracted Weber 2025 ("A Tunable Quantum Magnetometer …", MSc, Schätz group)
Table 3: **8 MEASURED** 25Mg+ ground-state |3,mF⟩↔|2,m′F⟩ qubit transition
frequencies (1775.6–1802.0 MHz) at the Zeeman-spectroscopy field
**B = 5.6454(9) G**, and it independently confirms A = 596.254376(54) MHz.

**The levels engine reproduces all 8 measured transitions to < 5.3 kHz** (clock
+0.31 kHz) and the listed responsivities (dω/dB) to ~0.0005 MHz/G — a real
MEASURED validation of `hyperfine_transitions` (the modes were calc cross-checks;
this is measurement). The clock responsivity +0.0248 MHz/G matches exactly.

Formalised the high-value, non-circular pieces:
- `weber2025` source (title added); `b_field_zeeman_weber_25mg` (precise MEASURED
  field 5.6454(9) G, input); `clock_transition_weber_25mg` (measured clock
  benchmark, 1788.8328(6) MHz).
- `_validate_weber_clock` (4 validations now): predicted 1788.8331 vs measured
  1788.8328 MHz, +0.31 kHz = 0.51σ. INDEPENDENT of Doerr's clock (1788.8322);
  the two measurements agree to ~0.6 kHz, consistent with the field difference.
- Full 8-transition comparison captured as a regression test
  (`test_reproduces_weber_2025_zeeman_manifold`).

Caveat: B = 5.6454 G was determined FROM these transitions (Zeeman
spectroscopy), so the field-SENSITIVE comparison is partly circular; the
field-insensitive clock and the responsivities are the independent tests. Tests
67 -> 68; substrate green (weber2025 honestly flagged unarchived).

---

## 2026-06-18 (later 9) — Radial modes, hyperfine spectrum, provenance fix

UW asked to add radial modes and the full ground-state hyperfine transitions.

**INTEGRITY FINDING.** Extracting Wittemer Table 3.2 revealed it is a CALCULATION
(Wübbena et al. 2012), NOT a measurement, and the experimental crystal is
mixed-isotope 25Mg+ + 26Mg+ (the table has a 25+25 equal-mass column too). My
earlier `omega_z_axial_stretch_2ion_25mg` benchmark was therefore MISLABELLED as
"measured". Decision (UW): relabel the mode references as calculated
cross-checks. Fixed: `observation_type: inferred -> simulated`, notes/conditions
clarified (Wübbena 2012, 25+25 column, not a measurement). Open schema tension
flagged: `measured_on` does not really apply to a `simulated` benchmark (kept as
the publication-context placeholder).

**Modes — radial.** Added `RadialModes` (transverse Hessian B_mm = (ω_r/ω_z)² −
Σ1/|u|³, B_mn = +1/|u|³; ω_radial,p = √μ_p·ω_z), with a zigzag-instability guard
(negative μ → raise). Reproduces the Wittemer 25+25 radial out-of-phase calc
EXACTLY: √(ω_r²−ω_z²) = 2.5699 MHz vs 2.57 (0.01σ). Verified vs numpy to ~1e-15
for N=2,3,4; N=6 at this ratio correctly flags the zigzag instability (real
physics). Added `omega_radial_rocking_2ion_25mg` benchmark.

**Levels — hyperfine spectrum.** Added `hyperfine_transitions(B)` (all 15
|F=3,mF⟩↔|F=2,mF′⟩ microwave transitions with |ΔmF|≤1; the clock is the (0,0)
entry) and `zeeman_splitting(F,mF,B)` (≈ g_F μ_B B; matches the Landé value to
~0.1%). The clock is now one line of the full spectrum.

**Runner.** Added `_validate_radial_rocking` (3 validations now); neutral
"benchmark" column (covers measured + calculated references). All consistent:
clock 1.09σ, axial stretch 2.17σ (the table's own √3-rounding), radial rocking
0.01σ. Tests 58 -> 67.

Follow-ups: a true MEASURED secular-frequency benchmark, and an unequal-mass
engine for the 25+26 experiment.

---

## 2026-06-18 (later 8) — Twin composition root

UW chose the composition root before the third (beams) engine: lock the
ledger-inputs → wall → benchmark → residual-in-sigma interface down with two
engines, before beams adds complexity (and while beams still lacks a measured
differential-AC-Stark benchmark to validate against).

**`spike/runner.py`** — registry of `Validation`s, a uniform `ValidationResult`
(predicted, measured, sigma_pred/meas, residual, n_sigma, status, consumed
inputs), `run_all` (captures per-validation errors as ERROR rows rather than
aborting), `uncovered_benchmarks` (measured benchmarks no engine covers), and one
rendered table. `validate_twin.py` is now a thin CLI shim → `runner.main`.
`check_clock`/`check_modes` became `_validate_*` returning `ValidationResult`.

```
benchmark                        engine  ...  residual/kHz  n_sigma  status
clock_transition_25mg            levels  ...  -2.65         1.09     ok
omega_z_axial_stretch_2ion_25mg  modes   ...  +21.67        2.17     ok
```

The runner exits nonzero on any tension (>3σ) or error. Tests prove it: detects
tension, surfaces a wall violation (engine consuming a benchmark) as ERROR, and
flags uncovered benchmarks. Adding the next engine = implement + one `_validate_*`
+ register.

**Cleanup (review):** refreshed `spike/README.md` (stale CLI output, wall section,
done follow-ups) and made the CI step labels engine-agnostic. Tests 50 -> 58.

---

## 2026-06-18 (later 7) — Second engine: axial normal modes

UW chose the self-contained "normal-mode (mode ratios)" shape for the second
engine (no trap-geometry extraction needed).

**`spike.engines.modes`** — axial normal modes of an N-ion chain. From the COM
secular frequency it computes mode frequencies sqrt(lambda_p)*omega_z, with
lambda_p the eigenvalues of the dimensionless axial Hessian at the ion
equilibrium positions (James 1998). Equilibrium via Newton, eigenvalues via a
pure-Python cyclic-Jacobi solver in `spike/linalg.py` (still no numpy).

Validation: COM (`input`, 1.30 MHz) -> predicted 2-ion stretch sqrt(3)*COM =
2.2517 MHz vs the measured stretch (`benchmark`, Wittemer Table 3.2, 2.23 MHz).
Added `omega_z_axial_stretch_2ion_25mg` as a benchmark. Residual +21.7 kHz =
2.17 sigma -> CONSISTENT; the sqrt(3) relation holds to ~1% in the 3-sig-fig data.

**Adversarial review** (3 agents). Physics verdict **correct** — equilibrium,
Hessian, and eigenvalues independently reproduced via scipy/numpy/sympy:
positions match to <6e-15, lambda_1=1 and lambda_2=3 PROVEN universal (exact
eigenvectors), N=3 -> {1,3,29/5}, N=4/5 match the literature, pure-Python eigvalsh
agrees with numpy to <5e-14, sqrt(3) confirmed for equal masses. Numerics
verdict correct_with_caveats; hardened the (engine-unreachable) edge cases:
- Newton now uses an N-independent max-component |g| criterion and RAISES on
  non-convergence (previously returned silently; N>=40 had a 2-norm noise floor —
  now converges).
- eigvalsh: scale-relative convergence + RAISES if max_sweeps is exhausted.
- axial_mode_eigenvalues guards against a negative eigenvalue before sqrt().
- solve() validates shape; docstrings corrected ("sqrt(3) analytically exact,
  ~1e-15 numerically"). +4 robustness tests.

Tests 36 -> 50; both engines CONSISTENT; substrate green.

---

## 2026-06-18 (later 6) — Full wall coverage: every engine input via the ledger

Two review-driven refinements so that *every* quantity the engine consumes —
constants included — comes through the wall.

1. **Closed the wall gap.** `validate_twin` had pulled the field
   (`b_field_quantization_freddy`) and the measured clock outside the wall, so
   only A and I were kind-checked. Added `Ledger.input_quantity()` (must be
   `kind:input`) and `benchmark_quantity()` (must be `kind:benchmark`); the
   runner now routes the field via `input_quantity` and the clock via
   `benchmark_quantity`, and `from_ledger` uses `input_quantity` for A/I. A field
   accidentally flipped to `benchmark` (or the clock to `input`) now fails loudly.

2. **g_J / g_I → ledger inputs.** The only remaining magic numbers in the levels
   engine. Added `g_factor_electron_2s12` (= free-electron g_e, sourced to
   `codata2018`) and `g_factor_nuclear_25mg` (= -0.34218 from the 25Mg moment,
   sourced to `stone2005`) as `input` records; `from_ledger` reads them through
   the wall, falling back to `constants.py` only if absent. Updated
   `qubit_quadratic_zeeman_coeff_25mg.derived_from` to include them (it genuinely
   depends on g_J/g_I). Values equal the old constants, so the prediction is
   unchanged (residual still -2651 Hz, 1.09 sigma).

Substrate 18 -> 20 records (+3 warnings: the two new sources are `verified:false`
pending a click-through, plus the known `doerr2024`). Tests 33 -> 36. The
pattern is now set for the next engine: every physical input is ledger-driven.

---

## 2026-06-18 (later 5) — First physics engine: the levels (Breit-Rabi) spike

Started the twin: a `spike/` package (UW chose in-repo, split later) with the
first engine. Physics/solver code stays OUT of the substrate (`records/`,
`validator/`); the spike imports the ledger and is checked against benchmarks.

**`spike.engines.levels`** — 25Mg+ ground-state (2S_1/2) hyperfine+Zeeman
structure, closed-form Breit-Rabi (exact for J=1/2, no numpy; sign-safe for the
inverted A<0 manifold). `from_ledger` consumes ONLY `input` records (refuses
benchmarks). `spike/validate_twin.py` predicts the clock from ledger inputs and
compares to the measured benchmark:

    predicted 1788.829549 MHz  vs  measured 1788.832200 MHz
    residual -2651 Hz = 1.09 sigma  -> CONSISTENT

So the hand-built forward validation (later 4) is now real engine code.

**Adversarial review** (3 agents: physics re-derivation, code, break-it). All
three: physics CORRECT — the closed form matches an independent 12x12
diagonalization to ~1e-7 Hz from 0 to 1e5 G; conservation, antisymmetry,
sign-handling, state counts all hold. Findings were robustness/docs, now fixed:
- **(medium) I=0 reachable crash** — `nuclear_spin_24mg/26mg` are I=0 `input`s;
  `from_ledger` built then `clock_transition` threw a bare `KeyError`. Added a
  domain guard (I >= 1/2; clock needs half-integer I) + tests.
- **(medium) runnability docstring** — dropped the false `python spike/...py`
  claim (relative imports); `python -m spike.validate_twin` from a source
  checkout is the supported path.
- **(low) sigma robustness** — a missing uncertainty no longer silently flips
  the run to "tension"; non-finite sigma now raises, 0/0 handled.
- **(low) missing-`kind`** — `quantity()` no longer defaults a kindless record
  to `input` (aligns with the wall).
- **(low/doc) g_I "independent"** — corrected: g_I enters the m_F=0 clock at
  ~12 Hz/5.5 G via (g_J mu_B - g_I mu_N), which the engine retains exactly
  (engine K 2195.77 vs the g_I-neglected ledger 2195.37 Hz/G^2; noted).

Tests 30 -> 33; substrate green; CI runs the spike tests + the clock
reproduction.

---

## 2026-06-18 (later 4) — Forward clock validation across the wall

Promoted the back-inferred field to a proper **forward validation** (UW: "go").
The wall is now exercised end-to-end on real data:

- **inputs** → `hyperfine_splitting_calc_25mg` (Itano nu_hf0), `qubit_quadratic_zeeman_coeff_25mg`
  (K), and **`b_field_quantization_freddy`** — a NEW independent field input
  (~5.5 G, calibrated from the ground-state Zeeman splitting, *not* from the
  clock, so no leakage);
- **prediction** → `clock_transition_predicted_25mg` = nu_hf0 + K·B² = 1788.8295
  MHz (`input`, all-input closure);
- **benchmark** → `clock_transition_25mg` = 1788.8322 MHz (measured);
- **residual** → `clock_transition_residual_25mg` = +2.7(24) kHz (`benchmark`:
  a transition-frequency residual, invariant 4, and benchmark by inheritance).

So the original ~69 kHz "tension" reduces to a **2.7 kHz residual, consistent
with zero** within the ~0.1 G field uncertainty (field stated as ~5.5 G; the
clock implies 5.61 G). `b_field_from_clock_25mg` is kept as the field-space
mirror of the same check. Numbers re-verified numerically; validator green
(18 records, 1 expected warning).

*Review follow-up:* `b_field_quantization_freddy` and `clock_transition_predicted_25mg`
had `configuration: null`, which asserts apparatus-INdependence — wrong for a
field and a value computed at it. Both now omit `configuration` (provisional),
consistent with the sibling field records; the remaining `null`s are all genuine
atomic constants (incl. K, which depends only on nu_hf0 + fundamental constants).
Also relabelled `clock_transition_residual_25mg` uncertainty `statistical ->
systematic` (the 2.4 kHz is the ~0.1 G field systematic propagated through the
prediction; the measured statistical part is ~0.2 kHz).

---

## 2026-06-18 (later 3) — Clock-transition tension resolved (quadratic Zeeman)

**Decision (UW):** dissolve the ~69 kHz tension by the magnetic field; Itano &
Wineland is the source of truth for the hyperfine constant "for now".

**Physics.** The |F=3,mF=0> <-> |F=2,mF=0> clock transition is only *first-order*
field-insensitive; it carries a second-order (quadratic) Zeeman shift
nu(B) = nu_hf0 + K·B², with K = g_J² mu_B² / (2 h² nu_hf0). Taking Itano's
zero-field nu_hf0 = 3|A| = 1788.763 MHz as truth, the measured offset
nu_clock - nu_hf0 = 69.1 kHz implies **B = 5.609(11) G** — consistent with
Doerr's independently stated ~5.5 G (at exactly 5.5 G the predicted clock is
1788.8295 MHz, a 2.7 kHz residual, within field calibration). Verified
numerically (K = 2195 Hz/G²). So the offset is physics, not an error.

**Also corrected an extraction misreading:** `omega_{0,0}` denotes the
mF=0<->mF=0 component, NOT "zero-field-extrapolated" — the Doerr value is the
finite-field transition. `clock_transition_25mg` extraction_note/conditions/
caveats updated accordingly.

**Added two derived records** (honoring the wall): `qubit_quadratic_zeeman_coeff_25mg`
(K; `input`, all-input closure) and `b_field_from_clock_25mg` (the inferred field;
`benchmark`, since it derives from the benchmark clock transition — a diagnostic,
NOT an independent field input). Validator green: 15 records, 1 expected warning.

---

## 2026-06-18 (later 2) — Validator hardening from code review

A code review (UW + external) raised three real issues; all fixed, with
regression tests.

- **High — validator could crash on malformed records.** `schema_check` called
  `rec.get()` before the `isinstance(rec, dict)` guard, and `graph_check` called
  `.get()` on `source` when it was a bare string (e.g. `source: bad`). Fixes:
  guard first in `schema_check`; type-guard nested mappings in `graph_check`; and
  **short-circuit the graph pass when any field-level (schema) error exists** —
  graph invariants assume schema-valid, well-typed data. Reproduced the exact
  crashing inputs end-to-end → now a clean error report (`source: 'bad' is not of
  type 'object'`), graph checks deferred, exit 1, no traceback.
- **Medium — `configuration: {}` slipped through.** A `confirmed`,
  apparatus-dependent record could pass with an empty configuration (the field
  merely had to exist). Added `minProperties: 1` to the configuration object in
  `record.schema.json`; `null` is still allowed for apparatus-independent
  quantities. *(Follow-up review:* `minProperties:1` still let an all-null object
  through, e.g. `{trap: null}`, which resolves nothing — added a graph error for
  dict configs where every slot is null, with a clear "use configuration: null"
  message.*)*
- **Medium — the wall lint was unverified.** The self-tests imported only
  `_check_inheritance_and_cycles`, so the AC-Stark hard rule and the
  input/benchmark default-kind lints were untested despite the README/CI claiming
  wall coverage. Extracted the lint into a testable `_check_kind_wall(...)` and
  added tests for: AC-Stark-as-input (hard error), AC-Stark-as-benchmark (clean),
  benchmark-by-default-marked-input (warn), input-by-default-marked-benchmark
  (warn), and a clean record.
- **Env caveat.** A system Python (Anaconda 3.9 + jsonschema 3.2.0) fails before
  collection. Improved the dependency import to emit a clear "jsonschema >= 4.18
  / Python >= 3.10" message, and added a venv note to the README.

Result: `pytest` 8 → **18 passing**; substrate validation still green (exit 0,
1 expected `doerr2024` warning).

---

## 2026-06-18 (later) — Source corpus expansion, consistent PDF naming, primary hyperfine reference

**Trigger.** UW added a reference PDF for the hyperfine constant (Itano &
Wineland 1981) plus a set of key journal papers, and asked to rename the PDFs
consistently.

**Identification (verified).** A multi-agent workflow (identify → adversarial
verify per PDF, 27 agents) identified all 13 added journal-paper PDFs by reading
each PDF first page AND cross-checking Crossref. Every verifier agreed. Results
registered in `sources.yaml` (11 new keys; `clos2016` already present):
`friedenauer2006/2008`, `schmitz2009`, `schneider2012` (review), `clos2014`,
`clos2016_suppmat`, `wittemer2018/2019_prl/2020`, `hasse2024`, `colla2025`. All
`verified: true` (PDF + DOI). Schema change: added `review` to the sources
`degree` enum.

**Naming convention adopted: `filename == citation key`.** All 25 local PDFs
renamed to `<key>.pdf`, so a record's `source.ref` maps 1:1 to its local file
(see `sources/pdf/README.md`). Notable: `wittemer2019` (PhD thesis) vs
`wittemer2019_prl` (PRL 123, 180502, same year) disambiguated; `enderlein2013`
corrects the old `..._2012` filename.

**Primary hyperfine reference wired in.** `hyperfine_a_constant_25mg` rewritten
from the earlier derived estimate (|A| = Δ/3 from the coarse Clos splitting) to
the PRIMARY value from Itano & Wineland: **A = -596.254376(54) MHz**, quoted from
the Abstract and Eq. (2); `observation_type: fitted`, `status: confirmed`. Added
`hyperfine_splitting_calc_25mg` = 3|A| (derived) to keep the derived/inheritance
demonstration on real numbers.

**Integrity catches.**
- The extraction agent's uncertainty was off by ×1000 (it reported 54 kHz); the
  printed "(54)" on `-596.254376` is **54 Hz**. Corrected on intake.
- **TENSION (flagged for reconciliation):** the derived zero-field splitting
  3|A| = 1788.763 MHz sits ~69 kHz BELOW the in-house measured
  `clock_transition_25mg` (Doerr, 1788.8322 MHz). Either the Doerr value is at
  finite field (not truly zero-field) or it is mis-extracted. This is precisely
  the input-vs-benchmark discrepancy the wall exists to surface — noted in both
  records' caveats.

**Validator.** Green: `pytest` 8/8; `python validator/validate.py` exit 0, 13
records, 1 warning (`doerr2024` still has no resolvable identifier). The Itano
warning cleared (now `verified: true`).

**Next steps (added).**
- [x] Reconcile the ~69 kHz clock-transition tension — RESOLVED (2026-06-18) as
  the second-order Zeeman shift at B ≈ 5.6 G; see the "later 3" entry above.
- [ ] The 11 new journal papers are registered but not yet cited by any record;
  mine them for further input/benchmark values (e.g. beam waists/powers,
  decoherence rates from wittemer2018/clos2016).

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
