# Source catalog

A curated catalog of every source registered in
[`registries/sources.yaml`](../registries/sources.yaml), with bibliographic
metadata, an own-words summary, and the **figures of merit (FOMs)** each source
contributes to the digital twin. It is a reading companion to the registry (the
machine-readable truth) and to [`docs/LOGBOOK.md`](LOGBOOK.md) (how each value
was extracted).

**How to read this.** Sources fall into three tiers:

| tier | meaning |
|------|---------|
| **Active** | currently cited by ≥1 ledger record — supplies inputs/benchmarks the twin consumes or reproduces |
| **Registered (lineage)** | in the corpus per the task card, resolvable, but not yet feeding a record |
| **Internal** | local PAULA analysis notebooks (own primary sources, no public identifier) |

**Provenance flags** (from the registry): `verified` = the DOI/permalink/URN was
confirmed to resolve; `archived: false` = no public copy found. The validator
WARNs when a referenced source lacks a resolvable identifier (`link`/`doi`/`urn`)
or is not marked `verified`. Most extracted values remain `status: provisional`;
confirmed records are marked explicitly in the ledger.

**Apparatus generations** (see [`registries/generations.yaml`](../registries/generations.yaml)):
`legacy` (Clos/Wittemer era) → `hasse` → `freddy` (Doerr/Weber era). A figure of
merit is only comparable within or across generations when the calibration is
known to carry over — flagged per record.

---

## Digital-twin figures of merit at a glance

The twin's headline FOMs are the **7 σ-validations** — engine predictions (from
`input` records) checked against `benchmark` records across the input/benchmark
wall. All currently agree within 3σ:

| benchmark (FOM) | engine | source | n_σ | what it tests |
|-----------------|--------|--------|-----|---------------|
| `clock_transition_25mg` = 1788.8322 MHz | levels | doerr2024 | 1.09 | Breit-Rabi clock from B(freddy) + A/I/g |
| `clock_transition_weber_25mg` = 1788.8328 MHz | levels | weber2025 | 0.51 | same, from the **measured** Weber field |
| `omega_z_axial_stretch_2ion` = 2.23 MHz | modes | wittemer2019 | 2.17 | √3·COM axial stretch (calc cross-check) |
| `omega_radial_rocking_2ion` = 2.57 MHz | modes | wittemer2019 | 0.01 | radial rocking from COM frequencies |
| `doppler_cooling_limit_25mg` = 1.00 mK | cooling | doerr2024 | 0.03 | T_D = ℏΓ/2k_B (theory consistency) |
| `doppler_cooled_occupation_25mg` = n̄ 10(1) | cooling | clos2017 | 0.42 | **measured** Doppler-cooled occupation |
| `bdd_ac_stark_shift_25mg` = ~10 MHz | acstark | hasse2025 | 0.15 | **measured** BDD light shift |

Three are independent **measurements** (Weber clock, Clos n̄, Hasse BDD); the rest
are theory/consistency cross-checks. The key **input anchors** are the hyperfine
constant A (itano_wineland_1981), the laser table (clos2017), the Raman geometry
(doerr2024), and the Lamb-Dicke parameter η=0.32 (wittemer2019).

---

## Master index

| key | author | year | type | resolves | records |
|-----|--------|------|------|:--------:|:-------:|
| **clos2017** | G. Clos | 2017 | PhD | ✅ | 23 |
| **doerr2024** | F. Dörr | 2024 | MSc | ⚠️ local | 20 |
| **itano_wineland_1981** | Itano & Wineland | 1981 | article | ✅ | 6 |
| **wittemer2019** | M. Wittemer | 2019 | PhD | ✅ | 5 |
| **kaufmann2022** | I. Kaufmann | 2022 | MSc | ⚠️ local | 5 |
| **weber2025** | A. Weber | 2025 | MSc | ⚠️ local | 2 |
| **thomm2021** | R. Thomm | 2021 | MSc | ⚠️ local | 8 |
| **friedenauer2010** | A. Friedenauer | 2010 | PhD (LMU) | ✅ | 1 |
| **hasse2025** | F. Haße | 2025 | PhD | ✅ | 1 |
| **codata2018** | Tiesinga et al. | 2021 | dataset | ⚠️ unverified | 1 |
| **stone2005** | N. J. Stone | 2005 | review | ⚠️ unverified | 1 |
| enderlein2013 | M. Enderlein | 2013 | PhD | ✅ | 0 |
| kalis2017 | H. Kalis | 2017 | PhD | ⚠️ local | 0 |
| **kalis2016** | Kalis et al. | 2016 | article | ✅ | 0 |
| clos2016 | Clos et al. | 2016 | article | ✅ | 0 |
| palani2023 | Palani et al. | 2023 | article | ✅ | 0 |
| friedenauer2006 | Friedenauer et al. | 2006 | article | ✅ | 0 |
| friedenauer2008 | Friedenauer et al. | 2008 | article | ✅ | 0 |
| schmitz2009 | Schmitz et al. | 2009 | article | ✅ | 0 |
| schneider2012 | Schneider et al. | 2012 | review | ✅ | 0 |
| clos2014 | Clos et al. | 2014 | article | ✅ | 0 |
| clos2016_suppmat | Clos et al. | 2016 | suppl. | ✅ | 0 |
| wittemer2018 | Wittemer et al. | 2018 | article | ✅ | 0 |
| wittemer2019_prl | Wittemer et al. | 2019 | article | ✅ | 0 |
| wittemer2020 | Wittemer et al. | 2020 | article | ✅ | 0 |
| hasse2024 | Hasse et al. | 2024 | article | ✅ | 0 |
| colla2025 | Colla et al. | 2025 | article | ✅ | 0 |
| schmitz2010 | H. Schmitz | 2010 | PhD (LMU) | ❌ | 0 |
| matjeschk2008 | R. Matjeschk | 2008 | Diplom | ❌ | 0 |
| pacher2014 | J. Pacher | 2014 | Diplom | ❌ | 0 |
| harlos2015 | J. Harlos | 2015 | MSc | ❌ | 0 |
| paula_mg_details_2014 | PAULA group | 2014 | notebook | 🔒 local | 0 |
| paula_mg_scatterrate_2015 | PAULA group | 2015 | notebook | 🔒 local | 0 |
| paula_bem_shims_2017 | PAULA group | 2017 | notebook | 🔒 local | 0 |
| paula_mode_orientations_2017 | PAULA group | 2017 | notebook | 🔒 local | 0 |
| paula_trapsim_2018 | PAULA group | 2018 | notebook | 🔒 local | 0 |

✅ public DOI/permalink resolves · ⚠️ local PDF only or identifier not re-verified ·
❌ no resolvable copy found · 🔒 internal/local-only

---

# Part A — Active sources (supply twin figures of merit)

### clos2017 — *Trapped atomic ions for fundamental studies of closed and open quantum systems*
- **Govinda Clos**, PhD, Albert-Ludwigs-Universität Freiburg, 2017
- DOI [10.6094/UNIFR/12400](https://doi.org/10.6094/UNIFR/12400) · [FreiDok](https://freidok.uni-freiburg.de/data/12400) · ✅ verified · thesis behind PRL 117, 170401 (clos2016)

**Summary.** The foundational characterization of the Freiburg single-/few-²⁵Mg⁺
apparatus (legacy generation). It carries the canonical laser-parameter table
(**Tab. 3.2**, "Specifications of all laser systems") — per beam: wavelength,
polarisation, transition, natural linewidth, detuning Δ/Γ, and saturation I/Isat —
which is the origin of the twin's entire cooling/detection/preparation beam layer
(BD, BDD [Clos labels it "BDdet"], BDX, RD, RP). It also develops the two-photon-Raman / Lamb-Dicke
formalism (§3.4), the AC-Stark light-shift theory (Eq. 2.2.24, δ_AC = ℏΩ²/4δ),
and reports the **measured Doppler-cooled mean phonon number n̄ = 10(1)**
(Fig. 3.14) at ω₁/2π = 1.915 MHz plus the ~1 mK Doppler limit.

**FOMs for the twin** (23 records — all `input` except where noted):
- Beam layer: `mg_p32_natural_linewidth` 41.8 MHz, `mg_p12_natural_linewidth` 41.3 MHz, `bd/rd_laser_wavelength` 279.635/280.353 nm, `bd_detection_waist` ~50 µm, and the detuning + saturation of all five beams (`bd/bdd/bdx_*`, `rd/rp_repump_*`).
- Field: `b_field_quantization` 5.85 G (legacy).
- Cooling: `omega_z_axial_clos_25mg` 1.915 MHz (input) → **`doppler_cooled_occupation_25mg` n̄=10(1)** *(benchmark)*.
- Internal structure: nuclear spins of ²⁴Mg/²⁶Mg (=0), hyperfine levels F=2/3, `hyperfine_splitting_25mg_f2_f3` 1.79 GHz.

### doerr2024 — *Advanced Interferometer Techniques for Measuring Near-Resonant Light Shifts and Superresolving Trapped-Ion Dynamics*
- **Frederike Dörr**, MSc, Albert-Ludwigs-Universität Freiburg, 2024
- ⚠️ not publicly archived; author/year confirmed via [RTG-DynCAM profile](https://rtg-dyncam.de/members-new/members/frederike-doerr/); title confirmed from the thesis title page; local PDF only · defines the **freddy** generation

**Summary.** Master thesis on interferometric measurement of near-resonant
(microwave AC-Zeeman) light shifts and superresolved trapped-ion dynamics on the
current apparatus. In characterizing the setup it defines the Raman two-photon
(TPSR) beam geometry — the four combinations **CC/OC/AC/ROC** and their effective
k-vectors — the **30° radial-mode tilt**, the single-ion mode spectrum (lf 1.3 /
mf 3.0 / hf 4.5 MHz), the 20 GHz Raman single-photon detuning, the freddy
quantization field (~5.5 G), the clock transition (1788.8322 MHz), and a set of
measured microwave Rabi rates on 8 ground-state hyperfine transitions.

**FOMs for the twin** (20 records):
- Raman geometry (`input`): `raman_detuning_from_p32` 20 GHz, the four `raman_*_combination_25mg` Δk/k vectors, `radial_mode_tilt_25mg` 30°, `omega_radial_mf/hf_25mg` 3.0/4.5 MHz.
- Field (`input`): `b_field_quantization_freddy` ~5.5 G → drives the clock prediction.
- Benchmarks: **`clock_transition_25mg`** 1788.8322 MHz with its `clock_transition_residual_25mg` (+2.7(24) kHz), `doppler_cooling_limit_25mg` ~1 mK, and the **8 `mw_rabi_*_doerr`** measured microwave Rabi rates (drive diagnostic).

### itano_wineland_1981 — *Precision measurement of the ground-state hyperfine constant of ²⁵Mg⁺*
- **W. M. Itano & D. J. Wineland**, Phys. Rev. A **24**, 1364 (1981), NBS Boulder
- DOI [10.1103/PhysRevA.24.1364](https://doi.org/10.1103/PhysRevA.24.1364) · ✅ verified

**Summary.** The primary precision measurement of the ²⁵Mg⁺ ground-state
magnetic-dipole hyperfine constant **A = −596.254376(54) MHz**. This single atomic
constant, with the nuclear spin I = 5/2, is the backbone of the twin's `levels`
(Breit-Rabi) engine — every clock-transition and Zeeman prediction descends from it.

**FOMs for the twin** (6 records): `hyperfine_a_constant_25mg` −596.254376 MHz,
`nuclear_spin_25mg` 5/2, the derived `hyperfine_splitting_calc_25mg` 3|A| =
1788.76 MHz, `qubit_quadratic_zeeman_coeff_25mg` 2.195×10¹¹ Hz/T², the predicted
clock, and `b_field_from_clock_25mg` (field back-inferred from the clock, benchmark).

### wittemer2019 — *Particle creation and memory effects in a trapped-ion quantum simulator*
- **Matthias Wittemer**, PhD, Albert-Ludwigs-Universität Freiburg, 2019
- DOI [10.6094/UNIFR/151582](https://doi.org/10.6094/UNIFR/151582) · ✅ verified

**Summary.** PhD on open-system / particle-creation analogues in the trapped-ion
simulator (legacy generation). It supplies the **axial Lamb-Dicke parameter
η = 0.32** at ω/2π = 1.920 MHz — the anchor for the twin's `sideband` engine — and
the single-ion secular frequencies (Tab. 3.2). Its mode table also lists the
calculated (Wübbena 2012) 2-ion stretch/rocking frequencies used as `modes`-engine
cross-checks. Also confirms the laser table independently (its Tab. 3.1).

**FOMs for the twin** (5 records): `raman_axial_lamb_dicke_25mg` 0.32, `omega_z_axial_com_25mg`
1.30 MHz, `omega_radial_com_25mg` 2.88 MHz (inputs); `omega_z_axial_stretch_2ion_25mg`
2.23 MHz and `omega_radial_rocking_2ion_25mg` 2.57 MHz (calculated benchmarks).

### kaufmann2022 — MSc thesis (microwave Rabi calibration)
- **Ingolf Kaufmann**, MSc, ~2022 · ⚠️ not publicly archived; name corroborated via a DPG-2022 (Erlangen) abstract; local PDF only

**Summary.** Master thesis from the Freiburg multi-ion-trap context providing a
second, independent set of measured microwave Rabi rates on the ²⁵Mg⁺
ground-state hyperfine transitions (legacy generation). Used alongside Doerr's set
in the twin's drive diagnostic to expose the apparatus-dominated (antenna) spread
that the atomic Clebsch-Gordan model alone cannot explain.

**FOMs for the twin** (5 records, all benchmarks): `mw_rabi_*_kaufmann` — five
measured microwave Rabi rates (3,3↔2,2; 3,1↔2,2; 3,−1↔2,−1; 3,−1↔2,−2; 3,−3↔2,−2).

### weber2025 — *A Tunable Quantum Magnetometer Based on Single Trapped Ions*
- **Andreas Weber**, MSc, Albert-Ludwigs-Universität Freiburg, 2025 (supervisor T. Schätz)
- ⚠️ not publicly archived; local PDF only · freddy generation

**Summary.** Master thesis using single ²⁵Mg⁺ ions as a tunable magnetometer.
Its **Table 3 (p. 29)** reports 8 *measured* ground-state |3,m_F⟩↔|2,m′_F⟩ qubit
transition frequencies, from which a Zeeman-spectroscopy field **B = 5.6454(9) G**
(Eq. 4.2) is determined and the hyperfine constant A is confirmed. This gives the
twin its second, *independent measured* clock benchmark and a precise operating field.

**FOMs for the twin** (2 records): `b_field_zeeman_weber_25mg` 5.6454 G (input) →
**`clock_transition_weber_25mg`** 1788.8328 MHz (benchmark, the levels engine's
independent measured test).

### kalis2016 — *Motional-mode analysis of trapped ions*
- **H. Kalis, F. Hakelberg, M. Wittemer, M. Mielenz, U. Warring, T. Schätz**, Phys. Rev. A **94**, 023401 (2016)
- DOI [10.1103/PhysRevA.94.023401](https://doi.org/10.1103/PhysRevA.94.023401) · ✅ verified · the group's own published method

**Summary.** The group's own published **tickle / secular-excitation** method for
measuring motional-mode frequencies (and orientations). A finite resonant excitation
pulse (duration t_exc) drives a mode as a classical oscillator → amplitude
A ∝ sin([ω_exc−ω_i]t_exc/2)/(ω_exc²−ω_i²) (Eq. 3, a **sinc** of FWHM ≈ 1/t_exc, *not*
a Lorentzian); the Doppler modulation (index β = ⟨u,k_BD⟩A) shifts carrier population
to motional sidebands ±vω_i with Bessel weights J_v(β)², giving the normalized
fluorescence F = ∏_i Σ_v J_v(β_i)²·g(Δ_BD+vω_i) (Eq. 2, |v|≤15). It is the canonical
reference behind [`spike/engines/tickle.py`](../spike/engines/tickle.py); the
[kalis2017](#kalis2017--initialization-of-quantum-states-in-a-two-dimensional-ion-trap-array)
thesis adds Gaussian mode-frequency noise (~1 kHz) for long t_exc.

**FOMs for the twin** (0 records). Applied to the PAULA `Tickle/PDQ_*_FScan` data
(`python -m spike.plot_tickle`, significance-gated over all scan files): axial
**1.299 MHz** (ledger nominal 1.30, 10/10 files), radial **3.224 MHz** (+7.5%, 3/4)
and **4.712 MHz** (+4.7%, 4/4) — the radials sit well above the nominal inputs
(an unresolved refinement; not yet recorded). Per-mode narrow "calibration" scans
that miss the mode are correctly rejected by the engine's F-test + edge guard.

### thomm2021 — *State Detection of Trapped Magnesium Ions*
- **Robin Thomm**, MSc, Albert-Ludwigs-Universität Freiburg, June 2021 (supervisor T. Schätz)
- ⚠️ not publicly archived; local PDF only · legacy generation

**Summary.** The Freiburg/PAULA **state-detection reference** for ²⁵Mg⁺ — exactly the
source for building out *electronic-state readout* and *motional effects*. It (a)
establishes the **electronic readout**: photon counts follow weighted Poissonians,
P↓ is estimated by **maximum likelihood** (2-Poissonian single ion, 3-Poissonian two
ions), with an optimised **t_det = 30 µs** and measured bright/dark fidelities; and
(b) demonstrates **motional-state engineering** on a single ion — resolved-sideband
ground-state cooling of all three modes, plus Fock, coherent-displaced and squeezed
states read out from carrier/blue-sideband flops and reconstructed as **Wigner
functions**. These are the *actual apparatus* numbers; they supersede the LMU
friedenauer2010 detection placeholders.

**FOMs for the twin** (8 records this pass; more catalogued below):
- *Electronic readout* (`detection`): **`mg_detection_time_25mg`** 30 µs (input);
  **`mg_bright_counts_25mg`** λ↓=2.682(9) / **`mg_dark_counts_25mg`** λ↑=0.036(3)
  counts per 30 µs (inputs, single ion → ~89 kHz bright rate);
  **`mg_readout_fidelity_bright_25mg`** 99.4(6)% / **`mg_readout_fidelity_dark_25mg`**
  97.4(6)% (benchmarks; two-ion 99(1)%/96(1)%). The fidelities are **ML
  state-probabilities**, not single-shot thresholds (a hard threshold on these means
  gives ~93% bright — recorded in the caveats).
- *Motional ground state* (`motion`): **`mg_rsb_cooled_nbar_axial_lf_25mg`** 0.07(2),
  **`mg_rsb_cooled_nbar_radial_mf_25mg`** 0.11(3), **`mg_rsb_cooled_nbar_radial_hf_25mg`**
  0.07(4) (benchmarks; P(n=0)=94/90/94%, 3D ground state 79(4)%, T≈26 µK) — the RSB
  counterpart to the Doppler `doppler_cooled_occupation_25mg` (n̄≈10).

**Further FOMs catalogued (not yet records)** — ready for the next build:
- Mode frequencies axial 1.3 / radial 2.7 / 4.4 MHz; RF drive 56 MHz; B≈5.5 G;
  hyperfine ω₀≈1774 MHz; BD cooling 20 MHz (=Γ/2), T_Doppler≈1 mK — corroborate
  existing `omega_*`, `clock_*`, `bd_*` records.
- **Lamb-Dicke η=0.39 at 1.3 MHz** — *consistent* with `raman_axial_lamb_dicke_25mg`
  0.32 at 1.92 MHz via η∝ω^(−1/2) (0.32·√(1.92/1.30)=0.389); noted as a cross-check,
  not a conflict.
- *Motional-state engineering* (for a future Wigner/tomography engine): Fock |1⟩
  89(3)% / |2⟩ 76(2)%; coherent displacement calib. 0.209(4)/µs (Raman),
  0.193(5)/µs (E-field), |α|=0.67–0.84; squeezing |ξ|=r=0.74(3) (~1.3 dB); Fock
  population fit Ωₙ,ₙ₊ₛ ∝ e^(−η²/2) Lₙ^|s|(η²), Hilbert truncation n_max=18; Fock
  |1⟩ Wigner negativity confirmed.
- *Raman/apparatus*: UV power 120–140 mW (post-BBO), 5–10 mW/beam; fine-structure
  splitting 2.7 THz; active MW↔Raman phase stabilisation (1 kHz servo, <1 ms
  settling, <5 ms coherence per shot).

**✅ Detuning history (resolved, UW 2026-06-19).** The Raman virtual-level detuning
from P3/2 |F=4,m_F=4⟩ evolved **50 → 200 GHz** (Thomm 2021, new BBO cavity) **→ 20 GHz**
(doerr2024 = current). The final reduction was deliberate — a new laser system + BBO
doubling-stage trouble, traded a smaller detuning for **larger Rabi rates** (Ω∝1/Δ_R).
So `raman_detuning_from_p32` = 20 GHz stands as the current value (Thomm's 50/200 are
earlier epochs); recorded in that record's caveats.

### friedenauer2010 — *Simulation of the Quantum Ising Model in an Ion Trap*
- **Axel Friedenauer**, PhD, **LMU München / MPQ Garching** (NOT Freiburg), 2010
- DOI [10.5282/edoc.11595](https://doi.org/10.5282/edoc.11595) · URN urn:nbn:de:bvb:19-115958 · ✅ verified

**Summary.** PhD behind Nature Physics 4, 757 (2008) (first quantum-magnet
simulation), on the earlier **LMU/MPQ** ²⁵Mg⁺ apparatus — a *different, predecessor*
setup. Only apparatus-independent atomic quantities transfer to Freiburg; here, the
²⁵Mg⁺ **saturation intensity Isat ≈ 255 mW/cm²** (an atomic constant). Its 30 µm BD
waist is LMU-specific and is *not* used as a Freiburg value.

**FOMs for the twin** (1 record): `mg_saturation_intensity` 2550 W/m² (= 255 mW/cm²).

### hasse2025 — *Observation of dynamic processes demonstrated in a trapped-ion quantum simulator*
- **Florian Haße**, PhD, Albert-Ludwigs-Universität Freiburg, 2025 (supervisor T. Schätz)
- DOI [10.6094/UNIFR/274764](https://doi.org/10.6094/UNIFR/274764) · ✅ verified · defines the **hasse** generation

**Summary.** Recent PhD on observing dynamic processes (phase-stable travelling
waves; superresolution). It quantifies the **measured AC-Stark shift of the
far-detuned BDD beam (~2π × 10 MHz)** and states that the resonant repumpers RD/RP
induce *no significant* shift — the experimental basis for the twin's far-detuned vs
near-resonant AC-Stark scoping. Its glossary independently corroborates Doerr's
CC/OC/AC/ROC TPSR beam-combination naming.

**FOMs for the twin** (1 record): **`bdd_ac_stark_shift_25mg`** ~10 MHz (benchmark,
the acstark engine's measured test).

### codata2018 — *CODATA Recommended Values of the Fundamental Physical Constants: 2018*
- **Tiesinga, Mohr, Newell, Taylor**, Rev. Mod. Phys. **93**, 025010 (2021) + [NIST portal](https://physics.nist.gov/cuu/Constants/)
- DOI [10.1103/RevModPhys.93.025010](https://doi.org/10.1103/RevModPhys.93.025010) · ⚠️ not re-resolved this session

**Summary.** The reference fundamental constants — Planck constant h, Bohr/nuclear
magnetons μ_B/μ_N, Boltzmann constant k_B, and the free-electron g-factor g_e
(used to approximate g_J of the ²S₁/₂ state). Hard-coded in `spike/constants.py`
and registered here as a record for the electronic g-factor.

**FOMs for the twin** (1 record): `g_factor_electron_2s12` 2.00232 (≈ g_e).

### stone2005 — *Table of nuclear magnetic dipole and electric quadrupole moments*
- **N. J. Stone**, At. Data Nucl. Data Tables **90**, 75 (2005)
- DOI [10.1016/j.adt.2005.04.001](https://doi.org/10.1016/j.adt.2005.04.001) · ⚠️ not re-resolved this session

**Summary.** Standard compilation of nuclear moments — supplies the ²⁵Mg nuclear
magnetic moment μ = −0.85545 μ_N, from which the nuclear g-factor g_I = −0.34218
is derived for the Breit-Rabi engine.

**FOMs for the twin** (1 record): `g_factor_nuclear_25mg` −0.34218.

---

# Part B — Registered lineage corpus (not yet feeding the twin)

Corpus members registered per the task card, available to source future records.
Freiburg-group lineage unless noted; several are deliberately flagged as having
no public copy.

| key | citation | role |
|-----|----------|------|
| **enderlein2013** | PhD, Freiburg, FreiDok 8886 (URN) | Optical ion trapping; lineage apparatus thesis |
| **kalis2017** | PhD, Freiburg, 2017 — ⚠️ local PDF | **DAQ system + `.dat` data-file FORMAT** only (see [DATA_FORMAT.md](DATA_FORMAT.md)); the example data are the group's own measurements (2026-06-15), *not* from this thesis |
| **clos2016** | PRL **117**, 170401 (2016) | Thermalization in an isolated quantum system (clos2017 paper) |
| **clos2016_suppmat** | Suppl. to PRL 117, 170401 | Supplemental methods/parameters |
| **clos2014** | PRL **112**, 113003 (2014) | Decoherence-assisted spectroscopy of a single Mg⁺ ion |
| **palani2023** | PRA **107**, L050601 (2023) | High-fidelity ion transport (multi-layer array) |
| **friedenauer2006** | Appl. Phys. B **84**, 371 (2006) | The 280 nm all-solid-state laser source for Mg⁺ |
| **friedenauer2008** | Nature Phys. **4**, 757 (2008) | First trapped-ion quantum-magnet simulation (LMU) |
| **schmitz2009** | PRL **103**, 090504 (2009) | Quantum walk of a trapped ion in phase space |
| **schneider2012** | Rep. Prog. Phys. **75**, 024401 (2012) | Review: quantum simulations with trapped ions |
| **wittemer2018** | PRA **97**, 020102(R) (2018) | Quantum memory effects and limitations |
| **wittemer2019_prl** | PRL **123**, 180502 (2019) | Phonon-pair creation by inflating fluctuations |
| **wittemer2020** | Phil. Trans. R. Soc. A **378**, 20190230 (2020) | Harmonic-oscillator toolkit under extreme conditions |
| **hasse2024** | PRA **109**, 053105 (2024) | Phase-stable travelling waves / superresolution |
| **colla2025** | Nat. Commun. **16**, 2502 (2025) | Time-dependent level renormalisation, ultrastrong coupling |
| schmitz2010 | PhD, LMU, ~2010 — ❌ no public copy | LMU predecessor thesis (print-only) |
| matjeschk2008 | Diplom, 2008 — ❌ no public copy | Co-author PRL 103, 090504 |
| pacher2014 | Diplom, 2014 — ❌ unverified (spelling) | Cited for a Doppler-cooling detail; needs confirmation |
| harlos2015 | MSc, 2015 — ❌ no public copy | Corpus member, title unconfirmed |

---

# Part C — Internal PAULA notebooks (own primary sources, local-only)

Mathematica notebooks kept under `sources/Mathematica/` (untracked). Citable as the
origin of simulated/derived quantities; **no public identifier** → the validator
WARNs when one is referenced (an accepted accessibility flag).

| key | notebook | likely role |
|-----|----------|-------------|
| paula_mg_details_2014 | `Notes_PAULA_Mg_details_2014_07_10.nb` | Mg level-structure details |
| paula_mg_scatterrate_2015 | `Notes_PAULA_Mg_streurate_2015_07_17.nb` | Mg scattering-rate (Streurate) notes |
| paula_bem_shims_2017 | `TrapSim_PAULA_BEM_opt_ExpZone_Shims_2017_05_30.nb` | BEM shim-field optimisation |
| paula_mode_orientations_2017 | `TrapSim_PAULA_Mode_orientations_2017_07_11.nb` | Motional-mode orientations (cf. the 30° tilt) |
| paula_trapsim_2018 | `Notes_PAULA_TrapSimulations_2018_10_31.nb` | Trap simulations |

---

## Open provenance items

- **doerr2024 / kaufmann2022 / weber2025** are not publicly archived (local PDF
  only). Their values carry the highest accessibility risk under FAIR; archiving or
  a citable secondary record would close it. The validator emits these as WARNs.
- **codata2018 / stone2005** DOIs were not re-resolved in the latest session
  (flagged `verified: false`); the values themselves are standard.
- **pacher2014** spelling/existence unconfirmed — do not cite a record to it
  without verifying the thesis.
- The Γ(P₃/₂) linewidth differs by source (Clos 41.8 MHz vs Hasse 42.7 MHz); the
  twin uses 41.8 MHz and notes the spread where it matters (BDD AC-Stark).
