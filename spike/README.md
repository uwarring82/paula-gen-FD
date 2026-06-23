# Twin spike

A throwaway-ish, end-to-end proof that the structured parameter layer can drive
**physics engines**. Each engine consumes `input` records from the
[ledger](../records/) and is checked against `benchmark` records — the
input/benchmark wall, executed.

> **Boundary.** Per the task card, solver/physics code is *out of scope* for the
> substrate (`records/` + `validator/`). This package is kept separate and may
> later graduate to its own repo(s) — `iontrap-levels` / `-modes` / `-optics`
> plus a twin composition root. Pure Python, no numpy: the `2S_1/2` Breit-Rabi
> problem is closed-form, and the small normal-mode linear algebra is in
> [`linalg.py`](linalg.py).

## Composition root

[`runner.py`](runner.py) is the twin's composition root: it runs every
registered engine validation against its benchmark and renders one table.

```
$ python -m spike.validate_twin
TWIN VALIDATION — engines reproduce benchmarks from ledger inputs

benchmark                        engine   subsystem       predicted        reference        residual    n_sigma  status
clock_transition_25mg            levels   internal_state  1788.829549 MHz  1788.832200 MHz  -2.65 kHz   1.09     ok
clock_transition_weber_25mg      levels   internal_state  1788.833107 MHz  1788.832800 MHz  +0.31 kHz   0.51     ok
omega_z_axial_stretch_2ion_25mg  modes    motion          2.251666 MHz     2.230000 MHz     +21.67 kHz  2.17     ok
omega_radial_rocking_2ion_25mg   modes    motion          2.569903 MHz     2.570000 MHz     -0.10 kHz   0.01     ok
doppler_cooling_limit_25mg       cooling  motion          1.0030 mK        1.0000 mK        +3.04 uK    0.03     ok
doppler_cooled_occupation_25mg   cooling  motion          10.42            10               +0.421      0.42     ok
bdd_ac_stark_shift_25mg          acstark  optics          10.450000 MHz    10.000000 MHz    +450.00 kHz 0.15     ok

7 validation(s): 7 ok, 0 not ok (threshold 3 sigma).
```

It exits nonzero if any result is in tension (> 3σ) or errors (e.g. an engine
tries to consume a non-input), and flags any measured benchmark that no engine
covers. **Adding an engine** is: implement it, write one `_validate_*` returning
a `ValidationResult`, and add it to `REGISTRY`.

## Engine: `levels` (Breit-Rabi)

[`engines/levels.py`](engines/levels.py) models the ²⁵Mg⁺ ground-state
(3s ²S₁/₂) hyperfine + Zeeman structure from the magnetic-dipole constant `A`,
nuclear spin `I`, electronic/nuclear g-factors, and the field `B`. Exact for
`J = 1/2` and **sign-safe**: `A < 0` (the inverted ²⁵Mg⁺ manifold, F=3 below F=2)
enters only through `A²` and diagonal `m_J m_I` terms. Reproduces the measured
`clock_transition_25mg` benchmark (residual −2.65 kHz = 1.09σ).

## Engine: `modes` (axial normal modes)

[`engines/modes.py`](engines/modes.py) computes the axial normal modes of an
N-ion chain from the COM secular frequency: `omega_p = sqrt(lambda_p) * omega_z`
with `lambda_p` the eigenvalues of the dimensionless axial Hessian at the ion
equilibrium positions (James 1998; `lambda_1 = 1` COM, `lambda_2 = 3` stretch).
Reproduces the 2-ion axial stretch (`omega_z_axial_stretch_2ion_25mg`):
predicted √3·COM = 2.2517 MHz vs measured 2.23 MHz, +21.7 kHz = 2.17σ (the √3
relation holds to ~1% in the 3-sig-fig data).

## Engine: `drive` (microwave Rabi couplings)

[`engines/drive.py`](engines/drive.py) predicts the *relative* magnetic-dipole
Rabi couplings of the ground-state |F,m_F⟩↔|F′,m_F′⟩ transitions from the
Clebsch-Gordan coefficients (pure-Python Racah formula; verified against sympy).
This is the **atomic** part only: the measured microwave Rabi rates
([`records/control.yaml`](../records/control.yaml), from Doerr/Kaufmann) are
*apparatus-dominated* — the MW antenna polarization + frequency response swing
them ~6× across the manifold (mirror pairs with equal |CG| differ 5×). So the
runner reports a **drive DIAGNOSTIC** (|CG| vs measured → apparatus factor), not
a σ-validation; an absolute-rate engine would need an antenna model.

## Engine: `projection` (Raman → motional mode)

[`engines/projection.py`](engines/projection.py) predicts which motional mode
each Raman (TPSR) beam combination addresses, from the effective k-vector
direction (Δk = k_B − k_R) and the radial-mode tilt — both consumed from the
ledger (`raman_*_combination_25mg` Δk vectors + `radial_mode_tilt_25mg`). The
coupling onto a mode is the direction cosine |Δk̂·ê_mode|; the axial (lf) mode is
along z, the two radial modes (mf/hf) lie in the x-y plane tilted 30°. The
geometry **reproduces Doerr 2024's documented addressing** (CC carrier-only, OC
axial, AC all three at 45°, ROC radial) — shown as a **projection DIAGNOSTIC**
(combination × mode, with a ✓/✗ against Doerr). It is the *geometric* part only;
absolute η (|Δk| + mode frequency) is a future sideband engine.

## Engine: `sideband` (absolute Lamb-Dicke + sideband Rabi)

[`engines/sideband.py`](engines/sideband.py) turns the `projection` direction
cosines into *absolute* Lamb-Dicke parameters: η(comb, mode, ω) = |Δk/k·ê_mode| ·
k·z̄(ω), with z̄ ∝ ω^(−1/2), anchored to the measured η = 0.32 (OC→lf at 1.92 MHz).
It gives the first-order sideband Rabi rates Ω_{n,n±1} = η√(n+1|n)·Ω₀ and the
**Raman differential AC-Stark shift** that moves the sideband resonance,
δ_AC,diff ≈ (ω_HF/Δ_R)·Ω₀ ≈ 0.09·Ω₀ (the |↓⟩,|↑⟩ light shifts differ by the
hyperfine splitting / Raman detuning). Capability + diagnostic (no independent
benchmark beyond the η anchor; the differential shift has no measured value). It
also carries the **carrier Debye-Waller thermal dephasing**: even a Δn=0 carrier is
motion-sensitive at finite η, with Ω_{n,n} = Ω₀·e^(−η²/2)·L_n(η²) (`carrier_rabi_
factor`, Laguerre L_n), so a thermal state (mean n̄, `thermal_pn`) is a spread of
Rabi frequencies and the flop **dephases** — `thermal_carrier_flip` is the exact
flip-probability sum and `thermal_coherence` its envelope |Σ P_n e^{iΩ_n t}| (the
cheap inversion handle). This is the leading *motional* contribution to a Raman
carrier-flop envelope, ledger-anchored from η + n̄.

## Engine: `acstark` (far-detuned light shift)

[`engines/acstark.py`](engines/acstark.py) is the far-detuned light shift
δ_AC = sΓ²/(8δ) = Ω²/(4δ) (Clos Eq. 2.2.24). It reproduces Hasse's **measured**
BDD shift (~2π×10 MHz → predicted 10.45 MHz, 0.15σ) — a real σ-validation. The
shift is only meaningful **far** from resonance: BDD (−10Γ) shifts coherently,
while the near-resonant BD/BDX/RD/RP *scatter* (cooling engine), which Hasse
confirms ("RD/RP induce no significant ac Stark shift"). The same validated
formula feeds the Raman differential shift in `sideband`.

## Engine: `scatter` (Raman off-resonant scattering + differential AC-Stark)

[`engines/scatter.py`](engines/scatter.py) is the loss budget of a two-photon
stimulated-Raman (TPSR) carrier drive, where both beams sit Δ_R = 20 GHz red of
`3P_3/2` (`raman_detuning_from_p32`). From that single detuning + the linewidth Γ
it gives the coherent rate `Ω = Ω_BΩ_R/(2Δ_R)`, the **off-resonant scattering
rate** `Γ_sc = Γ(Ω_B²+Ω_R²)/(4Δ_R²)` — which for balanced beams collapses to the
clean `Γ_sc = (Γ/Δ_R)·Ω`, so the **spontaneous-emission floor is detuning-only**:
`P_SE(π) = π·Γ/Δ_R ≈ 0.66 %` here — and the **differential AC-Stark shift**
`δ_AC ≈ (ω_HF/Δ_R)·Ω` (delegated to `sideband.raman_differential_stark_factor`).
`flip_probability(t, Ω, Γ_sc, δ_AC)` is the forward flop: AC-Stark speeds Ω_eff and
caps the amplitude at `Ω²/Ω_eff²`, scattering decoheres with the envelope
`e^{-(3/4)Γ_sc t}` (full-depolarisation leading order; `CONTRAST_DECAY_FACTOR`).
`RamanScatter.from_ledger` consumes Δ_R, Γ, ω_HF (all `input`). Capability +
diagnostic — there is **no measured Raman-scattering decoherence rate** to
σ-validate against (cf. the unanchored differential shift in `sideband`); the
tests pin the formula identities, the balanced collapse, and the `1/Δ_R` scaling.

## Engine: `raman_optical` (polarization+power-resolved light shifts & scattering)

[`engines/raman_optical.py`](engines/raman_optical.py) resolves the Raman AC-Stark
shift and off-resonant scattering by **each beam's polarization and power** — both of
which the scalar `scatter` engine ignored. It works in the **fine-structure
|P_J' mJ'> basis** for the excited manifold (at Δ_R = 20 GHz the 3P hyperfine/Zeeman
is unresolved, so F', mF' are *not* good quantum numbers for P), with the ground
|S₁/₂ F mF> decomposed into |mJ, mI> (nuclear spin mI a spectator of the optical
dipole). Each beam carries a (σ⁺,π,σ⁻) intensity decomposition + a relative power;
`light_shift`, `scatter_rate`, and the mI-conserving `two_photon_rabi` sum over
q, J'. The dimensionless ratios `differential_stark_per_rabi` / `scatter_per_rabi`
(δ_AC and Γ_sc per unit two-photon Rabi) are anchored to the **measured** flop, so
the absolute field scale cancels. The **scalar+vector polarizability**
(`scalar_vector_shift`; tensor = 0 for ²S₁/₂) is the transparent summary.

The old |F',mF'>-basis sum (`coupling`, via a Wigner 6j) is retained **only as a
degenerate-limit cross-check**: it equals the |mJ> sum to machine precision (basis
independence — tested). Validated against the cycling transition, the sublevel/F-
independent sum rule, and exact D2:D1 = 2:1 line strengths. Key result: since
Δ_R ≪ Δ_FS the laser sits essentially on P₃/₂ alone, so the **vector shift is not
fine-structure-suppressed** (~0.5× scalar for circular light) — a **10% circular
contamination of R2 changes the differential shift ~47%**. For the OC run it gives
δ_AC = −44.5 kHz (vs the scalar +18.6 kHz) and Γ_sc ~2× the scalar estimate. Records
(provisional): `mg_fine_structure_splitting_3p_25mg` (derived, 2.746 THz) +
`raman_{b1,b3,r1,r2}_polarization_25mg` (Clos Tab. 3.2; geometry/B3 **flagged**). See
[ADR-0008](../docs/decisions/0008-polarization-power-resolved-raman-optics.md).

## Engine: `raman_dephasing` (relative-phase noise of the two beams)

[`engines/raman_dephasing.py`](engines/raman_dephasing.py) is the decoherence of a
two-photon flop from **relative optical phase/frequency noise** between the two Raman
beams — a channel *independent* of the carrier Debye-Waller motional dephasing, the
off-resonant scattering, and the AC-Stark shift. The spin tracks the beat-note phase
φ_B − φ_R; for beams from one laser the common-mode phase cancels (the two-photon
frequency is RF-set), leaving path-imbalance (delay τ = ΔL/c high-passes the laser
noise), AOM/DDS, fibre, and pointing noise. The contrast envelope has two limits:
Lorentzian/white → `exp(−π·Δν·t)` (exponential, rate π·Δν), or quasi-static Gaussian
→ `exp(−(t/T₂)²)`. It's **capability + diagnostic** (no measured Raman mutual
linewidth to anchor): its job is to convert an observed residual decay rate into the
mutual linewidth Δν / coherence time T_φ it implies (`mutual_linewidth_from_rate`),
as the **alternative to a hot motional state** — the two are degenerate in one flop.
See [ADR-0007](../docs/decisions/0007-raman-scatter-vs-dephasing-in-flop-twin.md).

## Discriminator: sideband thermometry (`twin_sideband`)

[`twin_sideband.py`](twin_sideband.py) (`python -m spike.twin_sideband`) settles the
question `twin_oc_flop` left degenerate — was the OC carrier-flop contrast loss
**motional** or **Raman-beam dephasing**? — with the RSB+BSB sideband scan
([`OC_Axial/1_1R_LF_MA`](../sources/data/OC_Axial/1_1R_LF_MA/)). That file drives a
blue-sideband pulse (counter 0) AND a red-sideband pulse (counter 1) per shot, so it
holds **two** flop blocks (`DatFile.counter_blocks`). The blue (adds a phonon) flops
fully while the red (subtracts) is **near-constant** — the ground state can't subtract
— and the RSB/BSB peak ratio is the direct thermometer n̄/(n̄+1) (`sideband.thermal_
sideband_flip`, η√n vs η√(n+1)). Result: **n̄ = 0.27 ± 0.13 (COLD)**, far below the
carrier's apparent n̄_eff = 1.06. Decomposing the carrier decay at the *measured* n̄:
**~64% Raman-beam dephasing** (Δν ≈ 21 kHz, T_φ ≈ 15 µs) + ~36% motional — the
apparent "hot ion" was Raman dephasing posing as motion. Figure
[`../docs/figures/twin_sideband_thermometry.png`](../docs/figures/twin_sideband_thermometry.png).

## Twin: stroboscopic OC carrier flop / phase-grating baseline (`twin_strobo`)

[`twin_strobo.py`](twin_strobo.py) (`python -m spike.twin_strobo`) is the displ=0
baseline of the **"active phase grating"** (`Strobo2.0/1_FlopN_3p3_2p2_PDQ_displ_strobo`).
The script cools to |0⟩, then drives **N = 50 stroboscopic OC carrier pulses, one per
`DELTA_t` = 0.769 µs = the lf motional period** (fr_oc_strobo at the qubit carrier −40
kHz), scanning the per-cycle pulse width `delta_t`. With `u_displ = 0` the motion stays
in |0⟩, so it is a clean stroboscopic carrier flop: P_up = ½(1−cos(2π·N·Ω_strobo·`delta_t`)),
N amplifying the per-cycle Rabi. The twin fits **Ω_strobo = 499 ± 3 kHz/cycle** (bare
Ω₀ = 538 kHz after the |0⟩ Debye-Waller e^(−η²/2)). Key result: the train lasts a
*fixed* N·DELTA_t = 38.5 µs yet keeps **~61% contrast** — vs the ~8% a continuous
carrier flop (T_φ ≈ 15 µs, `twin_sideband`) would leave — so the **stroboscopic
structure decouples the slow Raman-beam dephasing** (effective T_φ ≈ 77 µs). With a
displacement the flop rate would be position-modulated (the grating → collapse/revival):
the displ≠0 follow-up. Figure
[`../docs/figures/twin_strobo_flop.png`](../docs/figures/twin_strobo_flop.png).

It also **forward-simulates a detuning scan** of the pulse train (`delta_t` fixed at
0.02 µs) via [`engines/strobo_sim.py`](engines/strobo_sim.py) — an exact small spin⊗Fock
stroboscopic propagator (displacement matrix D(iη), per-cycle U_pulse × U_free(δ),
N=50 cycles). The result is the **stroboscopic comb**: full-contrast resonances at the
**carrier (δ=0) and the ±k teeth at ±k·f_lf**, each ~26 kHz wide (1/N·DELTA_t) with
finite-train sinc side lobes. These are **symmetric Floquet sidebands of the pulsed
drive — not motional red/blue sidebands**: at the exact strobe they are full-contrast
and *independent of η* (present even at η=0; the motion wraps to identity each cycle).
The motional coupling only shows up when the strobe is detuned off the period or the
motion is displaced. Figure [`../docs/figures/twin_strobo_detuning_scan.png`](../docs/figures/twin_strobo_detuning_scan.png).

The same propagator demonstrates the train as a **sampling mixer / heterodyne receiver**
(`heterodyne_beat`, `engines.strobo_sim.strobo_population_vs_cycles`): with `delta_t`
fixed, detuning the strobe by f_IF off a tooth makes the **cycle-domain** population
P_flip(N) the down-converted output — homodyne (the π flop) on the tooth, and a nutation
at f_IF off it (first turning point ~1/(2·f_IF·DELTA_t); 50/100/200 kHz → ~13/7/3
cycles). The detuning-scan comb is the corresponding superheterodyne image spectrum.
Figure [`../docs/figures/twin_strobo_heterodyne_beat.png`](../docs/figures/twin_strobo_heterodyne_beat.png).

It also estimates the **AC-Stark systematics vs N** (`ac_stark_vs_N`). In the pulse train
**B1 stays on continuously and only R2 is pulsed**; for a π pulse the R2 width is
`delta_t = 1/(2N·Ω_strobo)` (= 0.020 µs at N=50). The differential shift is **negative**
(|2,2⟩ is closer to P₃/₂ → the qubit resonance moves *down* — the sign of the −40 kHz
`fr_oc_strobo` offset and of the earlier `scatter` estimate). Its absolute scale is
anchored to the **observed** flop rate rather than the nominal powers/waists, which
over-predict the 2γ Rabi (~1.2 MHz vs observed ~0.54 MHz → κ≈0.46): Rabi-anchored
**B1 ≈ −33 kHz** (continuous → a constant detuning floor) and **R2 ≈ −100 kHz** (only
while pulsing), bracketed below by the scalar (ω_HF/Δ)·Ω₀ ≈ −48 kHz and above by the
nominal upper bound (−72 / −214 kHz). Because R2's on-time `1/(2Ω_strobo)` is fixed by
the π condition, its AC-Stark phase is **N-independent**, while **B1's grows ∝ N** (it
acts for the whole `N·DELTA_t` train — ~1.3 cycles of phase already at N=50; the *shape*
is robust to the absolute anchor). Figure
[`../docs/figures/twin_strobo_acstark_vs_N.png`](../docs/figures/twin_strobo_acstark_vs_N.png).

## Integrated twin: OC axial carrier flop (`twin_oc_flop`)

[`twin_oc_flop.py`](twin_oc_flop.py) (`python -m spike.twin_oc_flop`) is the twin of
the **OC orthogonal-carrier axial Raman flop**
([`sources/data/OC_Axial/0_Car_Flop`](../sources/data/OC_Axial/0_Car_Flop/)). It
fits the measured flop (`rabi.fit_rabi` → Ω/2π ≈ 208 kHz, t_π ≈ 2.4 µs) and composes
**four ledger-anchored channels**: coherent Rabi, the differential AC-Stark shift
(`scatter`), off-resonant scattering (`scatter`), and the **carrier Debye-Waller
motional dephasing** (`sideband`, η = 0.389 OC→lf × the RSB-cooled n̄ = 0.07). The
honest decomposition: at the cooled n̄ the ledger floor (scattering ~2 % + motional
~9 %) explains only ~**11 %** of the observed ~60 % contrast loss. The **effective-n̄
inversion** attributes the whole decay to motion through the same model and recovers
**n̄_eff = 1.06 ± 0.27** (300-replica shot bootstrap, ~3.7σ above the cooled
benchmark) — ~15× the cooled benchmark — i.e. the flop is consistent with a
near-unity-n̄ state (sideband-cooling underperformance / heating this run and/or
technical Raman intensity-phase / B-field dephasing), **not** the cooled 0.07. The
figure overlays the per-shot cloud, the bare ledger floor (n̄ = 0.07, which barely
decays), and the n̄_eff curve that tracks the data
([`../docs/figures/twin_oc_axial_carrier_flop.png`](../docs/figures/twin_oc_axial_carrier_flop.png)).
See [ADR-0007](../docs/decisions/0007-raman-scatter-vs-dephasing-in-flop-twin.md).

## Engines: `rabi` + `detection` (raw-data ingestion)

These consume **raw `.dat` measurement files** (kalis2017 DAQ format, see
[`../docs/DATA_FORMAT.md`](../docs/DATA_FORMAT.md)) via [`datfile.py`](datfile.py),
rather than ledger records — so they live behind [`analyze_data.py`](analyze_data.py)
(`python -m spike.analyze_data`), not the ledger-based composition root.

- [`engines/rabi.py`](engines/rabi.py) fits a damped Rabi flop
  y(t)=c+e^(−γt)(a·cos2πft+b·sin2πft) by a grid scan with an exact weighted
  linear least-squares solve for (c,a,b) at each (f,γ) — pure Python, via
  `linalg.solve`. On the kalis2017 |3,+3⟩↔|2,+2⟩ duration scan → Ω/2π = **53.3 kHz**
  (t_π = 9.4 µs), ~10% below `mw_rabi_3p3_2p2_doerr` (59.45 kHz) — the expected
  MW-power/day dependence of the apparatus-limited rate.
- [`engines/detection.py`](engines/detection.py) does Poissonian bright/dark
  discrimination: optimal count threshold + readout **fidelity** from the two
  means, the Mandel Q (super-Poissonian broadening), and the empirical fidelity
  from the per-shot histograms. kalis2017: F_Poisson = 0.992 vs F_empirical ≈ 0.97,
  the (one-sided) gap being the bright histogram's low-count tail — bright-state
  loss/depumping during detection (the reference channel matches Poisson). It also
  models that tail directly: **`transition_count_pmf`** gives a Poisson core plus the
  tail from a stochastic switch (bright→dark **depumping** → low counts, à la Thomm
  2021; dark→bright **leaking** → high counts), and **`ml_estimate_p_down`** is the
  OPTIONAL maximum-likelihood state readout (infer P↓ from a set of counts under the
  realistic mixture; default = pure Poisson). The raw counts stay primary.

The frequency scan's resonance (1775.49 MHz) also cross-checks the `levels` engine's
field-sensitive (3,+3)↔(2,+2) prediction (1775.60 MHz at the Weber field).

## Engine: `cooling` (Doppler scattering)

[`engines/cooling.py`](engines/cooling.py) is clean textbook two-level physics
(no apparatus model): from the natural linewidth Γ, detuning Δ and saturation
parameter s = I/Isat (all in the ledger from the laser table), the scattering
rate is R = (Γ/2) s/(1 + s + (2Δ/Γ)²) and the Doppler limit is
T_D = ℏΓ/(2k_B). It **reproduces the ~1 mK Doppler limit** (1.0030 mK vs the
theses' 1.0000 mK) and notes the optimal detuning −Γ/2 = −20.9 MHz *is* the BD
cooling setting. It also predicts the **Doppler-cooled mean phonon number**
n̄ = 1/(exp(2ω/Γ)−1) and reproduces Clos's **measured** n̄ = 10(1) at ω₁/2π =
1.915 MHz (predicted 10.42, 0.42σ) — the independent, *measured* motional
benchmark (vs the theory-consistency Doppler limit). The runner adds a cooling
diagnostic (scatter rate per beam).

## Engines: `readout` + `sideband_cooling` (Thomm 2021 diagnostics)

Two engines connect the Thomm-2021 readout/motional benchmarks to physics. They are
**diagnostics, not σ-validations** — the measured values are preparation/protocol-
limited, not zero-parameter predictions, and forcing a σ-test would be circular or
show a false tension.

- [`engines/readout.py`](engines/readout.py) turns the detection count model into
  readout figures of merit: the **single-shot discrimination fidelity** (94.8 % at
  Thomm's λ↓=2.682/λ↑=0.036) and the **Fisher-information / Cramér-Rao precision** of a
  maximum-likelihood P↓ estimate over N shots (overhead ×1.1 at p=0.5, ×3 near the
  bright end vs ideal QPN). The runner's `readout` diagnostic decomposes the measured
  **99.4 %/97.4 %** ensemble fidelities *by channel*: the bright 0.6 % deficit is
  **preparation-dominated** (ML of a perfect bright state averages above the 95 %
  single-shot cap), while the dark 2.6 % deficit is **detection** (dark-state
  off-resonant scatter / depumping during t_det — a systematic ML averaging does *not*
  remove), so the readout limits the **dark** channel.
- [`engines/sideband_cooling.py`](engines/sideband_cooling.py) estimates the
  resolved-sideband cooling floor n̄_min ∼ α(κ/2ω)² (off-resonant carrier limit,
  LBMW03; α an O(1) factor) and **inverts** each achieved n̄ for the implied effective
  cooling rate κ. The `sideband_cooling` diagnostic confirms all κ/ω < 1 (resolved
  regime) but κ varies across modes → the achieved n̄ (0.07/0.11/0.07) are
  protocol/per-mode-limited, not a single common floor; it also cross-checks Thomm's
  P(n=0)↔n̄ (Bose-Einstein).

## Engine: `tickle` (Kalis 2016 secular-frequency spectroscopy)

[`engines/tickle.py`](engines/tickle.py) measures motional-mode frequencies from a
**tickle** (secular-excitation) scan, following the group's own method (Kalis et al.,
PRA **94**, 023401 (2016); thesis kalis2017). A finite resonant excitation pulse drives
a mode as a classical oscillator → amplitude A ∝ sin([ω_exc−ω_i]t_exc/2)/(ω_exc²−ω_i²)
— a **sinc** of FWHM ≈ 1/t_exc, *not* a Lorentzian — and the Doppler modulation
(index β=⟨u,k_BD⟩A) moves carrier population to motional sidebands with Bessel weights
J_v(β)² (pure-Python `besselj`, |v|≤15), so the fluorescence F = Σ_v J_v(β)²·g(Δ_BD+vω)
(Eq. 2) dips at resonance. `fit_tickle` fits the robust leading-order sinc² dip for the
secular frequency f₀ (the depth is left free — its absolute value depends on the
detection sensitivity Δ_BD/Γ_w/beam-waist, not pinned here), with an **F-test + edge
guard** that rejects narrow calibration scans that miss the mode. On the PAULA
`Tickle/PDQ_*_FScan` data (`python -m spike.plot_tickle`,
[figure](../docs/figures/tickle_modes.png), all files): axial **1.299 MHz** (ledger
1.30, 10/10 files), radial **3.224/4.712 MHz** (ledger 3.0/4.5 → **+7.5%/+4.7%**, 3/4
and 4/4 — the radial nominals need refining).

## Integrated twin: `spin` + `twin` (prepare → drive → detect)

The engines above are isolated calculators. The **integrated twin** weaves them
into one experiment cycle that transforms a single quantum **state** through
state-preparation → manipulation → detection, simulated **per shot** so the count
cloud emerges, *and* as the ensemble average.

- [`engines/spin.py`](engines/spin.py) is the qubit **state** — a Bloch vector
  (z=+1 → |↓⟩, P_up=(1−z)/2) — with coherent operations: `prepare(eps)` (residual
  |↑⟩ population), `pulse(Ω, δ, t, φ)` (Rodrigues rotation by Ω_eff·t about the
  generalised-Rabi axis), `free(δ, t)` (z-precession). A resonant π-pulse takes
  |↓⟩→|↑⟩, π/2 → equator; a detuned pulse caps the flip at Ω²/Ω_eff².
- [`twin.py`](twin.py) composes them. The **microwave drive carries an AC-Zeeman
  shift and (quasi-static) dephasing**. Key physics: the AC-Zeeman shift exists
  only while the MW is *on* (the pulses) — during a Ramsey free gap (MW off) the
  spin precesses at the bare detuning. So the Rabi (MW-on) resonance *includes*
  the shift while the Ramsey fringe (free precession) *reveals* it, and the
  fringe's Gaussian contrast decay e^(−(τ/T₂\*)²) gives the dephasing
  (σ_δ = √2/2πT₂\*). Per shot a quasi-static δ_noise ∼ 𝒩(0, σ_δ) is drawn, the
  state is evolved, **projectively measured** (QPN), and a **detection count**
  drawn from Poisson(μ_bright) for |↓⟩ / Poisson(μ_dark) for |↑⟩.
  `detection_levels` sets μ_bright = R_scatter·η·t_det from the **Friedenauer**
  collection efficiency η = 5.6×10⁻³ (`mg_detection_efficiency_25mg`).
- [`twin_demo.py`](twin_demo.py) (`python -m spike.twin_demo`) is the **test case**:
  inject (Ω=50 kHz, AC-Zeeman=3 kHz, T₂\*=800 µs), simulate a Rabi flop and a
  Ramsey fringe per shot, then **infer them back** over *N*=16 Monte-Carlo replicas
  (so the inference carries an uncertainty, not a single noisy draw) — the Rabi flop
  pins Ω = 50.15 ± 0.74 kHz (and is blind to the small shift / slow dephasing), the
  Ramsey fit (`fit_ramsey`, cosine × Gaussian) recovers AC-Zeeman = 3.00 ± 0.04 kHz
  and T₂\* = 791 ± 58 µs (all consistent with injected within the scatter). The
  detection levels (μ_b≈6.3) are themselves derived from the Friedenauer efficiency,
  not hand-set. Writes
  [`../docs/figures/twin_rabi_ramsey_inference.png`](../docs/figures/twin_rabi_ramsey_inference.png).
- [`twin_freqscan.py`](twin_freqscan.py) (`python -m spike.twin_freqscan`) is the
  **frequency-domain** version (`make_seq_ramsey_freq`): a Rabi spectroscopy scan
  (broad line, ±25 kHz) whose dip is **pulled to** f₀+δ_ACZ, vs Ramsey scans at
  τ=100/300/600 µs whose **fringe combs sit on the bare f₀** (free precession is
  MW-off, light-shift-free) with spacing 1/τ_eff (τ_eff = τ + 1/πΩ, the finite-pulse
  correction). The detuning window is **resolution-matched** (±2.5 fringe spacings),
  so longer τ zooms in — an RPE-style ladder. The Rabi dip vs Ramsey comb offset
  gives AC-Zeeman = 3.02 ± 0.41 kHz (injected 3.00; envelope-pull systematic removed,
  robust stats over N=100), the Ramsey comb centre being ~13× sharper than the broad
  Rabi dip. Writes
  [`../docs/figures/twin_freqscan_rabi_ramsey.png`](../docs/figures/twin_freqscan_rabi_ramsey.png).
- [`twin_detection.py`](twin_detection.py) (`python -m spike.twin_detection`) — the
  **realistic readout**: the detection step in `simulate_counts` now lets a bright ion
  **depump** (`MWModel.depump_bright`) and a dark ion **leak** (`leak_dark`), so the
  bright count histogram grows the Thomm low-count tail. The demo (Thomm levels
  λ↓=2.682/λ↑=0.036, depump Γt=0.3) shows the histogram tail and a Rabi flop read out
  by the optional **ML readout** (RMSE 0.028 vs truth) vs a fixed **threshold** (RMSE
  0.090, biased low by the tail). Writes
  [`../docs/figures/twin_detection_depumping.png`](../docs/figures/twin_detection_depumping.png).

## Layout

```
spike/
  constants.py      CODATA constants + 25Mg atomic g-factors (sourced)
  ledger.py         loads records/*.yaml into a queryable Ledger (pyyaml only)
  linalg.py         tiny pure-Python solve() + eigvalsh() (Jacobi)
  bootstrap.py      fit uncertainties: robust summary + gaussian (perturb-by-sigma) + shot (resample-counts) resamplers
  engines/
    levels.py       2S_1/2 hyperfine+Zeeman engine (closed-form Breit-Rabi)
    modes.py        axial + radial normal modes (equilibrium + Hessian)
    drive.py        relative microwave Rabi couplings (Clebsch-Gordan)
    cooling.py      two-level scattering rate + Doppler limit + occupation
    projection.py   Raman combination -> motional mode (Delta_k . mode axis)
    sideband.py     absolute Lamb-Dicke + sideband Rabi + Raman differential AC-Stark + carrier Debye-Waller + thermal RSB/BSB flops (nbar thermometry)
    acstark.py      far-detuned single-beam light shift (BDD vs Hasse)
    scatter.py      Raman off-resonant scattering (Gamma_sc, P_SE/pi) + differential AC-Stark + flip_probability
    raman_optical.py polarization+power-resolved light shifts + scattering (|J',mJ'> basis; scalar/vector; 6j cross-check; absolute differential_stark_hz)
    raman_dephasing.py relative-phase noise of the two beams: contrast envelopes + mutual-linewidth/T_phi readout of the residual
    strobo_sim.py    stroboscopic spin-motion propagator (displacement matrix + U_pulse/U_free) -> detuning-scan comb (carrier + sidebands)
    grating_tomography.py  char.-function/Wigner transfer-function kernels (chi, double/single-sum kernels, exact eta=0, chi->W reconstruction) + Ramsey 2-pulse chi-interferometer (population->chi(Delta beta), |Delta beta|<=2eta disk) + self-checks
    rabi.py         damped Rabi-flop fit -> Omega (raw .dat duration scans)
    detection.py    bright/dark discrimination: threshold + fidelity + depumping/leak PMF + ML readout
    readout.py      single-shot fidelity + Fisher/Cramer-Rao P_down precision (diagnostic)
    tickle.py       Kalis-2016 secular-freq spectroscopy: Bessel sinc lineshape + fit_tickle
    sideband_cooling.py  RSB cooling floor n_min ~ (kappa/2omega)^2 + per-mode kappa inversion (diagnostic)
    spin.py         qubit STATE: Bloch vector + prepare/pulse/free operations
  datfile.py        reader for the PAULA DAQ .dat files (kalis2017 format)
  runner.py         composition root: registry, ValidationResult, table, diagnostics
  twin.py           integrated cycle: MWModel, Rabi/Ramsey time+freq sequences, ensemble + per-shot MC, fit_ramsey
  twin_demo.py      test case: infer AC-Zeeman + dephasing from Rabi vs Ramsey TIME scans (-> figure)
  twin_freqscan.py  Rabi vs Ramsey FREQUENCY scans (tau=100/300/600us), resolution-matched windows (-> figure)
  twin_detection.py realistic detection: depumping count tail + optional ML readout (-> figure)
  twin_oc_flop.py   OC axial Raman carrier-flop twin: coherent x AC-Stark x scatter x motional (n_bar) vs .dat + n_bar_eff inversion (-> figure)
  twin_sideband.py  sideband thermometry discriminator: RSB/BSB ratio -> n_bar (cold) -> carrier loss is mostly Raman dephasing (-> figure)
  twin_strobo.py    stroboscopic OC carrier flop (phase-grating n=0 baseline): per-cycle Rabi + stroboscopic dephasing-decoupling (-> figure)
  analyze_data.py   raw-data analysis (rabi + detection on the .dat examples)
  validate_twin.py  CLI shim -> runner.main
  test_levels.py    levels physics + Weber/Doerr benchmarks + hyperfine spectrum
  test_modes.py     axial+radial eigenvalues + benchmarks + linalg robustness
  test_drive.py     Clebsch-Gordan vs sympy + symmetry + polarization
  test_cooling.py   scatter rate + Doppler limit + occupation benchmark
  test_projection.py  mode-axis geometry + addressed-modes vs Doerr + from_ledger
  test_sideband.py  absolute eta scaling + anchor + sideband Rabi + Raman stark
  test_acstark.py   light-shift formula + BDD vs Hasse + far-detuned predicate
  test_datfile.py   .dat parsing: scan def, signal vs reference, histograms
  test_rabi.py      damped-cosine fit recovery + kalis duration scan
  test_detection.py Poisson stats + threshold/fidelity + depump/leak/ML + kalis histograms
  test_readout.py   single-shot fidelity + Fisher info limits + CRB>=QPN + from_ledger
  test_sideband_cooling.py  RSB floor scaling/inversion + offres carrier + resolved regime
  test_tickle.py    Bessel values + sinc shape + fluorescence dip + f0 recovery + real LF data
  plot_tickle.py    tickle spectroscopy of the 3 modes (Kalis fit -> figure)
  test_runner.py    runner: result math, tension detection, wall refusal, coverage
  test_spin.py      Bloch operations: pi/pi2/2pi pulses, detuned cap, free precession
  test_twin.py      sequences, dephasing decay, fit_ramsey recovery, count cloud, detection levels
```

## Run

```bash
pytest spike/                  # all engine + runner tests
python -m spike.validate_twin  # the twin validation table
python -m spike.twin_demo      # integrated twin: Rabi-vs-Ramsey TIME-scan inference -> figure
python -m spike.twin_freqscan  # Rabi vs Ramsey FREQUENCY scans (tau=100/300/600us) -> figure
python -m spike.twin_detection # realistic detection: depumping tail + optional ML readout -> figure
python -m spike.twin_oc_flop   # OC axial Raman carrier-flop twin: AC-Stark + off-resonant scattering vs data -> figure
python -m spike.plot_tickle    # tickle secular-frequency spectroscopy of the 3 modes -> figure
```

## Wall discipline

The wall is enforced at the ledger boundary: `Ledger.input_quantity()` refuses a
`benchmark`, `benchmark_quantity()` refuses an `input`. Each engine's
`from_ledger(...)` consumes only inputs (`levels` reads A/I/g_J/g_I; `modes`
reads the COM), so it cannot "predict" a number it was handed. The runner reads
each benchmark once, at comparison time, via `benchmark_quantity()`. If a
validation violates this, the runner surfaces it as an `ERROR` row and exits
nonzero.

## Not yet (follow-ups)

- Anchor the Raman **differential** AC-Stark shift to a measured value (none in the
  theses yet) and refine it beyond the leading-order ω_HF/Δ_R (Clebsch-Gordan
  differences, P₁/₂ vs P₃/₂, beam imbalance).
- Extend `modes` to the full N-ion radial spectrum and connect it to `sideband`.
- Graduate the spike to its own repo(s) once the schema/engines stabilise.
