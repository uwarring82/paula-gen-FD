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
benchmark beyond the η anchor; the differential shift has no measured value).

## Engine: `acstark` (far-detuned light shift)

[`engines/acstark.py`](engines/acstark.py) is the far-detuned light shift
δ_AC = sΓ²/(8δ) = Ω²/(4δ) (Clos Eq. 2.2.24). It reproduces Hasse's **measured**
BDD shift (~2π×10 MHz → predicted 10.45 MHz, 0.15σ) — a real σ-validation. The
shift is only meaningful **far** from resonance: BDD (−10Γ) shifts coherently,
while the near-resonant BD/BDX/RD/RP *scatter* (cooling engine), which Hasse
confirms ("RD/RP induce no significant ac Stark shift"). The same validated
formula feeds the Raman differential shift in `sideband`.

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

## Layout

```
spike/
  constants.py      CODATA constants + 25Mg atomic g-factors (sourced)
  ledger.py         loads records/*.yaml into a queryable Ledger (pyyaml only)
  linalg.py         tiny pure-Python solve() + eigvalsh() (Jacobi)
  engines/
    levels.py       2S_1/2 hyperfine+Zeeman engine (closed-form Breit-Rabi)
    modes.py        axial + radial normal modes (equilibrium + Hessian)
    drive.py        relative microwave Rabi couplings (Clebsch-Gordan)
    cooling.py      two-level scattering rate + Doppler limit + occupation
    projection.py   Raman combination -> motional mode (Delta_k . mode axis)
    sideband.py     absolute Lamb-Dicke + sideband Rabi + Raman differential AC-Stark
    acstark.py      far-detuned single-beam light shift (BDD vs Hasse)
  runner.py         composition root: registry, ValidationResult, table, diagnostics
  validate_twin.py  CLI shim -> runner.main
  test_levels.py    levels physics + Weber/Doerr benchmarks + hyperfine spectrum
  test_modes.py     axial+radial eigenvalues + benchmarks + linalg robustness
  test_drive.py     Clebsch-Gordan vs sympy + symmetry + polarization
  test_cooling.py   scatter rate + Doppler limit + occupation benchmark
  test_projection.py  mode-axis geometry + addressed-modes vs Doerr + from_ledger
  test_sideband.py  absolute eta scaling + anchor + sideband Rabi + Raman stark
  test_acstark.py   light-shift formula + BDD vs Hasse + far-detuned predicate
  test_runner.py    runner: result math, tension detection, wall refusal, coverage
```

## Run

```bash
pytest spike/                  # all engine + runner tests
python -m spike.validate_twin  # the twin validation table
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
