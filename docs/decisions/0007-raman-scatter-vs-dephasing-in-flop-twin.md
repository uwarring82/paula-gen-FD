# 0007 — Raman flop twin separates the ledger-anchored loss floor from the empirical residual

- **Status:** Accepted
- **Date:** 2026-06-22
- **Deciders:** UW (decision), Claude (implementation)

## Context

We wanted a digital twin of the OC (Orthogonal-Carrier) axial Raman **carrier
flop** (`sources/data/OC_Axial/0_Car_Flop/09_46_23_12_06_2026.dat`) that includes
the **AC-Stark shift** and the **spontaneous-emission decoherence from off-resonant
scattering** — the two effects a two-photon stimulated-Raman (TPSR) drive carries
because both beams sit ~20 GHz red of `3P_3/2` (`raman_detuning_from_p32`).

The physics is unambiguous and ledger-derivable from a single detuning Δ_R + the
linewidth Γ (`mg_p32_natural_linewidth`):

- coherent carrier Rabi `Ω = Ω_B Ω_R / (2 Δ_R)`;
- differential AC-Stark shift `δ_AC ≈ (ω_HF/Δ_R)·Ω` (the |↓⟩,|↑⟩ light shifts differ
  by the qubit splitting `ω_HF`, `hyperfine_splitting_25mg_f2_f3`) — speeds Ω_eff and
  caps the flop amplitude at `Ω²/Ω_eff²`;
- off-resonant scattering rate `Γ_sc = Γ(Ω_B²+Ω_R²)/(4Δ_R²)`, which for balanced
  beams collapses to `Γ_sc = (Γ/Δ_R)·Ω`, so the **scattering-per-π floor is detuning
  only**: `P_SE(π) = π·Γ/Δ_R ≈ 0.66 %` at 20 GHz.

When fitted to the data, the observed flop **contrast** decays at ≈1.0×10⁵ /s
(τ ≈ 10 µs, ~60 % loss over the scan), while the scattering floor predicts only
≈2×10³ /s (τ ≈ 0.5 ms, ~2 % loss). **Spontaneous scattering explains only ~2 % of
the observed decay.** That is not a defect — it is the designed consequence of a
large Raman detuning (`Γ_sc/Ω ∝ Γ/Δ_R`; Ozeri 2007, Wineland): the detuning was
deliberately set so scattering is a sub-percent-per-π error. The dominant decay is
**motional / thermal dephasing** (the carrier Debye-Waller spread over the thermal
phonon distribution), for which we have no zero-parameter ledger prediction here.

The decision: how should the twin present a decay it can only *partly* derive from
first principles?

## Decision

The twin (`spike/twin_oc_flop.py`, on `spike/engines/scatter.py`) **separates the
ledger-anchored physical floor from an explicitly-empirical residual**, and never
folds one into the other:

- **AC-Stark** and **scattering** are computed only from ledger `input` records
  (Δ_R, Γ, ω_HF) via `RamanScatter.from_ledger` — wall-safe, zero free parameters.
- The **residual** dephasing is reported as `γ_resid = γ_observed − (3/4)Γ_sc`,
  **labelled as motional/thermal dephasing, not predicted**. It is the only fitted
  decay term and it is named as such in the report and the figure.
- The figure plots **two** curves: the full twin (with the residual) over the data,
  *and* the "ledger floor" (coherent × AC-Stark × scattering only) — which visibly
  barely decays — so the reader sees how little scattering removes.
- The decoherence model is documented as leading-order: each scattering event fully
  depolarises the qubit (T₁=T₂=1/Γ_sc → envelope `e^{-(3/4)Γ_sc t}`,
  `CONTRAST_DECAY_FACTOR`). Rayleigh/Raman branching and qubit-manifold leakage are
  named as refinements, not silently assumed away.

There is **no measured Raman-scattering decoherence rate** in the theses to anchor
against, so `scatter.py` ships **capability + diagnostic** tests (formula identities,
the balanced collapse `Γ_sc=(Γ/Δ_R)Ω`, the `1/Δ_R` scaling, the 0.66 %/π floor), not
a σ-validation — consistent with how `sideband.py` treats the (also unanchored)
differential AC-Stark ratio.

## Consequences

- **+** The twin reproduces the flop *and* quantifies each loss channel honestly:
  the user asked for AC-Stark + scattering, gets both, plus the result that they are
  small here and the decay is dephasing-dominated.
- **+** The physical channels stay zero-parameter and wall-enforced; only the
  residual is fitted, and it is visibly flagged (no fabricated "prediction").
- **+** `engines/scatter.py` is reusable for any TPSR flop (CC/OC/AC/ROC) and for a
  future Raman-gate error budget (`P_SE/π`).
- **−→+** The dominant decay (motional dephasing) was originally a phenomenological
  residual; the follow-up below makes it a third ledger-anchored channel.
- **−** Beam imbalance (the data's IR/green powers differ) only enters through the
  `balance` handle as a scalar; a full treatment needs the two single-beam Rabi
  rates, which are not recorded.

## Update — 2026-06-22: the motional channel is now predicted (not fitted)

The follow-up is implemented. `spike/engines/sideband.py` gained the **carrier
Debye-Waller thermal dephasing**: Ω_{n,n} = Ω₀·e^(−η²/2)·L_n(η²) (`carrier_rabi_
factor`), the thermal flip-probability sum (`thermal_carrier_flip`), its envelope
|Σ P_n e^{iΩ_n t}| (`thermal_coherence`), and the leading-order spread
σ_Ω = Ω₀ η² √(n̄(n̄+1)) (`thermal_dephasing_rate`). The twin now composes **four
ledger-anchored channels** (coherent, AC-Stark, scattering, motional), the motional
one from η (sideband, OC→lf at `omega_z_axial_com_25mg`) × the RSB-cooled n̄
(`mg_rsb_cooled_nbar_axial_lf_25mg` = 0.07). The empirical exponential residual is
**gone**.

The honest result stands and is now quantitative: at the cooled n̄ = 0.07 the ledger
floor (scattering ~2 % + motional ~9 %) explains only ~**11 %** of the observed ~60 %
contrast loss. The twin adds an **effective-n̄ inversion** — solve the same thermal
model for the n̄ that reproduces the observed decay — giving **n̄_eff ≈ 1.1**, ~16× the
cooled benchmark. So this OC flop is consistent with a **near-unity-n̄ motional state**
(sideband-cooling underperformance / heating in this run and/or technical Raman
intensity-phase / B-field dephasing), not the cooled 0.07. The twin separates and
quantifies; the residual is flagged for the experimentalist, not fabricated away.

Note n̄ = 0.07 is consumed via `benchmark_quantity` (it is a measured benchmark, from
Thomm 2021's 200 GHz epoch): the motional prediction is "anchored to a measured n̄",
honestly labelled, and the twin is a standalone driver — not a wall-enforced
σ-validation in the runner. The inversion uses the analytic `thermal_coherence`
envelope (not a per-n̄ curve fit), keeping `build()` ~instant.

Remaining: predict n̄ itself from first principles (Doppler/RSB cooling engines) to
turn the inversion into a closed cross-check, and confirm whether n̄_eff ≈ 1 is a
cooling/heating issue or technical dephasing by a dedicated measurement (a coherence
scan, or a flop right after Doppler vs RSB cooling).
