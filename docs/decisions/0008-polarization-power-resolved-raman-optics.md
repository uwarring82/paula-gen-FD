# 0008 — Polarization- and power-resolved Raman optics in the fine-structure basis

- **Status:** Accepted
- **Date:** 2026-06-22
- **Deciders:** UW (decision + physics steer), Claude (implementation)

## Context

The Raman light shifts (ADR-0007's AC-Stark channel) and the off-resonant scattering
were treated **polarization-blind and power-blind**: the differential AC-Stark shift
was a leading-order `omega_HF/Delta_R` scalar (`sideband.raman_differential_stark_
factor`) and the scattering carried a single scalar `balance` handle (`scatter.py`).
UW asked to **carefully distinguish the polarization and power of the individual
Raman beams**, because both significantly affect the shifts and the spontaneous
emission.

A first cut summed over the excited **|F', mF'>** hyperfine sublevels. UW flagged the
real problem: **F', mF' are not good quantum numbers for the 3P levels** — at
Delta_R = 20 GHz the excited hyperfine + Zeeman structure is completely unresolved,
so an F'-resolved sum is, at best, a basis-complete bookkeeping device. UW asked to
**set it up more versatile**.

## Decision

Build `spike/engines/raman_optical.py` in the **fine-structure |P_{J'} mJ'> basis**
for the excited manifold (no excited hyperfine), with the **ground** state kept in
resolved |S_{1/2} F mF> (there F *is* good: hyperfine 1.79 GHz >> Zeeman at 5.5 G)
and decomposed into electronic-spin components |mJ, mI> with the **nuclear spin mI a
spectator** of the optical dipole:

    |F mF> = sum_mJ <S_{1/2} mJ, I (mF-mJ) | F mF> |mJ, mI=mF-mJ>,
    <P_{J'} mJ'| d_q | S_{1/2} mJ> = CG(1/2 mJ; 1 q | J' mJ') * d0   (d0 common; the
                                     2:1 line strength is the (2J'+1) multiplicity).

Single-beam light shift / scattering of |F mF> are the mJ-weighted sums over q, J';
the two-photon Rabi conserves mI. Detunings are signed: the laser is ~20 GHz RED of
3P_3/2 and ~2.7 THz BLUE of 3P_1/2.

- **Versatility / honesty.** This is basis-independent in the excited manifold (the
  honest far-detuned limit). The old **F'-basis sum is kept only as a degenerate-
  limit CROSS-CHECK**: `coupling`/`line_coupling_sq` (via a Wigner 6j) equal the |mJ>
  sum to **machine precision** (basis independence, asserted in the tests). The
  transparent summary is the **scalar + vector polarizability** (`scalar_vector_
  shift`; tensor = 0 for 2S_1/2). An optional `p_hyperfine` hook is retained but
  defaults OFF (excited hyperfine degenerate).

- **Anchoring.** The absolute field scale is unknown but CANCELS in the dimensionless
  ratios `differential_stark_per_rabi` and `scatter_per_rabi` (delta_AC and Gamma_sc
  in units of the two-photon Rabi). The twin multiplies by the **measured** flop
  Rabi: `delta_AC[Hz] = ratio * rabi_hz`, `Gamma_sc[1/s] = ratio * 2pi * rabi_hz`
  (matching the scalar engine's conventions). So the engine refines the polarization/
  power/angular factors while staying pinned to the data.

- **Records (provisional, seeded + flagged; ADR-0001 pattern).**
  `mg_fine_structure_splitting_3p_25mg` = 2.7457 THz is **derived** from the two
  recorded D-line wavelengths (solid). The per-beam polarizations
  `raman_{b1,b3,r1,r2}_polarization_25mg` are seeded from Clos Tab. 3.2 (B1 pi; R1/R2
  balanced linear, C=0) as (sigma+, pi, sigma-) intensity fractions along B, **flagged
  PROVISIONAL**: the lab-frame -> B projection geometry is not fully tabulated and
  **B3's polarization is not stated anywhere** (placeholder, large uncertainty).
  POWER is a per-RUN setting (the .dat `pwr_b1/pwr_r2`), read by the twin.

## Consequences

- **+** Polarization and power now drive the shifts and scattering. Because
  Delta_R << Delta_FS the laser sits essentially on P_3/2 alone, so the **vector
  light shift is NOT fine-structure-suppressed** (~0.5x the scalar for circular
  light): a **10% circular contamination of R2 changes the differential shift by
  ~47%** — the large polarization sensitivity UW pointed to.
- **+** For this OC run (B1 pi + R2 balanced-linear, equal powers) the resolved
  differential is `delta_AC = -44.5 kHz` (vs the scalar +18.6 kHz — now with a sign
  and a ~4.4% amplitude cap), and `Gamma_sc` ~2x the scalar estimate. The crude
  factors were explicitly order-of-magnitude; this is the documented refinement.
- **+** The core is rigorously validated: Wigner-6j canonical values, the cycling
  transition, the sublevel/F-independent sum rule, exact D2:D1 = 2:1, and **basis
  independence to 1e-12**.
- **-** The dominant flop decay is still motional/technical (ADR-0007); a larger
  AC-Stark cap (4.4%) and 2x scattering raise the ledger floor only from ~11% to
  ~13% of the observed decay. The headline (n_bar_eff ~ 1) is unchanged.
- **- (needs UW)** The polarization records are provisional: the beam-frame ->
  quantization-axis geometry and B3's polarization need apparatus confirmation; any
  real circular contamination would change `delta_AC` substantially. The excited
  hyperfine is treated as degenerate (~1-2% detuning spread at 20 GHz).
