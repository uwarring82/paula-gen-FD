# Apparatus parameters

**Auto-generated from `records/*.yaml` — do not edit by hand.** Regenerate with
`python tools/gen_apparatus.py`. The records are the source of truth (each carries
its full provenance: source, uncertainty, conditions, caveats); this is a flat
index of every apparatus-specific value, grouped by subsystem.

83 records (79 provisional, 4 confirmed). `kind`: **input** = consumed by the engines; **benchmark** = a measured/inferred value the engines are checked against. See [STATE_OF_THE_TWIN.md](STATE_OF_THE_TWIN.md) and [SOURCES.md](SOURCES.md).

## Internal state (²⁵Mg⁺ hyperfine / Zeeman)

| name | kind | value | units | gen | status | source |
|------|------|-------|-------|-----|--------|--------|
| `clock_transition_25mg` | benchmark | 1.789e+09 | Hz | freddy | provisional | doerr2024 |
| `clock_transition_predicted_25mg` | input | 1.789e+09 | Hz | freddy | provisional | itano_wineland_1981 |
| `clock_transition_residual_25mg` | benchmark | 2660 | Hz | freddy | provisional | doerr2024 |
| `clock_transition_weber_25mg` | benchmark | 1.789e+09 | Hz | freddy | provisional | weber2025 |
| `g_factor_electron_2s12` | input | 2.00232 | dimensionless | legacy | provisional | codata2018 |
| `g_factor_nuclear_25mg` | input | -0.34218 | dimensionless | legacy | provisional | stone2005 |
| `hyperfine_a_constant_25mg` | input | -5.963e+08 | Hz | legacy | confirmed | itano_wineland_1981 |
| `hyperfine_level_f2_25mg` | input | 2 | dimensionless | legacy | confirmed | clos2017 |
| `hyperfine_level_f3_25mg` | input | 3 | dimensionless | legacy | confirmed | clos2017 |
| `hyperfine_splitting_25mg_f2_f3` | input | 1.79e+09 | Hz | legacy | provisional | clos2017 |
| `hyperfine_splitting_calc_25mg` | input | 1.789e+09 | Hz | legacy | provisional | itano_wineland_1981 |
| `nuclear_spin_24mg` | input | 0 | dimensionless | legacy | provisional | clos2017 |
| `nuclear_spin_25mg` | input | 2.5 | dimensionless | legacy | confirmed | itano_wineland_1981 |
| `nuclear_spin_26mg` | input | 0 | dimensionless | legacy | provisional | clos2017 |
| `qubit_quadratic_zeeman_coeff_25mg` | input | 2.195e+11 | Hz/T^2 | legacy | provisional | itano_wineland_1981 |

## Optics (laser + Raman beams)

| name | kind | value | units | gen | status | source |
|------|------|-------|-------|-----|--------|--------|
| `bd_cooling_detuning` | input | -2.09e+07 | Hz | legacy | provisional | clos2017 |
| `bd_cooling_saturation` | input | 0.5 | dimensionless | legacy | provisional | clos2017 |
| `bd_detection_waist` | input | 5e-05 | m | legacy | provisional | clos2017 |
| `bd_laser_wavelength` | input | 2.796e-07 | m | legacy | provisional | clos2017 |
| `bdd_ac_stark_shift_25mg` | benchmark | 1e+07 | Hz | hasse | provisional | hasse2025 |
| `bdd_far_cooling_detuning` | input | -4.18e+08 | Hz | legacy | provisional | clos2017 |
| `bdd_far_cooling_saturation` | input | 20 | dimensionless | legacy | provisional | clos2017 |
| `mg_fine_structure_splitting_3p_25mg` | input | 2.746e+12 | Hz | legacy | provisional | clos2017 |
| `mg_p12_natural_linewidth` | input | 4.13e+07 | Hz | legacy | provisional | clos2017 |
| `mg_p32_natural_linewidth` | input | 4.18e+07 | Hz | legacy | provisional | clos2017 |
| `mg_saturation_intensity` | input | 2550 | W/m^2 | legacy | provisional | friedenauer2010 |
| `raman_ac_combination_25mg` | input | [-1.41421, 0, -1.41421] | k-vector Delta_k/k (trap x,y,z) | freddy | provisional | doerr2024 |
| `raman_axial_lamb_dicke_25mg` | input | 0.32 | dimensionless | legacy | provisional | wittemer2019 |
| `raman_b1_polarization_25mg` | input | [0, 1, 0] | polarization intensity fractions (sigma+, pi, sigma-) along B | freddy | provisional | clos2017 |
| `raman_b3_polarization_25mg` | input | [0, 1, 0] | polarization intensity fractions (sigma+, pi, sigma-) along B | freddy | provisional | doerr2024 |
| `raman_beam_path_jitter_25mg` | benchmark | 0.0059 | m/s | freddy | provisional | paula_oc_axial_2026 |
| `raman_cc_combination_25mg` | input | [0, 0, 0] | k-vector Delta_k/k (trap x,y,z) | freddy | provisional | doerr2024 |
| `raman_detuning_from_p32` | input | 2e+10 | Hz | freddy | provisional | doerr2024 |
| `raman_mutual_linewidth_25mg` | benchmark | 21000 | Hz | freddy | provisional | paula_oc_axial_2026 |
| `raman_oc_combination_25mg` | input | [0, 0, 1.41421] | k-vector Delta_k/k (trap x,y,z) | freddy | provisional | doerr2024 |
| `raman_r1_polarization_25mg` | input | [0.5, 0, 0.5] | polarization intensity fractions (sigma+, pi, sigma-) along B | freddy | provisional | clos2017 |
| `raman_r2_polarization_25mg` | input | [0.5, 0, 0.5] | polarization intensity fractions (sigma+, pi, sigma-) along B | freddy | provisional | clos2017 |
| `raman_roc_combination_25mg` | input | [1.41421, 0, 0] | k-vector Delta_k/k (trap x,y,z) | freddy | provisional | doerr2024 |
| `rd_laser_wavelength` | input | 2.804e-07 | m | legacy | provisional | clos2017 |
| `rd_repump_detuning` | input | -2.065e+07 | Hz | legacy | provisional | clos2017 |
| `rd_repump_saturation` | input | 0.5 | dimensionless | legacy | provisional | clos2017 |
| `rp_repump_detuning` | input | -2.065e+07 | Hz | legacy | provisional | clos2017 |
| `rp_repump_saturation` | input | 1 | dimensionless | legacy | provisional | clos2017 |

## Detection / readout

| name | kind | value | units | gen | status | source |
|------|------|-------|-------|-----|--------|--------|
| `bdx_detection_detuning` | input | -4.18e+06 | Hz | legacy | provisional | clos2017 |
| `bdx_detection_saturation` | input | 0.2 | dimensionless | legacy | provisional | clos2017 |
| `mg_bright_counts_25mg` | input | 2.682 | dimensionless | legacy | provisional | thomm2021 |
| `mg_dark_counts_25mg` | input | 0.036 | dimensionless | legacy | provisional | thomm2021 |
| `mg_detection_efficiency_25mg` | input | 0.0056 | dimensionless | legacy | provisional | friedenauer2010 |
| `mg_detection_time_25mg` | input | 3e-05 | s | legacy | provisional | thomm2021 |
| `mg_fluorescence_count_rate_25mg` | input | 2.5e+05 | Hz | legacy | provisional | friedenauer2010 |
| `mg_readout_fidelity_bright_25mg` | benchmark | 0.994 | dimensionless | legacy | provisional | thomm2021 |
| `mg_readout_fidelity_dark_25mg` | benchmark | 0.974 | dimensionless | legacy | provisional | thomm2021 |

## Motion (modes, cooling)

| name | kind | value | units | gen | status | source |
|------|------|-------|-------|-----|--------|--------|
| `doppler_cooled_occupation_25mg` | benchmark | 10 | dimensionless | legacy | provisional | clos2017 |
| `doppler_cooling_limit_25mg` | benchmark | 0.001 | K | legacy | provisional | doerr2024 |
| `mg_rsb_cooled_nbar_axial_lf_25mg` | benchmark | 0.07 | dimensionless | legacy | provisional | thomm2021 |
| `mg_rsb_cooled_nbar_radial_hf_25mg` | benchmark | 0.07 | dimensionless | legacy | provisional | thomm2021 |
| `mg_rsb_cooled_nbar_radial_mf_25mg` | benchmark | 0.11 | dimensionless | legacy | provisional | thomm2021 |
| `mg_sideband_nbar_oc_lf_ma_25mg` | benchmark | 0.27 | dimensionless | freddy | provisional | paula_oc_axial_2026 |
| `omega_radial_com_25mg` | input | 2.88e+06 | Hz | legacy | provisional | wittemer2019 |
| `omega_radial_hf_25mg` | input | 4.5e+06 | Hz | freddy | provisional | doerr2024 |
| `omega_radial_mf_25mg` | input | 3e+06 | Hz | freddy | provisional | doerr2024 |
| `omega_radial_rocking_2ion_25mg` | benchmark | 2.57e+06 | Hz | legacy | provisional | wittemer2019 |
| `omega_z_axial_clos_25mg` | input | 1.915e+06 | Hz | legacy | provisional | clos2017 |
| `omega_z_axial_com_25mg` | input | 1.3e+06 | Hz | legacy | provisional | wittemer2019 |
| `omega_z_axial_stretch_2ion_25mg` | benchmark | 2.23e+06 | Hz | legacy | provisional | wittemer2019 |
| `radial_mode_tilt_25mg` | input | 30 | degree | freddy | provisional | doerr2024 |

## Fields (trap / magnetic)

| name | kind | value | units | gen | status | source |
|------|------|-------|-------|-----|--------|--------|
| `b_field_from_clock_25mg` | benchmark | 0.0005609 | T | freddy | provisional | itano_wineland_1981 |
| `b_field_quantization` | input | 0.000585 | T | legacy | provisional | clos2017 |
| `b_field_quantization_freddy` | input | 0.00055 | T | freddy | provisional | doerr2024 |
| `b_field_zeeman_weber_25mg` | input | 0.0005645 | T | freddy | provisional | weber2025 |

## Control (microwave)

| name | kind | value | units | gen | status | source |
|------|------|-------|-------|-----|--------|--------|
| `mw_rabi_3m1_2m1_doerr` | benchmark | 12376 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3m1_2m1_kaufmann` | benchmark | 6553 | Hz | legacy | provisional | kaufmann2022 |
| `mw_rabi_3m1_2m2_doerr` | benchmark | 7003 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3m1_2m2_kaufmann` | benchmark | 2075 | Hz | legacy | provisional | kaufmann2022 |
| `mw_rabi_3m3_2m2_doerr` | benchmark | 11574 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3m3_2m2_kaufmann` | benchmark | 3497 | Hz | legacy | provisional | kaufmann2022 |
| `mw_rabi_3p0_2m1_doerr` | benchmark | 24121 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3p0_2p0_doerr` | benchmark | 25355 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3p1_2p0_doerr` | benchmark | 45249 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3p1_2p2_doerr` | benchmark | 8591 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3p1_2p2_kaufmann` | benchmark | 4845 | Hz | legacy | provisional | kaufmann2022 |
| `mw_rabi_3p3_2p2_doerr` | benchmark | 59453 | Hz | freddy | provisional | doerr2024 |
| `mw_rabi_3p3_2p2_kaufmann` | benchmark | 23776 | Hz | legacy | provisional | kaufmann2022 |
