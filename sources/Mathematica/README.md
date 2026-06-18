# Internal PAULA analysis notebooks — local, registered, untracked

These Mathematica notebooks are **own primary-source analysis** for the PAULA
apparatus (trap simulations, BEM electrode optimisation, motional-mode
orientations, Mg level/scattering-rate notes, 2014–2018). They are kept **local
and untracked** (`.gitignore` excludes `sources/Mathematica/*` except this file)
because they are unpublished working analysis and `.nb` files diff poorly — but
each is **registered as a citable source** in
[`registries/sources.yaml`](../../registries/sources.yaml) (keys prefixed
`paula_`) so parameter records can cite them as the origin of `simulated` /
`derived` quantities.

| File | Registry key |
|------|--------------|
| `Notes_PAULA_Mg_details_2014_07_10.nb` | `paula_mg_details_2014` |
| `Notes_PAULA_Mg_streurate_2015_07_17.nb` | `paula_mg_scatterrate_2015` |
| `TrapSim_PAULA_BEM_opt_ExpZone_Shims_2017_05_30.nb` | `paula_bem_shims_2017` |
| `TrapSim_PAULA_Mode_orientations_2017_07_11.nb` | `paula_mode_orientations_2017` |
| `Notes_PAULA_TrapSimulations_2018_10_31.nb` | `paula_trapsim_2018` |

A record citing one of these uses a precise `source.loc` (notebook section /
cell / output label). Because these have no public identifier, the validator
WARNS when one is referenced — a visible, accepted accessibility flag.
