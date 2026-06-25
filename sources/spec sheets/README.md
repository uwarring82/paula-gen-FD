# Spec sheets — local only

Manufacturer datasheets / instruction manuals for apparatus hardware. Like the source
PDFs (`sources/pdf/`), these are **third-party copyrighted** and are **kept local, not
committed** (ADR-0002); `.gitignore` excludes everything here except this README.

Each is registered in [`registries/sources.yaml`](../../registries/sources.yaml) so records
and notes can cite it. Currently:

- `intraact_aom220uv_specs.pdf` — IntraAction **ASM-2202B3** acousto-optic modulator
  (UV fused silica, 220 MHz; the R2/B1 Raman-beam AOM). Registry key `intraaction_asm2202b3`.
  Source of the finite-sound-velocity Rabi model in
  [`docs/notes/aom_finite_sound_velocity_rabi.md`](../../docs/notes/aom_finite_sound_velocity_rabi.md).
