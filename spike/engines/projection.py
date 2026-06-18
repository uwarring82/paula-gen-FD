"""
Mode-projection engine — which motional modes each Raman (TPSR) beam combination
addresses, from the effective k-vector direction and the radial-mode tilt.

The Lamb-Dicke coupling of a combination to a mode is proportional to the
projection of the effective k-vector onto that mode's axis,

    eta_mode  ∝  |dk_hat . e_mode|   (a direction cosine in [0, 1]),

so a combination "addresses" a mode exactly when this projection is non-zero.
The axial low-frequency (lf) mode lies along z; the two radial modes (mf, hf)
lie in the trap x-y plane, tilted by `tilt` from the x / y axes (Doerr 2024,
Fig. 2.4). The effective k-vectors come from the ledger (the four
raman_*_combination_25mg records, normalised), the tilt from radial_mode_tilt_25mg.

Geometry reproduces Doerr's documented addressing:
  CC  (dk ~ 0)          -> carrier only (no mode)
  OC  (+z)              -> axial lf
  AC  (-(x+z)/sqrt2)    -> all three (lf at 45 deg; mf, hf via the tilted axes)
  ROC (+x)              -> radial mf, hf
This is the *geometric* part only: absolute eta needs |Delta_k| and the mode
frequency (a future sideband engine); here every projection is parameter-free
(directions + tilt).

Two modelling notes:
  * The single direct dot product equals Doerr's pedagogical "two-step" radial
    projection (onto x, then onto the tilted radial axis) ONLY because the
    recorded Delta_k all lie in the x-z plane (Delta_k_y = 0) and the radial
    modes lie in the x-y plane (e_z = 0), so the cross terms vanish. It is not a
    separate algorithm and would NOT coincide for an out-of-x-z-plane k-vector.
  * projection() returns |Delta_k_hat . e_mode|: the sign (the relative
    spin-dependent displacement phase across modes) is intentionally discarded.
    Correct for "is the mode addressed?" and the relative Lamb-Dicke MAGNITUDE;
    it must be restored if mode couplings are ever combined coherently.
"""
from __future__ import annotations

import math

MODES = ("lf", "mf", "hf")

# Canonical ledger record names for the four TPSR combinations.
_COMBO_RECORDS = {
    "CC": "raman_cc_combination_25mg",
    "OC": "raman_oc_combination_25mg",
    "AC": "raman_ac_combination_25mg",
    "ROC": "raman_roc_combination_25mg",
}


def mode_axes(tilt_deg: float) -> dict[str, tuple[float, float, float]]:
    """Unit vectors of the three modes in the trap (x, y, z) frame.

    lf is axial (along z); mf and hf are the two orthogonal radial modes in the
    x-y plane, mf at `tilt` from +x and hf at `tilt` from +y (i.e. perpendicular
    to mf). mf/hf label which physical frequency (3.0 / 4.5 MHz) sits on which
    tilted axis — the geometry only fixes the two orthogonal radial directions.
    """
    t = math.radians(tilt_deg)
    return {
        "lf": (0.0, 0.0, 1.0),
        "mf": (math.cos(t), math.sin(t), 0.0),
        "hf": (-math.sin(t), math.cos(t), 0.0),
    }


def _dot(a, b) -> float:
    return sum(x * y for x, y in zip(a, b))


def _unit(v) -> tuple[float, ...]:
    n = math.sqrt(_dot(v, v))
    return tuple(v) if n == 0.0 else tuple(x / n for x in v)


class ModeProjection:
    """Geometric projection of the four TPSR k-vectors onto the motional modes."""

    def __init__(self, tilt_deg: float, dk_directions: dict[str, tuple]):
        self.tilt_deg = float(tilt_deg)
        self.axes = mode_axes(tilt_deg)
        # normalise dk directions (the CC zero vector stays zero -> carrier only);
        # guard the dimensionality so a malformed record can't silently zip-truncate
        self.dk = {}
        for k, v in dk_directions.items():
            vec = [float(c) for c in v]
            if len(vec) != 3:
                raise ValueError(
                    f"Delta_k '{k}' must be a 3-vector (trap x,y,z), got {len(vec)} component(s)"
                )
            self.dk[k] = _unit(vec)

    @classmethod
    def from_ledger(cls, ledger, tilt_name: str = "radial_mode_tilt_25mg",
                    combos: dict[str, str] | None = None):
        """Build from the ledger: the radial tilt + the four combination Delta_k
        vectors, all `input` (wall-enforced via input_quantity)."""
        combos = combos or _COMBO_RECORDS
        tilt = ledger.input_quantity(tilt_name).value
        dk = {label: ledger.input_quantity(name).value for label, name in combos.items()}
        return cls(tilt_deg=tilt, dk_directions=dk)

    def projection(self, comb: str, mode: str) -> float:
        """Direction cosine |dk_hat . e_mode| in [0, 1] — the relative Lamb-Dicke
        coupling of combination `comb` onto motional mode `mode`."""
        return abs(_dot(self.dk[comb], self.axes[mode]))

    def projections(self, comb: str) -> dict[str, float]:
        return {m: self.projection(comb, m) for m in MODES}

    def addressed_modes(self, comb: str, threshold: float = 1e-9) -> tuple[str, ...]:
        """Modes with non-zero projection — those `comb` can drive sidebands on.
        threshold 1e-9 separates the clean-orthogonality FP residual (~1e-16) from
        any physical projection (the smallest here is 0.354)."""
        return tuple(m for m in MODES if self.projection(comb, m) > threshold)
