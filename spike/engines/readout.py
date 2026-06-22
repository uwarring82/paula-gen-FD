"""
Readout engine — turns the detection count model into state-readout figures of
merit: single-shot discrimination fidelity, the Fisher-information / Cramer-Rao
precision of a maximum-likelihood P_down estimate over N shots, and the optimal
detection time.

It builds on engines.detection.transition_count_pmf (Poisson core + depumping /
leaking tails). Two honest caveats about the Thomm 2021 numbers:

 * The reported state fidelities (99.4 % bright, 97.4 % dark) are ENSEMBLE
   state-PREPARATION+readout probabilities, dominated by preparation: a perfect
   single-shot detector cannot exceed the prepared-state purity. The DETECTION
   contribution alone (single-shot discrimination at lambda_down=2.682,
   lambda_up=0.036) is ~95 %, and averaging many shots drives the readout
   contribution far below the ~0.6 % preparation error. So this engine predicts the
   DETECTION-limited readout, not the prep-limited 99.4 % — they are different
   quantities (a diagnostic, not a sigma test of 99.4 %).
 * The optimal t_det (Thomm 30 us) is the point where the lesser of the rising
   bright fidelity and the leak-limited dark fidelity peaks; locating it needs the
   depumping/leak RATES, which are not yet ledger-calibrated.

Pure Python (math only).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .detection import optimal_threshold, transition_count_pmf


def single_shot_fidelity(lam_bright: float, lam_dark: float,
                         depump_bright: float = 0.0, leak_dark: float = 0.0,
                         threshold: int | None = None) -> float:
    """Single-shot discrimination fidelity F = 1 - 1/2[P(<n_c|bright)+P(>=n_c|dark)]
    using the REALISTIC bright/dark count distributions (transition_count_pmf), i.e.
    including depumping (bright low-count tail) and leaking (dark high-count tail).
    threshold None -> the Poisson-optimal n_c."""
    if threshold is None:
        threshold = optimal_threshold(lam_bright, lam_dark)[0]
    kmax = int(lam_bright + 10.0 * math.sqrt(lam_bright + 1.0)) + 5
    e_bright = sum(transition_count_pmf(k, lam_bright, lam_dark, depump_bright)
                   for k in range(0, threshold))                          # bright read as dark
    e_dark = sum(transition_count_pmf(k, lam_dark, lam_bright, leak_dark)
                 for k in range(threshold, kmax + 1))                     # dark read as bright
    return 1.0 - 0.5 * (e_bright + e_dark)


def fisher_information_p_down(p: float, lam_bright: float, lam_dark: float,
                             depump_bright: float = 0.0, leak_dark: float = 0.0) -> float:
    """Fisher information per shot for the mixture weight P_down, with the realistic
    component PMFs B, D: I(p) = sum_k (B(k)-D(k))^2 / (p B(k) + (1-p) D(k)). For
    perfectly separated B, D this reduces to 1/(p(1-p)) (the binomial/QPN limit);
    overlap (small lambda_bright) reduces it."""
    kmax = int(lam_bright + 10.0 * math.sqrt(lam_bright + 1.0)) + 5
    info = 0.0
    for k in range(0, kmax + 1):
        b = transition_count_pmf(k, lam_bright, lam_dark, depump_bright)
        d = transition_count_pmf(k, lam_dark, lam_bright, leak_dark)
        mix = p * b + (1.0 - p) * d
        if mix > 0.0:
            info += (b - d) ** 2 / mix
    return info


def p_down_uncertainty(p: float, n_shots: int, lam_bright: float, lam_dark: float,
                       depump_bright: float = 0.0, leak_dark: float = 0.0) -> float:
    """Cramer-Rao lower bound on the maximum-likelihood P_down estimate from n_shots
    fluorescence measurements: sigma(P_down) = 1/sqrt(N * I(p)). Compare the ideal
    quantum-projection noise sqrt(p(1-p)/N): the ratio is the readout overhead from
    finite photon counts."""
    info = fisher_information_p_down(p, lam_bright, lam_dark, depump_bright, leak_dark)
    if info <= 0.0 or n_shots <= 0:
        return float("nan")
    return 1.0 / math.sqrt(n_shots * info)


def qpn_uncertainty(p: float, n_shots: int) -> float:
    """Ideal quantum-projection noise sqrt(p(1-p)/N) (a perfect single-shot readout)."""
    if n_shots <= 0:
        return float("nan")
    return math.sqrt(p * (1.0 - p) / n_shots)


@dataclass
class ReadoutModel:
    """Readout figures of merit at a fixed detection time, from the ledger's measured
    bright/dark count levels (Thomm 2021)."""
    lam_bright: float
    lam_dark: float
    depump_bright: float = 0.0
    leak_dark: float = 0.0

    @classmethod
    def from_ledger(cls, ledger, depump_bright: float = 0.0, leak_dark: float = 0.0):
        # The bright/dark COUNT LEVELS are inputs (the readout FIDELITIES, compared
        # against separately, are the benchmarks). Use input_quantity so the wall is
        # explicit: this diagnostic consumes only inputs.
        return cls(lam_bright=ledger.input_quantity("mg_bright_counts_25mg").value,
                   lam_dark=ledger.input_quantity("mg_dark_counts_25mg").value,
                   depump_bright=depump_bright, leak_dark=leak_dark)

    def single_shot_fidelity(self) -> float:
        return single_shot_fidelity(self.lam_bright, self.lam_dark,
                                    self.depump_bright, self.leak_dark)

    def p_down_uncertainty(self, p: float, n_shots: int) -> float:
        return p_down_uncertainty(p, n_shots, self.lam_bright, self.lam_dark,
                                  self.depump_bright, self.leak_dark)
