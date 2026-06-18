"""
Fundamental and atomic constants for the twin spike (SI units).

Fundamental constants are CODATA 2018. The atomic g-factors are species
constants; they are NOT yet ledger records — TODO: graduate g_J / g_I (and the
nuclear magnetic moment) to `input` records so the wall covers them too. Until
then they live here, sourced, and the levels engine falls back to these when the
ledger does not provide them.
"""

# --- CODATA 2018 fundamental constants (SI) ---------------------------------
H_PLANCK = 6.62607015e-34        # J*s          (exact, SI definition)
MU_B = 9.2740100783e-24          # J/T          Bohr magneton
MU_N = 5.0507837461e-27          # J/T          nuclear magneton
G_E = 2.00231930436256           # —            free-electron g-factor

# Convenience: magnetons in frequency units
MU_B_OVER_H = MU_B / H_PLANCK     # Hz/T  (= 1.39962449e10; 1.39962 MHz/G)
MU_N_OVER_H = MU_N / H_PLANCK     # Hz/T

# --- 25Mg+ ground-state (3s 2S_1/2) atomic constants ------------------------
# Electronic g-factor of a 2S_1/2 state: L=0, so g_J = g_S = g_e to the level
# this spike needs (sub-ppm QED/relativistic corrections ignored).
G_J_2S12 = G_E

# 25Mg nuclear g-factor: g_I = mu_I / (I * mu_N) with the tabulated moment
# mu_I(25Mg) = -0.85545 mu_N and I = 5/2  ->  g_I = -0.34218.
# (It enters the m_F=0 clock only weakly, ~12 Hz at 5.5 G, via the combination
# g_J mu_B - g_I mu_N; it dominates the m_F != 0 Zeeman slopes.)
# Source: nuclear-moment tables (Stone 2005 / IAEA).
MU_I_25MG = -0.85545             # in units of mu_N
NUCLEAR_SPIN_25MG = 2.5          # I = 5/2 (also in the ledger: nuclear_spin_25mg)
G_I_25MG = MU_I_25MG / NUCLEAR_SPIN_25MG   # = -0.34218

# Unit helpers
GAUSS = 1.0e-4                    # 1 G in Tesla
MHZ = 1.0e6
