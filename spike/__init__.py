"""
iontrap-reference twin SPIKE.

A throwaway-ish end-to-end proof that the structured parameter layer can drive
physics engines: each engine consumes `input` records from the ledger and is
checked against `benchmark` records. Per the task card, solver/physics code is
out of scope for the substrate (records/ + validator/); this package is kept
separate and may later graduate to its own repo(s) (iontrap-levels / -fields /
-optics + a twin composition root).

First engine: `spike.engines.levels` — the 25Mg+ ground-state hyperfine+Zeeman
structure (Breit-Rabi), which reproduces the clock_transition_25mg benchmark.
"""
