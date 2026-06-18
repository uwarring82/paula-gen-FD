"""
CLI entry point for the twin validation.

The orchestration lives in spike/runner.py (the composition root); this module
is kept as the stable command:

    python -m spike.validate_twin

It loads the ledger, runs every registered engine validation against its
benchmark, prints one table, and exits nonzero if any result is in tension.
"""
from __future__ import annotations

import sys

from .runner import main

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
