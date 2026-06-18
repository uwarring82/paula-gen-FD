"""
Ledger loader for the twin spike.

Loads the YAML parameter records into a queryable in-memory ledger. The custom
loader mirrors validator/validate.py (it parses unsigned-exponent scientific
notation such as 2.1954e11 as a float and keeps ISO dates as strings) so the
spike stays decoupled from the validator (only pyyaml is required).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import yaml


# --- YAML loader: 2.1954e11 -> float; 2026-06-18 -> str ---------------------
class _RecordLoader(yaml.SafeLoader):
    pass


_RecordLoader.yaml_implicit_resolvers = {
    k: list(v) for k, v in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
for _ch, _mappers in list(_RecordLoader.yaml_implicit_resolvers.items()):
    _RecordLoader.yaml_implicit_resolvers[_ch] = [
        (tag, rx) for (tag, rx) in _mappers if tag != "tag:yaml.org,2002:timestamp"
    ]
_FLOAT_RX = re.compile(
    r"""^[-+]?(?:
            [0-9][0-9_]*\.[0-9_]*(?:[eE][-+]?[0-9]+)?
          | \.[0-9_]+(?:[eE][-+]?[0-9]+)?
          | [0-9][0-9_]*[eE][-+]?[0-9]+
        )$
        | ^[-+]?\.(?:inf|Inf|INF)$
        | ^\.(?:nan|NaN|NAN)$""",
    re.X,
)
_RecordLoader.add_implicit_resolver("tag:yaml.org,2002:float", _FLOAT_RX, list("-+0123456789."))


@dataclass(frozen=True)
class Quantity:
    """A resolved numeric value with its 1-sigma uncertainty (SI), in the units
    declared by the record. Uncertainty is symmetrised to a scalar here."""

    name: str
    value: float
    units: str
    sigma: float
    kind: str  # "input" | "benchmark"

    def __repr__(self) -> str:
        return f"{self.name} = {self.value:.6g} +/- {self.sigma:.2g} {self.units} [{self.kind}]"


def _sigma(unc) -> float:
    if not isinstance(unc, dict):
        return float("nan")
    v = unc.get("value")
    if isinstance(v, dict):  # asymmetric {lower, upper} -> max half-width
        return max(abs(float(v.get("lower", 0.0))), abs(float(v.get("upper", 0.0))))
    return abs(float(v)) if v is not None else float("nan")


class Ledger:
    """All parameter records, indexed by name."""

    def __init__(self, records: dict):
        self._records = records

    @classmethod
    def load(cls, root=None) -> "Ledger":
        root = Path(root) if root else Path(__file__).resolve().parent.parent
        records: dict = {}
        for path in sorted((root / "records").glob("*.yaml")):
            with path.open("r", encoding="utf-8") as fh:
                doc = yaml.load(fh, Loader=_RecordLoader) or []
            for rec in doc:
                if isinstance(rec, dict) and "name" in rec:
                    records[rec["name"]] = rec
        return cls(records)

    def __contains__(self, name) -> bool:
        return name in self._records

    def record(self, name: str) -> dict:
        return self._records[name]

    def quantity(self, name: str) -> Quantity:
        r = self._records[name]
        return Quantity(
            name=name,
            value=float(r["value"]),
            units=r.get("units", ""),
            sigma=_sigma(r.get("uncertainty")),
            kind=r.get("kind", ""),   # a missing kind is never silently treated as "input"
        )

    def value(self, name: str) -> float:
        return float(self._records[name]["value"])

    def input_quantity(self, name: str) -> Quantity:
        """A quantity that MUST be `kind: input` — the wall: engines consume only
        inputs. Raises if the record is a benchmark (or has no/!=input kind)."""
        q = self.quantity(name)
        if q.kind != "input":
            raise ValueError(
                f"'{name}' must be kind:input to be consumed by an engine, but kind={q.kind!r}"
            )
        return q

    def benchmark_quantity(self, name: str) -> Quantity:
        """A quantity that MUST be `kind: benchmark` — a validation target, never
        consumed. Raises if the record is not a benchmark."""
        q = self.quantity(name)
        if q.kind != "benchmark":
            raise ValueError(
                f"'{name}' must be kind:benchmark to be a validation target, but kind={q.kind!r}"
            )
        return q

    def names(self):
        return list(self._records)

    def by_kind(self, kind: str):
        return [n for n, r in self._records.items() if r.get("kind") == kind]
