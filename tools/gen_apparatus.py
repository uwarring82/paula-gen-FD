#!/usr/bin/env python3
"""
Generate docs/APPARATUS.md — a single, readable overview of EVERY apparatus parameter,
grouped by subsystem, from the records/*.yaml substrate. The records remain the source
of truth (each is self-describing: value, units, source, scope, generation, status);
this file is a convenience index so a reader does not have to grep four YAMLs.

    python tools/gen_apparatus.py        # rewrites docs/APPARATUS.md

It is checked in CI / by test_apparatus_doc.py to stay in sync with the records.
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
from validator.validate import load_yaml  # noqa: E402  (canonical record loader)

_OUT = _ROOT / "docs" / "APPARATUS.md"
# subsystem -> (display heading, order). Records carry scope.subsystem (ADR-0005).
_SUBSYS = [
    ("internal_state", "Internal state (²⁵Mg⁺ hyperfine / Zeeman)"),
    ("optics", "Optics (laser + Raman beams)"),
    ("detection", "Detection / readout"),
    ("motion", "Motion (modes, cooling)"),
    ("fields", "Fields (trap / magnetic)"),
    ("control", "Control (microwave)"),
]


def fmt_value(v) -> str:
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        a = abs(v)
        return f"{v:.4g}" if (a != 0 and (a < 1e-3 or a >= 1e5)) else f"{v:g}"
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(fmt_value(x) for x in v) + "]"
    return str(v).replace("|", "\\|")


def _records():
    out = []
    for path in sorted((_ROOT / "records").glob("*.yaml")):
        data = load_yaml(path) or []
        for r in data:
            if isinstance(r, dict) and "name" in r:
                r["_file"] = path.name
                out.append(r)
    return out


def render(records) -> str:
    by_sub = {}
    for r in records:
        sub = (r.get("scope") or {}).get("subsystem", "other")
        by_sub.setdefault(sub, []).append(r)

    n_prov = sum(1 for r in records if r.get("status") == "provisional")
    lines = [
        "# Apparatus parameters",
        "",
        "**Auto-generated from `records/*.yaml` — do not edit by hand.** Regenerate with",
        "`python tools/gen_apparatus.py`. The records are the source of truth (each carries",
        "its full provenance: source, uncertainty, conditions, caveats); this is a flat",
        "index of every apparatus-specific value, grouped by subsystem.",
        "",
        f"{len(records)} records ({n_prov} provisional, {len(records) - n_prov} confirmed). "
        "`kind`: **input** = consumed by the engines; **benchmark** = a measured/inferred "
        "value the engines are checked against. See "
        "[STATE_OF_THE_TWIN.md](STATE_OF_THE_TWIN.md) and "
        "[SOURCES.md](SOURCES.md).",
    ]
    seen = set()
    for sub, heading in _SUBSYS + [(s, s) for s in sorted(by_sub) if s not in {x[0] for x in _SUBSYS}]:
        recs = by_sub.get(sub)
        if not recs or sub in seen:
            continue
        seen.add(sub)
        lines += ["", f"## {heading}", "",
                  "| name | kind | value | units | gen | status | source |",
                  "|------|------|-------|-------|-----|--------|--------|"]
        for r in sorted(recs, key=lambda x: x["name"]):
            src = (r.get("source") or {}).get("ref", "—")
            lines.append("| `%s` | %s | %s | %s | %s | %s | %s |" % (
                r["name"], r.get("kind", "—"), fmt_value(r.get("value")),
                str(r.get("units", "—")).replace("|", "\\|"), r.get("generation", "—"),
                r.get("status", "—"), src))
    return "\n".join(lines) + "\n"


def main() -> int:
    text = render(_records())
    _OUT.write_text(text, encoding="utf-8")
    print(f"wrote {_OUT.relative_to(_ROOT)} ({text.count(chr(10))} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
