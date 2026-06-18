# 0002 — Source PDFs kept out of version control; cite resolvable links

- **Status:** Accepted
- **Date:** 2026-06-18
- **Deciders:** UW (decision), Claude (implementation)

## Context

The lineage corpus is ~148 MB of thesis/dissertation PDFs across 12 documents.
Three handling options were considered: commit them directly, track via Git LFS,
or keep them local and cite resolvable links. The task card is explicit:
*"register each in the sources registry with a stable citation key + resolvable
link, not just a dumped PDF."* Two further forces:

- **Copyright** — the theses are (in part) third-party copyrighted; the
  repository should not redistribute them.
- **FAIR** — *Findable* and *Accessible* are better served by a resolvable
  permalink/DOI than by a binary blob whose provenance is opaque.

## Decision

The PDFs stay **local and untracked**. `.gitignore` excludes `sources/pdf/*`
(keeping a tracked `sources/pdf/README.md` that explains the policy). The
[`registries/sources.yaml`](../../registries/sources.yaml) registry carries, for
each source, a stable citation key plus a **resolvable link or DOI** (preferring
institutional-repository permalinks, e.g. FreiDok). A pre-commit
`check-added-large-files` hook (2 MB cap) guards against the PDFs ever slipping
into history.

Numerical facts extracted from the sources are not themselves copyrightable; the
curation/structuring is what this repository licenses (see ADR-0003).

## Consequences

- **+** Lean git history; clean clones; no copyright redistribution.
- **+** Every number is traceable to a *resolvable* source, not a local file
  path that means nothing to anyone else.
- **−** A contributor must obtain the PDFs out of band to re-verify an
  extraction. Acceptable: the `source.loc` (page/table/eq) is precise enough to
  recover the number once the PDF is in hand.
- **−** Links rot. Mitigation: prefer DOIs/permalinks; the registry carries a
  `verified` flag per source so unverified links are visible.
