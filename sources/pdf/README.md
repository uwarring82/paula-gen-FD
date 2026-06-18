# Local source PDFs — deliberately untracked

The thesis/dissertation/paper PDFs live **here on your machine but are not
committed** to git (`.gitignore` excludes `sources/pdf/*` except this file).

Why (see [ADR-0002](../../docs/decisions/0002-source-pdfs-out-of-vcs.md)):

- they are large binaries (~148 MB) and would bloat history irreversibly;
- they are (in part) third-party copyrighted — the repository references them,
  it does not redistribute them;
- FAIR *Findability/Accessibility* is better served by a resolvable
  permalink/DOI than by an opaque blob.

The canonical, citable list of these sources — with stable citation keys and
resolvable links/DOIs — is [`registries/sources.yaml`](../../registries/sources.yaml).
Each parameter record's `source.loc` (page / table / equation / section) is
precise enough to recover the number once you have the corresponding PDF in hand.

To re-verify an extraction, obtain the PDF via its registry link and drop it in
this folder; tooling and the `source.loc` references will line up.
