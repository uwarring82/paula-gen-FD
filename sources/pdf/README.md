# Local source PDFs — deliberately untracked

The thesis/dissertation/paper PDFs live **here on your machine but are not
committed** to git (`.gitignore` excludes `sources/pdf/*` except this file).

Why (see [ADR-0002](../../docs/decisions/0002-source-pdfs-out-of-vcs.md)):

- they are large binaries and would bloat history irreversibly;
- many are third-party copyrighted — the repository references them, it does
  not redistribute them;
- FAIR *Findability/Accessibility* is better served by a resolvable
  permalink/DOI than by an opaque blob.

## Naming convention: `filename == citation key`

Every local PDF is named **`<citation-key>.pdf`**, where the citation key is its
entry in [`registries/sources.yaml`](../../registries/sources.yaml). So a record's
`source.ref` maps 1:1 to a local file: `source.ref: clos2017` ⇒ `clos2017.pdf`.
Supplementary material uses the `<parentkey>_suppmat.pdf` form. Degree/type,
title, authors, and the resolvable link/DOI all live in the registry, so the
filename carries no information that would be lost.

Current local PDFs (25), by domain:

- **Theses:** `matjeschk2008` · `schmitz2010` · `enderlein2013` · `pacher2014` ·
  `harlos2015` · `clos2017` · `wittemer2019` · `kaufmann2022` · `friedenauer2010`
  · `weber2025` · `doerr2024` · `hasse2025`
- **Journal papers / refs:** `friedenauer2006` · `friedenauer2008` · `schmitz2009`
  · `schneider2012` · `clos2014` · `clos2016` · `clos2016_suppmat` · `wittemer2018`
  · `wittemer2019_prl` · `wittemer2020` · `hasse2024` · `colla2025` ·
  `itano_wineland_1981`

Note: `wittemer2019` is the PhD thesis; the same-year PRL is `wittemer2019_prl`
(keys disambiguated). `enderlein2013` corrects the old `..._2012` filename
(actual year 2013).

To re-verify an extraction, open `<key>.pdf` here (or fetch it via its registry
link) and go to the record's `source.loc` (page / table / equation / section).
