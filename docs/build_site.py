"""
Build the static documentation site under docs/ (GitHub Pages friendly):

  docs/index.html               -- homepage: short repo intro + featured technical note(s)
  docs/notes/<name>.html        -- each technical note, rendered from its Markdown

Rendering uses **pandoc** (gfm + `$`-math via MathJax + implicit figure captions). Each note's
HTML is written **next to its .md** (in docs/notes/), so every relative link inside the note
(`../figures/*.png`, `../../registries/*`, sibling `*.md`, ...) stays valid with no rewriting.

Re-run after editing a note or regenerating figures (`python docs/figures/make_aom_rise_figs.py`):

    python docs/build_site.py

Requires pandoc on PATH (the generated HTML itself needs no tooling to view — only a browser,
which pulls MathJax from a CDN to typeset the equations).
"""
from __future__ import annotations

import html
import os
import re
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))            # .../docs

# --- the notes to publish ---------------------------------------------------
NOTES = [
    dict(
        md="notes/aom_finite_sound_velocity_rabi.md",
        out="notes/aom_finite_sound_velocity_rabi.html",
        title="Finite acoustic-transit effects on the stroboscopic R2 Rabi rate",
        blurb=("How the finite speed of sound in the single-pass Raman AOM shapes the "
               "stroboscopic drive: why switching only R2 preserves the spin-rotation area, "
               "why the ion still sees an effectively longer pulse with a hard floor, the "
               "single- vs double-pass comparison, and the detuning-comb fingerprint."),
        thumb="figures/aom_strobe_sequence.png",
        date="2026-06-25",
    ),
    dict(
        md="notes/strobo_grating_transfer_function.md",
        out="notes/strobo_grating_transfer_function.html",
        title="The stroboscopic phase grating as a phase-space probe — transfer function",
        blurb=("Self-contained derivation of the measurement transfer function of the Strobo2.0 "
               "active phase grating: which function of the motional state the spin-flip signal "
               "measures, and how a two-pulse Ramsey version reads the characteristic — and hence "
               "Wigner — function from qubit populations alone."),
        thumb="figures/grating_ramsey_phasespace.png",
        date="2026-06-23",
    ),
]

# --- styling (shared by every page) -----------------------------------------
CSS = """
:root{--ink:#1a1f29;--muted:#5b6675;--line:#e2e6ec;--bg:#f6f7f9;--card:#fff;
      --accent:#1c4f8b;--accent2:#0b7285;--code:#f2f4f7;}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);
     font:16px/1.62 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
a{color:var(--accent);text-decoration:none}
a:hover{text-decoration:underline}
.wrap{max-width:900px;margin:0 auto;padding:0 22px 72px}
.topbar{font-size:14px;color:var(--muted);padding:14px 0;border-bottom:1px solid var(--line);margin-bottom:8px}
.topbar a{color:var(--muted);font-weight:600}
/* hero */
.hero{padding:38px 0 8px}
.hero h1{font-size:2.0rem;line-height:1.18;margin:.1em 0 .15em}
.hero .sub{font-size:1.12rem;color:var(--muted);margin:.2em 0 .6em;max-width:46em}
.badges{font-size:.82rem;color:var(--muted)}
.badges code{background:var(--code);padding:1px 6px;border-radius:5px}
/* content */
.content h1{font-size:1.85rem;line-height:1.2;margin:.2em 0 .5em}
.content h2{font-size:1.32rem;margin:1.9em 0 .5em;padding-bottom:.22em;border-bottom:2px solid var(--line)}
.content h3{font-size:1.08rem;margin:1.5em 0 .4em;color:var(--accent)}
.content p,.content li{max-width:48em}
.content img{max-width:100%;height:auto}
figure{margin:1.6em 0;text-align:center}
figure img{border:1px solid var(--line);border-radius:8px;box-shadow:0 1px 4px rgba(20,30,50,.07);background:#fff}
figcaption{font-size:.86rem;color:var(--muted);font-style:italic;margin-top:.5em;max-width:46em;margin-left:auto;margin-right:auto}
blockquote{margin:1.2em 0;padding:.4em 1.1em;border-left:4px solid var(--accent2);background:#eef6f7;color:#2b3b40;border-radius:0 6px 6px 0}
code{background:var(--code);padding:1.5px 5px;border-radius:5px;font-size:.9em;
     font-family:"SF Mono",ui-monospace,Menlo,Consolas,monospace}
pre{background:#0f1729;color:#e6edf3;padding:14px 16px;border-radius:9px;overflow-x:auto;font-size:.86rem;line-height:1.5}
pre code{background:none;color:inherit;padding:0}
table{border-collapse:collapse;width:100%;margin:1.3em 0;font-size:.92rem;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:7px 11px;text-align:left;vertical-align:top}
thead th{background:#eef1f5}
tbody tr:nth-child(even){background:#fafbfc}
hr{border:0;border-top:1px solid var(--line);margin:2.2em 0}
.foot{margin-top:54px;padding-top:18px;border-top:1px solid var(--line);font-size:.84rem;color:var(--muted)}
/* homepage cards */
.lead{font-size:1.06rem}
.wall{margin:1.5em 0;padding:16px 20px;background:var(--card);border:1px solid var(--line);
      border-left:4px solid var(--accent);border-radius:0 9px 9px 0}
.wall b{color:var(--accent)}
.cards{display:grid;gap:20px;margin:1.2em 0}
.card{display:grid;grid-template-columns:230px 1fr;gap:18px;background:var(--card);border:1px solid var(--line);
      border-radius:11px;overflow:hidden;transition:box-shadow .15s,transform .15s}
.card:hover{box-shadow:0 6px 20px rgba(20,30,50,.10);transform:translateY(-1px)}
.card img{width:100%;height:100%;object-fit:cover;background:#fff;border-right:1px solid var(--line)}
.card .body{padding:16px 18px}
.card h3{margin:.1em 0 .3em;font-size:1.12rem}
.card h3 a{color:var(--ink)}
.card .blurb{font-size:.93rem;color:var(--muted);margin:.2em 0 .6em}
.card .meta{font-size:.8rem;color:var(--muted)}
.card .go{font-weight:600}
@media(max-width:620px){.card{grid-template-columns:1fr}.card img{height:170px;border-right:0;border-bottom:1px solid var(--line)}
  .hero h1{font-size:1.6rem}}
"""

MATHJAX = (
    "<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],"
    "displayMath:[['\\\\[','\\\\]']]},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};</script>\n"
    "<script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>"
)


def _page(title: str, body: str) -> str:
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>{html.escape(title)}</title>\n"
        f"<style>{CSS}</style>\n{MATHJAX}\n</head>\n<body>\n{body}\n</body>\n</html>\n"
    )


# A standalone image paragraph -> <figure> with the alt text as <figcaption>. We do this
# OURSELVES (instead of pandoc's +implicit_figures) so the output is identical regardless of
# the pandoc version: older pandocs ignore/!=handle implicit_figures, which silently dropped
# the captions on the CI runner. Without the extension every pandoc emits <p><img .../></p>.
_IMG_P = re.compile(r"<p>\s*(<img\b[^>]*>)\s*</p>")


def _wrap_figures(frag: str) -> str:
    def repl(m: "re.Match") -> str:
        img = m.group(1)
        alt = re.search(r'alt="([^"]*)"', img)
        cap = f"<figcaption>{alt.group(1)}</figcaption>" if (alt and alt.group(1)) else ""
        return f"<figure>{img}{cap}</figure>"
    return _IMG_P.sub(repl, frag)


def _pandoc(md_abs: str) -> str:
    out = subprocess.run(
        ["pandoc", md_abs, "-f", "gfm+tex_math_dollars",
         "-t", "html", "--mathjax", "--no-highlight"],
        capture_output=True, text=True, check=True,
    )
    return _wrap_figures(out.stdout)


def build_note(note: dict) -> None:
    md_abs = os.path.join(_HERE, note["md"])
    body_frag = _pandoc(md_abs)
    topbar = ('<div class="wrap"><div class="topbar">'
              '<a href="../index.html">← iontrap-reference</a> &nbsp;·&nbsp; technical note'
              '</div></div>')
    body = (topbar + '<div class="wrap"><article class="content">\n' + body_frag +
            '<div class="foot">Rendered from <code>' + html.escape(note["md"]) +
            '</code> by <code>docs/build_site.py</code> (pandoc + MathJax). '
            'Data &amp; docs CC-BY-4.0, code MIT. '
            '<a href="../index.html">← back to home</a></div>\n</article></div>')
    out_abs = os.path.join(_HERE, note["out"])
    with open(out_abs, "w", encoding="utf-8") as fh:
        fh.write(_page(note["title"], body))
    print("wrote", os.path.relpath(out_abs))


def build_home() -> None:
    cards = []
    for n in NOTES:
        cards.append(
            '<div class="card">'
            f'<a href="{html.escape(n["out"])}"><img src="{html.escape(n["thumb"])}" '
            f'alt="{html.escape(n["title"])}"></a>'
            '<div class="body">'
            f'<h3><a href="{html.escape(n["out"])}">{html.escape(n["title"])}</a></h3>'
            f'<p class="blurb">{html.escape(n["blurb"])}</p>'
            f'<p class="meta">{html.escape(n["date"])} &nbsp;·&nbsp; '
            f'<a class="go" href="{html.escape(n["out"])}">read the note →</a></p>'
            '</div></div>'
        )
    body = (
        '<div class="wrap">\n'
        '<header class="hero">\n'
        '<h1>iontrap-reference</h1>\n'
        '<p class="sub">A citable source-of-truth for the Freiburg ²⁵Mg⁺ '
        '(PAULA) trapped-ion apparatus — a schema-validated parameter layer plus a '
        'digital-twin spike that consumes it.</p>\n'
        '<p class="badges">Data &amp; docs <code>CC-BY-4.0</code> · code <code>MIT</code> '
        '· <a href="../README.md">README</a> · '
        '<a href="STATE_OF_THE_TWIN.md">State of the Twin</a> · '
        '<a href="LOGBOOK.md">Logbook</a></p>\n'
        '</header>\n'
        '<section class="content">\n'
        '<p class="lead">The lineage theses and papers are the human-readable source; on top '
        'sits a structured, schema-validated <b>parameter layer</b> that cites <em>into</em> them, '
        'so every number is traceable to where it was measured or derived. A digital-twin '
        '<b>spike</b> (physics engines + twin compositions) consumes that layer and is validated '
        'against held-out benchmarks.</p>\n'
        '<div class="wall"><b>The one idea — the input / benchmark wall.</b> Every record is '
        'either an <b>input</b> the twin consumes (beam powers, detunings, atomic constants …) '
        'or a <b>benchmark</b> it must reproduce and is <em>never</em> fed into fitting (Stark '
        'shifts, decoherence/heating rates, transition-frequency residuals …). The separation '
        'is mandatory and machine-enforced, so calibration data cannot leak into the '
        'parameterisation and the validation actually proves something.</div>\n'
        '<h2>Technical notes</h2>\n'
        '<p>Self-contained derivations and models, rendered with figures and equations in place.</p>\n'
        '<div class="cards">\n' + "\n".join(cards) + '\n</div>\n'
        '</section>\n'
        '<div class="foot">Generated by <code>docs/build_site.py</code>. The parameter layer carries no '
        'physics of its own; the engines, σ-validations and tomography live in the '
        '<a href="../spike/">spike/</a> sibling. See the '
        '<a href="STATE_OF_THE_TWIN.md">State of the Twin</a> for what is validated vs diagnostic.</div>\n'
        '</div>'
    )
    out_abs = os.path.join(_HERE, "index.html")
    with open(out_abs, "w", encoding="utf-8") as fh:
        fh.write(_page("iontrap-reference — ²⁵Mg⁺ (PAULA) trap", body))
    print("wrote", os.path.relpath(out_abs))


def main() -> int:
    if not shutil.which("pandoc"):
        print("ERROR: pandoc not found on PATH (needed to render Markdown notes).", file=sys.stderr)
        return 1
    for n in NOTES:
        build_note(n)
    build_home()
    print("\nOpen docs/index.html in a browser (MathJax typesets the equations from a CDN).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
