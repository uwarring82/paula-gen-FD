"""
Reader for the PAULA DAQ `.dat` files (docs/DATA_FORMAT.md; format documented in
kalis2017).

A file is an XML metadata header (every line `#`-commented) — parameter settings
(`<ionproperties>`), the scanned-variable definition (`<parameter>`), waveform
profiles, shim fields and the control script — followed by un-commented blocks:

  * the scan DATA (tab-separated: scanned_var, value, error, timestamp_ms),
    written as two consecutive COUNTER blocks (the scanned variable restarts at the
    block boundary): the ion-fluorescence signal, and a near-constant ~0.013
    reference/normalisation counter (NOT the dark-ion state);
  * per-scan-point photon-count HISTOGRAMS (tab-separated: count, #shots) — one
    group per (scan point, counter), i.e. signal histograms AND reference histograms.

This reader extracts the settings, the scan definition, the signal series (the
higher-variance counter block — the reference is the near-constant one), and the
histograms. The bright signal-histogram mean independently corroborates the signal.
"""
from __future__ import annotations

import re

_ITEM = re.compile(r'name="([^"]+)"><value>([^<]*)</value>')
_SCAN = re.compile(
    r'name="([^"]+)">'
    r'<value_l>([^<]*)</value_l><value_u>([^<]*)</value_u>'
    r'<int_points>([^<]*)</int_points><exp_point>([^<]*)</exp_point>'
)


def _maybe_float(s: str):
    try:
        return float(s)
    except ValueError:
        return s


class DatFile:
    """A parsed DAQ scan file."""

    def __init__(self, path: str):
        self.path = path
        self.settings: dict = {}
        self.scan: dict | None = None          # {name, lower, upper, points, shots}
        data_rows = []                          # (x, y, err, t)
        hist_rows = []                          # (count, occ)
        with open(path, encoding="latin-1") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    m = _ITEM.search(ln)
                    if m:
                        self.settings[m.group(1)] = _maybe_float(m.group(2))
                    s = _SCAN.search(ln)
                    if s:
                        self.scan = {
                            "name": s.group(1),
                            "lower": float(s.group(2)), "upper": float(s.group(3)),
                            "points": int(float(s.group(4))), "shots": int(float(s.group(5))),
                        }
                    continue
                parts = ln.split()
                if len(parts) == 4:
                    try:
                        vals = tuple(float(v) for v in parts)
                    except ValueError:
                        continue
                    if all(v == v and abs(v) != float("inf") for v in vals):   # reject NaN/Inf
                        data_rows.append(vals)
                elif len(parts) == 2:
                    try:
                        hist_rows.append((int(float(parts[0])), int(float(parts[1]))))
                    except ValueError:
                        pass
        self._data = data_rows
        self._hist_rows = hist_rows

    @property
    def timestamp(self):
        """Measurement date/time 'YYYY-MM-DD HH:MM:SS' parsed from the filename
        (HH_MM_SS_DD_MM_YYYY), corroborated by the <metadata label> in the header.
        The data are the group's OWN primary measurement (the file is self-describing);
        kalis2017 only documents the .dat FORMAT, it is not the data's origin."""
        import os
        stem = os.path.basename(str(self.path)).rsplit(".", 1)[0]
        m = re.match(r"(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{4})$", stem)
        if not m:
            return None
        hh, mm, ss, dd, mo, yy = m.groups()
        return f"{yy}-{mo}-{dd} {hh}:{mm}:{ss}"

    # --- scan signal --------------------------------------------------------
    def signal(self):
        """(x, y, sigma) of the ion-fluorescence counter, sorted by x.

        The two counters are written as consecutive blocks (the scanned variable
        restarts at the boundary); the fluorescence signal is the higher-variance
        block (the reference counter sits at a near-constant ~0.013 floor). This is
        a structural split, robust to a deep resonance dip (unlike a per-point
        magnitude pick). err (col 3) is the per-point standard error."""
        if not self._data:
            return [], [], []
        blocks = [[self._data[0]]]
        for row in self._data[1:]:
            if len(blocks[-1]) > 1 and row[0] <= blocks[-1][0][0]:
                blocks.append([row])      # scanned variable reset -> next counter block
            else:
                blocks[-1].append(row)

        def _var(blk):
            ys = [r[1] for r in blk]
            m = sum(ys) / len(ys)
            return sum((y - m) ** 2 for y in ys) / len(ys)

        sig = sorted(max(blocks, key=_var), key=lambda r: r[0])
        return [r[0] for r in sig], [r[1] for r in sig], [r[2] for r in sig]

    def counter_blocks(self):
        """ALL counter blocks in FILE ORDER, each as (x, y, sigma) sorted by x. Unlike
        signal() (which keeps only the highest-variance block and treats the rest as a
        reference), this returns every block -- needed when BOTH counters carry real
        signal, e.g. a sequence that drives a BSB pulse (counter 0) AND an RSB pulse
        (counter 1) and gates each into its own detection window (the 1R_LF_MA / BSB_RSB
        scans). The caller maps blocks to physical meaning by apparatus convention."""
        if not self._data:
            return []
        blocks = [[self._data[0]]]
        for row in self._data[1:]:
            if len(blocks[-1]) > 1 and row[0] <= blocks[-1][0][0]:
                blocks.append([row])
            else:
                blocks[-1].append(row)
        out = []
        for blk in blocks:
            srt = sorted(blk, key=lambda r: r[0])
            out.append(([r[0] for r in srt], [r[1] for r in srt], [r[2] for r in srt]))
        return out

    # --- per-shot count histograms -----------------------------------------
    def histograms(self):
        """List of per-scan-point count histograms, each a dict {count: n_shots},
        split on the blank lines between groups (a group restarts at count 0)."""
        groups = []
        cur: dict = {}
        for count, occ in self._hist_rows:
            if count == 0 and cur:
                groups.append(cur)
                cur = {}
            cur[count] = cur.get(count, 0) + occ
        if cur:
            groups.append(cur)
        return groups

    def point_shots(self):
        """[(x, [shot_count, ...])] for the ion-fluorescence counter: the individual
        per-shot photon counts at each scan point, reconstructed from the per-point
        histograms. histogram[i] pairs with the i-th data row (file order); the signal
        counter is the higher-variance of the two npts-row blocks (the reference is the
        near-constant ~0.013 floor)."""
        sc = self.scan
        if not sc:
            return []
        npts = sc["points"]
        hists = self.histograms()
        rows = self._data

        def _var(block):
            m = [r[1] for r in block]
            mu = sum(m) / len(m) if m else 0.0
            return sum((v - mu) ** 2 for v in m) / len(m) if m else 0.0

        if len(rows) >= 2 * npts and len(hists) >= 2 * npts and _var(rows[npts:2 * npts]) > _var(rows[:npts]):
            block_rows, block_h = rows[npts:2 * npts], hists[npts:2 * npts]
        else:
            n = min(npts, len(rows), len(hists))
            block_rows, block_h = rows[:n], hists[:n]
        out = []
        for r, hist in zip(block_rows, block_h):
            shots = [c for c, occ in sorted(hist.items()) for _ in range(occ)]
            out.append((r[0], shots))
        return out

    @staticmethod
    def hist_mean(hist: dict) -> float:
        n = sum(hist.values())
        return sum(c * occ for c, occ in hist.items()) / n if n else float("nan")

    @staticmethod
    def hist_variance(hist: dict) -> float:
        n = sum(hist.values())
        if not n:
            return float("nan")
        m = DatFile.hist_mean(hist)
        return sum(occ * (c - m) ** 2 for c, occ in hist.items()) / n
