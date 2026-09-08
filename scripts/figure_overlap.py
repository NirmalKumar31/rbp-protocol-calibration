"""Do any two pieces of text in a figure occupy the same space?

    python scripts/figure_overlap.py            # every figure under results/figures/
    python scripts/figure_overlap.py --frac 0.1 # stricter

Why this exists. `f10_three_protocols` carried a text collision in every build it ever had.
Panel a's left-aligned title measured 3.85in against a 3.10in axes at the rcParams default of
10.8pt, overflowed 0.75in to the right, and landed panel b's own label inside it, so page 17 of
the manuscript read "5.4bx above in 88/94". Fourteen external audit rounds looked at that page,
several of them reporting a page-by-page visual inspection, and none saw it.

Nothing could have. Every check in this repository reads extracted TEXT, and extraction returns
both strings happily: two words at the same coordinates are two words. Overlap is a property of
their POSITIONS, and no gate here had ever looked at a position. A figure is the one artefact
where being right about the numbers and wrong about the layout are independent failures.

So this reads word bounding boxes out of the PDF and reports any pair whose intersection covers
a third of the smaller box. The threshold is not zero because glyph boxes for adjacent characters
legitimately touch, and kerned pairs can share a fraction of a point; a third of a box is well
clear of that and still caught the real defect at 100%.

This is a LAYOUT check and it is not a substitute for looking. It cannot see a line that runs off
the canvas, a legend covering a data point, or an axis label clipped by the figure edge, because
none of those are two words in one place.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

FIGS = ROOT / "results" / "figures"
WORD = re.compile(r'<word xMin="([\d.]+)" yMin="([\d.]+)" '
                  r'xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>')


def words(pdf):
    """Word boxes, in PDF points. Raises if the extractor is absent or fails."""
    r = subprocess.run(["pdftotext", "-bbox", str(pdf), "-"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"pdftotext failed on {pdf.name}: {r.stderr.strip()[:200]}")
    return [(float(a), float(b), float(c), float(d), t)
            for a, b, c, d, t in WORD.findall(r.stdout)]


def overlapping(ws, frac):
    """Pairs sharing at least `frac` of the smaller box. O(n^2); n is small per figure."""
    out = []
    for i, a in enumerate(ws):
        for b in ws[i + 1:]:
            ix = min(a[2], b[2]) - max(a[0], b[0])
            iy = min(a[3], b[3]) - max(a[1], b[1])
            if ix <= 0 or iy <= 0:
                continue
            smaller = min((a[2] - a[0]) * (a[3] - a[1]), (b[2] - b[0]) * (b[3] - b[1]))
            if smaller > 0 and (ix * iy) / smaller >= frac:
                out.append((a[4], b[4], (ix * iy) / smaller))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frac", type=float, default=0.30,
                    help="fraction of the smaller box that must be covered to count")
    ap.add_argument("paths", nargs="*", help="defaults to every PDF under results/figures/")
    a = ap.parse_args()

    files = [Path(p) for p in a.paths] or sorted(FIGS.glob("*.pdf"))
    if not files:
        sys.exit(f"no figures found under {FIGS}; refusing to pass having checked nothing")

    bad = 0
    for f in files:
        hits = overlapping(words(f), a.frac)
        if hits:
            bad += 1
            log(f"  {f.name}: {len(hits)} overlapping pair(s)")
            for x, y, r in hits:
                log(f"      {x!r} over {y!r}, {r:.0%} of the smaller box")
    log(f"\n  {len(files)} figures checked, {bad} with overlapping text")
    if bad:
        sys.exit(f"{bad} figure(s) draw text on top of text")


if __name__ == "__main__":
    main()
