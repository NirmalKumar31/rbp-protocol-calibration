"""Stage 1b: count usable windows per protein across the full ENCODE panel.

Answers the one question the whole design rests on: how many of the 139 K562 eCLIP
RBPs clear the pre-registered inclusion filter? Runs the deterministic part of
preprocessing -- chromosome filter, window centring, region classification, dedup --
which needs only the peak BEDs and data/interim/regions.pkl. The genome is not
touched, so this is seconds per protein rather than minutes.

The genome-dependent drops (out-of-bounds, wrong length, ambiguous base) and the
negative-matching drops are excluded here. On the 19-protein development panel those
cost 0-8 windows out of thousands, and the calibration block below reports that gap
against the real prepare.py numbers rather than assuming it.

    python scripts/panel_census.py --panel config/panel_full.tsv
"""

import argparse
import pickle
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rbp.data import annotation as ann  # noqa: E402
from rbp.data import encode  # noqa: E402
from rbp.data import windows as win  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def census_one(protein, cfg, index):
    """Per-chromosome usable window counts for one protein."""
    size = cfg.windows["size"]
    keep = set(ann.MAIN_CHROMS) - set(cfg.encode.get("exclude_chroms", []))
    pth = encode.peak_path(ROOT, protein, cfg.encode["cell_line"])

    per_chrom = Counter()
    n_raw = n_offchrom = n_noregion = n_dup = 0
    seen = set()
    for chrom, start, end, _ in win.read_peaks(pth):
        n_raw += 1
        if chrom not in keep:
            n_offchrom += 1
            continue
        w0, w1 = win.window_bounds(start, end, size)
        if ann.classify(index, chrom, w0, w1) is None:
            n_noregion += 1
            continue
        if (chrom, w0) in seen:
            n_dup += 1
            continue
        seen.add((chrom, w0))
        per_chrom[chrom] += 1

    return {"protein": protein, "peaks_raw": n_raw, "windows": sum(per_chrom.values()),
            "off_chrom": n_offchrom, "no_region": n_noregion, "duplicate": n_dup,
            "per_chrom": dict(per_chrom)}


def calibrate(rows, cfg):
    """Compare projected test windows against real prepare.py output where it exists.

    The projection is only trustworthy if this gap is small and consistent; printing
    it is what makes the >=100 claim checkable instead of asserted.

    Keyed by cell line. Reading data/processed/<PROTEIN>/ without one compared HepG2
    counts against K562 prep reports and reported a +154 mean error that meant nothing.
    """
    test = set(cfg.split["test"])
    cell = cfg.encode["cell_line"]
    out = []
    for r in rows:
        rep = ROOT / "data/processed" / cell / r["protein"] / "prep_report.txt"
        if not rep.exists():
            continue
        actual = None
        for ln in rep.read_text().splitlines():
            if ln.startswith("test_pairs:"):
                actual = int(ln.split(":")[1])
        if actual is None:
            continue
        proj = sum(n for c, n in r["per_chrom"].items() if c in test)
        out.append({"protein": r["protein"], "projected": proj, "actual": actual,
                    "diff": proj - actual})
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--panel", default="config/panel_full.tsv")
    p.add_argument("--cell-line", default=None, help="override encode.cell_line")
    p.add_argument("--index", default=str(ROOT / "data/interim/regions.pkl"))
    p.add_argument("--out", default="config/panel_census.tsv")
    a = p.parse_args()

    cfg = cfgmod.load(a.config)
    if a.cell_line:
        cfg["encode"]["cell_line"] = a.cell_line
    names = [r["protein"] for r in cfgmod.proteins(panel=a.panel)]
    test = cfg.split["test"]
    min_test = cfg.panel["min_test_pairs"]

    print("loading region index ...", flush=True)
    index = pickle.loads(Path(a.index).read_bytes())
    print(f"census over {len(names)} proteins; test chroms {test}; "
          f"threshold {min_test} test pairs\n")

    rows = []
    for i, name in enumerate(names, 1):
        r = census_one(name, cfg, index)
        r["test_windows"] = sum(n for c, n in r["per_chrom"].items() if c in test)
        rows.append(r)
        print(f"  [{i:3d}/{len(names)}] {name:10} {r['windows']:7,} windows  "
              f"{r['test_windows']:6,} test  "
              f"{'PASS' if r['test_windows'] >= min_test else 'fail'}", flush=True)

    rows.sort(key=lambda r: -r["test_windows"])
    passing = [r for r in rows if r["test_windows"] >= min_test]

    chroms = [c for c in ann.MAIN_CHROMS if c not in set(cfg.encode.get("exclude_chroms", []))]
    out = ROOT / a.out
    cols = ["protein", "peaks_raw", "windows", "test_windows", "off_chrom",
            "no_region", "duplicate"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols + chroms) + "\n")
        for r in rows:
            fh.write("\t".join([str(r[c]) for c in cols]
                               + [str(r["per_chrom"].get(c, 0)) for c in chroms]) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")

    cal = calibrate(rows, cfg)
    if cal:
        worst = max(cal, key=lambda c: abs(c["diff"]))
        print(f"\ncalibration against real prepare.py output ({len(cal)} proteins):")
        print(f"  mean projected-minus-actual  {sum(c['diff'] for c in cal)/len(cal):+.1f}")
        print(f"  worst case                   {worst['protein']} {worst['diff']:+d} "
              f"({worst['projected']} vs {worst['actual']})")

    print(f"\n{len(passing)}/{len(rows)} proteins clear {min_test} test windows")
    if passing:
        med = passing[len(passing) // 2]["test_windows"]
        print(f"  median test windows among passing: {med:,}")
        print(f"  total windows among passing: {sum(r['windows'] for r in passing):,}")


if __name__ == "__main__":
    main()
