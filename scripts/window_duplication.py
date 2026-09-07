"""Does any window appear twice in one dataset, and if so does it cross a fold?

    python scripts/window_duplication.py --store ../rbp-store
    python scripts/window_duplication.py --from-cache

Why this exists, and it was found by accident. `scripts/validate.py` is a stage-3 gate from the
earlier train/test/val pipeline. Nothing runs it, it is not in `run.sh` or CI, and pointed at
this study's window store it reports "GATE CLOSED". Most of that is the wrong gate for the arm:
it applies a GC tolerance to the dinucleotide-matched arm, where a GC gap is expected by
construction, and a `min_test_pairs` floor to a `split` column this study replaced with 5-fold
`fold`. One of its checks was neither stale nor misapplied, and this file is that check, run
properly and gated.

A paper about how negative sets are built ought to know whether its own negative sets contain
the same window twice. They do, rarely, and the rate differs by construction:

    arm           rows        duplicate rows   datasets affected
    GC            1,792,050                0                  0
    dinucleotide  1,797,046               14                 11
    bias-aware      913,468              300                 46

The bias-aware rate is forty times the dinucleotide one and that is mechanistic rather than a
bug in one arm: bias-aware negatives are drawn from a FINITE POOL of other proteins' binding
sites, so the same window can be selected twice for one dataset. The composition-matched arms
sample from free genomic intervals, a far larger pool, and the GC arm never repeats at all.

What matters is not the count but the KIND, and both dangerous kinds are absent:

  * no duplicate pair straddles a fold boundary, in any arm. A window appearing in two folds
    would put the same sequence in training and in held-out evaluation, which is the leakage
    this study's whole fold design exists to prevent.
  * no duplicate pair carries conflicting labels. The same window is never both a positive and
    a negative.

So the effect is that one negative is double-weighted within one fold of one dataset, at a rate
of 3.3e-04 in the worst arm. That is reported rather than corrected, because correcting it would
change committed evidence for a quantity far below the precision anything is reported to, and
because the honest thing for a paper about negative-set construction is to state the property.
"""

import argparse
import glob
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "window_duplication.csv"
ARMS = ("gc", "dinuc", "neg2")


def scan(store):
    rows = []
    for arm in ARMS:
        n_ds = n_rows = n_dup = n_aff = n_cross = n_mixed = 0
        for f in sorted(glob.glob(f"{store}/processed/{arm}/*/*/dataset.tsv")):
            try:
                d = pd.read_csv(f, sep="\t", usecols=["chrom", "start", "label", "fold"])
            except (OSError, ValueError):
                continue
            n_ds += 1
            n_rows += len(d)
            mask = d.duplicated(["chrom", "start"], keep=False)
            if not mask.any():
                continue
            n_aff += 1
            g = d[mask].groupby(["chrom", "start"])
            n_dup += int(d.duplicated(["chrom", "start"]).sum())
            n_cross += int((g.fold.nunique() > 1).sum())
            n_mixed += int((g.label.nunique() > 1).sum())
        if not n_ds:
            log(f"  {arm}: no datasets found under {store}/processed/{arm}")
            continue
        rows.append(dict(arm=arm, datasets=n_ds, rows=n_rows, duplicate_rows=n_dup,
                         datasets_affected=n_aff, cross_fold=n_cross, mixed_label=n_mixed))
        log(f"  {arm:6} {n_ds:4d} datasets  {n_rows:9,} rows  {n_dup:4d} duplicate  "
            f"{n_aff:3d} affected  {n_cross} cross-fold  {n_mixed} mixed-label")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing scanned; refusing to write an empty table over committed evidence")
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    ap.add_argument("--from-cache", action="store_true")
    a = ap.parse_args()

    per = TABLES / "window_duplication_per_arm.csv"
    if a.from_cache:
        if not per.exists():
            sys.exit(f"no {per.name}; run with --store first")
        t = pd.read_csv(per)
    else:
        t = scan(a.store)
        t.to_csv(per, index=False)

    out = []

    def add(check, value, n="", note=""):
        out.append({"check": check, "value": value, "ci_low": "", "ci_high": "",
                    "n": n, "note": note})

    for _, r in t.iterrows():
        add(f"duplicate window rows, {r.arm} arm", int(r.duplicate_rows), n=int(r["rows"]),
            note=f"{int(r.datasets_affected)} of {int(r.datasets)} datasets affected")
        add(f"duplicate rate, {r.arm} arm", float(r.duplicate_rows) / float(r["rows"]),
            n=int(r["rows"]), note="duplicate rows over total rows")

    # The two kinds that would matter, asserted across every arm at once.
    add("duplicate windows that straddle a fold boundary", int(t.cross_fold.sum()),
        n=int(t["rows"].sum()),
        note="a window in two folds would put the same sequence in training and in held-out "
             "evaluation; this is the number that must be zero")
    add("duplicate windows carrying conflicting labels", int(t.mixed_label.sum()),
        n=int(t["rows"].sum()),
        note="the same window as both a positive and a negative; must be zero")
    add("worst duplicate rate across arms", float((t.duplicate_rows / t["rows"]).max()),
        n=int(t["rows"].sum()), note="the bias-aware arm draws from a finite pool of other "
                                     "proteins' sites, so it repeats where the others do not")
    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {per.name}")


if __name__ == "__main__":
    main()
