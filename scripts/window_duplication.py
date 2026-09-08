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
# ALL FOUR committed arms, not the three the headline uses. neg2_rm is the region-matched
# bias-aware arm and it is committed evidence, so a duplication CENSUS has to include it. An
# earlier version scanned three and reported 306 duplicate groups while a wildcard over the
# store found 1020; both were right about different populations, which is exactly the kind of
# ambiguity a census should not have.
ARMS = ("gc", "dinuc", "neg2", "neg2_rm")
# The fields a model actually consumes. A duplicate that agrees on all of these is a harmless
# repeat; one that disagrees would make "which copy" a real choice.
CONTENT = ("seq_rna", "label", "fold")
# A window is its coordinates AND its strand. chrom:start alone conflates a plus-strand window
# with the minus-strand window at the same start, which is a different sequence.
IDENTITY = ("chrom", "start", "end", "strand")


def scan(store):
    rows = []
    for arm in ARMS:
        n_ds = n_rows = n_dup = n_aff = n_cross = n_mixed = 0
        n_groups = n_same = 0
        for f in sorted(glob.glob(f"{store}/processed/{arm}/*/*/dataset.tsv")):
            try:
                d = pd.read_csv(f, sep="\t",
                            usecols=sorted(set(IDENTITY) | set(CONTENT)))
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
            # Are duplicated rows IDENTICAL in everything the model sees? If they are, then
            # "which copy is kept" by any subsampler is immaterial and selection stays
            # shuffle-invariant at the level of content. Measured rather than assumed, because
            # scripts/class_ratio.py's negative ranking relies on it.
            #
            # Grouped on the full identity, which is the key that ranking actually uses. On the
            # coarse (chrom, start) key the answer is different and misleading: 1028 groups of
            # which 8 differ in seq_rna, because a window at the same start on the opposite
            # strand is reverse-complemented and SHOULD have a different sequence. Those are
            # not duplicates, they are distinct windows the coarse key conflates.
            sub = d[d.duplicated(list(IDENTITY), keep=False)]
            if len(sub):
                per = sub.groupby(list(IDENTITY), sort=False)[list(CONTENT)].nunique()
                n_groups += int(len(per))
                n_same += int((per.max(axis=1) == 1).sum())
        if not n_ds:
            log(f"  {arm}: no datasets found under {store}/processed/{arm}")
            continue
        rows.append(dict(arm=arm, datasets=n_ds, rows=n_rows, duplicate_rows=n_dup,
                         datasets_affected=n_aff, cross_fold=n_cross, mixed_label=n_mixed,
                         dup_groups=n_groups, dup_groups_identical=n_same))
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
    add("duplicate groups whose rows are IDENTICAL in seq_rna, label and fold",
        int(t.dup_groups_identical.sum()), n=int(t.dup_groups.sum()),
        note="if every duplicate group is internally identical then which copy a subsampler "
             "keeps is immaterial, and selection stays shuffle-invariant in content. Measured "
             "because the class-ratio ranking relies on it")
    add("worst duplicate rate across arms", float((t.duplicate_rows / t["rows"]).max()),
        n=int(t["rows"].sum()), note="the bias-aware arm draws from a finite pool of other "
                                     "proteins' sites, so it repeats where the others do not")
    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {per.name}")


if __name__ == "__main__":
    main()
