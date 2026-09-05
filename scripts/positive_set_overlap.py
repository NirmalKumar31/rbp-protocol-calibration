"""How much do the two composition-matched arms' POSITIVE sets differ?

    python scripts/positive_set_overlap.py --store ../rbp-store
    python scripts/positive_set_overlap.py --from-cache

WHY IT MATTERS. The design claim is that only the negatives change: the model, the peak set, the
folds and the estimator are held fixed. That is very nearly true and not exactly true. A
positive is retained only when its matcher finds an acceptable negative, and the GC and
dinucleotide matchers fail on different windows, so the two arms end up with positive sets that
overlap heavily rather than coincide. This measures the gap, and the Limitations quote it.

WHY THIS FILE EXISTS AT ALL, WHICH IS THE UNCOMFORTABLE PART. `positive_set_overlap.csv` was
committed, is cited as Supplementary Table S8, and the Discussion quotes three numbers out of it
-- median Jaccard 0.9972, minimum 0.9237, identical in 10 of 94 datasets -- and NO SCRIPT IN THE
REPOSITORY PRODUCED IT. It was computed once, by hand, and the working was lost.

That is the same failure as the +0.0397 that motivated scripts/audit_manuscript.py: a number in
the paper that could not be reproduced and could not fail. The audit did not catch this one
because it closes the manuscript-to-table link and nothing closed the table-to-script link.
scripts/provenance.py closes it now, and found this along with two unquoted leftovers.

The values below reproduce the committed table exactly, which is the check that this
reconstruction is the original computation rather than a plausible substitute.
"""

import argparse
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "gc", "dn": "dinuc"}


KEY = ["chrom", "start", "end", "strand"]


def positives(path):
    """The positive windows of one arm, as a set of (chrom, start, end, strand).

    Keyed on genomic position and not on row order: the two arms retain different subsets, so
    the files have different lengths and a positional comparison would be meaningless. That is
    the bug the original implementation had -- it keyed on `id`, a per-arm row index -- and the
    resulting column was a count ratio rather than a set overlap.

    STRAND IS PART OF THE POSITION. Adding it changes nothing here, and that is measured rather
    than assumed: across all 914,732 positives in both arms there are zero pairs of positives
    sharing a coordinate on opposite strands. The assertion below keeps it that way.
    """
    d = pd.read_csv(path, sep="\t")
    p = d[d.label == 1]
    for c in KEY:
        if c not in p.columns:
            sys.exit(f"{path} has no {c} column; cannot identify a positive window by position")
    without_strand = p.groupby(KEY[:3]).ngroups
    with_strand = p.groupby(KEY).ngroups
    if without_strand != with_strand:
        sys.exit(f"{path}: {with_strand - without_strand} positives share a coordinate on "
                 f"opposite strands. Every quoted overlap figure assumed this could not happen; "
                 f"recheck them before trusting this run.")
    return set(zip(p.chrom, p.start, p.end, p.strand))


def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    rows = []
    for n, ds in enumerate(list(pub.dataset)[:limit or None], 1):
        protein, cell = ds.split(":")
        sets = {}
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                sets = None
                break
            sets[arm] = positives(f)
        if not sets:
            continue
        a, b = sets["gc"], sets["dn"]
        rows.append({"dataset": ds, "protein": protein, "cell": cell,
                     "n_pos_gc": len(a), "n_pos_dn": len(b),
                     "jaccard": len(a & b) / len(a | b)})
        log(f"[{n:3d}] {ds:18s} gc {len(a):6d}  dn {len(b):6d}  J {rows[-1]['jaccard']:.4f}")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    out = TABLES / "positive_set_overlap.csv"
    t = pd.read_csv(out) if a.from_cache else build(a.store, a.n)
    if not a.from_cache:
        if a.n:
            # A PARTIAL RUN MUST NOT OVERWRITE THE COMMITTED TABLE. `--n 5` is a smoke test,
            # and letting it write truncated the released 94-row table to 5 rows in place. The
            # same shape of bug -- a sampling flag that silently replaces a full artefact --
            # has now bitten this repository twice.
            out = out.with_suffix(".partial.csv")
            log(f"  --n {a.n} given: writing {out.name}, not the committed table")
        t.to_csv(out, index=False)

    log("")
    log(f"  datasets                              {len(t)}")
    log(f"  median Jaccard of the positive sets   {t.jaccard.median():.4f}")
    log(f"  minimum                               {t.jaccard.min():.4f}")
    log(f"  identical in                          {int((t.jaccard == 1).sum())} of {len(t)}")
    log("")
    log("  'only the negatives change' is therefore very nearly true and not exactly true.")


if __name__ == "__main__":
    main()
