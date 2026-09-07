"""Does the arm ordering survive a change in class balance?

    python scripts/class_ratio.py --store ../rbp-store
    python scripts/class_ratio.py --from-cache

Specified in `docs/SENSITIVITY_SPEC.md` section B, committed before this was run.

WHY IT MATTERS. Every nested fit in this paper is estimated at roughly 1:1 positives to
negatives, and a genome-wide screen faces a far smaller positive fraction. AUROC is
prevalence-invariant, so the naive expectation is that nothing moves. The fitted logistic is
not: the intercept shifts with balance, and at fixed penalty C the effective shrinkage on each
coefficient shifts with the number of rows. The composition block and the model score do not
have to absorb that shift equally, and the whole quantity here is a DIFFERENCE between two fits.
So the question is whether the arm ordering and the span survive, not whether AUROC does.

HOW THE RATIO IS REACHED, and it can only be one way. The committed windows are already about
1:1, so no ratio is reachable by adding rows. Rows are subsampled WITHIN FOLD, so the fold
partition and its blocking are untouched: for 1:2 and 1:4 every negative is kept and the
positives are cut, for 2:1 every positive is kept and the negatives are cut. That means the
smaller ratios also shrink the dataset, which is a confound this design cannot separate from
balance itself, and it is reported rather than hidden: the row count at each ratio is a column.

NO INTERVAL IS COMPUTED FOR THE SUBSAMPLED RATIOS. The subsample introduces a variance
component that a protein-clustered bootstrap over one draw does not represent, and reporting a
protein-clustered interval as though it covered both would be the same error this project made
once already with SD versus SD/sqrt(n). Point estimates and the ordering are what the analysis
is for.
"""

import argparse
import hashlib
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "class_ratio.csv"
PER = TABLES / "class_ratio_per_dataset.csv"

ARMS = {"dn": "dinuc", "gc": "gc", "neg2": "neg2"}
# (label, positives per negative). 1.0 is the published balance and is the control.
RATIOS = (("1:1", 1.0), ("1:2", 0.5), ("1:4", 0.25), ("2:1", 2.0))
SEED = 20260907
MIN_PER_FOLD = 20            # rows of EITHER class below which a fold is not analysable
K = 4


def subsample(y, folds, ratio, dataset, arm, label):
    """Row mask giving `ratio` positives per negative, within every fold.

    The generator is seeded from the constant plus a stable hash of (dataset, arm, ratio), so
    the draw does not depend on iteration order and rerunning reproduces it exactly. Python's
    own hash() is salted per process and would not.
    """
    key = f"{dataset}|{arm}|{label}".encode()
    rng = np.random.default_rng(SEED + int(hashlib.sha256(key).hexdigest()[:8], 16))
    y, folds = np.asarray(y, dtype=int), np.asarray(folds)
    keep = np.zeros(len(y), dtype=bool)
    for f in np.unique(folds):
        inf = folds == f
        pos = np.flatnonzero(inf & (y == 1))
        neg = np.flatnonzero(inf & (y == 0))
        if ratio >= 1.0:                        # cut negatives, keep every positive
            want_pos, want_neg = len(pos), int(round(len(pos) / ratio))
        else:                                   # cut positives, keep every negative
            want_pos, want_neg = int(round(len(neg) * ratio)), len(neg)
        want_pos, want_neg = min(want_pos, len(pos)), min(want_neg, len(neg))
        if want_pos < MIN_PER_FOLD or want_neg < MIN_PER_FOLD:
            return None
        keep[rng.choice(pos, want_pos, replace=False)] = True
        keep[rng.choice(neg, want_neg, replace=False)] = True
    return keep


def build(store, limit=0):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    datasets = list(pub.dataset)[:limit or None]
    rows, excluded = [], {lab: 0 for lab, _ in RATIOS}
    for n, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        ok = True
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            seqs, y, folds = d.seq_rna.values, d.label.values, d.fold.values
            for lab, ratio in RATIOS:
                m = subsample(y, folds, ratio, ds, arm, lab)
                if m is None:
                    excluded[lab] += 1
                    continue
                # The published estimator, unchanged: out-of-fold k-mer score, then the
                # nested comparison against the composition block on the same rows.
                sc, _, _ = kmer_oof(seqs[m], y[m], folds[m], k=K)
                good = ~np.isnan(sc)
                g = gain_over_composition(seqs[m][good], sc[good], y[m][good],
                                          folds[m][good])
                rec[f"gain_{arm}_{lab}"] = g.delta
                rec[f"comp_{arm}_{lab}"] = g.auroc_composition
                rec[f"rows_{arm}_{lab}"] = int(good.sum())
        if not ok:
            continue
        rows.append(rec)
        log(f"  [{n:3d}/{len(datasets)}] {ds:18s} " + "  ".join(
            f"{lab} dn {rec.get(f'gain_dn_{lab}', float('nan')):+.4f}" for lab, _ in RATIOS))
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to write an empty table over committed evidence")
    return t, excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    ap.add_argument("--n", type=int, default=0, help="smoke test; writes .partial.csv")
    ap.add_argument("--from-cache", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")

    excluded = None
    if a.from_cache:
        if not PER.exists():
            sys.exit(f"no {PER.name}; run with --store first")
        t = pd.read_csv(PER)
    else:
        t, excluded = build(a.store, a.n)
        if a.n:
            p = PER.with_suffix(".partial.csv")
            t.to_csv(p, index=False)
            log(f"\n  --n {a.n}: wrote {p.name}, not the committed table")
            return
        t.to_csv(PER, index=False)

    out = []

    def add(check, value, n="", note=""):
        out.append({"check": check, "value": value, "ci_low": "", "ci_high": "",
                    "n": n, "note": note})

    add("datasets", len(t), n=len(t), note="the study panel, all three arms present")
    log("")
    published = None
    for lab, _ in RATIOS:
        means = {}
        for arm in ARMS:
            col = f"gain_{arm}_{lab}"
            if col not in t.columns:
                continue
            v = t[col].dropna()
            if len(v):
                means[arm] = float(v.mean())
                add(f"panel-mean contribution, {arm} arm, {lab}", means[arm], n=len(v))
        if len(means) < 3:
            add(f"span, {lab}", np.nan, note="an arm did not survive this ratio")
            log(f"  {lab}: an arm did not survive")
            continue
        span = max(means.values()) / min(means.values()) if min(means.values()) > 0 else np.nan
        add(f"span, {lab}", span, n=len(t),
            note="max over min of the three panel means, as in the paper")
        holds = means["dn"] > means["gc"] > means["neg2"]
        add(f"ordering dn > gc > neg2 holds, {lab}", float(holds), n=len(t),
            note="the claim under test; reported whichever way it comes out")
        if lab == "1:1":
            published = span
        log(f"  {lab}: dn {means['dn']:+.4f}  gc {means['gc']:+.4f}  neg2 {means['neg2']:+.4f}"
            f"   span {span:.3f}   ordering {'HOLDS' if holds else 'BREAKS'}")

    if published is not None:
        add("span at the published 1:1 balance", published, n=len(t),
            note="the control: this is the arm ordering and span the paper reports")
    if excluded is not None:
        for lab, v in excluded.items():
            add(f"dataset-arms excluded, {lab}", v, n=len(t),
                note=f"a fold would hold under {MIN_PER_FOLD} rows of one class")

    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {PER.name}")


if __name__ == "__main__":
    main()
