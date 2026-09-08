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

How the ratio is reached, and it can only be one way. The committed windows are already about
1:1, so no ratio is reachable by adding rows. Rows are subsampled WITHIN FOLD, so the fold
partition and its blocking are untouched: for 1:2 and 1:4 every negative is kept and the
positives are cut, for 2:1 every positive is kept and the negatives are cut. That means the
smaller ratios also shrink the dataset, which is a confound this design cannot separate from
balance itself, and it is reported rather than hidden: the row count at each ratio is a column.

No interval is computed for the subsampled ratios. The subsample introduces a variance
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


def _rank(ids, dataset, label, salt):
    """A deterministic order over rows, keyed on the WINDOW rather than on its row index.

    This is what makes the positive subsample arm-invariant. Ranking by a hash of the window's
    own identity means a positive present in two arms gets the same rank in both, so the
    retained positives coincide wherever the positive sets do.
    """
    return np.argsort([hashlib.sha256(
        f"{SEED}|{dataset}|{label}|{salt}|{i}".encode()).hexdigest() for i in ids])


def subsample(ids, y, folds, ratio, dataset, arm, label):
    """Row mask giving `ratio` positives per negative, within every fold.

    The arm confound this fixes. The first version drew with an RNG seeded on
    (dataset, arm, ratio) and selected by row index, so each arm kept a DIFFERENT random subset
    of positives. For the ratios below one, which are reached by cutting positives, the arm
    comparison then varied both the negative protocol and the positive sample realisation, when
    the whole point is to vary only the first. An audit was right to call that out.

    Positives are now ranked by a hash of the window's own (chrom, start) identity, with the
    arm deliberately ABSENT from the key, so any positive shared between two arms is kept or
    dropped in both. Negatives are ranked with the arm in the key, because they are different
    windows by construction and there is nothing to hold fixed.
    """
    y, folds = np.asarray(y, dtype=int), np.asarray(folds)
    ids = np.asarray(ids)
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
        keep[pos[_rank(ids[pos], dataset, label, "pos")[:want_pos]]] = True
        keep[neg[_rank(ids[neg], dataset, label, f"neg|{arm}")[:want_neg]]] = True
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
            ids = (d.chrom.astype(str) + ":" + d.start.astype(str)).values
            for lab, ratio in RATIOS:
                m = subsample(ids, y, folds, ratio, ds, arm, lab)
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
        cols = [f"gain_{arm}_{lab}" for arm in ARMS]
        if not all(c in t.columns for c in cols):
            add(f"span, {lab}", np.nan, note="an arm produced no column at this ratio")
            log(f"  {lab}: an arm did not survive")
            continue
        # The common intersection, and the first version did not take it. Panel means were
        # computed over whatever each arm retained, then the span and the ordering were
        # reported with n = 94 whichever datasets had actually survived. At 1:4 five datasets
        # drop, so those rows claimed a panel they were not computed on. Every arm mean, the
        # span and the ordering now come from the SAME datasets, and n is that number.
        common = t.dropna(subset=cols)
        means = {arm: float(common[f"gain_{arm}_{lab}"].mean()) for arm in ARMS}
        for arm in ARMS:
            add(f"panel-mean contribution, {arm} arm, {lab}", means[arm], n=len(common),
                note="over the datasets retained in ALL THREE arms at this ratio")
        add(f"datasets in all three arms, {lab}", len(common), n=len(t),
            note=f"of {len(t)}; the denominator every row at this ratio is computed on")
        span = max(means.values()) / min(means.values()) if min(means.values()) > 0 else np.nan
        add(f"span, {lab}", span, n=len(common),
            note="max over min of the three panel means, on the common intersection")
        holds = means["dn"] > means["gc"] > means["neg2"]
        add(f"ordering dn > gc > neg2 holds, {lab}", float(holds), n=len(common),
            note="the claim under test; reported whichever way it comes out")
        if lab == "1:1":
            published = span
        log(f"  {lab}: n={len(common):3d}  dn {means['dn']:+.4f}  gc {means['gc']:+.4f}  "
            f"neg2 {means['neg2']:+.4f}   span {span:.3f}   "
            f"ordering {'HOLDS' if holds else 'BREAKS'}")

    if published is not None:
        add("span at the published 1:1 balance", published, n=len(t),
            note="the control: this is the arm ordering and span the paper reports")
    # --from-cache cannot recompute an exclusion count: the excluded cells are absent from the
    # per-dataset table it reads. Carried, so the documented command does not degrade a
    # committed table. See rbp.utils.carry.
    from rbp.utils.carry import emit as carry_emit
    for lab, _ in RATIOS:
        if excluded is not None:
            add(f"dataset-arms excluded, {lab}", excluded[lab], n=len(t) * len(ARMS),
                note=f"of {len(t) * len(ARMS)} dataset-arm cells at this ratio "
                     f"({len(t)} datasets x {len(ARMS)} arms); a fold would hold under "
                     f"{MIN_PER_FOLD} rows of one class")
        else:
            carry_emit(out, OUT, f"dataset-arms excluded, {lab}", None, recomputed=False)

    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {PER.name}")


if __name__ == "__main__":
    main()
