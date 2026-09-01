"""R1o: how much of "what the model contributes" is just one order of composition up?

    python scripts/baseline_order.py --n 30

THE OBJECTION, AND IT IS THE SHARPEST ONE ANYONE MADE ABOUT WHAT THE QUANTITY IS. A bag of
4-mer counts carries NO positional information. The composition baseline stops at order 2. So
the headline "nested contribution of a sequence model over composition" is, definitionally,
order-3-and-4 composition beyond order-1-and-2 composition -- and where the baseline stops is
an arbitrary choice doing all the work. `nested.py` defends stopping at order 2 as "exactly
what a dinucleotide shuffle would have preserved", which is a good reason to choose it and not
a reason to believe the choice does not matter.

Nobody had measured it. This does: refit the baseline at order 3 (mono + di + tri, 84
features) and recompute the same nested contribution.

WHAT IT COSTS THE PAPER. Roughly half of every arm's contribution and most of the R1 contrast.
The word "sequence model" has to be qualified: what is measured is short-range compositional
signal one order above wherever the matcher stopped, not motif recognition and not anything
positional.

WHAT IT DOES NOT COST. The FOLD RANGE, which is the paper's actual claim. Everything shrinks
roughly proportionally, so the protocol dependence is not an artefact of where the baseline
stops. That is the result worth having: the magnitude is a function of an analyst's choice,
the protocol dependence is not.

Subsampled by default because the order-3 fit is 84 features x 5 folds x 3 arms per dataset.
The subsample is size-stratified and its own R1 contrast is reported next to the published one
so the reader can see whether the subsample is representative.
"""

import argparse
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "processed/gc", "dn": "processed/dinuc", "neg2": "processed/neg2"}


def log(m):
    print(m, flush=True)


def composition(seqs, max_k):
    """Mono, di, ... up to max_k, as frequencies. Order 2 is the paper's baseline."""
    blocks = []
    for k in range(1, max_k + 1):
        vocab = {"".join(p): i for i, p in enumerate(product("ACGU", repeat=k))}
        M = np.zeros((len(seqs), len(vocab)))
        for r, s in enumerate(seqs):
            for i in range(len(s) - k + 1):
                j = vocab.get(s[i:i + k])
                if j is not None:
                    M[r, j] += 1.0
            M[r] /= max(len(s) - k + 1, 1)
        blocks.append(M)
    return np.hstack(blocks)


def oof_auc(X, y, folds):
    out = np.full(len(y), np.nan)
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        out[te] = LogisticRegression(max_iter=2000).fit(X[tr], y[tr]).decision_function(X[te])
    return float(roc_auc_score(y, out))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=30, help="size-stratified subsample")
    a = p.parse_args()
    warnings.filterwarnings("ignore")
    store = Path(a.store)

    panel = pd.read_csv(TABLES / "rehearsal_binding_gc.csv").sort_values("pairs")
    step = max(len(panel) // a.n, 1)
    sub = panel.iloc[::step].head(a.n)

    rows = []
    for i, r in enumerate(sub.itertuples(), 1):
        rec, ok = {"dataset": r.dataset, "protein": r.protein, "cell": r.cell}, True
        for arm, sd in ARMS.items():
            f = store / sd / r.cell / r.protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
            m = ~np.isnan(sc)
            seqs, y, fo = d.seq_rna.values[m], d.label.values[m], d.fold.values[m]
            z = standardise(sc[m]).reshape(-1, 1)
            for order in (2, 3):
                X = composition(seqs, order)
                rec[f"gain{order}_{arm}"] = (oof_auc(np.hstack([X, z]), y, fo)
                                             - oof_auc(X, y, fo))
        if ok:
            rows.append(rec)
            log(f"[{i:3d}/{len(sub)}] {r.dataset:18s} "
                + "  ".join(f"{k} {rec[f'gain2_{k}']:+.4f}->{rec[f'gain3_{k}']:+.4f}"
                            for k in ARMS))
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "baseline_order_per_dataset.csv", index=False)

    out = []
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(t), size=(2000, len(t)))

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in idx])
        lo, hi = np.percentile(b, [2.5, 97.5])
        out.append({"check": check, "value": float(v.mean()), "ci_low": float(lo),
                    "ci_high": float(hi), "n": len(t), "note": note})
        return float(v.mean())

    log(f"\n=== R1o: baseline order, n = {len(t)} size-stratified datasets ===\n")
    log(f"  {'arm':6s} {'over mono+di':>14s} {'over mono+di+tri':>18s} {'removed':>9s}")
    for arm in ARMS:
        g2 = add(f"gain over order-2 baseline, {arm} arm", t[f"gain2_{arm}"])
        g3 = add(f"gain over order-3 baseline, {arm} arm", t[f"gain3_{arm}"])
        frac = 1 - g3 / g2 if g2 else np.nan
        out.append({"check": f"fraction removed by order 3, {arm} arm", "value": float(frac),
                    "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
        log(f"  {arm:6s} {g2:+14.4f} {g3:+18.4f} {100 * frac:8.0f}%")

    c2 = add("R1 contrast (dn-gc), order-2 baseline", t.gain2_dn - t.gain2_gc)
    c3 = add("R1 contrast (dn-gc), order-3 baseline", t.gain3_dn - t.gain3_gc)
    out.append({"check": "fraction of the R1 contrast removed by order 3",
                "value": float(1 - c3 / c2), "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
    log(f"\n  R1 contrast  order-2 {c2:+.4f}   order-3 {c3:+.4f}   "
        f"{100 * (1 - c3 / c2):.0f}% removed")

    # THE POINT: the range is what the paper claims, and it must survive.
    for order in (2, 3):
        means = [t[f"gain{order}_{arm}"].mean() for arm in ARMS]
        fr = max(means) / min(means) if min(means) > 0 else np.nan
        out.append({"check": f"fold range across protocols, order-{order} baseline",
                    "value": float(fr), "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
        log(f"  fold range across protocols, order-{order} baseline: {fr:.2f}x")

    pd.DataFrame(out).to_csv(TABLES / "baseline_order.csv", index=False)
    log("\nwrote baseline_order.csv and baseline_order_per_dataset.csv")


if __name__ == "__main__":
    main()
