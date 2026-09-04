"""B9: two choices in how the AUROC is computed, not in what is modelled.

    python scripts/auroc_aggregation.py --store ../rbp-store

TWO CHOICES, BOTH MADE ONCE AND NEVER REVISITED.

1. POOLED VERSUS FOLD-AVERAGED. Every AUROC here pools all five folds' out-of-fold scores into
   one ranking and computes a single AUROC. The alternative computes an AUROC per fold and
   averages. They differ, and not only by noise: pooling compares scores ACROSS folds, so a
   fold whose predictor is shifted relative to the others contributes cross-fold comparisons
   that no within-fold analysis would make. Fold-averaging weights a small fold as heavily as a
   large one, which is why pooling was chosen; the cost of that choice is measured here.

2. RAW VERSUS RANK-NORMALISED WITHIN FOLD. If the five fold predictors are on slightly
   different scales, pooling them is comparing incommensurable numbers. Rank-normalising each
   fold's scores to [0, 1] before pooling removes any monotone per-fold shift by construction.
   If the pooled AUROC is an artefact of scale drift between folds, this is where it shows.

Both are computed on the nested contribution, the paper's actual estimand, not on the model
AUROC alone -- an increment can move even where neither of its two terms does much.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from deep_model_contrast import MIN_COVERAGE, MODELS, arm_roots, oof  # noqa: E402
from nested_scale import panel  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import _oof_scores, composition_features  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"


def log(m):
    print(m, flush=True)


def rank_within(s, folds):
    """Rank-normalise to [0, 1] inside each fold, so no per-fold shift survives pooling."""
    out = np.full(len(s), np.nan)
    for f in np.unique(folds):
        m = folds == f
        v = s[m]
        r = np.argsort(np.argsort(v)).astype(float)
        out[m] = r / (len(r) - 1) if len(r) > 1 else 0.5
    return out


def gains(y, folds, s_comp, s_full):
    """The nested increment under three aggregations of the same two predictors."""
    ok = np.isfinite(s_comp) & np.isfinite(s_full)
    y, folds = y[ok], folds[ok]
    c, f = s_comp[ok], s_full[ok]
    pooled = roc_auc_score(y, f) - roc_auc_score(y, c)

    per = []
    for fold in np.unique(folds):
        m = folds == fold
        if len(np.unique(y[m])) < 2:
            continue
        per.append(roc_auc_score(y[m], f[m]) - roc_auc_score(y[m], c[m]))
    averaged = float(np.mean(per)) if per else np.nan

    rc, rf = rank_within(c, folds), rank_within(f, folds)
    ranked = roc_auc_score(y, rf) - roc_auc_score(y, rc)
    return float(pooled), averaged, float(ranked), len(per)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--only", default="")
    a = p.parse_args()
    datasets = [x.strip() for x in a.only.split(",") if x.strip()] or panel(Path(a.store))
    roots = arm_roots(Path(a.store))
    log(f"=== B9: AUROC aggregation, {len(datasets)} datasets ===\n")

    rows = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        for arm, (dataroot, scoreroot) in roots.items():
            f = dataroot / cell / protein / "dataset.tsv"
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="\t")
            ids, got, bad = set(d.id), {}, False
            for model in MODELS:
                if model == "kmer":
                    continue
                s = oof(scoreroot, cell, protein, model)
                if s is None:
                    bad = True
                    break
                got[model] = s
                ids &= set(s.id)
            if bad or len(ids) / len(d) < MIN_COVERAGE:
                continue
            dd = d[d.id.isin(ids)].reset_index(drop=True)
            y, folds = dd.label.values, dd.fold.values
            sc, _, _ = kmer_oof(dd.seq_rna.values, y, folds, k=4)
            scores = {"kmer": sc}
            for model, s in got.items():
                scores[model] = dd[["id"]].merge(s, on="id", how="left").score.to_numpy()
            comp, _ = composition_features(dd.seq_rna.values, True)
            s_comp = _oof_scores(comp, y, folds, "l2")
            for model, raw in scores.items():
                s_full = _oof_scores(np.column_stack([comp, standardise(raw)]),
                                     y, folds, "l2")
                pooled, avg, ranked, nf = gains(y, folds, s_comp, s_full)
                rows.append({"dataset": ds, "arm": arm, "model": model, "n_folds": nf,
                             "pooled": pooled, "fold_averaged": avg, "rank_pooled": ranked})
        log(f"[{i:3d}/{len(datasets)}] {ds}")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no rows produced")
    t.to_csv(TABLES / "auroc_aggregation_per_dataset.csv", index=False)

    out = []
    log(f"\n  {'arm':6} {'model':11} {'pooled':>9} {'fold-avg':>9} {'rank-pooled':>12}")
    for arm in ("dn", "gc", "neg2"):
        for model in MODELS:
            g = t[(t.arm == arm) & (t.model == model)]
            if not len(g):
                continue
            p_, a_, r_ = g.pooled.mean(), g.fold_averaged.mean(), g.rank_pooled.mean()
            for name, v in (("pooled", p_), ("fold-averaged", a_), ("rank-pooled", r_)):
                out.append({"check": f"gain, {arm} arm, {model}, {name}",
                            "value": float(v), "n": len(g)})
            out += [{"check": f"fold-averaged minus pooled, {arm} arm, {model}",
                     "value": float(a_ - p_), "n": len(g),
                     "note": f"max |delta| {(g.fold_averaged - g.pooled).abs().max():.4f}"},
                    {"check": f"rank-pooled minus pooled, {arm} arm, {model}",
                     "value": float(r_ - p_), "n": len(g),
                     "note": f"max |delta| {(g.rank_pooled - g.pooled).abs().max():.4f}"}]
            log(f"  {arm:6} {model:11} {p_:9.4f} {a_:9.4f} {r_:12.4f}")

    # THE ORDERING IS THE CLAIM. It must survive both alternatives, or the paper's headline is a
    # property of one way of aggregating folds.
    held = 0
    for name in ("pooled", "fold_averaged", "rank_pooled"):
        ok = True
        for model in MODELS:
            v = [t[(t.arm == arm) & (t.model == model)][name].mean()
                 for arm in ("dn", "gc", "neg2")]
            if not (v[0] > v[1] > v[2]):
                ok = False
        held += int(ok)
        out.append({"check": f"protocol ordering holds for every model, {name}",
                    "value": int(ok), "n": len(MODELS)})
    out.append({"check": "aggregations on which the ordering holds for every model",
                "value": held, "n": 3})
    pd.DataFrame(out).to_csv(TABLES / "auroc_aggregation.csv", index=False)
    log(f"\n  ordering holds under {held}/3 aggregations")
    log("\nwrote auroc_aggregation.csv and auroc_aggregation_per_dataset.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
