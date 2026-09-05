"""How much of the nested contribution is the outer-fold channel the Methods disclose?

    python scripts/cross_fitting.py --store ../rbp-store
    python scripts/cross_fitting.py --from-cache

THE CHANNEL. Every published contribution is produced in two stages. A base model is fitted
under the study's folds and scores each window out of fold; that score then enters a second
cross-validation as a covariate beside the composition block. For outer test fold $i$ the
nested model is fitted on rows of the other folds, and those rows carry base scores from models
whose own training sets included fold $i$. Information about fold $i$ therefore reaches the
nested fit through the covariate, and the composition columns have no analogous route.

The Methods name this and do not remove it. An external review called that the release's
central validity issue, and it was right about one thing in particular: a disclosure is not a
measurement. The floor experiment bounds the channel only for a 2-mer, which is the least
overfitting model we could have chosen, so it says nothing about the transformer.

WHAT THIS SCRIPT DOES. Runs the estimator again with the channel closed, and reports the
difference. For outer test fold $i$:

  * the covariate on each outer-TRAINING row $j$ comes from a base model trained on the folds
    excluding BOTH $i$ and $j$, so no information from the outer test fold can reach the fitted
    coefficients;
  * the covariate on each outer-TEST row comes from the base model trained on everything except
    fold $i$, which is the existing published score and is already honest for that row;
  * the composition block is unchanged, so the two arms of the comparison differ in exactly the
    thing being measured.

With five folds that is ten extra base fits per dataset, one per unordered pair $\\{i, j\\}$,
on top of the five the published run already does.

WHY THE k-MER CLASSES AND NOT ALL THREE. A base fit here is a penalised logistic regression and
costs milliseconds, so the full panel runs on a laptop. Doing the same for the CNN and
SpliceBERT means four times the GPU sweep, and the honest statement is that it has not been
done; what the k-mer classes establish is the SIGN and the SCALE of the channel where it can be
measured exactly, and whether the published claim -- that the channel is one-directional --
survives contact with an exact computation. It bounds nothing for a fine-tuned transformer, by
the same argument the Methods already make in the other direction.

TWO MODELS, BECAUSE THEY ASK DIFFERENT QUESTIONS. The 4-mer is the published headline model, so
its cross-fitted contribution is directly comparable to a number in the paper. The 2-mer is the
floor model, whose true contribution is zero by construction, so its cross-fitted value
separates the two explanations the Results give for a positive floor: the conditioning effect,
which cross-fitting leaves untouched, and this channel, which cross-fitting removes.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.baseline import kmer_matrix  # noqa: E402
from rbp.eval.nested import composition_features  # noqa: E402
from rbp.stats import fit_full, standardise  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "gc", "dn": "dinuc", "neg2": "neg2"}
KS = (4, 2)
C_DEFAULT = 1.0


def _fit_logit(X, y, C=C_DEFAULT):
    from sklearn.linear_model import LogisticRegression
    return LogisticRegression(max_iter=3000, C=C).fit(X, y)


def base_scores(X, y, folds):
    """Base-model scores under both fitting regimes, in one pass over the fold pairs.

    Returns (published, crossfit) where

      published[r]      score for row r from the model trained on every fold but r's own.
                        This is exactly what rbp.eval.baseline.oof_scores produces and what
                        every published number used.
      crossfit[i][r]    score for row r used when fold i is the OUTER TEST fold. For a row in
                        fold i this is published[r], which is already clean. For a row in any
                        other fold j it comes from the model trained on the complement of
                        {i, j}, which is the whole point.

    The complement models are symmetric in i and j, so there are ten of them and not twenty.
    """
    y = np.asarray(y, dtype=int)
    folds = np.asarray(folds)
    uniq = np.unique(folds)
    published = np.full(len(y), np.nan)
    for f in uniq:
        tr = folds != f
        if len(np.unique(y[tr])) < 2:
            continue
        published[folds == f] = _fit_logit(X[tr], y[tr]).decision_function(X[folds == f])

    crossfit = {int(i): published.copy() for i in uniq}
    for a in range(len(uniq)):
        for b in range(a + 1, len(uniq)):
            i, j = int(uniq[a]), int(uniq[b])
            tr = (folds != i) & (folds != j)
            if len(np.unique(y[tr])) < 2:
                continue
            m = _fit_logit(X[tr], y[tr])
            # The model that saw neither fold. It supplies fold j's covariate when i is the
            # outer test fold, and fold i's when j is.
            for outer, scored in ((i, j), (j, i)):
                sel = folds == scored
                crossfit[outer][sel] = m.decision_function(X[sel])
    return published, crossfit


def _pooled(comp, col_by_fold, y, folds):
    """Pooled out-of-fold AUROC of [composition | score], score supplied per outer fold.

    `col_by_fold[i]` is the score column to use when fold i is held out: its outer-training
    entries are what the nested model is fitted on, its test entries what it is scored with.
    """
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y, dtype=int)
    folds = np.asarray(folds)
    out = np.full(len(y), np.nan)
    for f in np.unique(folds):
        te, tr = folds == f, folds != f
        if len(np.unique(y[tr])) < 2:
            continue
        col = col_by_fold[int(f)]
        X = np.column_stack([comp, standardise(col)])
        beta, _ = fit_full(X[tr], y[tr], "l2")
        out[te] = np.column_stack([np.ones(te.sum()), X[te]]) @ beta
    ok = np.isfinite(out)
    return float(roc_auc_score(y[ok], out[ok]))


def _pooled_comp(comp, y, folds):
    from sklearn.metrics import roc_auc_score
    from rbp.eval.nested import _oof_scores
    s = _oof_scores(comp, y, folds, "l2")
    ok = np.isfinite(s)
    return float(roc_auc_score(np.asarray(y)[ok], s[ok]))


def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    datasets = list(pub.dataset)[:limit or None]
    rows = []
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
            comp, _ = composition_features(seqs, True, standardise_cols=True)
            a_comp = _pooled_comp(comp, y, folds)
            rec[f"comp_{arm}"] = a_comp
            rec[f"n_{arm}"] = int(len(d))
            for k in KS:
                X, _ = kmer_matrix(list(seqs), k)
                pubs, cf = base_scores(X, y, folds)
                naive = {int(i): pubs for i in np.unique(folds)}
                rec[f"k{k}_pub_{arm}"] = _pooled(comp, naive, y, folds) - a_comp
                rec[f"k{k}_cf_{arm}"] = _pooled(comp, cf, y, folds) - a_comp
        if not ok:
            continue
        rows.append(rec)
        log(f"[{n:3d}/{len(datasets)}] {ds:18s} " + "  ".join(
            f"{a} k4 {rec[f'k4_pub_{a}']:+.4f}->{rec[f'k4_cf_{a}']:+.4f}" for a in ARMS))
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

    per = TABLES / "cross_fitting_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it without --from-cache")
    else:
        t = build(a.store, a.n)
        t.to_csv(per, index=False)

    # Protein-clustered, matching every other interval in the paper: the fifteen proteins
    # measured in both cell lines contribute two correlated rows each.
    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(4000)]

    out = []

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})

    add("datasets", np.full(len(t), len(t)))
    for k in KS:
        for arm in ARMS:
            pub = t[f"k{k}_pub_{arm}"].to_numpy(float)
            cf = t[f"k{k}_cf_{arm}"].to_numpy(float)
            add(f"{k}-mer contribution as published, {arm} arm", pub)
            add(f"{k}-mer contribution fully cross-fitted, {arm} arm", cf)
            add(f"{k}-mer outer-fold channel, {arm} arm", pub - cf,
                "published minus cross-fitted; positive means the channel inflated the "
                "published value")
            out.append({"check": f"{k}-mer datasets where cross-fitting LOWERS the "
                                 f"contribution, {arm} arm",
                        "value": float((pub > cf).sum()), "ci_low": "", "ci_high": "",
                        "n": len(t), "note": "the one-directional claim predicts all of them"})

    # The span is the headline, and it is the quantity a common bias cannot move. Recomputed
    # end to end under cross-fitting rather than inferred from the per-arm shifts.
    for k in KS:
        pubm = {a: t[f"k{k}_pub_{a}"].mean() for a in ARMS}
        cfm = {a: t[f"k{k}_cf_{a}"].mean() for a in ARMS}

        def span(m):
            """max/min, or NaN where the arms straddle zero.

            A ratio of panel means is only a span while they share a sign. Cross-fitting can
            take the 2-mer -- whose truth is zero -- to a small negative mean in one arm, and
            dividing by that produced a NEGATIVE "span" that would be read as a magnitude.
            """
            lo, hi = min(m.values()), max(m.values())
            return float(hi / lo) if lo > 0 else float("nan")

        out.append({"check": f"{k}-mer three-arm span, as published",
                    "value": span(pubm), "ci_low": "", "ci_high": "", "n": len(t), "note": ""})
        out.append({"check": f"{k}-mer three-arm span, fully cross-fitted",
                    "value": span(cfm), "ci_low": "", "ci_high": "", "n": len(t),
                    "note": "the headline quantity, with the channel closed; NaN where an "
                            "arm's mean is not positive, which makes a ratio meaningless"})

    r = pd.DataFrame(out)
    r.to_csv(TABLES / "cross_fitting.csv", index=False)
    log("")
    for _, x in r.iterrows():
        ci = (f"  [{x.ci_low:+.4f}, {x.ci_high:+.4f}]"
              if isinstance(x.ci_low, float) else "")
        log(f"  {x['check']:62s} {x['value']:+.4f}{ci}")
    log("\nwrote cross_fitting.csv and cross_fitting_per_dataset.csv")


if __name__ == "__main__":
    main()
