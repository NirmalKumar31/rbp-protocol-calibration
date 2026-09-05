"""Pre-registered robustness checks, run before any money is spent.

Four things that were promised and have to be verified rather than assumed:

  --what estimator   pooled out-of-fold AUROC vs the pair-weighted average of per-fold
                     AUROCs. `cv.estimator: pooled_oof` in params.yaml commits to the
                     first, with the second as a sensitivity check and a requirement that
                     they agree. Committing to a check and not running it is worse than
                     not committing.
  --what baseline    does the k-mer baseline get stronger with regularisation tuned per
                     dataset? It currently LOSES to 19 composition features on 6 small
                     datasets, which is a real weakness and an obvious reviewer attack.
  --what leakage     the sequence-similarity audit at panel scale. It was measured on 17
                     proteins; the claim is about 187.
  --what extremity   compositional extremity: how close the negatives could possibly be
                     in dinucleotide space, and what that predicts.

    python scripts/robustness.py --what all
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from rbp.eval import baseline  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
CELLS = ("K562", "HepG2")


def datasets(min_pairs, arm="dinuc"):
    """The arm selects both the panel file and the dataset directory; see
    rbp.utils.panel for why those two must not be chosen independently."""
    out, _missing = panelmod.datasets(CELLS, arm, min_pairs)
    return out


# ---------------------------------------------------------------------------------------

def estimator(paths, k):
    """Pooled out-of-fold AUROC against the pair-weighted per-fold average.

    They answer the same question differently. Pooled scores every pair once and is
    invariant to fold size, but mixes predictions from k different models. The per-fold
    average never mixes score scales but is noisy on small folds. If they disagree, the
    choice of estimator is doing work and has to be justified rather than declared.
    """
    rows = []
    for i, (name, path) in enumerate(sorted(paths.items()), 1):
        df = pd.read_csv(path, sep="\t")
        y, folds = df.label.to_numpy(), df.fold.to_numpy()
        s, _, _ = baseline.oof_scores(df.seq_rna.tolist(), y, folds, k=k)
        ok = np.isfinite(s)
        pooled = roc_auc_score(y[ok], s[ok])

        per, wts = [], []
        for f in np.unique(folds):
            sel = (folds == f) & ok
            if sel.sum() < 10 or len(np.unique(y[sel])) < 2:
                continue
            per.append(roc_auc_score(y[sel], s[sel]))
            wts.append(int(sel.sum()))
        if not per:
            continue
        per, wts = np.array(per), np.array(wts, dtype=float)
        rows.append({"dataset": name, "pooled": pooled,
                     "fold_mean_weighted": float(np.average(per, weights=wts)),
                     "fold_mean_unweighted": float(per.mean()),
                     "fold_sd": float(per.std(ddof=1)) if len(per) > 1 else np.nan,
                     "n_folds_used": len(per)})
        if i % 40 == 0:
            print(f"  [{i}/{len(paths)}]", flush=True)

    t = pd.DataFrame(rows)
    t["diff_weighted"] = t.pooled - t.fold_mean_weighted
    t["diff_unweighted"] = t.pooled - t.fold_mean_unweighted
    t.to_csv(TABLES / "robustness_estimator.csv", index=False)

    print(f"\n=== ESTIMATOR SENSITIVITY, {len(t)} datasets ===")
    print(f"pooled OOF AUROC          median {t.pooled.median():.4f}")
    print(f"pair-weighted fold mean   median {t.fold_mean_weighted.median():.4f}")
    print(f"unweighted fold mean      median {t.fold_mean_unweighted.median():.4f}")
    for c in ("diff_weighted", "diff_unweighted"):
        d = t[c]
        print(f"\n{c}: median {d.median():+.4f}  mean {d.mean():+.4f}  "
              f"max |{d.abs().max():.4f}|")
        print(f"  |diff| > 0.01: {int((d.abs() > 0.01).sum())}/{len(t)}   "
              f"> 0.02: {int((d.abs() > 0.02).sum())}/{len(t)}")
    print(f"\nrank agreement (Spearman): "
          f"{t[['pooled','fold_mean_weighted']].corr(method='spearman').iloc[0,1]:.4f}")
    worst = t.loc[t.diff_weighted.abs().idxmax()]
    print(f"worst disagreement: {worst.dataset} pooled {worst.pooled:.4f} vs "
          f"weighted {worst.fold_mean_weighted:.4f}")
    return t


# ---------------------------------------------------------------------------------------

def baseline_strength(paths, k, cs=(0.01, 0.1, 1.0, 10.0)):
    """Is the k-mer baseline as strong as it can be?

    It currently loses to 19 composition features on 6 small datasets, because 256 4-mer
    features overfit ~900 pairs at C=1.0. A baseline that is weaker than it needs to be
    makes every model look better than it is, and "your baseline was weak" is the first
    thing a reviewer says. Tuning C per dataset -- chosen on the VALIDATION fold, never on
    test -- fixes it or proves it cannot be fixed.
    """
    from sklearn.linear_model import LogisticRegression

    from rbp.data.splits import fold_roles

    rows = []
    for i, (name, path) in enumerate(sorted(paths.items()), 1):
        df = pd.read_csv(path, sep="\t")
        y, folds = df.label.to_numpy(), df.fold.to_numpy()
        X, _ = baseline.kmer_matrix(df.seq_rna.tolist(), k)
        n_folds = len(np.unique(folds))
        r = {"dataset": name, "pairs": int(len(df) // 2)}

        # fixed C: one value everywhere
        for C in cs:
            s, _, _ = baseline.oof_scores(df.seq_rna.tolist(), y, folds, k=k, C=C)
            ok = np.isfinite(s)
            r[f"C{C}"] = roc_auc_score(y[ok], s[ok])

        # tuned C: chosen on the VALIDATION fold, scored on the TEST fold. Selecting C by
        # out-of-fold AUROC -- which is what an earlier version of this function did -- is
        # selection on test and inflates the tuned number, exactly the flattery this check
        # exists to detect.
        tuned = np.full(len(y), np.nan)
        picked = []
        for f in range(n_folds):
            te_f, va_f, tr_f = fold_roles(f, n_folds)
            te, va = folds == te_f, folds == va_f
            tr = np.isin(folds, tr_f)
            if len(np.unique(y[tr])) < 2 or len(np.unique(y[va])) < 2 or va.sum() < 10:
                continue
            best, best_c = -1.0, cs[0]
            for C in cs:
                m = LogisticRegression(max_iter=3000, C=C).fit(X[tr], y[tr])
                a = roc_auc_score(y[va], m.decision_function(X[va]))
                if a > best:
                    best, best_c = a, C
            m = LogisticRegression(max_iter=3000, C=best_c).fit(X[tr], y[tr])
            tuned[te] = m.decision_function(X[te])
            picked.append(best_c)
        ok = np.isfinite(tuned)
        r["tuned_auroc"] = roc_auc_score(y[ok], tuned[ok]) if ok.sum() > 20 else np.nan
        r["best_C"] = max(set(picked), key=picked.count) if picked else np.nan
        r["default_auroc"] = r["C1.0"]
        r["oracle_auroc"] = max(r[f"C{C}"] for C in cs)   # selection-on-test, for contrast
        rows.append(r)
        if i % 40 == 0:
            print(f"  [{i}/{len(paths)}]", flush=True)

    t = pd.DataFrame(rows)
    t["improvement"] = t.tuned_auroc - t.default_auroc
    t["oracle_gap"] = t.oracle_auroc - t.tuned_auroc
    t.to_csv(TABLES / "robustness_baseline.csv", index=False)

    print(f"\n=== BASELINE STRENGTH, {len(t)} datasets ===")
    print(f"default C=1.0            median AUROC {t.default_auroc.median():.4f}")
    print(f"tuned on validation fold median AUROC {t.tuned_auroc.median():.4f}")
    print(f"oracle (selects on test) median AUROC {t.oracle_auroc.median():.4f}"
          f"   <- not usable, shown to size the selection bias")
    print(f"\nhonest improvement from tuning: median {t.improvement.median():+.4f}  "
          f"mean {t.improvement.mean():+.4f}  max {t.improvement.max():+.4f}")
    print(f"selection-on-test would have overstated it by "
          f"{t.oracle_gap.median():+.4f} (median)")
    print("\nby dataset size:")
    for lo, hi in ((0, 1000), (1000, 3000), (3000, 10000), (10000, 10 ** 9)):
        s = t[(t.pairs >= lo) & (t.pairs < hi)]
        if len(s):
            print(f"  {lo:5d}-{hi if hi < 10**9 else '+':>6} pairs: {len(s):3d} datasets, "
                  f"modal C={s.best_C.mode().iloc[0] if s.best_C.notna().any() else '-'}, "
                  f"median improvement {s.improvement.median():+.4f}")
    print("\nfixed-C comparison (all datasets, one C for everything):")
    for C in cs:
        print(f"  C={C:<6} median AUROC {t[f'C{C}'].median():.4f}")
    return t


# ---------------------------------------------------------------------------------------

def leakage(paths, kmer=32, sample=None):
    """Shared long k-mers between a dataset's folds.

    Chromosome-level folds guarantee no shared locus, but not that homologous sequence --
    repeats, paralogues, gene families -- is kept together. A 101-nt window holds 70
    distinct 32-mers, and two unrelated sequences sharing even one 32-mer has probability
    ~4^-32, so a shared 32-mer is proof of common origin rather than coincidence.

    Measured on 17 proteins earlier; the paper's claim covers 187.
    """
    from rbp.audit.leakage import build_reference, overlap_profile
    items = sorted(paths.items())
    if sample:
        rng = np.random.default_rng(7)
        items = [items[i] for i in rng.choice(len(items), min(sample, len(items)),
                                              replace=False)]
    rows = []
    for i, (name, path) in enumerate(items, 1):
        df = pd.read_csv(path, sep="\t")
        # fold 0 as the held-out set, the rest as reference: one representative iteration
        te = df[df.fold == 0].seq_rna.tolist()
        tr = df[df.fold != 0].seq_rna.tolist()
        if not te or not tr:
            continue
        ref = build_reference(tr, kmer)
        frac = overlap_profile(te, ref, kmer)
        rows.append({"dataset": name, "n_test": len(te), "n_train": len(tr),
                     "mean_shared": float(frac.mean()),
                     "frac_any": float((frac > 0).mean()),
                     "frac_gt50": float((frac > 0.5).mean()),
                     "frac_gt90": float((frac > 0.9).mean())})
        if i % 20 == 0:
            print(f"  [{i}/{len(items)}]", flush=True)

    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "robustness_leakage.csv", index=False)
    print(f"\n=== LEAKAGE AUDIT, {len(t)} datasets, {kmer}-mers, fold 0 held out ===")
    print(f"mean fraction of a held-out window's {kmer}-mers also in training: "
          f"{t.mean_shared.mean():.5f}")
    print(f"held-out windows sharing ANY {kmer}-mer:  "
          f"median {t.frac_any.median():.5f}  max {t.frac_any.max():.5f}")
    print(f"sharing >50% (near-duplicate):            "
          f"median {t.frac_gt50.median():.5f}  max {t.frac_gt50.max():.5f}")
    print(f"sharing >90% (effectively duplicate):     "
          f"median {t.frac_gt90.median():.5f}  max {t.frac_gt90.max():.5f}")
    print("\nworst datasets by near-duplicate fraction:")
    print(t.nlargest(5, "frac_gt50")[["dataset", "n_test", "frac_any",
                                      "frac_gt50"]].to_string(index=False))
    return t


# ---------------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--what", default="all",
                   choices=["estimator", "baseline", "leakage", "all"])
    p.add_argument("--k", type=int, default=4)
    # 0 MEANS ALL, AND IT IS THE DEFAULT. This defaulted to a 60-dataset random subsample, so
    # rerunning the script overwrote the committed 94-row table with a different and smaller
    # computation, silently. The manuscript quotes the full-panel figures.
    p.add_argument("--leakage-sample", type=int, default=0,
                   help="0 = every dataset in the panel; a positive N takes a random N")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    TABLES.mkdir(parents=True, exist_ok=True)

    paths = datasets(cfg.cv["min_pairs"])
    print(f"{len(paths)} datasets\n")

    if a.what in ("estimator", "all"):
        estimator(paths, a.k)
    if a.what in ("baseline", "all"):
        baseline_strength(paths, a.k)
    if a.what in ("leakage", "all"):
        leakage(paths, sample=a.leakage_sample)


if __name__ == "__main__":
    main()
