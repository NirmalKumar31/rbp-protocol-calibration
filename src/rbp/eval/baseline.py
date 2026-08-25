"""The k-mer baseline, fit under the study's own cross-validation folds.

Two jobs. First, it is the honest floor for the benchmark: if a bag of 4-mers already
reaches 0.96 on a protein, a transformer beating it by 0.01 is a small gain and the paper
has to say so. Second, and the reason this exists before any GPU work, it is a complete
stand-in for a trained model -- it produces out-of-fold binding scores AND variant delta
scores -- so the entire downstream pipeline can be exercised end to end for free.

The fold discipline is the part that is easy to get wrong. A variant sits at a genomic
position, that position is on a chromosome, and that chromosome belongs to exactly one
fold. Scoring the variant with a model that trained on its chromosome leaks. So variant
scoring uses the fold model that did NOT see that chromosome, which is the same rule as
out-of-fold prediction applied to a different kind of input.
"""

import numpy as np

C_DEFAULT = 1.0


def kmer_matrix(seqs, k):
    """Sparse k-mer count matrix and the fitted vectoriser.

    The vectoriser is returned because variant scoring has to transform new sequences
    with the SAME vocabulary. Refitting on the variant windows would produce a different
    feature space and the coefficients would refer to the wrong columns -- silently, since
    the shapes can still match.
    """
    from sklearn.feature_extraction.text import CountVectorizer
    vec = CountVectorizer(analyzer="char", ngram_range=(k, k), lowercase=False)
    return vec.fit_transform(seqs), vec


def fit_fold_models(seqs, y, folds, k=4, C=C_DEFAULT):
    """One model per fold, each trained on every fold except its own.

    Returns (models, vectoriser) where models[f] is the model that has NOT seen fold f,
    so it is the correct scorer for anything on fold f's chromosomes.
    """
    from sklearn.linear_model import LogisticRegression
    X, vec = kmer_matrix(list(seqs), k)
    y = np.asarray(y, dtype=int)
    folds = np.asarray(folds)
    models = {}
    for f in np.unique(folds):
        tr = folds != f
        if len(np.unique(y[tr])) < 2:
            continue
        models[int(f)] = LogisticRegression(max_iter=3000, C=C).fit(X[tr], y[tr])
    return models, vec


def oof_scores(seqs, y, folds, k=4, C=C_DEFAULT):
    """Out-of-fold decision values: every row scored exactly once, by a model that
    never saw its chromosome."""
    seqs = list(seqs)
    X, vec = kmer_matrix(seqs, k)
    models, _ = fit_fold_models(seqs, y, folds, k, C)
    folds = np.asarray(folds)
    out = np.full(len(seqs), np.nan)
    for f, m in models.items():
        te = folds == f
        out[te] = m.decision_function(X[te])
    return out, models, vec


def score_sequences(model, vec, seqs):
    """Decision values for new sequences under an existing model and vocabulary."""
    return model.decision_function(vec.transform(list(seqs)))


def variant_delta(models, vec, ref_seqs, alt_seqs, folds):
    """Predicted binding disruption per variant: score(reference) - score(alternate).

    Signed rather than absolute, because the downstream test standardises |delta| and
    keeping the sign here lets a caller check direction. Variants whose fold has no
    model, or that fall on an excluded chromosome, come back NaN rather than being
    quietly scored by the wrong model.
    """
    ref_seqs, alt_seqs = list(ref_seqs), list(alt_seqs)
    folds = np.asarray(folds, dtype=float)          # NaN marks "no fold for this variant"
    Xr, Xa = vec.transform(ref_seqs), vec.transform(alt_seqs)
    out = np.full(len(ref_seqs), np.nan)
    for f in np.unique(folds[np.isfinite(folds)]):
        m = models.get(int(f))
        if m is None:
            continue
        sel = folds == f
        out[sel] = m.decision_function(Xr[sel]) - m.decision_function(Xa[sel])
    return out


def evaluate(df, k=4, C=C_DEFAULT, seq_col="seq_rna", label_col="label",
             fold_col="fold"):
    """Out-of-fold AUROC with a DeLong interval, plus the scores themselves."""
    from .delong import auc_ci
    y = df[label_col].to_numpy()
    s, models, vec = oof_scores(df[seq_col].tolist(), y, df[fold_col].to_numpy(), k, C)
    ok = np.isfinite(s)
    auroc, lo, hi = auc_ci(s[ok], y[ok])
    return {"k": k, "n": int(ok.sum()), "auroc": auroc, "ci_low": lo, "ci_high": hi,
            "scores": s, "models": models, "vec": vec}


def sweep_k(df, ks=(3, 4, 5, 6), **kw):
    """Out-of-fold AUROC for several k, to see what motif width carries the signal."""
    out = {}
    for k in ks:
        r = evaluate(df, k=k, **kw)
        out[k] = {"auroc": r["auroc"], "ci_low": r["ci_low"], "ci_high": r["ci_high"],
                  "n": r["n"]}
    return out
