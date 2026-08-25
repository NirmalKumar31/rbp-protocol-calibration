"""DeLong's test for two AUROCs measured on the same data.

Why this and not a bootstrap or a t-test across proteins. When SpliceBERT scores 0.871
and RNA-FM scores 0.864 on the same test windows, the two AUROCs are *correlated*: they
are computed from the same positives and the same negatives, so a window that happens to
be easy makes both look good. Treating them as independent -- comparing their separate
confidence intervals, or running a two-sample test -- overstates the standard error of
the difference and turns real differences into "overlapping intervals, no conclusion".

DeLong's estimator computes the covariance between the two AUROCs directly from the
per-observation structural components, so the variance of the difference includes the
`-2 cov` term that the independent treatment throws away.

Implementation is the O(n log n) midrank form (Sun and Xu, 2014) rather than the
textbook O(n*m) double loop. At the panel level we compare models over ~900,000 pooled
pairs, where the double loop is 2e11 operations.

Scope: this compares models on ONE dataset. Comparing architectures across the whole
panel is a different question with a different unit of analysis (the protein, not the
window), and belongs in the hierarchical model.
"""

import numpy as np
from scipy.stats import norm


def midrank(x):
    """Ranks with ties averaged, 1-based, in the original order.

    Ties must get their average rank, not an arbitrary order. Model scores tie often --
    a saturated sigmoid returns exactly 1.0 for many windows -- and breaking those ties
    arbitrarily would make the AUROC depend on input ordering.
    """
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    sorted_x = x[order]
    n = len(x)
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and sorted_x[j] == sorted_x[i]:
            j += 1
        ranks[i:j] = 0.5 * (i + j - 1) + 1.0
        i = j
    out = np.empty(n, dtype=float)
    out[order] = ranks
    return out


def structural_components(scores, y):
    """(aucs, V10, V01) for one or more score vectors against shared labels.

    V10 has one entry per positive and V01 one per negative: each is that observation's
    contribution to the AUROC. Their covariance across models is what makes the paired
    comparison possible.
    """
    y = np.asarray(y)
    scores = np.atleast_2d(np.asarray(scores, dtype=float))
    pos = y == 1
    neg = ~pos
    m, n = int(pos.sum()), int(neg.sum())
    if m == 0 or n == 0:
        raise ValueError("need both classes present")

    k = scores.shape[0]
    aucs = np.empty(k)
    V10 = np.empty((k, m))
    V01 = np.empty((k, n))
    for i in range(k):
        px, nx = scores[i][pos], scores[i][neg]
        tz = midrank(np.concatenate([px, nx]))
        tx, ty = midrank(px), midrank(nx)
        aucs[i] = (tz[:m].sum() / m - (m + 1) / 2.0) / n
        V10[i] = (tz[:m] - tx) / n
        V01[i] = 1.0 - (tz[m:] - ty) / m
    return aucs, V10, V01


def auc_cov(scores, y):
    """(aucs, covariance matrix of the aucs)."""
    aucs, V10, V01 = structural_components(scores, y)
    m, n = V10.shape[1], V01.shape[1]
    # ddof=1: these are sample covariances of the structural components
    s10 = np.cov(V10, ddof=1) if V10.shape[0] > 1 else np.array([[np.var(V10[0], ddof=1)]])
    s01 = np.cov(V01, ddof=1) if V01.shape[0] > 1 else np.array([[np.var(V01[0], ddof=1)]])
    return aucs, np.atleast_2d(s10) / m + np.atleast_2d(s01) / n


def delong_test(score_a, score_b, y):
    """Compare two AUROCs on the same labels.

    Returns dict with both AUROCs, the difference, its standard error, z and a two-sided
    p-value. When the two score vectors are identical the variance of the difference is
    zero; that is reported as p=1 rather than a division by zero, because identical
    predictions are genuinely indistinguishable.
    """
    aucs, cov = auc_cov(np.vstack([score_a, score_b]), y)
    diff = float(aucs[0] - aucs[1])
    var = float(cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1])
    var = max(var, 0.0)                      # tiny negatives from floating point
    se = float(np.sqrt(var))
    if se == 0.0:
        return {"auc_a": float(aucs[0]), "auc_b": float(aucs[1]), "diff": diff,
                "se": 0.0, "z": 0.0, "p": 1.0}
    z = diff / se
    return {"auc_a": float(aucs[0]), "auc_b": float(aucs[1]), "diff": diff,
            "se": se, "z": float(z), "p": float(2.0 * norm.sf(abs(z)))}


def auc_ci(score, y, level=0.95):
    """AUROC with a DeLong confidence interval.

    Logit-transformed, so the interval cannot exceed [0, 1]. Near the ceiling -- and
    TARDBP reaches 0.988 -- a symmetric interval on the raw scale runs past 1.0 and has
    to be clipped, which quietly makes it too narrow on that side.
    """
    aucs, cov = auc_cov(np.atleast_2d(score), y)
    a = float(aucs[0])
    se = float(np.sqrt(max(cov[0, 0], 0.0)))
    z = norm.ppf(0.5 + level / 2.0)
    if se == 0.0 or not 0.0 < a < 1.0:
        return a, max(a - z * se, 0.0), min(a + z * se, 1.0)
    eta = np.log(a / (1.0 - a))
    se_eta = se / (a * (1.0 - a))
    lo, hi = eta - z * se_eta, eta + z * se_eta
    return a, float(1.0 / (1.0 + np.exp(-lo))), float(1.0 / (1.0 + np.exp(-hi)))


def pairwise(scores_by_model, y, correct=True):
    """All pairwise DeLong comparisons for a dict of model -> scores.

    With 5 architectures there are 10 comparisons per dataset, so the p-values are BH
    corrected within the family by default. Reporting 10 uncorrected p-values and quoting
    the smallest is the standard way to manufacture a significant result.
    """
    import pandas as pd

    from ..stats import benjamini_hochberg
    names = list(scores_by_model)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            r = delong_test(scores_by_model[a], scores_by_model[b], y)
            rows.append({"model_a": a, "model_b": b, **r})
    out = pd.DataFrame(rows)
    if correct and len(out):
        out["q"] = benjamini_hochberg(out.p.values)
    return out
