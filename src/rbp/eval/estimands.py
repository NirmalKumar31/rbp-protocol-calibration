"""The nested contribution, measured four other ways, because AUROC is a choice.

WHY. Every number in this study is an AUROC increment. Two reviewers made the same point
independently: "no protocol-free measure of contribution exists" is a claim about AUROC only,
and a rank-based summary has known blind spots -- it is insensitive to calibration, it weights
the whole ROC equally when only the top matters, and an increment in it is bounded above by
1 - baseline in a way that manufactures part of any comparison between unequal baselines.

So the same nested comparison is computed on four further estimands, all from the same two
out-of-fold linear predictors so nothing is refitted:

  delta deviance      2 * (LL_full - LL_comp). The likelihood-ratio statistic, on the scale the
                      model is actually fitted on. Unbounded above, so it cannot be squeezed
                      by a high baseline the way an AUROC increment is.
  McFadden increment  The same thing normalised by the null deviance, so it is comparable
                      across datasets of different size and prevalence.
  average precision   AP(full) - AP(comp). Prevalence-sensitive by design, and the summary a
                      reader who cares about the top of the ranking would use.
  IDI                 Integrated discrimination improvement: the change in
                      (mean predicted risk among positives) - (mean among negatives). A
                      calibration-sensitive summary from the risk-prediction literature, which
                      is where the nested-model question was formalised first.

AND ONE DIAGNOSTIC THAT IS NOT AN INCREMENT. `residual_auroc` regresses the model score on the
composition block and takes the AUROC of the residual. It answers a different and sharper
question: does the part of the score ORTHOGONAL to composition discriminate at all? An
increment can be small because the baseline is high; a residual AUROC at 0.5 means there is
nothing there.
"""

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

CLIP = 30.0


def _prob(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -CLIP, CLIP)))


def _loglik(y, z):
    """Bernoulli log likelihood of labels y under linear predictor z."""
    p = np.clip(_prob(z), 1e-12, 1.0 - 1e-12)
    return float(np.sum(y * np.log(p) + (1 - y) * np.log1p(-p)))


def alternatives(y, s_comp, s_full, comp=None, score=None):
    """Five summaries of the same nested comparison, from the same fitted predictors.

    `s_comp` and `s_full` are out-of-fold linear predictors. `comp` and `score` are needed
    only for the residualised diagnostic and may be omitted.
    """
    y = np.asarray(y, dtype=int)
    ll_c, ll_f = _loglik(y, s_comp), _loglik(y, s_full)
    base = float(y.mean())
    ll_0 = _loglik(y, np.full(len(y), np.log(base / (1 - base))))

    out = {
        "auroc_gain": float(roc_auc_score(y, s_full) - roc_auc_score(y, s_comp)),
        "delta_deviance": 2.0 * (ll_f - ll_c),
        "mcfadden_gain": (1.0 - ll_f / ll_0) - (1.0 - ll_c / ll_0) if ll_0 != 0 else np.nan,
        "ap_gain": float(average_precision_score(y, s_full)
                         - average_precision_score(y, s_comp)),
    }

    # IDI, on the probability scale both models are fitted to produce.
    pc, pf = _prob(s_comp), _prob(s_full)
    pos, neg = y == 1, y == 0
    if pos.any() and neg.any():
        out["idi"] = float((pf[pos].mean() - pf[neg].mean())
                           - (pc[pos].mean() - pc[neg].mean()))
    else:
        out["idi"] = np.nan

    # The residual diagnostic: strip everything the composition block can linearly explain.
    if comp is not None and score is not None:
        X = np.column_stack([np.ones(len(y)), np.asarray(comp, dtype=float)])
        sc = np.asarray(score, dtype=float)
        beta, *_ = np.linalg.lstsq(X, sc, rcond=None)
        resid = sc - X @ beta
        out["residual_auroc"] = (float(roc_auc_score(y, resid))
                                 if np.std(resid) > 0 else 0.5)
    else:
        out["residual_auroc"] = np.nan
    return out
