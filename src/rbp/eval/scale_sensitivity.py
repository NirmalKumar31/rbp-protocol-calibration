"""Two specification choices in the nested fit, measured rather than argued.

BOTH ARE REAL OBJECTIONS AND THEY ARE DIFFERENT IN KIND.

1. THE COVARIATE SCALE (an inconsistency, so it gets fixed). The nested model is a logistic
   regression on nineteen composition features plus one model score. Logistic regression is
   not invariant to a nonlinear transform of a covariate, and the three model classes arrive
   on different scales: the 4-mer's out-of-fold score is `decision_function`, already a log
   odds, while the CNN's and SpliceBERT's are sigmoid probabilities. So the model-class
   comparison put a log odds in one arm and a probability in the others. AUROC is invariant to
   monotone transforms, so no model-alone number is affected; the nested increment is not.

2. THE STANDARDISATION WINDOW (improper, so it gets bounded). Columns are centred and scaled
   over the whole dataset before the out-of-fold loop, so a test row's own mean and standard
   deviation enter its features. It is label-free and applied identically to every arm, so it
   cannot manufacture a contribution -- but "cannot manufacture" is an argument, and a
   reviewer is entitled to a number.

The 2x2 is computed on one pass so the composition baseline is fitted twice per (dataset, arm)
rather than once per cell, which is where the cost is.
"""

import numpy as np
from sklearn.metrics import roc_auc_score

from ..stats import standardise
from .delong import delong_test
from .nested import _oof_scores, composition_features, to_logit

# The 4-mer is native log odds, so a logit transform of it is meaningless and it carries one
# scale only. The neural models carry two, and their difference IS item 1.
SCALES = {"kmer": ("native",), "cnn": ("probability", "logit"),
          "splicebert": ("probability", "logit")}


def dataset_rows(dataset, arm, seqs, label, folds, scores, method="l2"):
    """One row per (model, scale, standardisation) for a single (dataset, arm).

    `scores` maps model name to an array aligned with `seqs`.
    """
    y = np.asarray(label, dtype=int)
    rows = []
    for within in (False, True):
        comp, _ = composition_features(seqs, True, standardise_cols=not within)
        s_comp = _oof_scores(comp, y, folds, method, within)
        okc = np.isfinite(s_comp)
        auc_comp = float(roc_auc_score(y[okc], s_comp[okc]))
        for model, raw in scores.items():
            for scale in SCALES[model]:
                v = to_logit(raw) if scale == "logit" else np.asarray(raw, dtype=float)
                col = v if within else standardise(v)
                s_full = _oof_scores(np.column_stack([comp, col]), y, folds, method, within)
                ok = okc & np.isfinite(s_full)
                r = delong_test(s_full[ok], s_comp[ok], y[ok])
                rows.append({
                    "dataset": dataset, "arm": arm, "model": model, "scale": scale,
                    "standardisation": "within_fold" if within else "whole_dataset",
                    "comp_auroc": auc_comp, "full_auroc": r["auc_a"], "gain": r["diff"],
                    "se": r["se"], "n": int(ok.sum()),
                })
    return rows
