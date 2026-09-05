"""The composition control: does a model's score predict binding beyond composition?

THE PROBLEM. Our negatives are matched on region, on GC content and on distance from
peaks, but GC constrains only G+C. It leaves the G/C split, the A/U split and every
dinucleotide frequency free, and we measured how much room that leaves: TARDBP positives
run GU at +1.94 log2 over negatives, PTBP1 runs AA at -1.51. A model can score well by
detecting "G-rich, C-poor" without learning anything sequence-specific, and the headline
AUROC cannot tell the two apart.

WHAT WE TRIED FIRST, AND WHY IT FAILED. The obvious control is to rebuild the negatives
with a dinucleotide-preserving shuffle, so composition is held fixed by construction and
any remaining signal must be positional. Two independent things broke it:

  * TARDBP's motif SURVIVES the shuffle. UGUGU appears in 92.2% of positives and still
    87.0% after shuffling, because the motif IS a dinucleotide repeat -- preserve UG and
    GU counts and you rebuild UGUGU. For repeat motifs the control cannot separate the
    two concepts even in principle.
  * Shuffled sequence is detectable per se. A classifier trained only to tell real
    negatives from shuffled negatives, with no positives involved at all, reaches AUROC
    ~0.69 (range 0.625-0.770). The shuffled arm therefore has a floor, and subtracting
    that floor from an AUROC is not valid arithmetic: AUROCs do not decompose additively
    and the floor was measured on a different discrimination.

WHAT WE DO INSTEAD. Ask the same question with regression, where "beyond" has an exact
meaning. Fit two nested models on the real data:

    reduced   bound ~ composition
    full      bound ~ composition + model_score

and test whether adding the score improves the fit. The coefficient on model_score is
the answer, and the likelihood-ratio statistic is its significance. Nothing is shuffled,
so there is no artefact to subtract.

The composition block is deliberately the mononucleotide and dinucleotide frequencies:
exactly what a dinucleotide shuffle would have preserved. So this asks the same
scientific question the shuffle was meant to ask, without the shuffle.

This is the same technique as the conservation control in variants/conservation.py --
fit with and without the nuisance variable, report the adjusted coefficient -- applied
to a different confounder. They share the estimator code in rbp.stats, which is the
point: one tested implementation, two uses.

WHAT THIS CONTROL DOES AND DOES NOT SAY. It measures information beyond mononucleotide
and dinucleotide composition. That is NOT the same as "is there a real motif", and
conflating the two would be a serious misreading. Measured on our own data: the single
UG dinucleotide frequency alone separates TARDBP's bound from unbound windows at AUROC
0.978, and PTBP1's CU frequency alone reaches 0.939. Those proteins bind (UG)n and
CU-rich tracts, so for them counting dinucleotides IS detecting the motif, and the
control correctly reports almost no additional information while the model is in fact
finding the real biology.

So the honest claim is about what a model adds over a composition baseline, and the
composition baseline has to be reported next to every headline AUROC rather than
mentioned in a limitations paragraph. Proteins whose motif is a short repeat are the
cases where the two concepts are genuinely inseparable, and they should be named.

A second caveat: adjusting for a confounder measured with error leaves residual
confounding. Mono+di frequencies are not the whole of "composition", so a surviving
coefficient is evidence of sequence-specific signal, not proof of a motif.
"""

from dataclasses import dataclass

import numpy as np

from ..stats import bootstrap_indices, coef_se, lr_test, percentile_ci, standardise, wald_p

ALPHABET = "ACGU"
DINUCS = [a + b for a in ALPHABET for b in ALPHABET]
BOOT = 2000


def _counts(seqs, k):
    """(n_seq x 4^k) k-mer count matrix, vectorised over sequences.

    Encodes each sequence as base-4 digits and reduces sliding windows to integers, which
    turns the whole panel into a couple of numpy passes. A per-sequence Python loop over
    900,000 windows is minutes; this is under a second.
    """
    lut = np.full(128, -1, dtype=np.int64)
    for i, c in enumerate(ALPHABET):
        lut[ord(c)] = i
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), -1)
    dig = lut[arr]
    n, L = dig.shape
    if L < k:
        return np.zeros((n, 4 ** k), dtype=np.int64)

    idx = np.zeros((n, L - k + 1), dtype=np.int64)
    valid = np.ones_like(idx, dtype=bool)
    for off in range(k):
        d = dig[:, off:off + L - k + 1]
        valid &= d >= 0
        idx = idx * 4 + np.where(d >= 0, d, 0)
    idx = np.where(valid, idx, 4 ** k)          # unknown bases go to a discard bin
    out = np.zeros((n, 4 ** k + 1), dtype=np.int64)
    np.add.at(out, (np.arange(n)[:, None], idx), 1)
    return out[:, :4 ** k]


def entropy(seqs):
    """Shannon entropy of each sequence's base distribution, in bits.

    Included because composition is not only frequencies: ELAVL1's signal is homopolymer
    runs (87.2% of positives contain a run of 5+ vs 59.3% of negatives) and TARDBP's is
    an alternating repeat, and both show up as low entropy at similar frequencies.
    """
    c = _counts(seqs, 1).astype(float)
    tot = c.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1
    p = c / tot
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.where(p > 0, p * np.log2(p), 0.0).sum(axis=1)
    return h


def composition_features(seqs, include_entropy=True, standardise_cols=True):
    """The nuisance block: mono and dinucleotide frequencies, plus entropy.

    One column from each frequency family is dropped. Frequencies sum to one, so keeping
    all of them makes the design matrix singular; the dropped level is absorbed into the
    intercept exactly as with any categorical reference level. This does not change what
    the block spans, so the test on model_score is unaffected -- but leaving the
    singularity in would make the pseudo-inverse silently pick one of infinitely many
    equivalent solutions, and the coefficients would not be reproducible.
    """
    seqs = list(seqs)
    mono = _counts(seqs, 1).astype(float)
    di = _counts(seqs, 2).astype(float)
    mono /= np.maximum(mono.sum(axis=1, keepdims=True), 1)
    di /= np.maximum(di.sum(axis=1, keepdims=True), 1)

    cols = [mono[:, :-1], di[:, :-1]]
    names = list(ALPHABET[:-1]) + DINUCS[:-1]
    if include_entropy:
        cols.append(entropy(seqs)[:, None])
        names.append("entropy")
    X = np.column_stack(cols)
    if not standardise_cols:
        return X, names
    return np.column_stack([standardise(X[:, j]) for j in range(X.shape[1])]), names


@dataclass
class Nested:
    coef: float
    ci_low: float
    ci_high: float
    se: float
    p_wald: float
    p_boot: float
    lr_stat: float
    lr_p: float
    n: int
    n_boot_ok: int

    @property
    def survives(self):
        """A positive coefficient whose interval excludes zero.

        The sign matters. A NEGATIVE coefficient that excludes zero means the model score
        anti-predicts binding once composition is accounted for, which is a failure, not
        a smaller success -- so a two-sided test alone would be the wrong summary.
        """
        return self.ci_low > 0


def test_score(seqs, score, label, blocks=None, n_boot=BOOT, seed=0, method="firth",
               include_entropy=True):
    """Does `score` predict `label` beyond the composition of `seqs`?

    `blocks` should be gene labels. Windows in the same gene are not independent -- our
    measured concentration reaches 8.2 windows per gene -- so the interval is a cluster
    bootstrap over genes. Bootstrapping rows instead would report an interval that is too
    narrow, in the direction that flatters the result.
    """
    y = np.asarray(label, dtype=int)
    comp, _ = composition_features(seqs, include_entropy)
    s = standardise(score)[:, None]
    full = np.column_stack([comp, s])

    coefs, ses = coef_se(full, y, method)
    coef, se = float(coefs[-1]), float(ses[-1])
    stat, _, lr_p = lr_test(comp, full, y, method)

    rng = np.random.default_rng(seed)
    boots = []
    for idx in bootstrap_indices(len(y), n_boot, rng, blocks):
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            b, _ = coef_se(full[idx], y[idx], method)
            boots.append(float(b[-1]))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            continue

    lo, hi, p_boot = percentile_ci(boots) if len(boots) >= 50 else (np.nan,) * 3
    return Nested(coef, lo, hi, se, wald_p(coef, se), p_boot,
                  stat, lr_p, len(y), len(boots))


def _fold_standardise(X, tr, te):
    """Centre and scale every column on the TRAINING rows only, applied to both sides.

    A column whose training rows are constant is left at zero rather than divided by zero.
    """
    mu = X[tr].mean(axis=0)
    sd = X[tr].std(axis=0)
    good = sd > 0
    out_tr = np.zeros((int(tr.sum()), X.shape[1]))
    out_te = np.zeros((int(te.sum()), X.shape[1]))
    out_tr[:, good] = (X[tr][:, good] - mu[good]) / sd[good]
    out_te[:, good] = (X[te][:, good] - mu[good]) / sd[good]
    return out_tr, out_te


def _oof_scores(X, y, folds, method="l2", within_fold=False):
    """Out-of-fold linear predictor for a design matrix, using the study's own folds.

    Every comparison in the paper has to be out-of-fold on both sides. An in-sample
    composition AUROC against an out-of-fold model AUROC flatters composition, and on
    TARDBP it did exactly that: in-sample composition reached 0.986 against the model's
    out-of-fold 0.984, which reads as "composition beats the model" and is an artefact.

    `within_fold` standardises each column on the training rows of each fold instead of over
    the whole dataset. The whole-dataset form is what every published number uses. It is
    label-free and applied identically to both arms, so it cannot manufacture a contribution,
    but it does let a test row's own mean and sd into its features, and a reviewer is entitled
    to ask what that costs. See scripts/nested_scale.py, which measures it.
    """
    from ..stats import fit_full
    y = np.asarray(y, dtype=int)
    folds = np.asarray(folds)
    out = np.full(len(y), np.nan)
    for f in np.unique(folds):
        te = folds == f
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        Xtr, Xte = _fold_standardise(X, tr, te) if within_fold else (X[tr], X[te])
        beta, _ = fit_full(Xtr, y[tr], method)
        out[te] = np.column_stack([np.ones(te.sum()), Xte]) @ beta
    return out


def composition_auroc(seqs, label, folds=None, method="l2", include_entropy=True):
    """AUROC of composition alone: the bar the model has to clear.

    Reported next to every model AUROC. Without it a model at 0.86 looks strong, and if
    composition alone reaches 0.84 the honest reading is completely different.

    Pass `folds` to get the out-of-fold value, which is the only one comparable to a
    model's out-of-fold AUROC. Omitting them returns the in-sample fit, which is useful
    only as a diagnostic upper bound.
    """
    from sklearn.metrics import roc_auc_score

    from ..stats import fit_full
    X, _ = composition_features(seqs, include_entropy)
    y = np.asarray(label, dtype=int)
    if folds is None:
        beta, _ = fit_full(X, y, method)
        s = np.column_stack([np.ones(len(X)), X]) @ beta
    else:
        s = _oof_scores(X, y, folds, method)
    ok = np.isfinite(s)
    return float(roc_auc_score(y[ok], s[ok]))


@dataclass
class Gain:
    """How much a model score adds to composition, as an effect size.

    The AUROC difference is the headline, not the p-value. With ~900,000 pooled pairs any
    non-zero contribution is "significant" -- EWSR1 alone gave lr_p = 2e-175 -- so a
    p-value says almost nothing about whether the contribution matters. A difference of
    +0.003 [+0.002, +0.004] and a difference of +0.140 [+0.130, +0.150] can carry the
    same p-value and mean entirely different things.
    """
    auroc_composition: float
    auroc_with_score: float
    delta: float
    delta_ci_low: float
    delta_ci_high: float
    delta_p: float
    n: int

    @property
    def helps(self):
        return self.delta_ci_low > 0


def to_logit(p, eps=1e-6):
    """log(p/(1-p)) for a probability, clipped off the endpoints.

    A neural model's sigmoid output saturates, and logit(0) is not a number. The clip is at
    1e-6, which bounds the transformed value at +/-13.8 and so cannot dominate a
    standardised column through a single saturated row.
    """
    q = np.clip(np.asarray(p, dtype=float), eps, 1.0 - eps)
    return np.log(q / (1.0 - q))


def gain_over_composition(seqs, score, label, folds, method="l2",
                          include_entropy=True, level=0.95,
                          score_scale="raw", within_fold=False):
    """Out-of-fold AUROC gain from adding `score` to a composition-only model.

    Both arms are fit and scored out-of-fold on the same folds, then compared with DeLong's
    PAIRED estimator. Pairing is required rather than optional: the two score vectors come from
    overlapping training data on identical test rows, are strongly correlated, and treating them
    as independent overstates the variance of their difference by about fourfold here.

    THAT IS NOT THE SAME AS DELONG BEING VALID FOR THIS COMPARISON, and this docstring used to
    say only the first thing while the Methods said only the second, so the code and the paper
    read as contradicting each other. Demler et al. (2012) show DeLong's test is not valid for
    strictly nested models, because under the null the added predictor's contribution
    degenerates. The exposure is bounded in two ways, both stated in the Methods: every
    panel-level interval -- which is where all the headlines are -- uses a protein-clustered
    bootstrap and not this standard error at all, and where it IS used, for the per-dataset
    significance counts, it is conservative against the recommended penalised likelihood-ratio
    alternative (80 significant against 89 in the GC arm).

    So: use the pairing, do not read the per-dataset p-value as a test of nested significance.

    TWO SPECIFICATION CHOICES, BOTH DEFAULTING TO WHAT THE PUBLISHED NUMBERS USED.

    `score_scale="logit"` transforms the score before it enters the design. This matters
    because logistic regression is NOT invariant to a nonlinear transform of a covariate,
    and the three model classes arrive on different scales: the 4-mer's out-of-fold score is
    already a log odds, while the CNN's and SpliceBERT's are sigmoid probabilities. Adding a
    probability where the natural scale is a log odds is a real inconsistency between arms of
    the model-class comparison, not a matter of taste. AUROC itself is invariant, so no
    model-alone number moves.

    `within_fold=True` standardises on each fold's training rows instead of over the whole
    dataset. See _oof_scores.
    """
    from .delong import delong_test

    y = np.asarray(label, dtype=int)
    comp, _ = composition_features(seqs, include_entropy,
                                   standardise_cols=not within_fold)
    if score_scale == "logit":
        score = to_logit(score)
    elif score_scale != "raw":
        raise ValueError(f"score_scale={score_scale!r}; expected 'raw' or 'logit'")
    col = np.asarray(score, dtype=float) if within_fold else standardise(score)
    full = np.column_stack([comp, col])

    s_comp = _oof_scores(comp, y, folds, method, within_fold)
    s_full = _oof_scores(full, y, folds, method, within_fold)
    ok = np.isfinite(s_comp) & np.isfinite(s_full)

    r = delong_test(s_full[ok], s_comp[ok], y[ok])
    from scipy.stats import norm
    z = norm.ppf(0.5 + level / 2.0)
    return Gain(r["auc_b"], r["auc_a"], r["diff"],
                r["diff"] - z * r["se"], r["diff"] + z * r["se"], r["p"], int(ok.sum()))


def run(datasets, model_cols, seq_col="seq_rna", label_col="label", fold_col="fold",
        gene_col=None, n_boot=BOOT, seed=0, method="firth"):
    """The composition control for every (dataset, model) pair.

    Reports BOTH views of the same question, deliberately: the out-of-fold AUROC gain
    (the effect size, and the headline) and the regression coefficient with its
    likelihood-ratio p-value (the parametric test, FDR corrected). They can disagree, and
    where they do the gain is the number to trust -- a tiny gain with an astronomical
    p-value means the sample is large, not that the model matters.
    """
    import pandas as pd

    from ..stats import benjamini_hochberg
    rows = []
    for name, df in datasets.items():
        seqs = df[seq_col].tolist()
        y = df[label_col].to_numpy()
        folds = df[fold_col].to_numpy()
        blocks = df[gene_col].to_numpy() if gene_col and gene_col in df else None
        for m in model_cols:
            s = df[m].to_numpy()
            g = gain_over_composition(seqs, s, y, folds)
            r = test_score(seqs, s, y, blocks=blocks, n_boot=n_boot, seed=seed,
                           method=method)
            rows.append({"dataset": name, "model": m,
                         "auroc_composition": g.auroc_composition,
                         "auroc_with_score": g.auroc_with_score,
                         "delta_auroc": g.delta, "delta_ci_low": g.delta_ci_low,
                         "delta_ci_high": g.delta_ci_high, "delta_p": g.delta_p,
                         "helps": g.helps, "coef": r.coef, "coef_ci_low": r.ci_low,
                         "coef_ci_high": r.ci_high, "lr_stat": r.lr_stat,
                         "lr_p": r.lr_p, "p_boot": r.p_boot, "n": r.n,
                         "survives": r.survives})
    out = pd.DataFrame(rows)
    if len(out):
        out["q_lr"] = benjamini_hochberg(out.lr_p.values)
        out["q_delta"] = benjamini_hochberg(out.delta_p.values)
        out["survives_fdr"] = (out.q_lr < 0.05) & (out.coef > 0)
    return out
