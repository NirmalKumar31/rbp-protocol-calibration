"""Estimators and tests shared by the variant analysis and the binding analysis.

These started life inside `variants/conservation.py`, which is where the
conservation-controlled test lives. The composition-controlled binding test turned out
to need exactly the same machinery -- fit two nested models, compare them, bootstrap
the interval by cluster, correct across the family of tests -- so they were lifted here
rather than copied. That is not a tidying exercise: the two controls being literally the
same code is the argument that they are the same technique applied to different
confounders, and a second copy would be free to drift away from the tested one.

`conservation.py` re-exports what it always exported, so nothing that imported from
there had to change.
"""

import numpy as np
from scipy.stats import chi2, norm
from sklearn.linear_model import LogisticRegression

CLIP = 30.0        # logit clip; exp(30) is already 1e13, so this only guards overflow
METHODS = ("firth", "l2", "none")


def standardise(x):
    """Standardise, so coefficients are comparable across variables of different scale."""
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return np.zeros_like(x) if sd == 0 else (x - x.mean()) / sd


def _design(X):
    return np.column_stack([np.ones(len(X)), np.asarray(X, dtype=float)])


def _prob(Xd, beta):
    return 1.0 / (1.0 + np.exp(-np.clip(Xd @ beta, -CLIP, CLIP)))


# ---------------------------------------------------------------------------------------
# Estimators
# ---------------------------------------------------------------------------------------

def firth_fit(X, y, max_iter=200, tol=1e-8):
    """Firth penalised logistic regression. Returns (beta_with_intercept, se).

    Newton-Raphson on the score equation modified by a Jeffreys prior:

        U*(b) = X' (y - p + h * (0.5 - p))

    where h is the diagonal of the weighted hat matrix. The `h * (0.5 - p)` term nudges
    each observation slightly toward 0.5, which is what keeps coefficients finite when
    the classes are separable. This is the standard fix for rare-event logistic
    regression, and our groups can have as few as 8 events.
    """
    Xd = _design(X)
    y = np.asarray(y, dtype=float)
    beta = np.zeros(Xd.shape[1])

    for _ in range(max_iter):
        p = _prob(Xd, beta)
        w = p * (1.0 - p)
        cov = np.linalg.pinv(Xd.T @ (Xd * w[:, None]))
        h = np.einsum("ij,jk,ik->i", Xd, cov, Xd) * w      # hat-matrix diagonal
        step = cov @ (Xd.T @ (y - p + h * (0.5 - p)))
        beta = beta + step
        if not np.all(np.isfinite(beta)):
            raise FloatingPointError("Firth diverged")
        if np.max(np.abs(step)) < tol:
            break

    p = _prob(Xd, beta)
    w = p * (1.0 - p)
    cov = np.linalg.pinv(Xd.T @ (Xd * w[:, None]))
    return beta, np.sqrt(np.clip(np.diag(cov), 0.0, None))


def coef_se(X, y, method):
    """Coefficient vector (no intercept) and standard errors, for the chosen estimator."""
    if method == "firth":
        beta, se = firth_fit(X, y)
        return beta[1:], se[1:]
    # sklearn >=1.8 deprecated `penalty`; C=inf is unpenalised, C=1.0 is the old L2 default
    C = np.inf if method == "none" else 1.0
    m = LogisticRegression(C=C, max_iter=3000).fit(X, y)
    return m.coef_[0], np.full(X.shape[1], np.nan)


def fit_full(X, y, method="firth"):
    """Fit and return (beta_with_intercept, se_with_intercept), any estimator."""
    if method == "firth":
        return firth_fit(X, y)
    C = np.inf if method == "none" else 1.0
    m = LogisticRegression(C=C, max_iter=3000).fit(X, y)
    beta = np.concatenate([m.intercept_, m.coef_[0]])
    return beta, np.full(len(beta), np.nan)


# ---------------------------------------------------------------------------------------
# Likelihood and the nested-model comparison
# ---------------------------------------------------------------------------------------

def loglik(X, y, beta):
    """Binomial log-likelihood at `beta`."""
    Xd = _design(X)
    y = np.asarray(y, dtype=float)
    eta = np.clip(Xd @ beta, -CLIP, CLIP)
    # log(1+exp(eta)) computed stably
    return float(np.sum(y * eta - np.logaddexp(0.0, eta)))


def penalised_loglik(X, y, beta):
    """Firth's penalised log-likelihood: l(b) + 0.5 log|I(b)|.

    Only comparable between fits that share the SAME design matrix. The log-determinant
    of the information matrix grows with the number of columns, so evaluating this on a
    2-column model and a 3-column model and subtracting measures the extra column's
    existence rather than its usefulness. See `lr_test`.
    """
    Xd = _design(X)
    p = _prob(Xd, beta)
    w = p * (1.0 - p)
    info = Xd.T @ (Xd * w[:, None])
    sign, logdet = np.linalg.slogdet(info)
    if sign <= 0:                      # singular information matrix
        logdet = -np.inf
    return loglik(X, y, beta) + 0.5 * logdet


def firth_fit_constrained(X, y, zero_cols, max_iter=200, tol=1e-8):
    """Firth fit with some coefficients held at zero, using the FULL design's penalty.

    This is what makes a penalised likelihood-ratio test valid. The naive approach fits
    the reduced model on its own smaller design, but then the two penalties are
    log-determinants of differently sized matrices and the difference is dominated by
    dimension rather than fit. Measured: that version rejected a pure-noise column in
    100% of null simulations at alpha=0.05.

    Holding the coefficient at zero inside the full design keeps the penalty functional
    identical on both sides, so it cancels. `zero_cols` indexes the columns of X, not of
    the design matrix.
    """
    Xd = _design(X)
    free = np.ones(Xd.shape[1], dtype=bool)
    free[[c + 1 for c in np.atleast_1d(zero_cols)]] = False
    y = np.asarray(y, dtype=float)
    beta = np.zeros(Xd.shape[1])
    ix = np.ix_(free, free)

    for _ in range(max_iter):
        p = _prob(Xd, beta)
        w = p * (1.0 - p)
        info = Xd.T @ (Xd * w[:, None])
        # the hat matrix uses the full design, so the penalty is the full model's
        h = np.einsum("ij,jk,ik->i", Xd, np.linalg.pinv(info), Xd) * w
        score = Xd.T @ (y - p + h * (0.5 - p))
        step = np.linalg.pinv(info[ix]) @ score[free]
        beta[free] += step
        if not np.all(np.isfinite(beta)):
            raise FloatingPointError("constrained Firth diverged")
        if np.max(np.abs(step)) < tol:
            break
    return beta


def lr_test(X_reduced, X_full, y, method="firth"):
    """Likelihood-ratio test of the extra columns in X_full. Returns (stat, df, p).

    `X_reduced` must be the leading columns of `X_full`; the test is on the trailing ones.

    For method="firth" this is the PROFILE penalised likelihood ratio: the null fit holds
    the tested coefficients at zero inside the full design rather than refitting a smaller
    model, so the penalty cancels instead of biasing the statistic. Verified against the
    uniform null in tests/unit/test_nested.py.

    NON-CONVERGENCE IS REPORTED IN BOTH DIRECTIONS, and it was not always. A negative
    statistic cannot happen for nested models sharing an objective, so it signals a
    convergence failure and is clipped to zero, i.e. reported as no evidence. An INFINITE
    statistic is the same failure with the opposite sign, and it used to pass straight
    through: chi2.sf(inf) is 0.0, so a diverged fit was reported as overwhelming evidence.
    The guard was therefore asymmetric in the direction that flatters a result.

    A non-finite statistic now returns nan for the p-value, so the failure is visible
    instead of being read as either extreme. Callers that treat nan as "no evidence" behave
    conservatively; callers that propagate it make the failure obvious.

    This matters on near-separable designs, which is exactly where the control is most
    interesting: a score that separates the classes almost perfectly drives the penalised
    likelihood toward its bound, and whether the fit diverges then depends on BLAS
    summation order and so on the machine. tests/unit/test_nested.py's pure-composition
    fixture sits on that boundary and failed intermittently in CI for that reason.
    """
    df = X_full.shape[1] - X_reduced.shape[1]
    if df <= 0:
        raise ValueError("X_full must have more columns than X_reduced")
    k = X_reduced.shape[1]
    if not np.allclose(X_full[:, :k], X_reduced):
        raise ValueError("X_reduced must be the leading columns of X_full")

    if method == "firth":
        b_f, _ = firth_fit(X_full, y)
        b_r = firth_fit_constrained(X_full, y, list(range(k, X_full.shape[1])))
        stat = 2.0 * (penalised_loglik(X_full, y, b_f)
                      - penalised_loglik(X_full, y, b_r))
    else:
        b_r, _ = fit_full(X_reduced, y, method)
        b_f, _ = fit_full(X_full, y, method)
        stat = 2.0 * (loglik(X_full, y, b_f) - loglik(X_reduced, y, b_r))

    if not np.isfinite(stat):
        return float("nan"), df, float("nan")
    stat = float(max(stat, 0.0))
    return stat, df, float(chi2.sf(stat, df))


def wald_p(coef, se):
    if se is None or not np.isfinite(se) or se == 0:
        return np.nan
    return float(2.0 * norm.sf(abs(coef / se)))


# ---------------------------------------------------------------------------------------
# Resampling
# ---------------------------------------------------------------------------------------

def bootstrap_indices(n, n_boot, rng, blocks=None):
    """Yield resampled row indices, by row or by whole block.

    Block resampling is the honest choice when rows are not independent -- windows in the
    same gene share expression, local sequence and a conservation baseline. Resampling
    rows understates the uncertainty; our measured concentration reaches 8.2 windows per
    gene, so the understatement is not small.
    """
    if blocks is None:
        for _ in range(n_boot):
            yield rng.integers(0, n, n)
        return
    blocks = np.asarray(blocks)
    uniq = np.unique(blocks)
    member = {b: np.flatnonzero(blocks == b) for b in uniq}
    for _ in range(n_boot):
        picked = rng.choice(uniq, size=len(uniq), replace=True)
        yield np.concatenate([member[b] for b in picked])


def percentile_ci(values, ci=(2.5, 97.5)):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return np.nan, np.nan, np.nan
    lo, hi = np.percentile(v, ci)
    # two-sided bootstrap p, floored at 1/n since we cannot resolve below that
    tail = min((v <= 0).mean(), (v >= 0).mean())
    return float(lo), float(hi), float(max(2.0 * tail, 1.0 / len(v)))


# ---------------------------------------------------------------------------------------
# Multiple comparisons
# ---------------------------------------------------------------------------------------

def benjamini_hochberg(pvals):
    """BH-adjusted p-values (q-values), input order preserved, NaNs passed through.

    With 85 (protein, model) tests at a 95% interval, about four will exclude zero by
    chance. BH controls the expected proportion of false discoveries among those called
    significant, which is the right family-wise notion here: we are screening many pairs
    and want the surviving set to be mostly real, not to guarantee zero errors.
    """
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]     # enforce monotonicity
    out = np.empty(m)
    out[order] = np.clip(adj, 0.0, 1.0)
    q[ok] = out
    return q
