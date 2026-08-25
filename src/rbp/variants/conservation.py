"""The conservation-controlled test: the scientific core of the project.

Question: does a model's predicted binding disruption (delta) predict ClinVar
pathogenicity beyond what sequence conservation already explains?

It has to be asked this way because conservation is a confounder. Disruptive variants sit
at conserved positions, and conserved positions are more likely pathogenic -- for reasons
that may have nothing to do with this protein. In our data conservation ALONE separates
pathogenic from benign noncoding variants at AUROC 0.954, so a delta score can look
predictive purely by tracking it.

So we fit two models per group and compare delta's coefficient:

    alone       label ~ |delta|
    controlled  label ~ |delta| + conservation

A signal "survives" when the controlled coefficient's interval excludes zero with the
expected positive sign.

Three deliberate departures from a naive implementation:

  * Firth penalised likelihood by default. After restriction some groups have very few
    pathogenic variants, where maximum likelihood is biased and can diverge outright
    under separation. Firth keeps estimates finite. sklearn's L2 default is also
    available, but it shrinks coefficients toward zero -- which for a paper arguing
    "the signal mostly is not there" would bias in our own favour.
  * 2000 bootstrap resamples, not 200. At 200 the interval bounds moved 5-10% between
    runs, which was enough for a quoted interval to disagree with its own figure.
  * Benjamini-Hochberg FDR across the whole results table. With 85 tests at a 95%
    interval, roughly four will exclude zero by chance.
"""

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

# The estimators, the bootstrap and BH live in rbp.stats because the composition control
# in eval/nested.py needs exactly the same ones. Re-exported here under their original
# names so callers and tests are unaffected by the move.
from ..stats import METHODS, benjamini_hochberg, firth_fit  # noqa: F401
from ..stats import bootstrap_indices, coef_se as _coef_se, standardise as _z, wald_p

BOOT = 2000
CI = (2.5, 97.5)


# ---------------------------------------------------------------------------------------
# The test
# ---------------------------------------------------------------------------------------

@dataclass
class Fit:
    coef: float
    ci_low: float
    ci_high: float
    se: float
    p_wald: float
    p_boot: float
    n_boot_ok: int

    @property
    def survives(self):
        """Positive sign and an interval that excludes zero."""
        return self.ci_low > 0

    @property
    def significant(self):
        return self.ci_low > 0 or self.ci_high < 0


def fit_delta_coef(delta, label, conservation=None, n_boot=BOOT, seed=0,
                   method="firth", blocks=None):
    """Coefficient of |delta| predicting label, optionally controlling for conservation.

    Returns a Fit with the point estimate, a percentile bootstrap interval, and two
    p-values: a Wald p from the Firth standard error, and a bootstrap p from the
    proportion of resamples on the wrong side of zero. They answer the same question
    differently, so agreement between them is reassuring and disagreement is a warning.

    `blocks` makes the bootstrap CLUSTER-AWARE. Variants in the same gene are not
    independent -- they share expression, local sequence and a conservation baseline --
    so resampling individual rows understates the uncertainty. Passing gene labels
    resamples whole genes instead, which is the standard cluster bootstrap and is the
    only honest choice given our measured concentration (up to 8.2 windows per gene).
    """
    y = np.asarray(label, dtype=int)
    cols = [_z(np.abs(delta))]
    if conservation is not None:
        cols.append(_z(conservation))
    X = np.column_stack(cols)

    coef, se = _coef_se(X, y, method)
    coef, se = float(coef[0]), float(se[0])
    p_wald = wald_p(coef, se)

    rng = np.random.default_rng(seed)
    boots = []
    for idx in bootstrap_indices(len(y), n_boot, rng, blocks):
        if len(np.unique(y[idx])) < 2:            # a one-class resample cannot be fit
            continue
        try:
            b, _ = _coef_se(X[idx], y[idx], method)
            boots.append(float(b[0]))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            continue

    if len(boots) < 50:
        return Fit(coef, np.nan, np.nan, se, p_wald, np.nan, len(boots))

    b = np.asarray(boots)
    lo, hi = np.percentile(b, CI)
    # two-sided bootstrap p, floored at 1/n_boot since we cannot resolve below that
    tail = min((b <= 0).mean(), (b >= 0).mean())
    p_boot = float(max(2.0 * tail, 1.0 / len(b)))
    return Fit(coef, float(lo), float(hi), se, p_wald, p_boot, len(b))


def test_group(delta, label, conservation, seed=0, n_boot=BOOT, method="firth",
               min_n=20):
    """Both fits plus descriptive stats for one (group, model) pair."""
    delta = np.asarray(delta, dtype=float)
    label = np.asarray(label, dtype=int)
    conservation = np.asarray(conservation, dtype=float)

    keep = ~(np.isnan(delta) | np.isnan(conservation))
    delta, label, conservation = delta[keep], label[keep], conservation[keep]

    out = {"n": int(len(label)), "n_pathogenic": int(label.sum()),
           "n_benign": int((label == 0).sum()), "n_dropped_nan": int((~keep).sum()),
           "corr_delta_conservation": np.nan, "conservation_coef": np.nan,
           "conservation_auroc": np.nan, "note": ""}

    if len(label) < min_n or len(np.unique(label)) < 2:
        out["note"] = "insufficient data"
        return out

    out["corr_delta_conservation"] = float(
        np.corrcoef(np.abs(delta), conservation)[0, 1])

    # how well conservation alone does, so the model has a stated benchmark to beat
    from sklearn.metrics import roc_auc_score
    out["conservation_auroc"] = float(roc_auc_score(label, conservation))
    out["delta_auroc"] = float(roc_auc_score(label, np.abs(delta)))

    alone = fit_delta_coef(delta, label, None, n_boot, seed, method)
    ctrl = fit_delta_coef(delta, label, conservation, n_boot, seed, method)

    X = np.column_stack([_z(np.abs(delta)), _z(conservation)])
    try:
        c2, _ = _coef_se(X, label, method)
        out["conservation_coef"] = float(c2[1])
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        pass

    for name, f in (("alone", alone), ("controlled", ctrl)):
        out.update({f"{name}_{k}": v for k, v in asdict(f).items()})
        out[f"{name}_survives"] = bool(f.survives)
        out[f"{name}_significant"] = bool(f.significant)

    # attenuation: how much of the naive effect conservation explains away
    if alone.coef != 0:
        out["attenuation"] = round(1.0 - ctrl.coef / alone.coef, 4)
    return out


def run(df, model_cols, label_col="label", conservation_col="conservation",
        group_col=None, seed=0, n_boot=BOOT, method="firth", min_n=20):
    """Run the test for every (group, model) pair. `group_col=None` pools everything."""
    groups = [(None, df)] if group_col is None else list(df.groupby(group_col, sort=True))
    rows = []
    for key, g in groups:
        for m in model_cols:
            r = {"group": key if key is not None else "ALL", "model": m}
            r.update(test_group(g[m].to_numpy(), g[label_col].to_numpy(),
                                g[conservation_col].to_numpy(),
                                seed=seed, n_boot=n_boot, method=method, min_n=min_n))
            rows.append(r)
    return add_fdr(pd.DataFrame(rows))


# ---------------------------------------------------------------------------------------
# Multiple comparisons
# ---------------------------------------------------------------------------------------

def add_fdr(res, alpha=0.05):
    """Add BH q-values and an FDR-corrected survival flag."""
    res = res.copy()
    for name in ("alone", "controlled"):
        col = f"{name}_p_wald"
        if col not in res:
            continue
        res[f"{name}_q"] = benjamini_hochberg(res[col].values)
        res[f"{name}_survives_fdr"] = (
            (res[f"{name}_q"] < alpha) & (res[f"{name}_coef"] > 0)).fillna(False)
    return res


def summarise(res):
    """Survival counts per model: the headline table."""
    ok = res[res.note == ""]
    agg = {"groups": ("group", "size"),
           "survives_alone": ("alone_survives", "sum"),
           "survives_controlled": ("controlled_survives", "sum"),
           "median_corr": ("corr_delta_conservation", "median"),
           "median_attenuation": ("attenuation", "median"),
           "total_n": ("n", "sum"),
           "total_pathogenic": ("n_pathogenic", "sum")}
    if "controlled_survives_fdr" in ok:
        agg["survives_controlled_fdr"] = ("controlled_survives_fdr", "sum")
    return (ok.groupby("model").agg(**agg).reset_index()
              .sort_values("survives_controlled", ascending=False))


# ---------------------------------------------------------------------------------------
# Power
# ---------------------------------------------------------------------------------------

def min_detectable_effect(n, prevalence, conservation_corr=0.0, conservation_effect=1.0,
                          target_power=0.8, effects=None, n_sim=200, n_boot=400,
                          seed=0, method="firth"):
    """Smallest true delta coefficient detected at `target_power`, by simulation.

    This is what turns "we found nothing" into "we could have found effects above X and
    did not". Without it a null result is a statement about sample size, not about
    biology.

    Simulates conservation-correlated delta, injects a known effect, and measures how
    often the controlled test recovers it with an interval excluding zero.
    """
    if effects is None:
        effects = [0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    rng = np.random.default_rng(seed)
    r = float(np.clip(conservation_corr, -0.99, 0.99))
    prevalence = float(np.clip(prevalence, 0.01, 0.99))
    intercept = np.log(prevalence / (1 - prevalence))

    for beta in effects:
        hits = tried = 0
        for _ in range(n_sim):
            cons = rng.standard_normal(n)
            # Simulated delta must be NON-NEGATIVE, because the real delta is a
            # magnitude and `fit_delta_coef` standardises |delta|. An earlier version
            # drew a signed normal and injected the effect on the signed value, then
            # fitted on its absolute value -- which destroys the relationship for a
            # symmetric variable and made this function badly understate power.
            latent = r * cons + np.sqrt(max(0.0, 1 - r * r)) * rng.standard_normal(n)
            delta = np.abs(latent)
            # the effect enters through exactly the transform the fit uses, so `beta` is
            # directly the standardised coefficient we are trying to recover
            logit = (intercept + beta * _z(delta) + conservation_effect * _z(cons))
            y = rng.binomial(1, 1.0 / (1.0 + np.exp(-logit)))
            if len(np.unique(y)) < 2:
                continue
            tried += 1
            try:
                f = fit_delta_coef(delta, y, cons, n_boot=n_boot,
                                   seed=int(rng.integers(1 << 30)), method=method)
            except Exception:
                continue
            if f.ci_low > 0:
                hits += 1
        if tried and hits / tried >= target_power:
            return beta
    return float("inf")


def power_table(res, **kw):
    """Minimum detectable effect for every group in a results table."""
    rows = []
    for _, r in res[res.note == ""].iterrows():
        n, npath = int(r["n"]), int(r["n_pathogenic"])
        mde = min_detectable_effect(n, npath / n if n else 0.01,
                                    conservation_corr=float(
                                        r.get("corr_delta_conservation") or 0.0), **kw)
        rows.append({"group": r["group"], "model": r["model"], "n": n,
                     "n_pathogenic": npath, "min_detectable_effect": mde,
                     "controlled_coef": r.get("controlled_coef"),
                     "underpowered": bool(mde == float("inf")
                                          or abs(r.get("controlled_coef") or 0) < mde)})
    return pd.DataFrame(rows)
