"""The genuinely unconditional refit, and why "attenuation" needs a calibrated null.

WHAT WENT WRONG THAT THIS FIXES. `cloud_analysis.py:388` builds every row of
`variant_specificity_refit.csv` by calling `fit_delta_coef(..., u.conservation, ...)`. The
`conservation` argument is what adds the phyloP column, so the row the paper labels
"controls = none" was ALREADY adjusted for conservation. Its published companion, "0.21%
attenuation after conditioning on phyloP", is therefore the difference between two
conservation-adjusted fits on slightly different row counts (18,998 against 18,830): it
measures dropping 168 rows. `fit_delta_coef` has always accepted `conservation=None`; nothing
ever passed it.

WHAT THIS SCRIPT IS. The same estimator, done correctly, on the k-mer arm --
`variant_scores.csv`, the one variant table committed to the repo. It is a METHODS
DEMONSTRATION, not R4: R4's claim is about SpliceBERT, whose per-variant scores live only in
GCS. What transfers is the procedure and the interpretation, both of which were wrong.

Two things are enforced here that were not there:

  identical rows   Both fits run on ONE row set, built once. The attenuation cannot absorb a
                   change in n, which is what the retracted 0.21% actually was.

  a calibrated     This is the part that matters, and its absence is why ~0% attenuation read
  null             as a triumph. Logistic regression is NOT COLLAPSIBLE: adding a strong
                   predictor of the outcome inflates the other coefficients even when the two
                   predictors are perfectly independent, because the latent residual scale is
                   fixed. So the null for "phyloP carries no information about the model's
                   signal" is not 0% attenuation -- it is a large NEGATIVE attenuation, i.e.
                   amplification. Without that number, an observed 0% cannot be told apart
                   from strong sharing.

  The null is constructed rather than argued, by FORWARD simulation: draw a score and a
  covariate with a known correlation, generate labels from them with coefficients matched to
  the observed conditional fit, then fit with and without the covariate exactly as the real
  analysis does. At correlation zero this gives the non-collapsibility null; sweeping the
  correlation gives the value that would reproduce the observed attenuation, which can then be
  compared against the correlation actually measured.

  A FIRST ATTEMPT AT THIS NULL WAS WRONG AND IS RECORDED HERE SO IT IS NOT REPEATED. It drew
  the covariate FROM the label, C = mu*y + noise, on the grounds that this makes C independent
  of the score. It does not make it a valid null: C is then a descendant of the outcome, so
  conditioning on it is conditioning on a collider and the amplification largely vanishes. It
  reported a null of -1.5% where the standard approximation for a covariate this strong,
  1 - sqrt(1 + 0.346*c^2), gives about -60%. Getting the null wrong in the anti-conservative
  direction would have turned "the signal is substantially shared with conservation" into "the
  signal is essentially independent of it" -- the same false conclusion the retracted version
  reached, arrived at a different way.

  A SECOND CORRECTION, ALSO RECORDED. The null was first simulated with a NORMAL covariate,
  and the closed form 1 - sqrt(1 + 0.346*c^2) was gated as a cross-check that "must agree"
  with it. Both were wrong. Non-collapsibility depends on the omitted covariate's whole
  distribution, not just its variance, and phyloP is skewed (+1.04); the covariate is now drawn
  through a Gaussian copula onto phyloP's own empirical marginal. And the closed form is a
  small-sigma probit approximation which at c = 2.12 is far outside its range and
  ANTI-CONSERVATIVE, so gating agreement with it gated the wrong direction. Exact
  Gauss-Hermite marginalisation says the simulation is correct. It is now reported as an
  order-of-magnitude reference and nothing is asserted about the gap. Seeds raised 12 -> 40,
  since the earlier estimate had a seed-to-seed sd of about 0.054.

Reference for the estimator's limits, which applies to the original framing too: Pepe,
Janes, Longton, Leisenring & Newcomb 2004, Am J Epidemiol -- a large, tightly bounded odds
ratio implies very little about added discrimination. A coefficient was the wrong summary for
"does the model add anything" regardless of what it was conditioned on.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
TABLES = ROOT / "results" / "tables"

from rbp.stats import coef_se                                          # noqa: E402
from rbp.variants import conservation as cons                          # noqa: E402
from rbp.utils.log import log  # noqa: E402

N_BOOT = 2000
SEED = 0
BLOCK = 1_000_000
CAL_SEEDS = 40
CAL_RHO = np.linspace(0.0, 0.40, 21)



def z(a):
    a = np.asarray(a, dtype=float)
    s = a.std()
    return (a - a.mean()) / (s if s else 1.0)


def rows():
    """One row set, built once, used by every fit below."""
    sc = pd.read_csv(TABLES / "variant_scores.csv")
    cv = pd.read_csv(TABLES / "variant_conservation.csv")[["vid", "conservation"]]
    d = sc.merge(cv, on="vid", how="left").dropna(subset=["delta", "conservation"]).copy()
    d["dataset"] = d.protein + ":" + d.cell
    d["ad"] = d.delta.abs()
    g = d.groupby("dataset").ad
    d["ad_within"] = (d.ad - g.transform("mean")) / g.transform("std").replace(0, np.nan)
    return d.dropna(subset=["ad_within"])


def dedupe(d, col):
    """Deduplicate on the column being analysed, as the original does."""
    u = d.sort_values(col, key=abs, ascending=False).drop_duplicates("vid")
    return u.assign(block=u.vid.str.split(":").str[0] + "_" +
                    (u.vid.str.split(":").str[1].astype(int) // BLOCK).astype(str))


def simulate_attenuation(rho, b, c, prev, n, marginal=None, method="firth"):
    """Attenuation under a KNOWN correlation between score and covariate.

    Forward simulation: draw the score and the covariate with correlation `rho`, generate
    labels from both with coefficients b and c, and run the same two fits the real analysis
    runs. rho = 0 is the non-collapsibility null.

    THE COVARIATE USES phyloP's OWN MARGINAL, NOT A NORMAL ONE. Non-collapsibility depends on
    the whole distribution of the omitted covariate and not merely its variance, and real
    phyloP is skewed (+1.04, sd 3.19, range -20 to +10). Drawing normals understated the null
    by roughly ten points. Correlation is induced with a Gaussian copula so `rho` still means
    what it says while the marginal is empirical.
    """
    a = np.log(prev / (1.0 - prev))
    vals = []
    for s in range(CAL_SEEDS):
        rng = np.random.default_rng(2000 + s)
        x = rng.standard_normal(n)
        u = rho * x + np.sqrt(max(1.0 - rho * rho, 0.0)) * rng.standard_normal(n)
        if marginal is None:
            cc = u
        else:                                   # copula: rank-map onto the empirical marginal
            cc = np.quantile(marginal, np.clip(norm.cdf(u), 1e-9, 1 - 1e-9))
        p = 1.0 / (1.0 + np.exp(-(a + b * z(x) + c * z(cc))))
        y = (rng.random(n) < p).astype(int)
        if len(np.unique(y)) < 2:
            continue
        bm, _ = coef_se(z(x)[:, None], y, method)
        bc, _ = coef_se(np.column_stack([z(x), z(cc)]), y, method)
        if float(bm[0]) != 0.0:
            vals.append(1.0 - float(bc[0]) / float(bm[0]))
    return float(np.mean(vals)) if vals else np.nan


def calibrate(observed_att, b, c, prev, n, marginal):
    """The null, an order-of-magnitude reference, and the rho reproducing what we observed."""
    curve = np.array([simulate_attenuation(r, b, c, prev, n, marginal) for r in CAL_RHO])
    null = float(curve[0])
    # NOT A CROSS-CHECK. 1 - sqrt(1 + 0.346 c^2) is a small-sigma probit approximation to the
    # logistic-normal attenuation factor. At c = 2.12 it is far outside the range where it
    # holds and it is ANTI-CONSERVATIVE: it understates the amplification, which biases the
    # excess-over-null upward, i.e. in the direction that flatters the result. An earlier
    # version of this script asserted that the two "must agree" to within 0.15 and gated it.
    # Exact Gauss-Hermite marginalisation says the simulation is the correct one. It is kept
    # only as an order-of-magnitude reference, and it is no longer gated as agreement.
    analytic = 1.0 - np.sqrt(1.0 + 0.346 * c * c)
    o = np.argsort(curve)
    implied = float(np.interp(observed_att, curve[o], CAL_RHO[o]))
    return null, float(analytic), implied


def main():
    d = rows()
    out = []

    for tag, col in (("pooled", "ad"), ("within_dataset", "ad_within")):
        u = dedupe(d, col)
        x, y = u[col].to_numpy(), u.label.to_numpy()
        cv, blk = u.conservation.to_numpy(), u.block.to_numpy()
        ta = (col == "ad")

        # THE TWO FITS, on identical rows. conservation=None is the whole correction.
        unc = cons.fit_delta_coef(x, y, None, n_boot=N_BOOT, seed=SEED, blocks=blk,
                                  take_abs=ta)
        con = cons.fit_delta_coef(x, y, cv, n_boot=N_BOOT, seed=SEED, blocks=blk,
                                  take_abs=ta)
        att = 1.0 - con.coef / unc.coef

        # phyloP's own coefficient, and the score's, are the structural parameters the
        # simulation is calibrated to
        xs = np.abs(x) if ta else x
        b, _ = coef_se(np.column_stack([z(xs), z(cv)]), y, "firth")
        phylop_coef = float(b[1])
        null_att, analytic, implied_rho = calibrate(att, float(b[0]), phylop_coef,
                                                    float(y.mean()), len(y), cv)

        rho = spearmanr(xs, cv)

        for check, val, f in (("coef, TRULY unconditional", unc.coef, unc),
                              ("coef, conditional on phyloP", con.coef, con)):
            out.append({"standardisation": tag, "check": check, "value": val,
                        "ci_low": f.ci_low, "ci_high": f.ci_high, "n": len(u), "note": ""})
        out += [
            {"standardisation": tag, "check": "attenuation fraction (identical rows)",
             "value": att, "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
             "note": "positive = conservation explains part of it"},
            {"standardisation": tag, "check": "phyloP coefficient in the joint fit",
             "value": phylop_coef, "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
             "note": "the strength the null is calibrated to"},
            {"standardisation": tag,
             "check": "NULL attenuation at rho=0 (simulated)",
             "value": null_att, "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
             "note": "non-collapsibility; this is what 'no sharing' predicts"},
            {"standardisation": tag,
             "check": "NULL attenuation at rho=0 (analytic reference only)",
             "value": analytic, "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
             "note": "small-sigma approximation, ANTI-CONSERVATIVE at this c; not a check"},
            {"standardisation": tag, "check": "excess attenuation over the null",
             "value": att - null_att, "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
             "note": "THE INTERPRETABLE QUANTITY: >0 means signal IS shared"},
            {"standardisation": tag,
             "check": "correlation implied by the observed attenuation",
             "value": implied_rho, "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
             "note": "what rho would have to be to produce what we saw"},
            {"standardisation": tag, "check": "spearman(|delta|, phyloP), MEASURED",
             "value": float(rho.statistic), "ci_low": np.nan, "ci_high": np.nan,
             "n": len(u), "note": f"p={rho.pvalue:.3g}; compare with the implied value"},
        ]

    res = pd.DataFrame(out)
    res.to_csv(TABLES / "unconditional_refit.csv", index=False)
    for tag in res.standardisation.unique():
        log(f"\n  [{tag}]")
        for _, x in res[res.standardisation == tag].iterrows():
            ci = f" [{x.ci_low:+.3f}, {x.ci_high:+.3f}]" if pd.notna(x.ci_low) else ""
            log(f"    {x.check:56} {x.value:+.4f}{ci}")
    log(f"\n  wrote {TABLES / 'unconditional_refit.csv'}")


if __name__ == "__main__":
    main()
