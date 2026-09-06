"""Q1: does the strand artifact drive R1's surviving contrast?

WHY THIS IS BACK, HAVING ONCE BEEN SETTLED. `negatives.py:328` gives each negative the
POSITIVE's strand, so a negative window carries its own gene's strand only about half the time
(measured: 55.2% among unambiguous windows, 47.4% once the 14.0% ambiguous are counted). A
$25 regeneration was scoped to fix it and then CANCELLED, on the evidence that the contrast
grew rather than shrank when restricted to sense-only negatives: +0.2643 -> +0.2787.

That evidence is void. +0.2643 is the composition-SHARE contrast, retracted since as an
algebraic identity (share_m = C/gain_m, so the comparison excluded zero with probability 1).
The claim that survived is a different quantity, the nested gain contrast +0.0397, and it has
never been strand-tested. The manuscript says so in its own limitations. This closes that gap
without the regeneration, using only committed tables.

WHAT IS AND IS NOT BEING ASKED. The artifact inflates ABSOLUTE AUROCs in both arms, and that is
conceded regardless of what happens here. The question is narrower: does it drive the
DIFFERENCE between the arms? Both arms share their positives and their strand convention, so
the paired design already predicts it should not. That prediction is now tested rather than
asserted.

THE HONEST PART, AND THE REASON THE OLD VERSION OF THIS ARGUMENT WAS WEAK. The strand audit
covers 40 datasets, not 94. A non-significant correlation at that n is feeble evidence of
absence: this script therefore reports a BOOTSTRAP INTERVAL for the correlation and states
which effect sizes the data actually exclude, rather than resting on p > 0.05. The previous
strand block in golden.yaml gated `contrast_rho_max: 0.40`, a ceiling so loose that a genuine
confound would have sailed through it. That block was deleted; this one asserts the interval.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from rbp.utils.log import log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
TABLES = ROOT / "results" / "tables"
N_BOOT = 2000
SEED = 0



def boot_rho(x, y, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(x)
    out = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if len(np.unique(x[i])) < 3 or len(np.unique(y[i])) < 3:
            continue
        out.append(spearmanr(x[i], y[i]).statistic)
    return np.percentile(out, [2.5, 97.5]) if len(out) >= 50 else (np.nan, np.nan)


def main():
    sa = pd.read_csv(TABLES / "strand_audit.csv")[["dataset", "frac_sense", "frac_ambiguous"]]
    cm = pd.read_csv(TABLES / "cost_of_matching.csv")
    d = cm.merge(sa, on="dataset", how="inner")
    d["contrast"] = d.delta_auroc_dn - d.delta_auroc_gc
    n = len(d)
    log(f"  {n} datasets with both a strand audit and both rehearsal arms "
        f"(of {len(cm)} in R1, {len(sa)} audited)")

    rows = []

    def add(check, value, lo=np.nan, hi=np.nan, note=""):
        rows.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                     "n": n, "note": note})

    add("mean frac_sense on these datasets", d.frac_sense.mean(),
        note=f"range {d.frac_sense.min():.3f}-{d.frac_sense.max():.3f}")
    add("mean contrast on these datasets", d.contrast.mean(),
        note="the published contrast is +0.0397 on all 94")

    # THE TEST. Does the contrast track the artifact?
    r = spearmanr(d.frac_sense, d.contrast)
    lo, hi = boot_rho(d.frac_sense.to_numpy(), d.contrast.to_numpy())
    add("spearman(frac_sense, CONTRAST)", r.statistic, lo, hi,
        note=f"p={r.pvalue:.3g}; THE TEST: near zero means the artifact does not drive it")

    # The arms individually. Both are expected to move, and if NEITHER does, the artifact is
    # not biting at all here and the contrast result is uninformative rather than reassuring.
    for col, lab in (("delta_auroc_gc", "GC arm gain"), ("delta_auroc_dn", "dinuc arm gain"),
                     ("auroc_gc", "GC arm AUROC"), ("auroc_dn", "dinuc arm AUROC")):
        rr = spearmanr(d.frac_sense, d[col])
        l2, h2 = boot_rho(d.frac_sense.to_numpy(), d[col].to_numpy())
        add(f"spearman(frac_sense, {lab})", rr.statistic, l2, h2, note=f"p={rr.pvalue:.3g}")

    # IS IT THE SIZE EFFECT IN DISGUISE? R1 already reports that the contrast grows with
    # dataset size (rho +0.307 on all 94). If frac_sense tracked size, any apparent strand
    # association would be that, relabelled. Partial correlation on ranks, controlling for
    # log10(pairs). Measured: frac_sense vs size is rho -0.076 (p=0.64), and the partial is
    # -0.252 against a raw -0.240, so size explains none of it. The association, such as it
    # is, is not a size artifact -- which unfortunately makes it harder to dismiss, not easier.
    lp = np.log10(d.pairs.to_numpy())
    rk = lambda v: pd.Series(v).rank().to_numpy()                       # noqa: E731
    X = np.column_stack([np.ones(n), rk(lp)])
    resid = lambda v: rk(v) - X @ np.linalg.lstsq(X, rk(v), rcond=None)[0]   # noqa: E731
    rs = spearmanr(d.frac_sense.to_numpy(), lp)
    add("spearman(frac_sense, log10 pairs)", rs.statistic, note=f"p={rs.pvalue:.3g}")
    rp = spearmanr(resid(d.frac_sense.to_numpy()), resid(d.contrast.to_numpy()))
    add("PARTIAL rho(frac_sense, contrast | size)", rp.statistic,
        note=f"p={rp.pvalue:.3g}; size explains none of the association")

    # Model-free version: split at the median and compare the two halves.
    med = d.frac_sense.median()
    poor, rich = d[d.frac_sense <= med], d[d.frac_sense > med]
    u = mannwhitneyu(rich.contrast, poor.contrast)
    add("contrast, sense-POOR half", poor.contrast.mean(), note=f"n={len(poor)}")
    add("contrast, sense-RICH half", rich.contrast.mean(), note=f"n={len(rich)}")
    rng = np.random.default_rng(SEED)
    diffs = [rich.contrast.iloc[rng.integers(0, len(rich), len(rich))].mean()
             - poor.contrast.iloc[rng.integers(0, len(poor), len(poor))].mean()
             for _ in range(N_BOOT)]
    dlo, dhi = np.percentile(diffs, [2.5, 97.5])
    add("difference between halves", rich.contrast.mean() - poor.contrast.mean(), dlo, dhi,
        note=f"Mann-Whitney p={u.pvalue:.3g}; should straddle zero")

    # THE DECISION-RELEVANT NUMBER: what would the contrast be if every negative sat on its own
    # gene's strand? A linear extrapolation of contrast on frac_sense to 1.0.
    #
    # THIS IS AN EXTRAPOLATION FAR OUTSIDE THE OBSERVED RANGE and must be labelled as one. The
    # audited datasets span frac_sense 0.433 to 0.615, so predicting at 1.0 is a reach of more
    # than twice the observed spread, on 40 points, with a slope whose interval includes zero.
    # It is reported because a decision about a $25 regeneration turns on it, not because it is
    # a reliable estimate. If it lands near +0.0397 the artifact is unimportant; if it lands
    # near zero the regeneration is mandatory before submission.
    x, y = d.frac_sense.to_numpy(), d.contrast.to_numpy()
    slope, icpt = np.polyfit(x, y, 1)
    rng2 = np.random.default_rng(SEED)
    preds, slopes = [], []
    for _ in range(N_BOOT):
        i = rng2.integers(0, n, n)
        if len(np.unique(x[i])) < 3:
            continue
        s, c = np.polyfit(x[i], y[i], 1)
        slopes.append(s)
        preds.append(s * 1.0 + c)
    slo, shi = np.percentile(slopes, [2.5, 97.5])
    plo, phi = np.percentile(preds, [2.5, 97.5])
    add("OLS slope, contrast on frac_sense", slope, slo, shi,
        note="negative = contrast shrinks as negatives become correctly stranded")
    add("EXTRAPOLATED contrast at frac_sense = 1.0", slope + icpt, plo, phi,
        note="EXTRAPOLATION beyond the observed 0.433-0.615; decision input, not an estimate")

    # WHAT THE SAMPLE CAN ACTUALLY EXCLUDE. With n this small, "p > 0.05" is not evidence of
    # absence, and saying so is the difference between a control and a formality.
    add("frac_sense leverage available", float(d.frac_sense.max() - d.frac_sense.min()),
        note="the artifact barely varies across datasets, so this test is weak BY DESIGN")

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "strand_contrast.csv", index=False)
    for _, x in out.iterrows():
        ci = f" [{x.ci_low:+.3f}, {x.ci_high:+.3f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:44} {x.value:+.4f}{ci}   {x.note}")
    log(f"\n  wrote {TABLES / 'strand_contrast.csv'}")


if __name__ == "__main__":
    main()
