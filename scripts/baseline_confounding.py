"""R1l: the protocol and the baseline it induces are not separable, and that IS the finding.

    python scripts/baseline_confounding.py

THE OBJECTION THIS ANSWERS, which is the strongest one left. R1k concedes that the nested
contribution tracks the composition baseline at Spearman -0.60. A referee will say: then "the
protocol changes the measured contribution" is just "the protocol changes the baseline, and the
baseline changes the contribution", which is trivially true and not worth a paper.

THE TEST, RUN NAIVELY, APPEARS TO AGREE WITH THE REFEREE. Regress gain on a flexible function
of the composition baseline, pooled over all 282 protocol-dataset cells, then add protocol
dummies: R2 rises 0.3955 -> 0.4065, F(2,276) = 2.6, p = 0.08. With dataset fixed effects,
0.8243 -> 0.8301, p = 0.045. Knowing which protocol produced a baseline adds almost nothing to
knowing the baseline.

WHY THAT TEST IS UNINFORMATIVE, AND THIS IS THE POINT. The three protocols do not overlap in
baseline. Their 10th-90th percentile ranges intersect in a window 0.005 AUROC wide containing
one cell per arm. There is essentially no region of the covariate space where two protocols can
be compared at the same baseline, so the regression cannot separate them ON THIS DATA and
neither can anything else. Protocol and baseline are not two variables that happen to be
correlated; the protocol IS the operation that sets the baseline.

WHAT SURVIVES, AND IT IS NOT NOTHING. If the whole pattern were compression, transplanting a
model's own d' increment across baselines would reproduce the other arm's gain exactly. It does
not: the transplant explains 21% of the gc->dn difference, 40% of gc->neg2 and 52% of dn->neg2.
And a single constant d' increment applied to each cell's own baseline predicts the 282 gains at
only R2 = 0.18, with residuals that are systematically ordered by protocol (dn +0.0247,
gc -0.0038, neg2 -0.0143). So the arms differ by more than the arithmetic of the ceiling -- but
how much more cannot be estimated, because the estimate would require overlap that does not
exist.

THE CONCLUSION THE PAPER SHOULD DRAW. "Is it the protocol or the baseline?" is a malformed
question, and the practical consequence is exactly the paper's recommendation: report the
composition baseline measured under the same protocol, because it is the only summary of what
the protocol did that a reader can act on, and never compare contributions across protocols
because there is no common scale on which to do so.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f as fdist
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
ARMS = ("gc", "dn", "neg2")
R2C = np.sqrt(2.0)


def dprime(a):
    return R2C * norm.ppf(np.clip(a, 1e-6, 1.0 - 1e-6))


def auroc(z):
    return norm.cdf(z / R2C)


def ols(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    return b, float(r @ r)


def main():
    d = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    long = pd.concat([pd.DataFrame({"ds": d.dataset, "arm": a, "comp": d[f"comp_{a}"],
                                    "gain": d[f"gain_{a}"]}) for a in ARMS],
                     ignore_index=True)
    rows = []

    print("=== 1. Does protocol add anything beyond the baseline? (the naive test) ===")
    c, y = long.comp.values, long.gain.values
    Xb = np.column_stack([np.ones(len(c)), c, c ** 2, c ** 3])
    _, ssb = ols(Xb, y)
    tss = float((y - y.mean()) @ (y - y.mean()))
    D = pd.get_dummies(long.arm, drop_first=True).values.astype(float)
    _, ssa = ols(np.column_stack([Xb, D]), y)
    df1, df2 = D.shape[1], len(y) - Xb.shape[1] - D.shape[1]
    F = ((ssb - ssa) / df1) / (ssa / df2)
    p = float(1 - fdist.cdf(F, df1, df2))
    # dataset fixed effects too, so every number quoted in the manuscript has a source
    DS = pd.get_dummies(long.ds, drop_first=True).values.astype(float)
    _, ssf = ols(np.column_stack([Xb, DS]), y)
    _, ssfa = ols(np.column_stack([Xb, DS, D]), y)
    rows += [{"check": "R2, baseline plus dataset fixed effects", "value": 1 - ssf / tss},
             {"check": "R2, baseline plus dataset FE plus protocol", "value": 1 - ssfa / tss},
             {"check": "R2, cubic in composition baseline alone", "value": 1 - ssb / tss},
             {"check": "R2, plus protocol dummies", "value": 1 - ssa / tss},
             {"check": "F statistic, protocol beyond baseline", "value": float(F)},
             {"check": "p value, protocol beyond baseline", "value": p}]
    print(f"  R2 baseline alone {1 - ssb / tss:.4f} -> with protocol {1 - ssa / tss:.4f}"
          f"   F({df1},{df2}) = {F:.2f}, p = {p:.3f}")

    print("\n=== 2. Why that test cannot answer the question: there is no overlap ===")
    q = {a: long.comp[long.arm == a] for a in ARMS}
    for a in ARMS:
        print(f"  {a:5s} median {q[a].median():.3f}  10th-90th "
              f"{q[a].quantile(.10):.3f}-{q[a].quantile(.90):.3f}")
        rows.append({"check": f"composition baseline median, {a} arm",
                     "value": float(q[a].median())})
        rows.append({"check": f"composition baseline 10th percentile, {a} arm",
                     "value": float(q[a].quantile(.10))})
        rows.append({"check": f"composition baseline 90th percentile, {a} arm",
                     "value": float(q[a].quantile(.90))})
    lo = max(q[a].quantile(.10) for a in ARMS)
    hi = min(q[a].quantile(.90) for a in ARMS)
    width = max(hi - lo, 0.0)
    n_common = int(((long.comp >= lo) & (long.comp <= hi)).sum())
    rows += [{"check": "common support width (AUROC)", "value": width},
             {"check": "cells inside the common support", "value": n_common}]
    print(f"  common support {lo:.3f}-{hi:.3f}, width {width:.4f} AUROC, "
          f"{n_common} of {len(long)} cells")
    print("  -> protocol and baseline are confounded BY CONSTRUCTION, not by accident")

    print("\n=== 3. What survives: compression does NOT explain the differences ===")
    for src, tgt in (("gc", "dn"), ("gc", "neg2"), ("dn", "neg2")):
        inc = dprime(d[f"full_{src}"]) - dprime(d[f"comp_{src}"])
        pred = (auroc(dprime(d[f"comp_{tgt}"]) + inc) - d[f"comp_{tgt}"]).mean()
        obs = d[f"gain_{tgt}"].mean() - d[f"gain_{src}"].mean()
        frac = (pred - d[f"gain_{src}"].mean()) / obs if obs else np.nan
        rows.append({"check": f"fraction of {src}->{tgt} difference explained by compression",
                     "value": float(frac)})
        print(f"  {src:4s} -> {tgt:4s}  observed {obs:+.4f}   compression explains {frac:5.0%}")

    kbar = float(np.mean(np.concatenate(
        [dprime(d[f"full_{a}"].values) - dprime(d[f"comp_{a}"].values) for a in ARMS])))
    pred = auroc(dprime(long.comp) + kbar) - long.comp
    r2_one = 1 - float((long.gain - pred) @ (long.gain - pred)) / tss
    rows += [{"check": "mean d' increment across all cells", "value": kbar},
             {"check": "R2, one constant increment on each cell's own baseline",
              "value": r2_one}]
    print(f"\n  one constant d' increment ({kbar:.3f}) + each cell's own baseline: "
          f"R2 = {r2_one:.4f}")
    for a in ARMS:
        m = (long.arm == a).values
        resid = float((long.gain[m] - pred[m]).mean())
        rows.append({"check": f"residual after the constant-increment model, {a} arm",
                     "value": resid})
        print(f"    {a:5s} residual {resid:+.4f}")

    # THE SPEARMAN IS PARTLY BETWEEN-ARM, and a statistician will decompose it, so do it here.
    from scipy.stats import spearmanr
    rho_all, _ = spearmanr(long.comp, long.gain)
    cen = long.copy()
    cen["comp"] = cen.comp - cen.groupby("arm").comp.transform("mean")
    cen["gain"] = cen.gain - cen.groupby("arm").gain.transform("mean")
    rho_c, _ = spearmanr(cen.comp, cen.gain)
    rows += [{"check": "spearman(baseline, gain), pooled", "value": float(rho_all)},
             {"check": "spearman(baseline, gain), arm-centred", "value": float(rho_c)}]
    print(f"\n  Spearman pooled {rho_all:+.3f}   arm-centred {rho_c:+.3f}")
    for a in ARMS:
        m = long.arm == a
        r, pv = spearmanr(long.comp[m], long.gain[m])
        rows.append({"check": f"spearman(baseline, gain), within {a} arm", "value": float(r),
                     "note": f"p={pv:.3f}"})
        print(f"    within {a:5s} {r:+.3f}  p={pv:.3f}")

    pd.DataFrame(rows).to_csv(TABLES / "baseline_confounding.csv", index=False)
    print("\nwrote baseline_confounding.csv")


if __name__ == "__main__":
    main()
