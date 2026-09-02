"""R1s: whose property is the protocol multiplier? The protein's, and not by as much as claimed.

    python scripts/multiplier_variance.py

THE CLAIM THIS EXISTS TO CHECK. README.md and docs/65 assert "the multiplier is a property of
the PROTEIN (68.9% of variance) not the model (1.5%)". That decomposition was run by a referee
agent during review and committed to no script and no table: the paper's most quotable
one-liner about R1g was, until now, sourced to a conversation. Either it reproduces and gets
gated, or it comes out of the README.

WHAT IS DECOMPOSED. One cell per (dataset, model): the LOG multiplier log(gain_dn / gain_gc),
defined only where both arms are positive, because a ratio across zero is not a multiplier.
Logs because the quantity is a ratio and its arithmetic mean is not the thing anyone means.

WHY A PERMUTATION NULL IS NOT OPTIONAL HERE. Protein has 79 levels on ~230 cells. Model class
has 3 and cell line has 2. A factor with 79 levels absorbs a large share of any variance,
including variance in pure noise, so "protein explains 68.9%" and "model explains 1.5%" are not
comparable numbers as stated -- most of the gap is degrees of freedom. Each factor's share is
therefore reported next to the share the SAME factor gets on relabelled data, and the honest
statistic is the excess over that null. This is the project's own lesson: a control is not
evidence until something has tried to break it.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
MODELS = ("kmer", "cnn", "splicebert")
N_PERM = 2000


def shares(y, factors):
    """Each factor's share of total sum of squares, from an OLS with all factors present.

    Type-II style: the drop in residual SS when a factor is added to the model containing the
    others, over the total SS. Shares need not sum to 1 under collinearity, and the residual is
    reported so the reader can see how much they miss by.
    """
    tss = float(((y - y.mean()) ** 2).sum())

    def rss(keys):
        if not keys:
            return tss
        X = np.column_stack([np.ones(len(y))] + [
            pd.get_dummies(factors[k], drop_first=True).to_numpy(dtype=float) for k in keys])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ b
        return float(r @ r)

    allk = list(factors)
    full = rss(allk)
    out = {k: (rss([j for j in allk if j != k]) - full) / tss for k in allk}
    out["residual"] = full / tss
    return out


def main():
    d = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv")
    rows = []
    for _, r in d.iterrows():
        for m in MODELS:
            g, h = r[f"{m}_gain_gc"], r[f"{m}_gain_dn"]
            if g > 0 and h > 0:
                rows.append({"protein": r.protein, "cell": r.cell, "model": m,
                             "log_mult": float(np.log(h / g))})
    t = pd.DataFrame(rows)
    y = t.log_mult.to_numpy()
    factors = {"protein": t.protein, "cell": t.cell, "model": t.model}

    print(f"n = {len(t)} (dataset x model) cells, both arms positive, "
          f"{t.protein.nunique()} proteins")
    print(f"panel multiplier, exp(mean log) = {np.exp(y.mean()):.2f}x   "
          f"ratio of means would be a different number and is not this one\n")

    obs = shares(y, factors)
    rng = np.random.default_rng(0)
    null = {k: [] for k in factors}
    for _ in range(N_PERM):
        p = rng.permutation(len(y))
        null_shares = shares(y[p], factors)
        for k in factors:
            null[k].append(null_shares[k])

    out = []
    print(f"  {'factor':10s} {'levels':>7s} {'share':>8s} {'null (relabelled)':>20s} "
          f"{'excess':>9s}")
    for k in factors:
        nl = np.array(null[k])
        exc = obs[k] - nl.mean()
        p = float((nl >= obs[k]).mean())
        out += [{"check": f"share of log-multiplier variance, {k}", "value": obs[k],
                 "n": len(t), "note": f"{t[k].nunique()} levels"},
                {"check": f"permutation null share, {k}", "value": float(nl.mean()),
                 "ci_low": float(np.percentile(nl, 2.5)),
                 "ci_high": float(np.percentile(nl, 97.5)), "n": len(t),
                 "note": f"p = {p:.4f}"},
                {"check": f"excess over the permutation null, {k}", "value": float(exc),
                 "n": len(t), "note": f"p = {p:.4f}"}]
        print(f"  {k:10s} {t[k].nunique():7d} {100 * obs[k]:7.1f}% "
              f"{100 * nl.mean():13.1f}% [{100 * np.percentile(nl, 97.5):.1f}] "
              f"{100 * exc:8.1f}%  p={p:.4f}")
    out.append({"check": "residual share", "value": obs["residual"], "n": len(t)})
    out.append({"check": "panel multiplier, exp(mean log)", "value": float(np.exp(y.mean())),
                "n": len(t)})
    print(f"  {'residual':10s} {'':7s} {100 * obs['residual']:7.1f}%")

    # THE CLAIM, IN THE FORM IT CAN BE MADE. Not "protein explains 46x more than model class"
    # -- that comparison is mostly degrees of freedom -- but that protein beats its own null
    # and model class does not beat its own.
    pe = obs["protein"] - np.mean(null["protein"])
    me = obs["model"] - np.mean(null["model"])
    out.append({"check": "protein excess minus model-class excess", "value": float(pe - me),
                "n": len(t)})
    print(f"\n  protein excess {100 * pe:+.1f}%  vs model-class excess {100 * me:+.1f}%")
    print("  -> the comparable statistic is the excess over each factor's OWN null,")
    print("     because 79 levels absorb variance that 3 levels cannot")

    pd.DataFrame(out).to_csv(TABLES / "multiplier_variance.csv", index=False)
    print("\nwrote multiplier_variance.csv")


if __name__ == "__main__":
    main()
