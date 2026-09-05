"""R1s: whose property is the protocol multiplier? The protein's, and not by as much as claimed.

    python scripts/multiplier_variance.py

PURPOSE. The claim that the protocol multiplier is a property of the protein rather than of the
model rests on a variance decomposition. This computes it, so that the claim is sourced to code
rather than asserted, and gates the result.

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

AND WHY ONE NULL IS NOT ENOUGH, which an adversarial statistician caught in the first version
of this script. Permuting y wholesale destroys the protein effect AND the within-dataset
correlation among the two or three model cells that share the same windows, labels and folds.
Because 79 proteins sit on 94 datasets, the protein factor is very nearly the DATASET factor
(using dataset as the factor gives 68.0%, slightly MORE than protein's 64.8%), so the wholesale
null is far too permissive and the "excess" it reports is inflated about five-fold.

The right null for the question "is this about the PROTEIN or about the individual experiment?"
permutes protein labels BETWEEN datasets while keeping each (dataset x model) block intact. On
that null the excess falls from +35.1 to +7.3 points, still significant. Both are reported,
because they answer different questions and quoting only the first was the error.

The strongest available evidence is neither: it is that the SAME protein's log multiplier
agrees across the two cell lines (r ~ 0.59 over the 15 proteins assayed in both), which is a
direct test and is what the section should lead with.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

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
                rows.append({"dataset": r.dataset, "protein": r.protein, "cell": r.cell,
                             "model": m,
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

    # THE SECOND NULL, and it is the one that answers the question actually being asked.
    # Wholesale permutation of y breaks the protein effect AND the within-dataset correlation
    # among the model cells that share windows, labels and folds. Since 79 proteins sit on 94
    # datasets, protein is very nearly the DATASET factor -- so that null is far too permissive.
    # This one permutes which PROTEIN each dataset belongs to, preserving the group-size
    # pattern and leaving every (dataset x model) block intact, so it isolates protein identity
    # from experiment identity.
    ds_of = t.dataset.to_numpy()
    uniq_ds = list(pd.unique(ds_of))
    prot_of_ds = t.drop_duplicates("dataset").set_index("dataset").protein.to_dict()
    sizes = pd.Series(prot_of_ds).value_counts().to_numpy()  # e.g. fifteen 2s, then 1s
    null_block = []
    for _ in range(N_PERM):
        order = list(rng.permutation(len(uniq_ds)))
        lab, at = {}, 0
        for j, sz in enumerate(sizes):
            for _ in range(sz):
                lab[uniq_ds[order[at]]] = f"P{j}"
                at += 1
        permuted = pd.Series([lab[x] for x in ds_of], index=t.index)
        null_block.append(shares(y, {"protein": permuted, "cell": t.cell,
                                     "model": t.model})["protein"])
    null_block = np.array(null_block)

    # And the DATASET factor itself, to show how little separates it from protein.
    ds_share = shares(y, {"dataset": t.dataset, "model": t.model})["dataset"]

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

    # PROTEIN IS NEARLY THE DATASET FACTOR, and the excess over the wholesale null was
    # therefore about five-fold too generous. Both numbers are reported.
    pb = obs["protein"] - null_block.mean()
    p_block = float((null_block >= obs["protein"]).mean())
    out += [{"check": "share of log-multiplier variance, dataset", "value": float(ds_share),
             "n": len(t), "note": "94 levels; protein's 79 is nearly this factor"},
            {"check": "block-preserving null share, protein", "value": float(null_block.mean()),
             "ci_low": float(np.percentile(null_block, 2.5)),
             "ci_high": float(np.percentile(null_block, 97.5)), "n": len(t),
             "note": f"p = {p_block:.4f}; permutes protein labels between datasets"},
            {"check": "excess over the block-preserving null, protein", "value": float(pb),
             "n": len(t), "note": f"p = {p_block:.4f}"}]
    print(f"\n  dataset as the factor: {100 * ds_share:.1f}%  "
          f"vs protein's {100 * obs['protein']:.1f}%")
    print("  -> protein (79 levels) is nearly the dataset factor (94 levels), so the")
    print("     wholesale null is too permissive. Against a null that permutes protein")
    print("     labels BETWEEN datasets and keeps each (dataset x model) block intact:")
    print(f"     null {100 * null_block.mean():.1f}% [{100 * np.percentile(null_block, 2.5):.1f}, "
          f"{100 * np.percentile(null_block, 97.5):.1f}]   "
          f"excess {100 * pb:+.1f}%  p={p_block:.4f}")

    # THE DIRECT TEST, which is stronger than any share: does the SAME protein get the same
    # multiplier in the other cell line? Fifteen proteins are assayed in both.
    both = t.pivot_table(index=["protein", "model"], columns="cell", values="log_mult")
    both = both.dropna()
    if len(both) >= 6:
        cells = list(both.columns)
        r = float(np.corrcoef(both[cells[0]], both[cells[1]])[0, 1])
        rho, pv = spearmanr(both[cells[0]], both[cells[1]])
        out += [{"check": "cross-cell-line correlation of the log multiplier", "value": r,
                 "n": len(both), "note": f"{both.index.get_level_values(0).nunique()} proteins "
                                         f"in both lines, x "
                                         f"{both.index.get_level_values(1).nunique()} models"},
                {"check": "cross-cell-line spearman of the log multiplier", "value": float(rho),
                 "n": len(both), "note": f"p = {pv:.4f}"}]

        # AND THE SAME TEST AT THE RIGHT CLUSTER LEVEL. Those rows are fifteen proteins by
        # three models sharing windows, labels and folds, so the p-value over 40 rows is
        # anti-conservative and contradicts this project's own rule of resampling proteins.
        # Collapsing over models costs an order of magnitude in p and is reported as primary.
        col = both.groupby(level=0).mean()
        cr, cp = pearsonr(col[cells[0]], col[cells[1]])
        csr, csp = spearmanr(col[cells[0]], col[cells[1]])

        # B13. AN INTERVAL, BECAUSE n IS 15. A correlation quoted with a p and no interval at
        # this sample size invites the reader to treat the point estimate as the finding. The
        # Fisher z interval is wide enough to make the honest reading obvious, and a
        # bias-corrected bootstrap over the same fifteen proteins is reported beside it so the
        # width is not an artefact of the normal approximation at small n.
        z = np.arctanh(cr)
        se = 1.0 / np.sqrt(len(col) - 3)
        lo_f, hi_f = np.tanh(z - 1.959964 * se), np.tanh(z + 1.959964 * se)
        rng = np.random.default_rng(7)
        a1, a2 = col[cells[0]].to_numpy(), col[cells[1]].to_numpy()
        draws = []
        for _ in range(4000):
            i = rng.integers(0, len(a1), len(a1))
            if np.std(a1[i]) > 0 and np.std(a2[i]) > 0:
                draws.append(float(np.corrcoef(a1[i], a2[i])[0, 1]))
        lo_b, hi_b = np.percentile(draws, [2.5, 97.5])
        out += [{"check": "cross-cell-line correlation, collapsed over models",
                 "value": float(cr), "n": len(col), "note": f"p = {cp:.4f}",
                 "ci_low": float(lo_f), "ci_high": float(hi_f)},
                {"check": "cross-cell-line correlation, collapsed, bootstrap CI",
                 "value": float(np.median(draws)), "n": len(col),
                 "note": f"{len(draws)} protein draws",
                 "ci_low": float(lo_b), "ci_high": float(hi_b)},
                {"check": "cross-cell-line spearman, collapsed over models",
                 "value": float(csr), "n": len(col), "note": f"p = {csp:.4f}"}]
        print(f"    Fisher z 95% CI [{lo_f:+.3f}, {hi_f:+.3f}]   "
              f"bootstrap [{lo_b:+.3f}, {hi_b:+.3f}] over {len(col)} proteins")
        print(f"    collapsed over models: r = {cr:+.3f} (p={cp:.4f})  "
              f"spearman {csr:+.3f} (p={csp:.4f})  over {len(col)} proteins")
        for m in both.index.get_level_values(1).unique():
            s = both.xs(m, level=1)
            if len(s) < 4:
                continue
            mr, mp = pearsonr(s[cells[0]], s[cells[1]])
            out.append({"check": f"cross-cell-line correlation, {m} only", "value": float(mr),
                        "n": len(s), "note": f"p = {mp:.4f}"})
            print(f"      {m:11s} r = {mr:+.3f} (p={mp:.3f}) over {len(s)} proteins")
        print("\n  THE DIRECT TEST: the same protein's log multiplier across cell lines")
        print(f"    r = {r:+.3f}  spearman {rho:+.3f} (p={pv:.4f})  over {len(both)} "
              f"protein x model pairs")
        print("    This is what the section should lead with; a variance share is weaker.")

    pd.DataFrame(out).to_csv(TABLES / "multiplier_variance.csv", index=False)
    print("\nwrote multiplier_variance.csv")


if __name__ == "__main__":
    main()
