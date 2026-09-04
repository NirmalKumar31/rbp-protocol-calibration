"""Where does the design effect of 1.35 come from? Until this script, nowhere.

    python scripts/design_effect.py                 # needs the window store
    python scripts/design_effect.py --from-cache    # re-gate the committed table

THE GAP. `cluster_intervals.py` hard-codes `DESIGN_EFFECT = 1.35` with the comment "Measured
by a referee on this data", and the manuscript describes it as "the product of two quantities
measured on these data": 1.10 for DeLong conditioning on fitted score vectors as fixed, and
1.23 for treating spatially clustered windows as independent. Neither had any computing code
here, and the committed per-window scores carry `id, label, fold, score` and no coordinates,
so the block bootstrap was not reproducible from the release either. A published count rested
on a constant the harness could not check -- the same class of defect as the fold-partition
one, and found the same way, by looking rather than by a gate.

WHAT IS MEASURED, per dataset, on the GC arm:

  clustering  SD of the nested gain when whole genomic blocks are resampled, over the DeLong
              SE. Swept over block length, because the right cluster is arguable: windows
              concentrate at ~8 per gene and genes routinely exceed 10 kb.
  fitting     SD of the gain when the model score is permuted within fold and the nested
              model refit, over the mean DeLong SE under that same null. This is the
              variability that fitting two score vectors induces and that DeLong conditions
              away.

THE ANSWER, and it runs the author's way: clustering is 1.112 at 10 kb, 1.120 at 100 kb and
1.134 at 1 Mb, so block length barely matters; fitting is 0.996, i.e. DeLong's SE is about
right under the null. The product is ~1.11 against the applied 1.35, which moves the count of
datasets where the 4-mer significantly helps from 80/94 to 77/94 rather than to 72/94. The
published figure is therefore CONSERVATIVE, and it is kept for that reason -- but it is now a
measurement with a script behind it rather than a number on trust.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
ARM = ("processed/gc", "gc")
BLOCKS_NT = (10_000, 100_000, 1_000_000)
N_BLOCK_BOOT = 400
N_PERM = 20
SEED = 7
APPLIED = 1.35


def log(m):
    print(m, flush=True)


def block_sd(chrom, start, s_full, s_comp, y, block_nt, rng):
    """SD of the AUROC difference when whole `block_nt` genomic blocks are resampled."""
    codes = pd.Series([f"{c}:{p // block_nt}" for c, p in zip(chrom, start)]) \
        .astype("category").cat.codes.to_numpy()
    groups = [np.flatnonzero(codes == k) for k in range(codes.max() + 1)]
    out = []
    for _ in range(N_BLOCK_BOOT):
        idx = np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))])
        yy = y[idx]
        if yy.min() == yy.max():
            continue
        out.append(roc_auc_score(yy, s_full[idx]) - roc_auc_score(yy, s_comp[idx]))
    return float(np.std(out, ddof=1)) if len(out) > 2 else np.nan


def measure(store, panel, rng):
    from rbp.eval.baseline import oof_scores as kmer_oof
    from rbp.eval.delong import delong_test
    from rbp.eval.nested import _oof_scores, composition_features, standardise

    rows = []
    for i, r in enumerate(panel.itertuples(), 1):
        f = store / ARM[0] / r.cell / r.protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t", usecols=["seq_rna", "label", "fold", "chrom", "start"])
        y = d.label.to_numpy(dtype=int)
        folds = d.fold.to_numpy()
        sc, _, _ = kmer_oof(d.seq_rna.values, y, folds, k=4)
        comp, _ = composition_features(d.seq_rna.values)
        s_comp = _oof_scores(comp, y, folds)
        s_full = _oof_scores(np.column_stack([comp, standardise(sc)]), y, folds)
        ok = np.isfinite(s_comp) & np.isfinite(s_full) & np.isfinite(sc)
        se = float(delong_test(s_full[ok], s_comp[ok], y[ok])["se"])

        row = {"dataset": r.dataset, "n": int(ok.sum()), "se_delong": se}
        for b in BLOCKS_NT:
            sd = block_sd(d.chrom.values[ok], d.start.values[ok], s_full[ok], s_comp[ok],
                          y[ok], b, rng)
            row[f"ratio_clustering_{b // 1000}kb"] = sd / se if se else np.nan

        gains, ses = [], []
        for _ in range(N_PERM):
            p = sc.copy()
            for k in np.unique(folds):
                m = folds == k
                p[m] = rng.permutation(p[m])
            sf = _oof_scores(np.column_stack([comp, standardise(p)]), y, folds)
            o = np.isfinite(sf) & np.isfinite(s_comp)
            rr = delong_test(sf[o], s_comp[o], y[o])
            gains.append(float(rr["diff"]))
            ses.append(float(rr["se"]))
        null_se = float(np.mean(ses))
        row["ratio_fitting"] = np.std(gains, ddof=1) / null_se if null_se else np.nan
        rows.append(row)
        if i % 10 == 0:
            log(f"  [{i}/{len(panel)}]")
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    per = TABLES / "design_effect_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it with --store")
    else:
        store = Path(a.store)
        if not (store / ARM[0]).exists():
            sys.exit(f"no window store at {store}. It is not published because the windows "
                     f"carry genomic sequence; use --from-cache to re-gate the table.")
        t = measure(store, pd.read_csv(TABLES / "three_arm_per_dataset.csv"),
                    np.random.default_rng(SEED))
        if t.empty:
            sys.exit("nothing measured; refusing to overwrite the committed table")
        t.to_csv(per, index=False)

    cols = [f"ratio_clustering_{b // 1000}kb" for b in BLOCKS_NT] + ["ratio_fitting"]
    out = []
    log(f"\n=== design-effect components, {ARM[1]} arm, {len(t)} datasets ===\n")
    for c in cols:
        v = t[c].dropna()
        out.append({"check": c, "value": float(v.median()), "n": len(v),
                    "note": f"mean {v.mean():.3f}, p90 {v.quantile(0.90):.3f}"})
        log(f"  {c:28s} median {v.median():.3f}   mean {v.mean():.3f}   "
            f"p90 {v.quantile(0.90):.3f}")

    product = float(t[cols[0]].median() * t.ratio_fitting.median())
    out.append({"check": "measured design effect, 10 kb blocks", "value": product,
                "n": len(t), "note": f"applied in cluster_intervals.py: {APPLIED}"})
    log(f"\n  measured product {product:.3f} against the applied {APPLIED}, so the applied "
        f"figure is {'CONSERVATIVE' if APPLIED > product else 'ANTI-CONSERVATIVE'}")

    # The count the paper reports, at each factor, so the choice is visible rather than
    # asserted. The gain and its DeLong SE are per dataset in the deep-contrast table.
    d = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv")
    log("\n  datasets where the 4-mer significantly helps, GC arm:")
    for f in (1.00, round(product, 2), APPLIED):
        n = int(((d.kmer_gain_gc - 1.96 * d.kmer_se_gc * f) > 0).sum())
        out.append({"check": f"datasets significant at design effect {f:.2f}", "value": n,
                    "n": len(d), "note": ""})
        log(f"    factor {f:.2f}  ->  {n}/{len(d)}")

    pd.DataFrame(out).to_csv(TABLES / "design_effect.csv", index=False)
    log("\nwrote design_effect.csv and design_effect_per_dataset.csv")


if __name__ == "__main__":
    main()
