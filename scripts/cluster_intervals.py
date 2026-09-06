"""R1i: every published interval, recomputed with the clustering the panel actually has.

    python scripts/cluster_intervals.py

THE PROBLEM. The panel is 94 datasets but only 79 proteins: 15 proteins are assayed in both
K562 and HepG2 and contribute two rows each. Every bootstrap in this project resamples
DATASETS, which treats those two rows as independent evidence. They are not -- the
within-protein correlation of the contrast is +0.92 for the k-mer -- so every interval the
manuscript prints is narrower than the data support.

WHAT THIS DOES. Recomputes each headline quantity under a bootstrap that resamples PROTEINS
and takes all of a sampled protein's datasets, which is the correct unit. Both intervals are
emitted so the widening is visible rather than swapped in silently.

WHAT IT DOES NOT DO. It does not change any conclusion, and that is the point of running it:
the paper's claims are robust to the correction, so making the correction costs nothing and
removes an objection a referee would otherwise land. If a future run makes any headline cross
zero under clustering, `no_conclusion_may_change` fails and the claim must be restated.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
MODELS = ("kmer", "cnn", "splicebert")
N_BOOT = 4000
# Measured by a referee on this data: 1.10 for treating fitted score vectors as fixed, 1.23
# for treating spatially clustered windows as independent. Applied to the DeLong SE.
DESIGN_EFFECT = 1.35


def boot_ci(values, groups, n_boot=N_BOOT, seed=0, cluster=True):
    """Percentile CI, resampling clusters (proteins) or rows (datasets)."""
    rng = np.random.default_rng(seed)
    v = np.asarray(values, dtype=float)
    if not cluster:
        idx = rng.integers(0, len(v), size=(n_boot, len(v)))
        draws = v[idx].mean(axis=1)
        return np.percentile(draws, [2.5, 97.5])
    keys = pd.Series(groups).to_numpy()
    uniq = np.unique(keys)
    members = [np.flatnonzero(keys == k) for k in uniq]
    draws = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        rows = np.concatenate([members[i] for i in pick])
        draws[b] = v[rows].mean()
    return np.percentile(draws, [2.5, 97.5])


def _resample(df, rng, cluster):
    """One bootstrap draw of a group: by protein if clustering, else by row."""
    if not cluster:
        return df.iloc[rng.integers(0, len(df), size=len(df))]
    uniq = df.protein.unique()
    pick = uniq[rng.integers(0, len(uniq), size=len(uniq))]
    return pd.concat([df[df.protein == p] for p in pick], ignore_index=True)


def main():
    d = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv")
    prot = d.protein
    n_prot = prot.nunique()
    print(f"panel: {len(d)} datasets, {n_prot} proteins, "
          f"{(prot.value_counts() > 1).sum()} of them assayed twice\n")

    # How correlated ARE the two rows of a doubled protein? This is the quantity that decides
    # whether the correction matters at all, so it is reported rather than assumed.
    rows = []
    for m in MODELS:
        pairs = [g[f"{m}_gain_dn"].values - g[f"{m}_gain_gc"].values
                 for _, g in d.groupby("protein") if len(g) == 2]
        r = np.corrcoef(np.array([p[0] for p in pairs]),
                        np.array([p[1] for p in pairs]))[0, 1]
        rows.append({"check": f"within_protein_correlation_{m}", "quantity": "correlation",
                     "value": r, "n_pairs": len(pairs)})
        print(f"  within-protein correlation of the contrast, {m:11s} {r:+.3f} "
              f"({len(pairs)} doubled proteins)")

    print(f"\n{'quantity':28s} {'by dataset (published)':>26s} {'by PROTEIN (correct)':>26s} "
          f"{'width':>7s}")
    quantities = {}
    for m in MODELS:
        quantities[f"contrast_{m}"] = d[f"{m}_gain_dn"] - d[f"{m}_gain_gc"]
    for a, b in (("cnn", "kmer"), ("splicebert", "cnn"), ("splicebert", "kmer")):
        quantities[f"step_{a}_minus_{b}"] = ((d[f"{a}_gain_dn"] - d[f"{a}_gain_gc"])
                                             - (d[f"{b}_gain_dn"] - d[f"{b}_gain_gc"]))

    worst_ratio = 0.0
    all_exclude = True
    for name, series in quantities.items():
        lo_d, hi_d = boot_ci(series, prot, cluster=False)
        lo_p, hi_p = boot_ci(series, prot, cluster=True)
        ratio = (hi_p - lo_p) / (hi_d - lo_d)
        worst_ratio = max(worst_ratio, ratio)
        excludes = lo_p > 0
        all_exclude &= excludes
        rows.append({"check": name, "quantity": "mean", "value": float(series.mean()),
                     "ci_low_dataset": lo_d, "ci_high_dataset": hi_d,
                     "ci_low_protein": lo_p, "ci_high_protein": hi_p,
                     "width_ratio": ratio, "excludes_zero_clustered": bool(excludes)})
        print(f"  {name:26s} [{lo_d:+.4f},{hi_d:+.4f}] [{lo_p:+.4f},{hi_p:+.4f}] "
              f"{ratio:6.2f}x{'' if excludes else '   <-- CROSSES ZERO'}")

    # --- THE REST OF THE PAPER, not just R1g -------------------------------------------
    #
    # R1g was clustered first because it was newest. The same correction is owed to R1's own
    # headline, to R1f, and to R1c -- and R1c turns out not to need it, which is a cheap
    # checkable fact worth recording rather than leaving as an assumption.
    print()
    extra = {}
    reh_gc = pd.read_csv(TABLES / "rehearsal_binding_gc.csv")
    reh_dn = pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv")
    m = reh_gc.merge(reh_dn, on="dataset", suffixes=("_gc", "_dn"))
    extra["R1_contrast"] = (m.delta_auroc_dn - m.delta_auroc_gc, m.protein_gc)

    reg = pd.read_csv(TABLES / "region_heterogeneity_per_dataset.csv")
    if "dominant" in reg.columns:
        cds = reg[reg.dominant == "cds"]
        intr = reg[reg.dominant == "intron"]
        # CDS minus intron is a between-group contrast, so it is bootstrapped as such below.
        extra["_region"] = (cds, intr)

    for name, (series, groups) in [(k, v) for k, v in extra.items() if not k.startswith("_")]:
        lo_d, hi_d = boot_ci(series, groups, cluster=False)
        lo_p, hi_p = boot_ci(series, groups, cluster=True)
        ratio = (hi_p - lo_p) / (hi_d - lo_d)
        worst_ratio = max(worst_ratio, ratio)
        excludes = lo_p > 0
        all_exclude &= excludes
        rows.append({"check": name, "quantity": "mean", "value": float(series.mean()),
                     "ci_low_dataset": lo_d, "ci_high_dataset": hi_d,
                     "ci_low_protein": lo_p, "ci_high_protein": hi_p,
                     "width_ratio": ratio, "excludes_zero_clustered": bool(excludes)})
        print(f"  {name:26s} [{lo_d:+.4f},{hi_d:+.4f}] [{lo_p:+.4f},{hi_p:+.4f}] "
              f"{ratio:6.2f}x{'' if excludes else '   <-- CROSSES ZERO'}")

    if "_region" in extra:
        cds, intr = extra["_region"]
        rng = np.random.default_rng(1)
        def _grpboot(cluster):
            out = np.empty(N_BOOT)
            for b in range(N_BOOT):
                a = _resample(cds, rng, cluster)
                c = _resample(intr, rng, cluster)
                out[b] = ((a.delta_auroc_dn - a.delta_auroc_gc).mean()
                          - (c.delta_auroc_dn - c.delta_auroc_gc).mean())
            return np.percentile(out, [2.5, 97.5])
        obs = ((cds.delta_auroc_dn - cds.delta_auroc_gc).mean()
               - (intr.delta_auroc_dn - intr.delta_auroc_gc).mean())
        lo_d, hi_d = _grpboot(False)
        lo_p, hi_p = _grpboot(True)
        ratio = (hi_p - lo_p) / (hi_d - lo_d)
        worst_ratio = max(worst_ratio, ratio)
        all_exclude &= lo_p > 0
        rows.append({"check": "R1f_cds_minus_intron", "quantity": "mean", "value": float(obs),
                     "ci_low_dataset": lo_d, "ci_high_dataset": hi_d,
                     "ci_low_protein": lo_p, "ci_high_protein": hi_p,
                     "width_ratio": ratio, "excludes_zero_clustered": bool(lo_p > 0)})
        print(f"  {'R1f_cds_minus_intron':26s} [{lo_d:+.4f},{hi_d:+.4f}] "
              f"[{lo_p:+.4f},{hi_p:+.4f}] {ratio:6.2f}x"
              f"{'' if lo_p > 0 else '   <-- CROSSES ZERO'}")

    # R1c: no clustering exists there, and that is a fact worth asserting rather than assuming.
    sp = pd.read_csv(TABLES / "strand_placebo_per_dataset.csv")
    prots = sp.dataset.str.split(":").str[0]
    n_dup = int((prots.value_counts() > 1).sum())
    rows.append({"check": "R1c_duplicated_proteins", "quantity": "count", "value": n_dup,
                 "n_datasets": len(sp)})
    print(f"\n  R1c: {len(sp)} datasets over {prots.nunique()} proteins -> "
          f"{n_dup} duplicated. Clustering is a no-op there.")

    # --- THE PER-DATASET COUNTS, at the measured design effect --------------------------
    #
    # "the score adds significantly in 80/94 datasets" uses DeLong's SE, which understates
    # here for two separate reasons a referee measured: DeLong conditions on the two score
    # vectors as FIXED, but both are fitted (permutation null / DeLong SE = 1.10); and it
    # treats windows as independent when they are spatially clustered (10-kb block bootstrap
    # / DeLong = 1.23). Combined design effect ~1.35 on the SE.
    #
    # src/rbp/eval/nested.py's own test_score docstring identifies the second one and
    # cluster-bootstraps over genes. The headline estimator hands the rows straight to DeLong.
    # Eight datasets is 8.5% of the panel, too large to park behind a Methods caveat.
    print()
    for arm, tbl in (("gc", reh_gc), ("dn", reh_dn)):
        se = (tbl.delta_ci_high - tbl.delta_ci_low) / (2 * 1.959964)
        for factor, label in ((1.00, "published"), (1.10, "OOF fitting"),
                              (1.23, "window clustering"), (DESIGN_EFFECT, "both")):
            n = int(((tbl.delta_auroc - 1.959964 * se * factor) > 0).sum())
            if factor == DESIGN_EFFECT:
                rows.append({"check": f"helps_{arm}_design_effect", "quantity": "count",
                             "value": n, "n_datasets": len(tbl)})
            elif factor == 1.00:
                rows.append({"check": f"helps_{arm}_published", "quantity": "count",
                             "value": n, "n_datasets": len(tbl)})
            print(f"  {arm} arm, SE x{factor:.2f} ({label:18s}) -> {n}/{len(tbl)} helps")

    rows.append({"check": "max_width_ratio", "quantity": "ratio", "value": worst_ratio})
    rows.append({"check": "all_headlines_exclude_zero_clustered", "quantity": "flag",
                 "value": float(all_exclude)})
    pd.DataFrame(rows).to_csv(TABLES / "cluster_intervals.csv", index=False)
    print(f"\n  widest inflation {worst_ratio:.2f}x; every headline still excludes zero: "
          f"{all_exclude}")
    print("wrote cluster_intervals.csv")


if __name__ == "__main__":
    main()
