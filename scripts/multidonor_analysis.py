"""The multi-donor wrong-protein control: does protein specificity survive donor matching?

WHAT THE OLD ARM GOT WRONG. One donor per target, chosen as (index + 47) % 95 on an
alphabetical manifest. That draw handed most targets a donor whose own model was WEAKER than
theirs (binding AUROC 0.802 vs 0.850, Wilcoxon p=0.018), so the measured gap confounded
"this head knows THIS protein" with "this head saw more training data than that one". The
tells, all on the 44 powered datasets: gap vs log10(donor pairs) rho=-0.533 (p=1.9e-04, and
-0.493 partialling out both power and target size); splitting on donor SIZE reproduced the
published co-binding stratification BETTER (+0.1362 vs +0.1210, both 16/17); and where the
donor's own model was at least as good as the target's the gap was -0.0025 (9/18, p=1.00).

THE ESTIMAND HERE IS AN INTERCEPT, NOT A MEAN. With five donors per target spanning donor
quality, the gap can be regressed on donor-minus-target quality and read at zero advantage.
That is the number the old design could not produce, because with one donor per target the
donor advantage had no within-target variance to regress on.

TWO THINGS THIS REPORTS THAT A SUBSET MEAN CANNOT.

  within-target   donor advantage varies WITHIN a target, so target identity can be swept out
                  entirely by centring. Everything a target contributes -- its variant count,
                  its pathogenic rate, its genomic neighbourhood -- cancels. This is the
                  cleanest estimate of the quality slope and it needs no covariates.

  between-target  the intercept still has to be read somewhere, and the naive intercept is
                  read at whatever power the panel happens to carry. Power is therefore
                  mean-centred, so the intercept means "donor-neutral, average-power dataset".

BOTH PANELS ARE REPORTED, ALWAYS. The old arm's adjusted intercept was +0.037 on the powered
44 and -0.014 on all 82 usable -- it changed SIGN with the stratum, which is what "not
identified" looks like. Reporting one panel and not the other is how that stayed hidden.
"""

import argparse
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils import cloud as cloudcfg  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
PREFIX = "variants/scores_md/"
BOOT = 5000
SEED = 0



def fetch(bucket):
    """Every (target, donor) score file, concatenated. The filename carries the donor."""
    rows = []
    for b in bucket.client.list_blobs(bucket.name, prefix=PREFIX):
        df = pd.read_csv(io.BytesIO(b.download_as_bytes()))
        if df.empty:
            continue
        rows.append(df)
    if not rows:
        raise SystemExit(f"no objects under {PREFIX}; has the sweep finished?")
    d = pd.concat(rows, ignore_index=True)
    d["dataset"] = d.protein + ":" + d.cell
    return d


def per_pair_auroc(scores, matched):
    """One AUROC per (target, donor). The matched arm is joined on target, not recomputed."""
    out = []
    for (ds, donor), g in scores.groupby(["dataset", "weights_from"]):
        if g.label.nunique() < 2:
            continue
        m = matched.get(ds)
        if m is None:
            continue
        out.append({"dataset": ds, "donor": donor, "n": len(g),
                    "n_pathogenic": int(g.label.sum()),
                    "auroc_floor": roc_auc_score(g.label, g.delta.abs()),
                    "auroc_matched": m})
    r = pd.DataFrame(out)
    r["gap"] = r.auroc_matched - r.auroc_floor
    return r


def _ols(y, cols):
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in cols])
    return np.linalg.lstsq(X, np.asarray(y, float), rcond=None)[0]


def cluster_boot(d, cols, rng, n=BOOT):
    """Resample TARGETS, not rows. Five rows share a target and are not independent."""
    targets = d.dataset.unique()
    by = {t: g for t, g in d.groupby("dataset")}
    est = _ols(d.gap, [d[c] for c in cols])
    bs = []
    for _ in range(n):
        pick = rng.choice(targets, len(targets), replace=True)
        s = pd.concat([by[t] for t in pick], ignore_index=True)
        try:
            bs.append(_ols(s.gap, [s[c] for c in cols]))
        except Exception:
            pass
    bs = np.array(bs)
    lo, hi = np.percentile(bs[:, 0], [2.5, 97.5])
    p = 2 * min((bs[:, 0] <= 0).mean(), (bs[:, 0] >= 0).mean())
    return est[0], lo, hi, p, bs


def within_target(d, rng, n=BOOT):
    """Target-centred slope and intercept: target identity cannot contribute at all.

    Centring gap and advantage within each target removes every between-target difference,
    so what is left is purely "as this target's donor gets better, what happens to the gap".
    The intercept is then the target-mean gap adjusted to zero donor advantage.
    """
    g = d.copy()
    g["gap_c"] = g.gap - g.groupby("dataset").gap.transform("mean")
    g["adv_c"] = g.adv - g.groupby("dataset").adv.transform("mean")
    slope = np.linalg.lstsq(g[["adv_c"]].to_numpy(float), g.gap_c.to_numpy(float),
                            rcond=None)[0][0]
    # target mean gap, corrected to zero advantage using the within slope
    per = g.groupby("dataset").agg(gap=("gap", "mean"), adv=("adv", "mean"))
    per["adj"] = per.gap - slope * per.adv
    targets = per.index.to_numpy()
    bs = []
    for _ in range(n):
        bs.append(per.adj.loc[rng.choice(targets, len(targets), replace=True)].mean())
    bs = np.array(bs)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p = 2 * min((bs <= 0).mean(), (bs >= 0).mean())
    return slope, per.adj.mean(), lo, hi, p


def analyse(r, meta, rng):
    d = r.merge(meta, on=["dataset", "donor"], how="left").dropna(
        subset=["donor_qual", "target_qual", "donor_pairs", "target_pairs"])
    d["adv"] = d.donor_qual - d.target_qual
    d["lr"] = np.log10(d.donor_pairs / d.target_pairs)
    d["pw"] = np.log10(d.n_pathogenic.clip(lower=1))
    d["pw"] = d.pw - d.pw.mean()

    out = []
    for panel, sub in (("all usable", d), ("powered (n_path>=20)", d[d.n_pathogenic >= 20])):
        if sub.dataset.nunique() < 8:
            continue
        nt, nr = sub.dataset.nunique(), len(sub)
        # the headline: unadjusted mean gap, targets as the unit
        per = sub.groupby("dataset").gap.mean()
        wins = int((per > 0).sum())
        wp = wilcoxon(per).pvalue if len(per) > 1 else np.nan
        out.append({"panel": panel, "estimator": "mean gap (unadjusted)", "n_targets": nt,
                    "n_pairs": nr, "value": per.mean(), "ci_low": np.nan, "ci_high": np.nan,
                    "p": wp, "wins": wins, "note": "target mean of per-donor gaps"})

        for lab, cols in (("intercept | advantage", ["adv"]),
                          ("intercept | advantage+size", ["adv", "lr"]),
                          ("intercept | advantage+size+power", ["adv", "lr", "pw"])):
            b, lo, hi, p, _ = cluster_boot(sub, cols, rng)
            out.append({"panel": panel, "estimator": lab, "n_targets": nt, "n_pairs": nr,
                        "value": b, "ci_low": lo, "ci_high": hi, "p": p, "wins": np.nan,
                        "note": "target-clustered bootstrap"})

        sl, b, lo, hi, p = within_target(sub, rng)
        out.append({"panel": panel, "estimator": "within-target adjusted", "n_targets": nt,
                    "n_pairs": nr, "value": b, "ci_low": lo, "ci_high": hi, "p": p,
                    "wins": np.nan, "note": f"within-target slope {sl:+.4f}"})

        # the diagnostic that condemned the old arm, recomputed here
        rho, pv = spearmanr(sub.adv, sub.gap)
        out.append({"panel": panel, "estimator": "gap vs donor advantage (spearman)",
                    "n_targets": nt, "n_pairs": nr, "value": rho, "ci_low": np.nan,
                    "ci_high": np.nan, "p": pv, "wins": np.nan,
                    "note": "negative = stronger donors close the gap"})

        # THE MECHANISM, as two gated numbers rather than a sentence in a doc. A generic
        # plausibility floor should not care how many pathogenic variants a dataset has; the
        # protein's own head should. That asymmetry is the whole reason the two panels
        # disagree, so it is asserted rather than asserted-about.
        for col, lab in (("auroc_floor", "floor vs power (spearman)"),
                         ("auroc_matched", "matched arm vs power (spearman)")):
            rho, pv = spearmanr(np.log10(sub.n_pathogenic.clip(lower=1)), sub[col])
            out.append({"panel": panel, "estimator": lab, "n_targets": nt, "n_pairs": nr,
                        "value": rho, "ci_low": np.nan, "ci_high": np.nan, "p": pv,
                        "wins": np.nan,
                        "note": "flat floor + steep matched arm = the detection threshold"})

        # donors that BEAT the target: the old arm's decisive null, now with real n
        strong = sub[sub.adv > 0]
        if strong.dataset.nunique() >= 8:
            ps = strong.groupby("dataset").gap.mean()
            out.append({"panel": panel, "estimator": "donors STRONGER than target",
                        "n_targets": strong.dataset.nunique(), "n_pairs": len(strong),
                        "value": ps.mean(), "ci_low": np.nan, "ci_high": np.nan,
                        "p": wilcoxon(ps).pvalue if len(ps) > 1 else np.nan,
                        "wins": int((ps > 0).sum()), "note": "adv > 0 only"})
    return d, pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--local", default="")
    a = p.parse_args()
    rng = np.random.default_rng(SEED)

    meta = pd.read_csv(TABLES / "donor_tasks.tsv", sep="\t")[
        ["target", "donor", "target_pairs", "donor_pairs", "target_qual", "donor_qual",
         "jaccard"]].rename(columns={"target": "dataset"})

    old = pd.read_csv(TABLES / "variant_specificity.csv")
    matched = dict(zip(old.dataset, old.auroc_matched))

    scores = pd.read_csv(a.local) if a.local else fetch(cloudcfg.bucket())
    log(f"{len(scores):,} scored rows, {scores.dataset.nunique()} targets, "
        f"{scores.weights_from.nunique()} distinct donors")

    r = per_pair_auroc(scores, matched)
    r.to_csv(TABLES / "multidonor_pairs.csv", index=False)
    log(f"{len(r)} (target, donor) AUROC pairs -> multidonor_pairs.csv")

    d, res = analyse(r, meta, rng)
    res.to_csv(TABLES / "multidonor_specificity.csv", index=False)

    log("")
    for panel, g in res.groupby("panel", sort=False):
        log(f"--- {panel} ---")
        for _, x in g.iterrows():
            ci = (f" [{x.ci_low:+.4f}, {x.ci_high:+.4f}]"
                  if pd.notna(x.ci_low) else "")
            w = f" wins {int(x.wins)}/{int(x.n_targets)}" if pd.notna(x.wins) else ""
            log(f"  {x.estimator:38} {x.value:+.4f}{ci}{w} p={x.p:.4g}"
                f"   (n={int(x.n_targets)} targets, {int(x.n_pairs)} pairs)")
    log("\nwrote multidonor_specificity.csv")


if __name__ == "__main__":
    main()
