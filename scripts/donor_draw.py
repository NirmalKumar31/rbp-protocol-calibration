"""Build the multi-donor control manifest.

WHY THIS EXISTS. The original wrong-protein control used ONE donor per target, picked as
(index + 47) % 95 on an alphabetical manifest. That is a single deterministic draw, so the
confidence interval carries no uncertainty about donor choice -- and it turned out the draw
handed most targets a donor whose own model was WEAKER than theirs (binding AUROC 0.802 vs
0.850, Wilcoxon p=0.018). The measured gap therefore confounds "the head knows this protein"
with "my head was trained on more data than the donor's": gap vs log10(donor pairs) is
rho=-0.533 (p=1.9e-04), and where the donor is the stronger model the gap is -0.0025 (9/18,
p=1.00).

WHAT THIS FIXES. Five donors per target, drawn to SPAN the donor-quality range rather than
matched to it. Matching was the obvious fix and it is the wrong one twice over: it destroys
the variance needed to estimate the quality slope, and the pool is thin exactly where it
matters -- only 36/44 powered targets have five donors within |dlog10 pairs|<=0.15 AND
|dAUROC|<=0.05, and three have none, because targets are systematically the strong datasets.

Spanning instead gives a regression with support on both sides of zero, so the estimand is
an INTERCEPT (the gap at zero donor advantage) rather than a subset mean. The hard
constraint is >=2 donors stronger than the target where the pool allows it, which is what
puts mass on the far side and identifies the intercept.

Donors are screened on JACCARD of the assigned ClinVar variant sets, not on the shared_frac
the first version used. shared_frac = |t&d|/|t| is normalised by the target alone, so it
cannot be large when the donor is small -- it correlates with donor set size at rho=+0.696
and was measuring donor size wearing a contamination costume.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils import cloud as cloudcfg  # noqa: E402
from rbp.utils.log import log  # noqa: E402

MANIFEST = "variants/variant_tasks.tsv"
DONORS = "variants/donor_tasks.tsv"
N_DONORS = 5
MIN_STRONGER = 2          # donors that must beat the target, where the pool allows
JACCARD_MAX = 0.02        # variant-set overlap ceiling; see module docstring
SEED = 20260826



def build(tables: Path, seed=SEED):
    man = pd.read_csv(tables / "variant_tasks.tsv", sep="\t")
    man["dataset"] = man.protein + ":" + man.cell

    ps = pd.read_csv(tables / "panel_summary.csv")
    fm = pd.read_csv(tables / "matched_four_models.csv")
    va = pd.read_csv(tables / "variant_assignments.csv")
    va["dataset"] = va.protein + ":" + va.cell

    pairs = dict(zip(ps.dataset, ps.pairs))
    qual = dict(zip(fm.dataset, fm.splicebert))
    vids = {d: set(g.vid) for d, g in va.groupby("dataset")}

    man["pairs"] = man.dataset.map(pairs)
    man["qual"] = man.dataset.map(qual)
    missing = man[man.qual.isna() | man.pairs.isna()]
    if len(missing):
        log(f"WARNING {len(missing)} datasets lack pairs/quality: {list(missing.dataset)}")

    rng = np.random.default_rng(seed)
    rows, short = [], []
    for _, t in man.iterrows():
        tv = vids.get(t.dataset, set())
        pool = []
        for _, d in man.iterrows():
            if d.protein == t.protein:            # never the same protein, either cell line
                continue
            dv = vids.get(d.dataset, set())
            u = tv | dv
            j = len(tv & dv) / len(u) if u else 0.0
            if j > JACCARD_MAX:
                continue
            pool.append((int(d.idx), d.dataset, d.pairs, d.qual, j))

        p = pd.DataFrame(pool, columns=["idx", "dataset", "pairs", "qual", "jaccard"]).dropna()
        if len(p) < N_DONORS:
            short.append((t.dataset, len(p)))
            take = p
        else:
            stronger = p[p.qual > t.qual]
            weaker = p[p.qual <= t.qual]
            # >=2 from the stronger side where it exists: that is what puts support on the
            # far side of zero and makes the intercept identifiable rather than extrapolated.
            n_s = min(MIN_STRONGER, len(stronger))
            picks = [stronger.sample(n_s, random_state=int(rng.integers(1 << 31)))] if n_s else []
            rest = pd.concat([stronger.drop(picks[0].index) if picks else stronger, weaker])
            # spread the remainder across quality quantiles instead of sampling flat, so the
            # slope is estimated over the whole range and not from a clump.
            need = N_DONORS - n_s
            rest = rest.sort_values("qual").reset_index(drop=True)
            if len(rest) <= need:
                picks.append(rest)
            else:
                cuts = np.linspace(0, len(rest) - 1, need).round().astype(int)
                picks.append(rest.iloc[sorted(set(cuts))])
                while sum(len(x) for x in picks) < N_DONORS:
                    left = rest.drop(pd.concat(picks).index, errors="ignore")
                    if left.empty:
                        break
                    picks.append(left.sample(1, random_state=int(rng.integers(1 << 31))))
            take = pd.concat(picks).drop_duplicates("idx").head(N_DONORS)

        for _, d in take.iterrows():
            rows.append({"target_idx": int(t.idx), "target": t.dataset,
                         "donor_idx": int(d.idx), "donor": d.dataset,
                         "target_pairs": t.pairs, "donor_pairs": d.pairs,
                         "target_qual": t.qual, "donor_qual": d.qual,
                         "jaccard": d.jaccard})

    out = pd.DataFrame(rows).reset_index(drop=True)
    out.insert(0, "task_idx", out.index)
    return out, short


def report(d):
    per = d.groupby("target").size()
    adv = d.donor_qual - d.target_qual
    log(f"{len(d)} tasks, {d.target.nunique()} targets, {per.min()}-{per.max()} donors each")
    log(f"donor advantage (donor - target binding AUROC): "
        f"min {adv.min():+.3f} median {adv.median():+.3f} max {adv.max():+.3f}")
    log(f"  donors STRONGER than target: {(adv > 0).sum()}/{len(d)} "
        f"({d[adv > 0].target.nunique()} targets have >=1)")
    lr = np.log10(d.donor_pairs / d.target_pairs)
    log(f"log10 donor/target pairs: min {lr.min():+.2f} "
        f"median {lr.median():+.2f} max {lr.max():+.2f}")
    log(f"jaccard: median {d.jaccard.median():.4f} max {d.jaccard.max():.4f}")
    n_bal = sum(1 for _, g in d.groupby("target")
                if (g.donor_qual > g.target_qual).sum() >= MIN_STRONGER)
    log(f"targets with >={MIN_STRONGER} stronger donors: {n_bal}/{d.target.nunique()}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables", default=str(ROOT / "results" / "tables"))
    p.add_argument("--upload", action="store_true")
    a = p.parse_args()

    tables = Path(a.tables)
    d, short = build(tables)
    report(d)
    if short:
        log(f"thin pools ({len(short)}): {short[:6]}")

    local = tables / "donor_tasks.tsv"
    d.to_csv(local, sep="\t", index=False)
    log(f"wrote {local}")

    if a.upload:
        bucket = cloudcfg.bucket()
        bucket.blob(DONORS).upload_from_string(
            d.to_csv(sep="\t", index=False), content_type="text/tab-separated-values")
        log(f"uploaded gs://{bucket.name}/{DONORS}")


if __name__ == "__main__":
    main()
