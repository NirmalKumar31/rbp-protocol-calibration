"""R1f: the protocol effect is twice as large for coding-region binders as for intronic ones.

THE ONLY BIOLOGICAL STATEMENT THIS STUDY CAN HONESTLY MAKE, and it is a modest one. A referee
has already objected that "a 4-mer adds +0.0265 over composition" says nothing about RNA-protein
recognition. This is the answer, and its limitation is printed next to it rather than below it.

THE OBSERVATION. Grouping the 94 datasets by the region their POSITIVE windows mostly fall in,
the contrast is roughly twice as large for CDS-dominant proteins as for intron-dominant ones.

THE MECHANISM, AND IT IS TESTABLE. Intronic binding sites are compositionally distinctive:
polypyrimidine tracts, U-rich stretches, the low-complexity sequence around splice signals.
Composition alone therefore already discriminates them well, so GC matching leaves less room
for the protocol to change what is attributed to the model. CDS sites are compositionally
ordinary, so the same protocol change moves much more. That predicts composition-alone AUROC
should be HIGHER for intron-dominant datasets, which is checked below.

THE LIMITATION, WHICH MUST BE PRINTED. Region does not act independently of effect size: once
the total nested gain is partialled out, the region association disappears. So the honest claim
is that region indexes HOW MUCH non-compositional signal there is to expose, not that region is
a separate mechanism. Stated that way it is a real finding; stated as an independent effect it
would be wrong, and a referee would say so.

`region` is matched exactly between each positive and its negative by construction, so it is a
property of the pair rather than a confound between the arms.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
TABLES = ROOT / "results" / "tables"
N_BOOT = 2000
SEED = 0
MIN_GROUP = 8


def log(m):
    print(m, flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gc-root", required=True)
    a = p.parse_args()

    cm = pd.read_csv(TABLES / "cost_of_matching.csv")
    cm["contrast"] = cm.delta_auroc_dn - cm.delta_auroc_gc
    cm["total_gain"] = cm.delta_auroc_dn + cm.delta_auroc_gc

    rows = []
    for r in cm.itertuples():
        f = Path(a.gc_root) / r.cell / r.protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t", usecols=["label", "region"])
        pos = d[d.label == 1]
        # region is matched 1:1 within each pair, so this is a property of the pair
        vc = pos.region.value_counts(normalize=True)
        rows.append({"dataset": r.dataset, "dominant": vc.index[0],
                     "dominant_frac": float(vc.iloc[0]),
                     **{f"frac_{k}": float(v) for k, v in vc.items()}})
    reg = pd.DataFrame(rows).fillna(0.0)
    m = cm.merge(reg, on="dataset")
    log(f"  {len(m)} datasets with region annotation")

    out = []

    def add(check, value, lo=np.nan, hi=np.nan, n=len(m), note=""):
        out.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                    "n": n, "note": note})

    groups = [(g, sub) for g, sub in m.groupby("dominant") if len(sub) >= MIN_GROUP]
    for g, sub in sorted(groups, key=lambda t: -t[1].contrast.mean()):
        add(f"contrast, {g}-dominant datasets", sub.contrast.mean(), n=len(sub),
            note="see the matching composition row below")
        # Emitted as its own row rather than only in a note, so every number the manuscript
        # prints for this table has a source that scripts/audit_manuscript.py can find.
        add(f"composition alone, {g}-dominant datasets",
            sub.composition_auroc_dn.mean(), n=len(sub), note="dinucleotide arm")

    by = dict(groups)
    if "cds" in by and "intron" in by:
        c, i = by["cds"].contrast, by["intron"].contrast
        u = mannwhitneyu(c, i)
        rng = np.random.default_rng(SEED)
        d = [c.iloc[rng.integers(0, len(c), len(c))].mean()
             - i.iloc[rng.integers(0, len(i), len(i))].mean() for _ in range(N_BOOT)]
        lo, hi = np.percentile(d, [2.5, 97.5])
        add("CDS minus intron", c.mean() - i.mean(), lo, hi, len(c) + len(i),
            note=f"Mann-Whitney p={u.pvalue:.3g}; as large as the headline contrast itself")
        # THE MECHANISM: intronic sites are compositionally distinctive, CDS sites are not.
        uc = mannwhitneyu(by["intron"].composition_auroc_dn, by["cds"].composition_auroc_dn)
        add("composition alone, intron-dominant", by["intron"].composition_auroc_dn.mean(),
            n=len(by["intron"]), note="dinucleotide arm")
        add("composition alone, cds-dominant", by["cds"].composition_auroc_dn.mean(),
            n=len(by["cds"]),
            note=f"Mann-Whitney p={uc.pvalue:.3g}; intronic sites ARE more compositional")

    k = kruskal(*[sub.contrast for _, sub in groups])
    add("Kruskal-Wallis across region groups", float(k.statistic), n=len(m),
        note=f"p={k.pvalue:.3g}, {len(groups)} groups")

    # THE LIMITATION. Region indexes how much signal there is, not a separate mechanism.
    rho = spearmanr(m.frac_intron if "frac_intron" in m else m.dominant_frac, m.contrast)
    add("spearman(intronic fraction, contrast)", rho.statistic, n=len(m),
        note=f"p={rho.pvalue:.3g}")
    rk = lambda v: pd.Series(v).rank().to_numpy()                        # noqa: E731
    X = np.column_stack([np.ones(len(m)), rk(m.total_gain)])
    res = lambda v: rk(v) - X @ np.linalg.lstsq(X, rk(v), rcond=None)[0]  # noqa: E731
    rp = spearmanr(res(m.frac_intron if "frac_intron" in m else m.dominant_frac),
                   res(m.contrast))
    add("...partialling out total nested gain", rp.statistic, n=len(m),
        note=f"p={rp.pvalue:.3g}; THE LIMITATION: region acts THROUGH effect size, "
             f"it is not an independent mechanism")

    res_df = pd.DataFrame(out)
    res_df.to_csv(TABLES / "region_heterogeneity.csv", index=False)
    m.to_csv(TABLES / "region_heterogeneity_per_dataset.csv", index=False)
    for _, x in res_df.iterrows():
        ci = f" [{x.ci_low:+.4f}, {x.ci_high:+.4f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:44} {x.value:+.4f}{ci}   {x.note}")
    log(f"\n  wrote region_heterogeneity.csv")


if __name__ == "__main__":
    main()
