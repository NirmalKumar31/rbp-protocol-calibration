"""R4: does the model add anything conservation and position do not?

WHY THIS IS THE PRIMARY ClinVar ANALYSIS AND THE AUROC LADDER IS NOT.

The obvious framing -- rank the model against phyloP and against a positional rule, and see
who wins -- is both a replication (Grimm 2015; Schreiber 2020) and circular. ClinVar
pathogenic assertions for noncoding SNVs rest heavily on PP3, which is operationalised through
CADD/REVEL/GERP, all conservation-derived; benign assertions rest on allele frequency. So
phyloP is partly a re-read of the labels: benign variants here have mean phyloP 0.022 against
pathogenic 5.06. A horse race between a model and a partial proxy for the answer key measures
the answer key.

Conditioning is immune to that. Leakage into the labels inflates the competitor's AUROC; it
does not manufacture INCREMENTAL value for the model. So the question becomes "how much of the
model's variant signal survives controlling for conservation and position", and the answer is
essentially all of it.

TWO THINGS THIS SCRIPT ADDS THAT THE COEFFICIENT TABLE ALONE CANNOT.

  attenuation   As a FRACTION, not a difference. Both coefficients drifting together would
                satisfy an absolute gate and would not satisfy this one.

  decile check  Whether the 1-Mb positional rule is merely conservation in disguise. If the
                rule only works because conserved variants cluster, the paper has one baseline
                and not two, and the whole "two independent leaks" framing collapses. Measured
                WITHIN each phyloP decile, so conservation is held approximately fixed.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
TABLES = ROOT / "results" / "tables"
N_DECILES = 10
BLOCK = 1_000_000


def log(m):
    print(m, flush=True)


def main():
    coef = pd.read_csv(TABLES / "variant_coefficients.csv")
    refit = pd.read_csv(TABLES / "variant_specificity_refit.csv")

    unc = coef[(coef.arm == "matched") & (coef.standardisation == "within_dataset")].iloc[0]
    rows = [{"check": "coef, matched, unconditional", "value": float(unc.coef),
             "ci_low": float(unc.ci_low), "ci_high": float(unc.ci_high), "n": int(unc.n_variants),
             "note": "within-dataset standardised |delta|, 1-Mb block clustered"}]
    for _, r in refit.iterrows():
        rows.append({"check": f"coef, {r.arm}, controls = {r.controls}", "value": float(r.coef),
                     "ci_low": float(r.ci_low), "ci_high": float(r.ci_high), "n": int(r.n),
                     "note": ""})

    m_phy = refit[(refit.arm == "matched") & (refit.controls == "conservation only")].iloc[0]
    att = 1.0 - float(m_phy.coef) / float(unc.coef)
    rows.append({"check": "attenuation fraction after conditioning on phyloP", "value": att,
                 "ci_low": np.nan, "ci_high": np.nan, "n": int(m_phy.n),
                 "note": "THE CLAIM: conservation does not explain the model away"})

    mm_phy = refit[(refit.arm == "mismatched") & (refit.controls == "conservation only")].iloc[0]
    sep = float(m_phy.ci_low) > float(mm_phy.ci_high)
    rows.append({"check": "arms separated under control", "value": float(sep),
                 "ci_low": np.nan, "ci_high": np.nan, "n": int(m_phy.n),
                 "note": f"matched CI low {m_phy.ci_low:.3f} > mismatched CI high "
                         f"{mm_phy.ci_high:.3f}"})

    # --- is the positional rule conservation in disguise? --------------------------------
    cv = pd.read_csv(TABLES / "variant_conservation.csv")[["vid", "conservation"]]
    va = pd.read_csv(TABLES / "variant_assignments.csv")
    # DEDUPLICATE BEFORE THE BLOCK STATISTIC, NOT AFTER. variant_assignments.csv carries 2.40
    # rows per variant (one per dataset the variant was scored in), and the leave-one-out
    # subtraction removes exactly ONE of them. Computing prevalence on the duplicated table
    # therefore left ~1.4 copies of a variant's own label inside its own 1-Mb block, which is
    # own-label leakage into the baseline the paper calls trivial. Verified constant per vid:
    # 0/27,492 variants disagree on label or position across copies, so the dedup is lossless
    # and which copy survives cannot matter.
    d = (va.merge(cv, on="vid", how="left")
           .dropna(subset=["conservation"])
           .drop_duplicates("vid"))
    blk = d.chrom + "_" + (d.pos // BLOCK).astype(str)
    gb = d.assign(_b=blk).groupby("_b").label
    tot, cnt = gb.transform("sum"), gb.transform("size")
    # leave-one-out: a variant never sees its own label
    d["prev"] = ((tot - d.label) / (cnt - 1)).where(cnt > 1, np.nan)
    d = d.dropna(subset=["prev"])

    rho, pv = spearmanr(d.prev, d.conservation)
    rows.append({"check": "spearman(positional rule, phyloP)", "value": float(rho),
                 "ci_low": np.nan, "ci_high": np.nan, "n": len(d), "note": f"p={pv:.3g}"})
    rows.append({"check": "positional rule, unstratified",
                 "value": float(roc_auc_score(d.label, d.prev)), "ci_low": np.nan,
                 "ci_high": np.nan, "n": len(d), "note": ""})

    d["dec"] = pd.qcut(d.conservation, N_DECILES, labels=False, duplicates="drop")
    per = []
    for k, g in d.groupby("dec"):
        if g.label.nunique() == 2:
            per.append(roc_auc_score(g.label, g.prev))
    rows.append({"check": "positional rule, MIN AUROC within a phyloP decile",
                 "value": float(min(per)), "ci_low": np.nan, "ci_high": np.nan, "n": len(per),
                 "note": "if this is near 0.5 the rule IS conservation and there is one "
                         "baseline, not two"})
    rows.append({"check": "positional rule, MAX AUROC within a phyloP decile",
                 "value": float(max(per)), "ci_low": np.nan, "ci_high": np.nan, "n": len(per),
                 "note": ""})

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "incremental_value.csv", index=False)
    for _, x in out.iterrows():
        ci = f" [{x.ci_low:+.3f}, {x.ci_high:+.3f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:52} {x.value:+.4f}{ci}")
    log(f"\n  wrote {TABLES / 'incremental_value.csv'}")


if __name__ == "__main__":
    main()
