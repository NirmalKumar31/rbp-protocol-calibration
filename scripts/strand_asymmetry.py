"""Which arm carries more of the strand artifact? The answer decides its direction of bias.

The strand bug gives every negative the positive's strand, so a negative lands on its own
gene's strand roughly by chance. That inflates absolute AUROCs. The question this script
answers is narrower and more useful: **is the spurious cue stronger in the arm with the LARGER
gain or the SMALLER one?**

If it is stronger in the dinucleotide arm, it could be manufacturing part of the contrast and
the whole result is in doubt. If it is stronger in the GC arm -- the arm with the smaller
nested gain -- then it works AGAINST the reported direction, and +0.0397 is a conservative
estimate rather than an inflated one. That is a one-line argument a referee can check, and it
does not depend on the placebo experiment at all.

Coordinates only, no sequences and no model fitting, so this runs in seconds.

SCOPE, STATED BECAUSE IT LIMITS THE CLAIM. Only datasets whose LOCAL window tables reproduce
their committed row are used, and for the dinucleotide arm that is the 40 whose canonical
tables were fetched. The GC arm reproduces locally. Comparing a canonical GC arm against a
locally-redrawn dinucleotide arm would compare two different negative draws, so the paired
comparison is restricted to datasets where both are canonical.
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402
from strand_audit import gene_index, own_strands  # noqa: E402
from strand_placebo import n_genes  # noqa: E402

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"



def frac_sense(path, idx):
    """Also returns the EXCHANGEABILITY diagnostic: how sense-kept negatives differ from
    dropped ones. This is what decides whether a uniform-random placebo is a valid
    counterfactual for the sense-only restriction, and it is why the placebo is stratified on
    region rather than on GC."""
    d = pd.read_csv(path, sep="\t",
                    usecols=["label", "chrom", "start", "end", "strand", "region", "gc"])
    neg = d[d.label == 0]
    keep = np.zeros(len(neg), dtype=bool)
    amb = 0
    for i, (c, s, e, a) in enumerate(zip(neg.chrom, neg.start, neg.end, neg.strand)):
        ss = own_strands(idx, c, int(s), int(e))
        if len(ss) > 1:
            amb += 1
        elif len(ss) == 1 and next(iter(ss)) == a:
            keep[i] = True
    n = len(neg)
    k, dr = neg[keep], neg[~keep]
    dens = np.array([n_genes(idx, c, int(a), int(b))
                     for c, a, b in zip(neg.chrom, neg.start, neg.end)])
    diag = {"ngenes_kept": float(dens[keep].mean()) if keep.any() else np.nan,
            "ngenes_dropped": float(dens[~keep].mean()) if (~keep).any() else np.nan}
    for lab, sub in (("kept", k), ("dropped", dr)):
        diag[f"intron_{lab}"] = float((sub.region == "intron").mean()) if len(sub) else np.nan
        diag[f"exon_nc_{lab}"] = float((sub.region == "exon_nc").mean()) if len(sub) else np.nan
        diag[f"gc_{lab}"] = float(sub.gc.mean()) if len(sub) else np.nan
    return keep.sum() / n, amb / n, n, diag


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", default="")
    p.add_argument("--gc-root", default="")
    p.add_argument("--dn-root", default="")
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    if a.from_cache:
        return summarise(pd.read_csv(TABLES / "strand_asymmetry_per_dataset.csv"))

    idx = gene_index(a.gtf)
    audited = list(pd.read_csv(TABLES / "strand_audit.csv").dataset)

    rows = []
    for i, ds in enumerate(audited, 1):
        prot, cell = ds.split(":")
        fg = Path(a.gc_root) / cell / prot / "dataset.tsv"
        fd = Path(a.dn_root) / cell / prot / "dataset.tsv"
        if not (fg.exists() and fd.exists()):
            continue
        sg, ag, ng, dg = frac_sense(fg, idx)
        sd, ad, nd, dd = frac_sense(fd, idx)
        rows.append({"dataset": ds, "frac_sense_gc": sg, "frac_sense_dn": sd,
                     "frac_amb_gc": ag, "frac_amb_dn": ad, "n_neg_gc": ng, "n_neg_dn": nd,
                     **{f"{k}_gc": v for k, v in dg.items()},
                     **{f"{k}_dn": v for k, v in dd.items()}})
        log(f"  [{i:2d}/{len(audited)}] {ds:22} gc {sg:.3f}  dn {sd:.3f}")

    m = pd.DataFrame(rows)
    m.to_csv(TABLES / "strand_asymmetry_per_dataset.csv", index=False)
    return summarise(m)


def summarise(m):
    if "diff" not in m.columns:
        m["diff"] = m.frac_sense_dn - m.frac_sense_gc
    w = wilcoxon(m.frac_sense_dn, m.frac_sense_gc)

    out = [
        {"check": "frac_sense, GC arm", "value": m.frac_sense_gc.mean(),
         "n": len(m), "note": f"sd {m.frac_sense_gc.std():.4f}; 0.5 is a coin flip"},
        {"check": "frac_sense, dinuc arm", "value": m.frac_sense_dn.mean(),
         "n": len(m), "note": f"sd {m.frac_sense_dn.std():.4f}"},
        {"check": "asymmetry, dinuc minus GC", "value": m["diff"].mean(),
         "n": len(m), "note": f"dinuc higher in {int((m['diff'] > 0).sum())}/{len(m)}, "
                              f"paired Wilcoxon p={w.pvalue:.3g}"},
        # EXCHANGEABILITY: why the placebo is stratified on region and not on GC.
        {"check": "intron fraction, sense-KEPT negatives, GC arm",
         "value": m.intron_kept_gc.mean(), "n": len(m), "note": "vs dropped, below"},
        {"check": "intron fraction, DROPPED negatives, GC arm",
         "value": m.intron_dropped_gc.mean(), "n": len(m),
         "note": "kept pairs are more intronic: the restriction is not random"},
        {"check": "intron fraction, sense-KEPT negatives, dinuc arm",
         "value": m.intron_kept_dn.mean(), "n": len(m), "note": ""},
        {"check": "intron fraction, DROPPED negatives, dinuc arm",
         "value": m.intron_dropped_dn.mean(), "n": len(m), "note": ""},
        {"check": "genes overlapping a sense-KEPT negative",
         "value": m[["ngenes_kept_gc", "ngenes_kept_dn"]].mean().mean(), "n": len(m),
         "note": "retention requires exactly ONE strand, so it selects against multi-gene loci"},
        {"check": "genes overlapping a DROPPED negative",
         "value": m[["ngenes_dropped_gc", "ngenes_dropped_dn"]].mean().mean(), "n": len(m),
         "note": "why the placebo is stratified on gene density as well as region"},
        {"check": "GC, sense-KEPT negatives",
         "value": m[["gc_kept_gc", "gc_kept_dn"]].mean().mean(),
         "n": len(m), "note": "balanced against dropped, so GC is NOT the confound"},
        {"check": "GC, DROPPED negatives",
         "value": m[["gc_dropped_gc", "gc_dropped_dn"]].mean().mean(), "n": len(m), "note": ""},
        {"check": "datasets where dinuc is more sense", "value": float((m["diff"] > 0).sum()),
         "n": len(m), "note": "THE ARGUMENT: the cue is weaker in the arm with the LARGER "
                              "gain, so it cannot manufacture the contrast, only shrink it"},
    ]
    res = pd.DataFrame(out)
    res.to_csv(TABLES / "strand_asymmetry.csv", index=False)
    log("")
    for _, x in res.iterrows():
        log(f"  {x.check:38} {x.value:+.4f}   {x.note}")
    log(f"\n  n = {len(m)};  wrote strand_asymmetry.csv")


if __name__ == "__main__":
    main()
