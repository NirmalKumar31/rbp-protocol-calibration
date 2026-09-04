"""B4: is the bias-aware arm's deficit label noise from co-binding?

    python scripts/cobinding_noise.py --store ../rbp-store

THE ALTERNATIVE EXPLANATION THIS TESTS, and it is the most serious one left for the bias-aware
arm. Its negatives are other RBPs' binding sites. Many RBPs co-bind the same transcripts, so
some of those "negatives" may be sites the TARGET also binds -- which is not a hard negative, it
is a mislabelled positive. Label noise of that kind depresses any model's measured contribution,
and it would do so without any of the compositional story this paper tells.

The construction already excludes a donor window within `negatives.min_peak_distance` of the
target's own positive WINDOWS. That is not the same as excluding it from the target's PEAKS: the
positives are one 101 nt window per peak midpoint, so a wide peak is represented by a window
that covers a fraction of it, and a donor site elsewhere in the same peak survives the filter.
This measures exactly that residual.

THE TEST. Per dataset, what fraction of bias-aware negatives fall inside a called peak of the
target protein? Then: does the contribution deficit, gain_neg2 - gain_gc, track that fraction
across datasets? If co-binding drives the deficit, datasets with more residual overlap should
show a larger one. If the deficit is compositional, it should not.

Stratifying rather than adjusting, because the fraction is not randomised and a regression
coefficient on it would invite a causal reading it cannot support.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data import encode  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

TABLES = ROOT / "results" / "tables"
PEAKROOT = ROOT.parent / "rna-binding-proteins"


def log(m):
    print(m, flush=True)


def peak_intervals(protein, cell):
    """Merged peak intervals per chromosome, as sorted (start, end) arrays."""
    try:
        path = encode.peak_path(PEAKROOT, protein, cell)
    except (FileNotFoundError, IndexError):
        return None
    if path is None or not Path(path).exists():
        return None
    p = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"], compression="gzip")
    out = {}
    for chrom, g in p.groupby("chrom"):
        iv = g[["start", "end"]].to_numpy()
        iv = iv[np.argsort(iv[:, 0])]
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        m = np.asarray(merged)
        out[chrom] = (m[:, 0], m[:, 1])
    return out


def inside_fraction(windows, peaks):
    """Fraction of window midpoints that fall inside a merged peak interval."""
    if peaks is None or not len(windows):
        return np.nan
    hit = 0
    for chrom, g in windows.groupby("chrom"):
        if chrom not in peaks:
            continue
        starts, ends = peaks[chrom]
        mid = ((g.start + g.end) // 2).to_numpy()
        i = np.searchsorted(starts, mid, side="right") - 1
        ok = (i >= 0) & (ends[np.clip(i, 0, len(ends) - 1)] > mid)
        hit += int(ok.sum())
    return hit / len(windows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    a = p.parse_args()
    cfg = cfgmod.load()
    margin = int(cfg.negatives["min_peak_distance"])
    three = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    root = Path(a.store) / "processed" / "neg2"

    rows = []
    for i, r in enumerate(three.itertuples(), 1):
        f = root / r.cell / r.protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t", usecols=["chrom", "start", "end", "label"])
        neg = d[d.label == 0]
        pk = peak_intervals(r.protein, r.cell)
        frac = inside_fraction(neg, pk)
        rows.append({"dataset": r.dataset, "protein": r.protein, "cell": r.cell,
                     "n_neg": len(neg), "frac_in_target_peak": frac,
                     "deficit": float(r.gain_neg2 - r.gain_gc),
                     "gain_neg2": float(r.gain_neg2), "gain_gc": float(r.gain_gc),
                     "comp_neg2": float(r.comp_neg2)})
        if i % 20 == 0:
            log(f"  [{i}/{len(three)}]")
    t = pd.DataFrame(rows).dropna(subset=["frac_in_target_peak"])
    if t.empty:
        sys.exit("no datasets measured")
    t.to_csv(TABLES / "cobinding_noise_per_dataset.csv", index=False)

    out = [{"check": "window exclusion margin applied at construction (nt)",
            "value": margin, "n": len(t)},
           {"check": "median fraction of bias-aware negatives inside a target peak",
            "value": float(t.frac_in_target_peak.median()), "n": len(t)},
           {"check": "mean fraction of bias-aware negatives inside a target peak",
            "value": float(t.frac_in_target_peak.mean()), "n": len(t)},
           {"check": "max fraction of bias-aware negatives inside a target peak",
            "value": float(t.frac_in_target_peak.max()), "n": len(t)}]

    rho, pv = spearmanr(t.frac_in_target_peak, t.deficit)
    out.append({"check": "spearman(residual co-binding, bias-aware deficit)",
                "value": float(rho), "n": len(t), "note": f"p={pv:.3f}"})

    # STRATIFY. If co-binding drives the deficit, the high-overlap half should show the larger
    # one. Terciles, so a monotone trend is visible rather than assumed linear.
    t = t.sort_values("frac_in_target_peak").reset_index(drop=True)
    k = len(t) // 3
    strata = [("lowest third", t.iloc[:k]), ("middle third", t.iloc[k:2 * k]),
              ("highest third", t.iloc[2 * k:])]
    log(f"\n=== B4: residual co-binding in the bias-aware arm, n={len(t)} datasets ===\n")
    log(f"  exclusion applied at construction: {margin} nt from the target's own windows")
    log(f"  fraction of negatives inside a called target peak: "
        f"median {t.frac_in_target_peak.median():.4f}, "
        f"mean {t.frac_in_target_peak.mean():.4f}, max {t.frac_in_target_peak.max():.4f}")
    log(f"\n  spearman(overlap, deficit) = {rho:+.3f}  (p={pv:.3f})\n")
    log(f"  {'stratum':16} {'n':>4} {'overlap':>9} {'deficit':>10}")
    for name, g in strata:
        out.append({"check": f"bias-aware deficit, {name} by residual co-binding",
                    "value": float(g.deficit.mean()), "n": len(g),
                    "note": f"mean overlap {g.frac_in_target_peak.mean():.4f}"})
        log(f"  {name:16} {len(g):4d} {g.frac_in_target_peak.mean():9.4f} "
            f"{g.deficit.mean():+10.4f}")
    lo, hi = strata[0][1].deficit.mean(), strata[-1][1].deficit.mean()
    out.append({"check": "deficit in the highest overlap third minus the lowest",
                "value": float(hi - lo), "n": len(t)})
    log(f"\n  highest third minus lowest third: {hi - lo:+.4f}")
    log("  a co-binding explanation predicts this is strongly NEGATIVE; "
        "a compositional one predicts near zero.")
    pd.DataFrame(out).to_csv(TABLES / "cobinding_noise.csv", index=False)
    log("\nwrote cobinding_noise.csv and cobinding_noise_per_dataset.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
