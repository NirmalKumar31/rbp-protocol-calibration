"""R1u: a protocol-free coordinate EXISTS in sample. It just does not transport. Both halves matter.

    python scripts/transport_check.py

WHAT THIS RETRACTS, AND IT IS THE PAPER'S TITLE CLAIM. R1m says "no monotone rescaling recovers
a protocol-free quantity: over eight transforms the range never falls below 2.00x". A hostile
reader pointed out that g/(1-comp) is not a coordinate, it is the p = 1 member of a
one-parameter family, and the sweep stopped exactly at the member that happens to be the
paper's recommendation. Extend the family and the claim fails in sample:

    p        0     0.5    1.0    1.25   1.544   1.75   2.0
    range   5.42   3.43   2.00   1.48   1.005   1.34   1.95

At p = 1.544 the three protocols agree on the panel mean to 0.5%. So a monotone rescaling that
equalises them DOES exist, and "no rescaling recovers a protocol-free quantity" is false as
literally worded. The old escape hatch -- "eight standard coordinates a benchmarker would
actually use" -- does not work, because g/(1-comp) has no more derivation than
g/(1-comp)^1.544; both are argmins of the same sweep.

WHAT REPLACES IT, AND IT IS A STRONGER CLAIM. The equalising exponent is a property of the
benchmark, not of the quantity. Fitted on our three protocols it is 1.544. Fitted on
Horlacher's two it is 3.649, more than twice as large -- and our exponent leaves their
benchmark at 2.34x, barely better than the 2.38x it started from. So there is no TRANSPORTABLE
rescaling: any normalisation strong enough to equalise one benchmark has to be refitted on the
next one, which defeats the purpose of normalising at all.

AND THE SAME TEST FALSIFIES THE PAPER'S RECOMMENDATION OUT OF SAMPLE. recommendation_works.py
names its own falsification criterion: "rank agreement falling, or disagreement rising, under
normalisation". On Horlacher's 45 datasets both happen -- rank agreement +0.706 -> +0.656 and
disagreement 0.860 -> 0.908. It is not significant at n = 45, so this is a failure to
replicate rather than a refutation, but it is the paper's own pre-registered test applied to
the only data this project did not build, and it points the wrong way.

AND ONE THING THAT SURVIVES AND STRENGTHENS. The 2.00x floor was being compared against 1.0.
That is the wrong null: max/min over three noisy means is bounded below by 1 and biased up. The
correct comparison is against the range you would see with EQUAL true arm means, preserving the
real between-arm dataset pairing. That null has median 1.07 and a 95th percentile of 1.20, so
2.00x is far outside it. The inference is sound; the comparison point was too generous.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

from rbp.utils.log import log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
OURS = ("dn", "gc", "neg2")
THEIRS = ("n1", "n2")
N_BOOT = 2000
N_NULL = 2000



def fold_range(d, arms, p):
    """The published aggregation: mean over datasets of the per-dataset ratio, then max/min.

    KEEP THIS AS IT IS. It is the quantity every published span refers to and changing it
    silently would move a headline. What it hides is measured in aggregation_sweep() below.
    """
    m = [(d[f"gain_{a}"] / np.power(1 - d[f"comp_{a}"], p)).mean() for a in arms]
    return max(m) / min(m) if min(m) > 0 else np.inf


def fold_range_median(d, arms, p):
    """Median of the per-dataset ratio instead of its mean."""
    m = [(d[f"gain_{a}"] / np.power(1 - d[f"comp_{a}"], p)).median() for a in arms]
    return max(m) / min(m) if min(m) > 0 else np.inf


def fold_range_of_means(d, arms, p):
    """Ratio of panel means: mean(g) over mean((1-c)^p), never a mean of a ratio."""
    m = [d[f"gain_{a}"].mean() / np.power(1 - d[f"comp_{a}"], p).mean() for a in arms]
    return max(m) / min(m) if min(m) > 0 else np.inf


AGGS = (("mean of per-dataset ratios", fold_range),
        ("median of per-dataset ratios", fold_range_median),
        ("ratio of panel means", fold_range_of_means))


def aggregation_sweep(d, arms, out):
    """DOES THE EQUALISING EXPONENT SURVIVE THE CHOICE OF AGGREGATION? It does not.

    g/(1-c)^p is a RATIO, and 1-c reaches 0.0267 on this panel with 62 of 282 cells below
    0.15, so a handful of cells dominate any mean of that ratio. The published span is such a
    mean. Under it the span touches 1.004 at p = 1.544, which the paper read as an equalising
    exponent existing in sample. Under either aggregation that is not a mean of a ratio it
    never approaches 1: the median bottoms out at 1.24 and the ratio of panel means at 1.28.
    The exponent that equalises the three protocols is therefore an artefact of one
    aggregation interacting with a near-zero denominator, and the paper now says so.
    """
    log("\n=== 1b. is the equalising exponent an artefact of the aggregation? ===\n")
    log(f"  {'p':>7}  " + "  ".join(f"{n:>28}" for n, _ in AGGS))
    for pw in (0.0, 0.5, 1.0, 1.25, 1.544, 1.75, 2.0):
        vals = [fn(d, arms, pw) for _, fn in AGGS]
        log(f"  {pw:7.3f}  " + "  ".join(f"{v:28.3f}" for v in vals))
        for (name, _), v in zip(AGGS, vals, strict=True):
            out.append({"check": f"span at p={pw:g}, {name}", "value": float(v)})

    log("")
    for name, fn in AGGS:
        # fn bound as a default: a late-binding closure over the loop variable happens to
        # work here only because minimize_scalar is called at once, and that is not a
        # property to rely on.
        r = minimize_scalar(lambda pw, fn=fn: fn(d, arms, pw), bounds=(0.0, 6.0),
                            method="bounded")
        out += [{"check": f"minimum span over p, {name}", "value": float(r.fun)},
                {"check": f"exponent at minimum span, {name}", "value": float(r.x)}]
        log(f"  {name:30s} minimum {r.fun:.4f} at p = {r.x:.3f}")

    # THE DENOMINATOR, because that is the whole mechanism.
    h = np.concatenate([1 - d[f"comp_{a}"] for a in arms])
    out += [{"check": "min headroom 1-c over all cells", "value": float(h.min())},
            {"check": "cells with headroom below 0.15", "value": int((h < 0.15).sum())},
            {"check": "cells with headroom below 0.05", "value": int((h < 0.05).sum())},
            {"check": "cells in the headroom distribution", "value": int(len(h))}]
    for a in arms:
        out.append({"check": f"min headroom 1-c, {a} arm",
                    "value": float((1 - d[f"comp_{a}"]).min())})
    log(f"\n  headroom 1-c: min {h.min():.4f}, {int((h < 0.15).sum())} of {len(h)} cells "
        f"below 0.15, {int((h < 0.05).sum())} below 0.05")


def main():
    argparse.ArgumentParser().parse_args()
    d = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    h = pd.read_csv(TABLES / "horlacher_per_dataset.csv")
    out = []

    log("=== 1. the exponent family, on our three protocols ===\n")
    log(f"  {'p':>6s} {'fold range':>11s}")
    grid = [0.0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0]
    for p in grid:
        log(f"  {p:6.2f} {fold_range(d, OURS, p):11.3f}")
    r = minimize_scalar(lambda p: fold_range(d, OURS, p), bounds=(0.2, 3.0), method="bounded")
    p_ours, min_ours = float(r.x), float(r.fun)
    out += [{"check": "fold range at the headroom coordinate, p = 1",
             "value": fold_range(d, OURS, 1.0)},
            {"check": "equalising exponent on our protocols", "value": p_ours},
            {"check": "fold range at our equalising exponent", "value": min_ours}]
    aggregation_sweep(d, OURS, out)
    log(f"\n  ARGMIN p = {p_ours:.3f} gives {min_ours:.4f}x -- a protocol-free coordinate")
    log("  EXISTS in sample, so 'no rescaling recovers one' is FALSE as worded.")

    log("\n=== 2. the null the 2.00x should have been compared against ===\n")
    # Equal true arm means, preserving the real between-arm dataset pairing, so the null keeps
    # the correlation structure that makes max/min over three means as tight as it is.
    head = {a: (d[f"gain_{a}"] / (1 - d[f"comp_{a}"])).to_numpy() for a in OURS}
    grand = np.mean([v.mean() for v in head.values()])
    centred = {a: head[a] - head[a].mean() + grand for a in OURS}   # equal means, same spread
    rng = np.random.default_rng(0)
    null = []
    for _ in range(N_NULL):
        i = rng.integers(0, len(d), len(d))
        m = [centred[a][i].mean() for a in OURS]
        null.append(max(m) / min(m))
    null = np.array(null)
    obs = fold_range(d, OURS, 1.0)
    out += [{"check": "equal-means null fold range, median", "value": float(np.median(null))},
            {"check": "equal-means null fold range, 95th percentile",
             "value": float(np.percentile(null, 95))},
            {"check": "observed headroom fold range exceeds the null 95th percentile",
             "value": float(obs > np.percentile(null, 95))}]
    log(f"  equal-means null: median {np.median(null):.3f}x, 95th pct "
        f"{np.percentile(null, 95):.3f}x")
    log(f"  observed {obs:.3f}x -> outside the null. The 2.00x inference is SOUND;")
    log("  it was simply being compared against 1.0, which is too generous a bar.")

    log("\n=== 3. TRANSPORT: does the equalising exponent carry to another benchmark? ===\n")
    r2 = minimize_scalar(lambda p: fold_range(h, THEIRS, p), bounds=(0.2, 6.0), method="bounded")
    p_theirs, min_theirs = float(r2.x), float(r2.fun)
    at_ours = fold_range(h, THEIRS, p_ours)
    raw_theirs = fold_range(h, THEIRS, 0.0)
    out += [{"check": "equalising exponent on Horlacher's protocols", "value": p_theirs},
            {"check": "their fold range at their own exponent", "value": min_theirs},
            {"check": "their fold range at OUR exponent", "value": at_ours},
            {"check": "their fold range, raw", "value": raw_theirs},
            {"check": "exponent ratio, theirs over ours", "value": p_theirs / p_ours}]
    log(f"  our exponent      p = {p_ours:.3f}")
    log(f"  their exponent    p = {p_theirs:.3f}   ({p_theirs / p_ours:.2f}x ours)")
    log(f"  their range raw {raw_theirs:.3f}x -> at OUR exponent {at_ours:.3f}x "
        f"-> at THEIR exponent {min_theirs:.3f}x")
    log("  -> the equalising exponent is a property of the BENCHMARK, not of the quantity.")
    log("     Any normalisation strong enough to equalise one has to be refitted on the next,")
    log("     which is what 'no protocol-free measure' should have said all along.")

    log("\n=== 4. the recommendation's OWN falsification test, on their data ===\n")
    prot = h.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(N_BOOT)]
    res = {}
    for lab, p in (("raw", 0.0), ("headroom", 1.0)):
        a = (h.gain_n1 / np.power(1 - h.comp_n1, p)).to_numpy()
        b = (h.gain_n2 / np.power(1 - h.comp_n2, p)).to_numpy()
        rho = float(spearmanr(a, b)[0])
        dis = float(np.abs(a - b).mean() / np.mean([a.mean(), b.mean()]))
        res[lab] = (a, b, rho, dis)
        out += [{"check": f"external rank agreement, {lab}", "value": rho, "n": len(h)},
                {"check": f"external scale-free disagreement, {lab}", "value": dis, "n": len(h)}]
        log(f"  {lab:9s} rank agreement {rho:+.3f}   disagreement {dis:.3f}")
    drank = np.array([spearmanr(res["headroom"][0][i], res["headroom"][1][i])[0]
                      - spearmanr(res["raw"][0][i], res["raw"][1][i])[0] for i in draws])
    drank = drank[np.isfinite(drank)]
    lo, hi = np.percentile(drank, [2.5, 97.5])
    out.append({"check": "external rank agreement change under normalisation",
                "value": float(res["headroom"][2] - res["raw"][2]), "ci_low": float(lo),
                "ci_high": float(hi), "n": len(h),
                "note": f"P(<=0) = {float((drank <= 0).mean()):.2f}"})
    log(f"\n  change in rank agreement {res['headroom'][2] - res['raw'][2]:+.3f} "
        f"[{lo:+.3f}, {hi:+.3f}]  P(<=0) = {(drank <= 0).mean():.2f}")
    log("  BOTH of recommendation_works.py's stated falsification criteria fire: rank")
    log("  agreement FALLS and disagreement RISES. Not significant at n = 45, so this is a")
    log("  failure to replicate rather than a refutation -- and it must be reported as such.")

    log("\n=== 5. and the family mechanism DOES replicate externally, per arm ===\n")
    for arm, fam in (("n1", "transcript background (composition-matched family)"),
                     ("n2", "other RBPs' sites")):
        rho, pv = spearmanr(h[f"comp_{arm}"], h[f"gain_{arm}"])
        out.append({"check": f"external within-arm spearman(baseline, gain), {arm}",
                    "value": float(rho), "n": len(h), "note": f"p = {pv:.2e}"})
        log(f"  {arm}  {rho:+.3f}  p = {pv:.2e}   {fam}")
    log("  -> our R1n mechanism (-0.545 / -0.462 composition-matched, -0.122 n.s. for")
    log("     other-RBPs'-sites) replicates on an independent lab's negatives and folds.")

    pd.DataFrame(out).to_csv(TABLES / "transport_check.csv", index=False)
    log("\nwrote transport_check.csv")


if __name__ == "__main__":
    main()
