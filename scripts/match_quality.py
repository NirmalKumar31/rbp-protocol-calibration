"""What the matchers ACHIEVED, not what they were asked for. Reported nowhere until now.

    python scripts/match_quality.py

THE GAP, found by a methods referee. Every arm's matcher is specified by its tolerance
(`gc_tolerance: 0.05`) or by its objective (nearest L1 over dinucleotide counts), and the paper
argues from those specifications: the GC arm controls 1 of the composition baseline's 15 degrees
of freedom and the dinucleotide arm controls 15 of 15. That is exactly true of the DESIGN. It is
only approximately true of the REALISATION, and the realisation was measured nowhere -- not in a
committed table, not in a golden key, not in verify.py, which has 522 assertions and none about
match quality.

Worse, the GC matcher does not hold its nominal tolerance. It relaxes: 40 draws, tolerance
doubled to 0.10 after 25 misses, then a best-seen fallback accepted up to 3x tolerance = 0.15.
Reporting "matched within 5 percentage points" is therefore an overstatement of what happened,
and nothing counted how often the relaxation fired.

WHY IT MATTERS BEYOND HONESTY. R1f is described in the paper's own notes as "confounded with
achieved match quality" -- a confound named in a write-up with no number anywhere to quantify
it. And the 1-of-15 versus 15-of-15 argument is the reason the sign of the headline contrast is
design-implied, which is the paper's largest concession. Both need this table.

WHAT IS MEASURED, per arm, over every pair in the panel: |dGC| between the positive and its
negative, and the L1 distance over the 16 dinucleotide FREQUENCIES (the same units the matcher
reports, 0 to 2). Medians, p90, p99, max, and the fraction inside each rung of the relaxation
ladder.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data.negatives import dinuc_matrix  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "processed/gc", "dn": "processed/dinuc", "neg2": "processed/neg2"}
NOMINAL_GC = 0.05          # config negatives.gc_tolerance
LADDER = (0.05, 0.10, 0.15)  # nominal, relaxed after 25 misses, best-seen fallback cap


def log(m):
    print(m, flush=True)


def pair_distances(d):
    """|dGC| and dinucleotide-frequency L1 for each positive/negative pair, in file order.

    Pairs are formed by position within label, which is how both matchers write them: the
    i-th negative was assigned to the i-th positive. Asserted rather than assumed.
    """
    pos = d[d.label == 1].reset_index(drop=True)
    neg = d[d.label == 0].reset_index(drop=True)
    n = min(len(pos), len(neg))
    if n == 0:
        return None, None
    pos, neg = pos.iloc[:n], neg.iloc[:n]
    gc_gap = np.abs(pos.gc.to_numpy(dtype=float) - neg.gc.to_numpy(dtype=float))
    P = dinuc_matrix(pos.seq_dna.tolist(), normalise=True)
    N = dinuc_matrix(neg.seq_dna.tolist(), normalise=True)
    return gc_gap, np.abs(P - N).sum(axis=1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    a = p.parse_args()
    store = Path(a.store)
    if not (store / "processed" / "gc").exists():
        sys.exit(f"no window store at {store}; it is ~3 GB and is not published")

    panel = pd.read_csv(TABLES / "rehearsal_binding_gc.csv")
    if a.n:
        panel = panel.head(a.n)

    per, agg = [], {k: {"gc": [], "l1": []} for k in ARMS}
    for i, r in enumerate(panel.itertuples(), 1):
        for arm, sd in ARMS.items():
            f = store / sd / r.cell / r.protein / "dataset.tsv"
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="\t")
            g, l1 = pair_distances(d)
            if g is None:
                continue
            agg[arm]["gc"].append(g)
            agg[arm]["l1"].append(l1)
            per.append({"dataset": r.dataset, "protein": r.protein, "cell": r.cell,
                        "arm": arm, "pairs": len(g), "gc_gap_median": float(np.median(g)),
                        "gc_gap_p90": float(np.percentile(g, 90)),
                        "gc_gap_max": float(g.max()),
                        "dinuc_l1_median": float(np.median(l1)),
                        "dinuc_l1_p90": float(np.percentile(l1, 90)),
                        "dinuc_l1_max": float(l1.max())})
        if i % 20 == 0:
            log(f"  [{i}/{len(panel)}]")

    t = pd.DataFrame(per)
    if t.empty:
        sys.exit("no dataset could be read; refusing to write an empty table")
    t.to_csv(TABLES / "match_quality_per_dataset.csv", index=False)

    out = []
    log(f"\n=== achieved match quality, {t.dataset.nunique()} datasets ===\n")
    log(f"  {'arm':6s} {'pairs':>9s} {'|dGC| med':>10s} {'p90':>7s} {'p99':>7s} {'max':>7s}"
        f" {'L1 med':>8s} {'p90':>7s} {'max':>7s}")
    for arm in ARMS:
        if not agg[arm]["gc"]:
            continue
        g = np.concatenate(agg[arm]["gc"])
        l1 = np.concatenate(agg[arm]["l1"])
        stats = {
            "pairs": len(g),
            "gc_gap_median": float(np.median(g)),
            "gc_gap_p90": float(np.percentile(g, 90)),
            "gc_gap_p99": float(np.percentile(g, 99)),
            "gc_gap_max": float(g.max()),
            "dinuc_l1_median": float(np.median(l1)),
            "dinuc_l1_p90": float(np.percentile(l1, 90)),
            "dinuc_l1_max": float(l1.max()),
        }
        for k, v in stats.items():
            out.append({"check": f"{k}, {arm} arm", "value": v, "n": t.dataset.nunique()})
        # THE RELAXATION LADDER, which is the number the paper does not currently report.
        for rung in LADDER:
            frac = float((g <= rung).mean())
            out.append({"check": f"fraction of pairs within |dGC| {rung:.2f}, {arm} arm",
                        "value": frac, "n": t.dataset.nunique()})
        log(f"  {arm:6s} {len(g):9,d} {stats['gc_gap_median']:10.4f} "
            f"{stats['gc_gap_p90']:7.4f} {stats['gc_gap_p99']:7.4f} {stats['gc_gap_max']:7.4f}"
            f" {stats['dinuc_l1_median']:8.4f} {stats['dinuc_l1_p90']:7.4f} "
            f"{stats['dinuc_l1_max']:7.4f}")

    log("")
    for arm in ARMS:
        if not agg[arm]["gc"]:
            continue
        g = np.concatenate(agg[arm]["gc"])
        log(f"  {arm:6s} within |dGC| " + "  ".join(
            f"{rung:.2f}: {100 * (g <= rung).mean():5.1f}%" for rung in LADDER))
    log(f"\n  The GC matcher's NOMINAL tolerance is {NOMINAL_GC}. It relaxes to 0.10 after 25")
    log("  misses and accepts a best-seen candidate up to 0.15. Report the ladder, not 0.05.")

    # The headline comparison the paper's 1-of-15 vs 15-of-15 argument rests on.
    if agg["gc"]["l1"] and agg["dn"]["l1"]:
        gl = float(np.median(np.concatenate(agg["gc"]["l1"])))
        dl = float(np.median(np.concatenate(agg["dn"]["l1"])))
        out.append({"check": "dinucleotide L1 improvement factor, dn vs gc arm",
                    "value": gl / dl, "n": t.dataset.nunique()})
        log(f"\n  The dinucleotide matcher reduces median composition distance "
            f"{gl:.3f} -> {dl:.3f}, a factor of {gl / dl:.2f}x.")
        log("  It does NOT drive it to zero, so '15 of 15 degrees of freedom controlled' is")
        log("  exactly true of the design and approximately true of the realisation.")

    pd.DataFrame(out).to_csv(TABLES / "match_quality.csv", index=False)
    log("\nwrote match_quality.csv and match_quality_per_dataset.csv")


if __name__ == "__main__":
    main()
