"""B14: where along the window does the discriminative signal sit?

    python scripts/positional_signal.py --store ../rbp-store

THE ASSERTION THIS TESTS. Methods justifies the CNN's global max-pool with "exploratory
analysis placed the discriminative signal approximately 15 nt off centre, varying by protein".
That claim carried no table, no figure and no script -- it is exactly the shape of the unsourced
number this project's manuscript audit exists to catch, and the audit could not catch it because
"15" is a bare integer inside a prose clause about design rationale.

THE MEASUREMENT. For each dataset and each of the 101 positions, the mutual information between
the base at that position and the label, in bits, over all windows. Mutual information rather
than a per-position AUROC because the base is a four-level categorical and a single position has
no natural ordering; and it needs no model, so nothing here depends on a fit.

Positives are 101 nt windows on the peak MIDPOINT, so position 50 is the midpoint and the offset
of the peak-information position from 50 is the quantity the Methods claim is about. Reported as
the median absolute offset over datasets, with its spread, and the fraction of datasets whose
peak sits within a few nucleotides of the centre.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from rbp.utils.log import log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
BASES = ("A", "C", "G", "U", "T")



def positional_mi(seqs, y):
    """Mutual information in bits between the base at each position and the label."""
    n = len(y)
    L = len(seqs[0])
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(n, L)
    y = np.asarray(y, dtype=bool)
    py = np.array([(~y).mean(), y.mean()])
    hy = -np.sum(py[py > 0] * np.log2(py[py > 0]))
    codes = {ord(b) for b in BASES}
    mi = np.zeros(L)
    for j in range(L):
        col = arr[:, j]
        tot = 0.0
        for b in codes:
            m = col == b
            k = int(m.sum())
            if not k:
                continue
            p1 = y[m].mean()
            pv = np.array([1 - p1, p1])
            h = -np.sum(pv[pv > 0] * np.log2(pv[pv > 0]))
            tot += (k / n) * h
        mi[j] = max(hy - tot, 0.0)
    return mi


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--arm", default="gc", choices=["gc", "dinuc", "neg2"])
    p.add_argument("--near", type=int, default=5, help="'near the centre' half-width, nt")
    a = p.parse_args()

    three = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    root = Path(a.store) / "processed" / a.arm
    rows, profiles = [], []
    for i, r in enumerate(three.itertuples(), 1):
        f = root / r.cell / r.protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t", usecols=["seq_rna", "label"])
        seqs = d.seq_rna.astype(str).tolist()
        if not seqs or len({len(s) for s in seqs}) != 1:
            continue
        mi = positional_mi(seqs, d.label.values)
        centre = len(mi) // 2
        peak = int(np.argmax(mi))
        rows.append({"dataset": r.dataset, "protein": r.protein, "cell": r.cell,
                     "length": len(mi), "centre": centre, "peak_position": peak,
                     "offset": peak - centre, "abs_offset": abs(peak - centre),
                     "mi_peak": float(mi[peak]), "mi_centre": float(mi[centre]),
                     "mi_mean": float(mi.mean())})
        profiles.append(pd.DataFrame({"dataset": r.dataset, "position": np.arange(len(mi)),
                                      "mi": mi}))
        if i % 20 == 0:
            log(f"  [{i}/{len(three)}]")
    if not rows:
        sys.exit("no datasets measured")
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "positional_signal_per_dataset.csv", index=False)
    pd.concat(profiles, ignore_index=True).to_csv(
        TABLES / "positional_signal_profile.csv", index=False)

    med = float(t.abs_offset.median())
    near = int((t.abs_offset <= a.near).sum())
    out = [{"check": "datasets measured for positional signal", "value": len(t), "n": len(t)},
           {"check": "median absolute offset of peak information from the window centre (nt)",
            "value": med, "n": len(t)},
           {"check": "mean absolute offset of peak information from the centre (nt)",
            "value": float(t.abs_offset.mean()), "n": len(t)},
           {"check": "IQR of the absolute offset (nt)",
            "value": float(t.abs_offset.quantile(0.75) - t.abs_offset.quantile(0.25)),
            "n": len(t)},
           {"check": f"datasets whose peak sits within {a.near} nt of the centre",
            "value": near, "n": len(t)},
           {"check": "median peak positional mutual information (bits)",
            "value": float(t.mi_peak.median()), "n": len(t)},
           {"check": "median positional mutual information at the centre (bits)",
            "value": float(t.mi_centre.median()), "n": len(t)}]
    pd.DataFrame(out).to_csv(TABLES / "positional_signal.csv", index=False)

    log(f"\n=== B14: positional signal, {a.arm} arm, n={len(t)} datasets ===\n")
    log("  absolute offset of the peak-information position from the centre:")
    log(f"    median {med:.1f} nt   mean {t.abs_offset.mean():.1f} nt   "
        f"IQR {t.abs_offset.quantile(0.25):.0f}-{t.abs_offset.quantile(0.75):.0f} nt   "
        f"max {t.abs_offset.max():.0f} nt")
    log(f"    within {a.near} nt of the centre in {near}/{len(t)} datasets")
    log(f"  peak MI {t.mi_peak.median():.5f} bits vs {t.mi_centre.median():.5f} at the centre")
    log("\n  the Methods claim is 'approximately 15 nt off centre, varying by protein'.")
    log("\nwrote positional_signal.csv, _per_dataset.csv and _profile.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
