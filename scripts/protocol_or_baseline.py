"""R1n: is it the protocol, or the baseline it leaves? A referee found a way to ask. It is the baseline.

    python scripts/protocol_or_baseline.py

R1l argued that protocol and baseline cannot be separated because the three protocols' baseline
distributions barely overlap. That argument was too comfortable. A hostile referee found two
places where they DO overlap, asked the question there, and got an answer.

THE NATURAL EXPERIMENT. The neg2 protocol usually raises the composition baseline relative to
GC matching, but not always: in 27 of 94 datasets it LOWERS it. If the protocol label carried
information, neg2's deficit would persist in those 27. It does not -- it reverses. Whichever
protocol leaves the lower baseline gets the higher contribution, regardless of which protocol
that is.

THE MATCHED COMPARISON. Pair each dinucleotide-arm dataset with the GC-arm dataset whose
composition baseline is closest, keep pairs within 0.02 AUROC, and compare. The raw contrast is
+0.0398; matched on baseline it is indistinguishable from zero.

WHY THE dn-vs-gc CONTRAST CANNOT BE ASKED THIS WAY DIRECTLY. comp_dn < comp_gc in 94 of 94
datasets, so within a dataset the two arms are perfectly rank-confounded with the baseline by
construction. The matching above borrows across datasets, which is why it is a weaker design
than the neg2 discordance and why both are reported.

WHAT THIS DOES TO THE PAPER. It converts the thesis from "the protocol determines the measured
contribution" -- which invites "yes, via the baseline, so what?" -- into the sharper and more
useful "the composition baseline is what determines the measurable contribution, and the
protocol label carries essentially no information beyond it." That is a statement a benchmark
builder can act on: report the baseline, because it is the whole story.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
N_BOOT = 4000
MATCH_TOL = 0.02


def ci(v, seed=0, n_boot=N_BOOT):
    v = np.asarray(v, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), size=(n_boot, len(v)))
    return np.percentile(v[idx].mean(axis=1), [2.5, 97.5])


def main():
    d = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    rows = []

    print("=== 1. Does protocol add anything once the baseline is known? ===")
    long = pd.concat([pd.DataFrame({"arm": a, "comp": d[f"comp_{a}"], "gain": d[f"gain_{a}"]})
                      for a in ("gc", "dn", "neg2")], ignore_index=True)
    c, y = long.comp.values, long.gain.values
    def ols(X):
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ b
        return float(r @ r)
    tss = float((y - y.mean()) @ (y - y.mean()))
    Xb = np.column_stack([np.ones(len(c)), c, c ** 2])
    D = pd.get_dummies(long.arm, drop_first=True).values.astype(float)
    inc_protocol = (ols(Xb) - ols(np.column_stack([Xb, D]))) / tss
    # and the reverse, for contrast
    Xd = np.column_stack([np.ones(len(c)), D])
    inc_baseline = (ols(Xd) - ols(np.column_stack([Xd, c[:, None], (c ** 2)[:, None]]))) / tss
    rows += [{"check": "incremental R2 of the protocol label, given the baseline",
              "value": float(inc_protocol)},
             {"check": "incremental R2 of the baseline, given the protocol label",
              "value": float(inc_baseline)}]
    print(f"  protocol given baseline: {100 * inc_protocol:5.2f}% of variance")
    print(f"  baseline given protocol: {100 * inc_baseline:5.2f}%   <- an order of magnitude more")

    print("\n=== 2. The natural experiment: 27 datasets where neg2 LOWERS the baseline ===")
    hi = (d.comp_neg2 > d.comp_gc).values
    rows.append({"check": "datasets where neg2 raises the composition baseline",
                 "value": int(hi.sum())})
    for lab, m, key in (("neg2 baseline HIGHER", hi, "concordant"),
                        ("neg2 baseline LOWER", ~hi, "discordant")):
        diff = (d.gain_neg2 - d.gain_gc).values[m]
        lo, up = ci(diff)
        rows += [{"check": f"neg2 minus gc gain, {key} datasets", "value": float(diff.mean()),
                  "ci_low": float(lo), "ci_high": float(up), "n": int(m.sum())},
                 {"check": f"neg2 higher in, {key} datasets", "value": int((diff > 0).sum()),
                  "n": int(m.sum())}]
        print(f"  {lab:22s} n={m.sum():3d}  {diff.mean():+.4f} [{lo:+.4f}, {up:+.4f}]"
              f"  neg2 higher in {int((diff > 0).sum())}/{m.sum()}")
    r, p = spearmanr(d.comp_neg2 - d.comp_gc, d.gain_neg2 - d.gain_gc)
    rows.append({"check": "within-dataset spearman(delta baseline, delta gain)",
                 "value": float(r), "note": f"p={p:.1e}"})
    print(f"  within-dataset spearman(delta baseline, delta gain) = {r:+.3f}  p={p:.1e}")
    print("  -> whichever protocol leaves the LOWER baseline gets the higher contribution,")
    print("     regardless of which protocol that is")

    print("\n=== 3. dn vs gc, matched on the composition baseline ===")
    rank = int((d.comp_dn < d.comp_gc).sum())
    rows.append({"check": "datasets where the dinuc baseline is lower than the GC baseline",
                 "value": rank})
    print(f"  comp_dn < comp_gc in {rank}/94 -- perfectly rank-confounded within dataset,")
    print(f"  so the comparison has to borrow across datasets:")
    g = d[["comp_gc", "gain_gc"]].values
    matched = []
    for cb, gg in d[["comp_dn", "gain_dn"]].values:
        j = int(np.argmin(np.abs(g[:, 0] - cb)))
        if abs(g[j, 0] - cb) < MATCH_TOL:
            matched.append(gg - g[j, 1])
    matched = np.array(matched)
    lo, up = ci(matched)
    raw = float((d.gain_dn - d.gain_gc).mean())
    rows += [{"check": "dn minus gc, matched on baseline", "value": float(matched.mean()),
              "ci_low": float(lo), "ci_high": float(up), "n": len(matched)},
             {"check": "dn minus gc, unmatched (the published contrast)", "value": raw}]
    print(f"  matched within {MATCH_TOL} AUROC: n={len(matched)}  "
          f"{matched.mean():+.4f} [{lo:+.4f}, {up:+.4f}]")
    print(f"  unmatched: {raw:+.4f}  -> the contrast does not survive matching")

    pd.DataFrame(rows).to_csv(TABLES / "protocol_or_baseline.csv", index=False)
    print("\nwrote protocol_or_baseline.csv")


if __name__ == "__main__":
    main()
