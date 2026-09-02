"""R1q: does the paper's own recommendation actually help? Measured, not asserted.

    python scripts/recommendation_works.py

THE GAP THIS FILLS. The paper's deliverable is "report the composition-only AUROC under the
same protocol alongside every headline AUROC". Every result up to here shows that NOT doing so
is a problem -- a 5.4-fold range, a floor of 2.00x under any rescaling, a protocol label that
adds 1% once the baseline is known. None of them shows that DOING it helps. A benchmarking
paper whose recommendation is untested is a complaint, not a contribution.

THE TEST. Put a reader in the position the recommendation is meant to rescue: they have two
papers using different negative-set protocols and want to compare what the models contributed.
Compare two ways of doing that.

  raw        compare the reported nested contributions directly
  headroom   compare gain / (1 - composition AUROC), which is what the two-number report buys

Two measures, because they answer different halves of the question:

  RANK AGREEMENT (Spearman across datasets between the two protocols). Scale-free, so no
  normalisation can flatter it. Answers "would a reader rank proteins the same way?"

  DISAGREEMENT ON A COMMON SCALE (mean |A - B| divided by the pooled panel mean of that
  coordinate). Raw and headroom-normalised numbers have different units, so the raw absolute
  difference is not comparable between them; dividing each by its own panel mean is what makes
  the comparison fair, and doing it any other way would rig the result.

WHAT WOULD FALSIFY THE RECOMMENDATION. Rank agreement falling, or disagreement rising, under
normalisation. Then the honest paper says "we can show the problem and we cannot offer a fix",
which is still publishable and considerably weaker. It does not fall.

R1m already showed the headroom coordinate is the least protocol-sensitive of eight, at 2.00x
against 5.42x. This is the same fact from the reader's side, and it is the one that belongs in
the Discussion.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
ARMS = ("gc", "dn", "neg2")
PAIRS = (("gc", "dn"), ("gc", "neg2"), ("dn", "neg2"))
N_BOOT = 4000


def main():
    d = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    prot = d.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == p) for p in uniq]
    rng = np.random.default_rng(0)
    draws = [np.concatenate([members[i] for i in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(N_BOOT)]

    raw = {a: d[f"gain_{a}"].to_numpy() for a in ARMS}
    head = {a: (d[f"gain_{a}"] / (1 - d[f"comp_{a}"])).to_numpy() for a in ARMS}
    m_raw = float(np.mean([raw[a].mean() for a in ARMS]))
    m_head = float(np.mean([head[a].mean() for a in ARMS]))

    rows = []
    print(f"n = {len(d)} datasets over {len(uniq)} proteins; intervals bootstrap the protein\n")
    print("A reader holding two papers that used different protocols, comparing what the")
    print("models contributed. Does the two-number report help?\n")
    print(f"  {'pair':16s} {'rank agreement':>28s} {'scale-free disagreement':>30s}")
    print(f"  {'':16s} {'raw -> headroom':>28s} {'raw -> headroom':>30s}")
    better_rank = better_dis = 0
    for a, b in PAIRS:
        r1, _ = spearmanr(raw[a], raw[b])
        r2, _ = spearmanr(head[a], head[b])
        dr = float(np.abs(raw[a] - raw[b]).mean() / m_raw)
        dh = float(np.abs(head[a] - head[b]).mean() / m_head)
        # protein-clustered interval on the IMPROVEMENT, which is the claim
        drank = np.array([spearmanr(head[a][i], head[b][i])[0]
                          - spearmanr(raw[a][i], raw[b][i])[0] for i in draws])
        lo, hi = np.percentile(drank[np.isfinite(drank)], [2.5, 97.5])
        better_rank += int(r2 > r1)
        better_dis += int(dh < dr)
        rows += [
            {"check": f"rank agreement, raw, {a} vs {b}", "value": float(r1)},
            {"check": f"rank agreement, headroom, {a} vs {b}", "value": float(r2)},
            {"check": f"rank agreement gain, {a} vs {b}", "value": float(r2 - r1),
             "ci_low": float(lo), "ci_high": float(hi)},
            {"check": f"scale-free disagreement, raw, {a} vs {b}", "value": dr},
            {"check": f"scale-free disagreement, headroom, {a} vs {b}", "value": dh},
        ]
        print(f"  {a + ' vs ' + b:16s} {r1:11.3f} -> {r2:.3f}  [{lo:+.3f},{hi:+.3f}]"
              f" {dr:16.3f} -> {dh:.3f}")

    rows += [{"check": "protocol pairs where rank agreement improves", "value": better_rank},
             {"check": "protocol pairs where disagreement shrinks", "value": better_dis},
             {"check": "number of protocol pairs", "value": len(PAIRS)}]
    print(f"\n  rank agreement improves in {better_rank}/{len(PAIRS)} pairs; "
          f"disagreement shrinks in {better_dis}/{len(PAIRS)}")
    print("  -> the recommendation is tested, not merely asserted")

    pd.DataFrame(rows).to_csv(TABLES / "recommendation_works.csv", index=False)
    print("\nwrote recommendation_works.csv")


if __name__ == "__main__":
    main()
