"""Does the composition finding survive stronger negatives?

The objection to every composition result in this project is: "your negatives are matched
on GC only, so of course composition discriminates. Match on dinucleotides and it goes
away." This script answers it by building both arms and comparing them directly.

THE PREDICTION I REGISTERED HERE BEFORE RUNNING IT, AND IT FAILED. I predicted that
datasets with extremity <= 0 would barely change, because I believed GC matching had already
reached the best achievable dinucleotide match for them. It had not: the "floor" I was
comparing against was the distance between two RANDOM windows (0.540), and a targeted
nearest-neighbour search reaches 0.220.

WHAT ACTUALLY HAPPENS. Dinucleotide matching drops composition-only AUROC from 0.793 to
0.609 and the model's AUROC from 0.819 to 0.702, so roughly 0.099 of a reported AUROC is
nucleotide composition. The model's real contribution RISES from +0.019 to +0.048, because
GC matching was crediting the model with work composition was doing.

So the central claim ("most of a reported AUROC is composition") is confirmed and now
demonstrated rather than inferred, and the secondary one ("the model adds almost nothing")
was too strong.

    python scripts/compare_negatives.py
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.eval import baseline, nested  # noqa: E402
from rbp.eval import extremity as ex

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def arm(path, k):
    """AUROC, composition baseline and gain for one prepared dataset."""
    df = pd.read_csv(path, sep="\t")
    y, folds, seqs = df.label.to_numpy(), df.fold.to_numpy(), df.seq_rna.tolist()
    res = baseline.evaluate(df, k=k)
    g = nested.gain_over_composition(seqs, res["scores"], y, folds)
    e = ex.from_dataset(df)
    return {"pairs": int(len(df) // 2), "auroc": res["auroc"],
            "composition": g.auroc_composition, "gain": g.delta,
            "gain_lo": g.delta_ci_low, "gain_hi": g.delta_ci_high,
            "extremity": e["extremity"], "l1_pos_neg": e["l1_pos_neg"],
            "l1_floor": e["l1_floor"]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--cell", default="K562")
    p.add_argument("--gc-dir", default="data/processed")
    p.add_argument("--dinuc-dir", default="data/processed_dinucmatch")
    a = p.parse_args()

    dinuc_root = ROOT / a.dinuc_dir / a.cell
    if not dinuc_root.exists():
        raise SystemExit(f"no dinucleotide-matched datasets in {dinuc_root}")
    names = sorted(d.name for d in dinuc_root.iterdir()
                   if (d / "dataset.tsv").exists())
    print(f"{len(names)} datasets with both arms\n")

    rows = []
    for i, name in enumerate(names, 1):
        gc_path = ROOT / a.gc_dir / a.cell / name / "dataset.tsv"
        if not gc_path.exists():
            continue
        g = arm(gc_path, a.k)
        d = arm(dinuc_root / name / "dataset.tsv", a.k)
        rows.append({"protein": name,
                     **{f"gc_{key}": v for key, v in g.items()},
                     **{f"dn_{key}": v for key, v in d.items()}})
        print(f"  [{i:2d}/{len(names)}] {name:9} "
              f"comp {g['composition']:.3f}->{d['composition']:.3f}  "
              f"gain {g['gain']:+.3f}->{d['gain']:+.3f}", flush=True)

    t = pd.DataFrame(rows)
    t["comp_change"] = t.dn_composition - t.gc_composition
    t["gain_change"] = t.dn_gain - t.gc_gain
    t["matchable"] = t.gc_extremity <= 0
    # Keyed by cell line. A single filename would have HepG2 overwrite K562 and silently
    # halve the comparison -- the same failure the flat peak directory had.
    t.insert(1, "cell", a.cell)
    t.to_csv(TABLES / f"compare_negatives_{a.cell}.csv", index=False)

    print(f"\n{'':=<74}")
    print(f"GC-MATCHED vs DINUCLEOTIDE-MATCHED NEGATIVES, {len(t)} datasets")
    print(f"{'':=<74}\n")
    print(f"{'':26} {'GC-matched':>12} {'dinuc-matched':>14} {'change':>9}")
    for c, lab in (("auroc", "model AUROC"), ("composition", "composition alone"),
                   ("gain", "gain over composition")):
        print(f"{lab:26} {t[f'gc_{c}'].median():12.4f} {t[f'dn_{c}'].median():14.4f} "
              f"{(t[f'dn_{c}'] - t[f'gc_{c}']).median():+9.4f}")

    print("\nachieved match quality (L1 over 16 dinucleotide frequencies, 0-2):")
    print(f"  GC-matched arm:    median {t.gc_l1_pos_neg.median():.3f}")
    print(f"  dinuc-matched arm: median {t.dn_l1_pos_neg.median():.3f}")
    print(f"  random-pair reference: median {t.gc_l1_floor.median():.3f}"
          f"   <- NOT a floor; the targeted match beats it")

    print("\n--- split by compositional extremity ---")
    for flag, lab in ((True, "extremity <= 0 (positives near their negatives)"),
                      (False, "extremity  > 0 (compositionally extreme positives)")):
        s = t[t.matchable == flag]
        if not len(s):
            continue
        print(f"\n{lab}  n={len(s)}")
        print(f"  composition alone {s.gc_composition.median():.4f} -> "
              f"{s.dn_composition.median():.4f}  "
              f"({s.comp_change.median():+.4f})")
        print(f"  gain              {s.gc_gain.median():+.4f} -> "
              f"{s.dn_gain.median():+.4f}  ({s.gain_change.median():+.4f})")

    print("\nDOES COMPOSITION STILL DISCRIMINATE UNDER THE STRONGER CONTROL?")
    print(f"  datasets where composition alone stays above 0.60: "
          f"{int((t.dn_composition > 0.60).sum())}/{len(t)}")
    print(f"  datasets where composition alone stays above 0.70: "
          f"{int((t.dn_composition > 0.70).sum())}/{len(t)}")
    print(f"  median composition AUROC under dinucleotide matching: "
          f"{t.dn_composition.median():.4f}")
    print(f"\n  corr(extremity, composition drop) = "
          f"{np.corrcoef(t.gc_extremity, t.comp_change)[0,1]:+.3f}"
          f"   <- POSITIVE means extreme datasets drop LESS, i.e. they resist matching")
    print(f"\nwrote results/tables/compare_negatives_{a.cell}.csv")


if __name__ == "__main__":
    main()
