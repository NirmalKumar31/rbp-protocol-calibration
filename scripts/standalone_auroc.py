"""The 4-mer's own out-of-fold AUROC in each arm, which Table 1 needs and no table held.

    python scripts/standalone_auroc.py              # needs the window store
    python scripts/standalone_auroc.py --from-cache

WHY. Table 1's middle column used to be labelled "apparent AUROC" while holding the
composition-PLUS-score AUROC (0.8092 gc, 0.6937 dn), and the prose and abstract quoted a drop
of 0.1095 -- which is the drop in the model's OWN AUROC (0.7981 to 0.6879), a quantity that
appeared nowhere. Subtracting the two printed columns gave 0.1155 instead, so the table and the
abstract could not be reconciled by a reader. `audit_manuscript.py` passed it because 0.1095
does exist in `cost_of_matching.csv`; matching a value is not matching a quantity.

This emits all three arms on the paired 94-dataset panel, so the table carries composition
alone, model alone, composition plus model, and the contribution, and every difference in the
prose can be read off it.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
DIRS = {"dn": "dinuc", "gc": "gc", "neg2": "neg2"}
N_BOOT = 4000


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    per = TABLES / "standalone_auroc_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
    else:
        from rbp.eval.baseline import oof_scores as kmer_oof

        store = Path(a.store) / "processed"
        if not (store / "gc").exists():
            sys.exit(f"no window store at {store}; use --from-cache to re-gate the table")
        panel = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
        rows = []
        for i, r in enumerate(panel.itertuples(), 1):
            row = {"dataset": r.dataset, "protein": r.protein, "cell": r.cell}
            ok = True
            for arm, sub in DIRS.items():
                f = store / sub / r.cell / r.protein / "dataset.tsv"
                if not f.exists():
                    ok = False
                    break
                d = pd.read_csv(f, sep="\t", usecols=["seq_rna", "label", "fold"])
                sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
                m = np.isfinite(sc)
                row[f"auroc_{arm}"] = float(roc_auc_score(d.label.values[m], sc[m]))
                row[f"n_{arm}"] = int(m.sum())
            if ok:
                rows.append(row)
            if i % 20 == 0:
                print(f"  [{i}/{len(panel)}]", flush=True)
        t = pd.DataFrame(rows)
        if t.empty:
            sys.exit("nothing computed; refusing to overwrite the committed table")
        t.to_csv(per, index=False)

    three = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    m = t.merge(three, on=["dataset", "protein", "cell"])
    out = []
    print(f"\n=== 4-mer standalone out-of-fold AUROC, n={len(t)} datasets ===\n")
    print(f"  {'arm':6s} {'composition':>12s} {'model alone':>12s} {'comp+model':>11s}"
          f" {'contribution':>13s} {'comp beats model':>17s}")
    for arm in DIRS:
        k = m[f"auroc_{arm}"]
        comp, full, gain = m[f"comp_{arm}"], m[f"full_{arm}"], m[f"gain_{arm}"]
        beats = int((comp > k).sum())
        out += [{"check": f"model alone, {arm} arm", "value": float(k.mean()), "n": len(m)},
                {"check": f"datasets where composition beats the model, {arm} arm",
                 "value": beats, "n": len(m)}]
        print(f"  {arm:6s} {comp.mean():12.4f} {k.mean():12.4f} {full.mean():11.4f}"
              f" {gain.mean():+13.4f} {beats:14d}/{len(m)}")
    # HOW HARD A PROTOCOL LOOKS AND HOW HARD IT IS FOR COMPOSITION ARE NEARLY THE SAME THING,
    # and this is the number that says so. Section 3.4 argues that a protocol acts on measured
    # contribution through the composition baseline it leaves; the premise is that the model's
    # own AUROC and the composition-only AUROC move together, and it was never quantified.
    # Protein-clustered because 94 datasets are 79 proteins.
    #
    # CORRELATE THE MODEL ALONE, NOT comp+model. The nested column is a superset of the
    # composition column by construction and correlating them gives +0.94 for free, which is
    # the same wrong-quantity trap this script exists to fix.
    long = pd.concat([pd.DataFrame({"arm": a, "alone": m[f"auroc_{a}"],
                                    "comp": m[f"comp_{a}"], "protein": m.protein})
                      for a in DIRS], ignore_index=True)
    r_p = float(pearsonr(long.alone, long.comp)[0])
    r_s = float(spearmanr(long.alone, long.comp)[0])
    rng = np.random.default_rng(7)
    by = {q: gg for q, gg in long.groupby("protein")}
    names = long.protein.unique()
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        bb = pd.concat([by[q] for q in rng.choice(names, len(names), replace=True)],
                       ignore_index=True)
        draws[b] = pearsonr(bb.alone, bb.comp)[0]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    out.append({"check": "pearson(model alone, composition alone), pooled over arms",
                "value": r_p, "ci_low": float(lo), "ci_high": float(hi), "n": len(long),
                "note": f"spearman {r_s:+.4f}; {N_BOOT} protein draws"})
    print(f"\n  pearson(model alone, composition alone), pooled: {r_p:+.4f} "
          f"[{lo:+.4f}, {hi:+.4f}]   spearman {r_s:+.4f}")
    for arm in DIRS:
        k = long.arm == arm
        rr = float(pearsonr(long.alone[k], long.comp[k])[0])
        ss = float(spearmanr(long.alone[k], long.comp[k])[0])
        out.append({"check": f"pearson(model alone, composition alone), {arm} arm",
                    "value": rr, "n": int(k.sum()), "note": f"spearman {ss:+.4f}"})
        print(f"    within {arm:5s} {rr:+.4f}  (spearman {ss:+.4f})")

    drop = float(m.auroc_gc.mean() - m.auroc_dn.mean())
    out.append({"check": "model-alone AUROC drop, gc to dn", "value": drop, "n": len(m),
                "note": f"lower in {int((m.auroc_dn < m.auroc_gc).sum())}/{len(m)} datasets"})
    print(f"\n  model-alone drop, GC to dinucleotide: {drop:.4f} "
          f"(lower in {int((m.auroc_dn < m.auroc_gc).sum())}/{len(m)} datasets)")
    pd.DataFrame(out).to_csv(TABLES / "standalone_auroc.csv", index=False)
    print("\nwrote standalone_auroc.csv and standalone_auroc_per_dataset.csv")


if __name__ == "__main__":
    main()
