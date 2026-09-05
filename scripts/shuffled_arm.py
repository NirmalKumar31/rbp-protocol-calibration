"""B5: the shuffled fourth arm, where difficulty and contribution move TOGETHER.

    python scripts/shuffled_arm.py --n 6      # smoke
    python scripts/shuffled_arm.py           # the panel

WHY THIS ARM EXISTS. Dinucleotide-preserving shuffling is what GraphProt, iDeepS and RBPsuite
actually do, and it is the protocol \\citet{tourne2026} indicts. It is not a matching procedure
at all: each negative is a permutation of its own positive with the dinucleotide counts held
exactly, so mononucleotide frequencies, dinucleotide frequencies, GC and sequence entropy are
IDENTICAL between a positive and its negative, row by row. The nineteen-feature composition
baseline therefore takes the same value on both members of every pair.

THE PREDICTION, and it has no free parameters. The composition baseline must sit at exactly
0.5, not approximately: every pair is a tie, so no threshold separates them. Whatever a model
scores is then credited entirely as contribution, because the baseline it is measured against
is uninformative by construction.

WHY IT MATTERS TO THIS PAPER RATHER THAN TO THAT LITERATURE. Every arm in the main analysis
shows difficulty and contribution moving in OPPOSITE directions: the protocol that makes
discrimination harder yields a smaller increment over composition. The shuffled arm is the one
place in our own data where they move together, and it is not a counterexample but the
mechanism: matching raises the baseline by making negatives compositionally similar in
DISTRIBUTION, while shuffling makes them identical in composition PAIRWISE and so removes the
baseline entirely. The two operations look alike and do opposite things to the measurement.

4-mer only. The neural arms were trained on the three matched protocols; scoring them here
would need a GPU retrain, and the point does not need one -- it is about the baseline, which is
model-free.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data import shuffle as sh  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import composition_features  # noqa: E402
from rbp.eval.nested import _oof_scores  # noqa: E402
from rbp.stats import standardise  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
SEED = 7



def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    key = "dataset" if "dataset" in pub.columns else pub.columns[0]
    datasets = list(pub[key])[:limit or None]
    root = Path(store) / "processed" / "gc"
    rows = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        f = root / cell / protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        pos = d[d.label == 1].reset_index(drop=True)
        seqs, dropped = sh.shuffled_negatives(pos.seq_rna.tolist(), seed=SEED, method="dinuc")
        keep = [j for j, s in enumerate(seqs) if s is not None]
        if len(keep) < 200:
            continue
        p = pos.loc[keep].reset_index(drop=True)
        neg = p.copy()
        neg["seq_rna"] = [seqs[j] for j in keep]
        neg["label"] = 0
        dd = pd.concat([p, neg], ignore_index=True)

        # THE CONSTRUCTION CHECK, per dataset, before anything is fitted. verify() compares
        # dinucleotide counts; if a single pair fails, the arm is not the arm it claims to be
        # and its baseline of 0.5 would be an accident rather than a consequence.
        bad = sum(0 if sh.verify(a, b, k=2) else 1
                  for a, b in zip(p.seq_rna, neg.seq_rna))
        y, fo = dd.label.values, dd.fold.values
        X, _ = composition_features(dd.seq_rna.values)
        s_comp = _oof_scores(X, y, fo)
        sc, _, _ = kmer_oof(dd.seq_rna.values, y, fo, k=4)
        s_full = _oof_scores(np.column_stack([X, standardise(sc)]), y, fo)
        good = np.isfinite(s_comp) & np.isfinite(s_full)
        r = delong_test(s_full[good], s_comp[good], y[good])
        r_alone = delong_test(standardise(sc)[good], s_comp[good], y[good])
        # HOW MANY COMPOSITION ROWS ARE LITERALLY TIED. This is the claim in its most direct
        # form and it does not go through an AUROC at all.
        Xr = np.round(X, 9)
        tied = int(sum(np.array_equal(Xr[j], Xr[j + len(p)]) for j in range(len(p))))
        rows.append({"dataset": ds, "protein": protein, "cell": cell, "pairs": len(p),
                     "dropped_failed": dropped["failed"],
                     "dropped_similar": dropped["too_similar"],
                     "dinuc_violations": bad, "tied_composition_rows": tied,
                     "comp_auroc": float(r["auc_b"]), "full_auroc": float(r["auc_a"]),
                     "kmer_gain": float(r["diff"]), "kmer_alone": float(r_alone["auc_a"])})
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} pairs {len(p):6d}  comp {r['auc_b']:.4f}  "
            f"kmer alone {r_alone['auc_a']:.4f}  gain {r['diff']:+.4f}  "
            f"tied {tied}/{len(p)}  dinuc violations {bad}")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no dataset could be built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "shuffled_arm_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it without --from-cache")
    else:
        t = build(a.store, a.n)
        t.to_csv(per, index=False)

    out = []
    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(2000)]

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})
        return float(v.mean())

    log(f"\n=== B5: the dinucleotide-shuffled arm, n = {len(t)}, {len(uniq)} proteins ===\n")

    # THE CONSTRUCTION, FIRST AND AS A HARD STOP. Everything below is a consequence of the
    # negatives being exact dinucleotide permutations of the positives. If they are not, the
    # 0.5 baseline is a coincidence and the section says nothing.
    v = int(t.dinuc_violations.sum())
    out.append({"check": "pairs whose dinucleotide counts differ from their source",
                "value": v, "n": len(t)})
    if v:
        sys.exit(f"{v} shuffled negatives do not preserve their source's dinucleotide counts")
    frac_tied = float((t.tied_composition_rows / t.pairs).mean())
    out.append({"check": "fraction of pairs with an identical composition feature vector",
                "value": frac_tied, "n": len(t),
                "note": "identical to 9 decimals, so the baseline cannot separate them"})
    log(f"  dinucleotide counts preserved in every pair; composition feature vectors "
        f"identical in {100 * frac_tied:.2f}% of pairs")

    c = add("composition AUROC, shuffled arm", t.comp_auroc)
    out.append({"check": "max |composition AUROC - 0.5|, shuffled arm",
                "value": float((t.comp_auroc - 0.5).abs().max()), "n": len(t)})
    g = add("4-mer nested contribution, shuffled arm", t.kmer_gain)
    al = add("4-mer standalone AUROC, shuffled arm", t.kmer_alone)
    log(f"  composition {c:.4f}   4-mer alone {al:.4f}   contribution {g:+.4f}")

    # THE WHOLE POINT: the contribution equals the model's own AUROC minus a half, because the
    # baseline contributes nothing. Stating the residual makes that an assertion rather than an
    # observation about two numbers that happen to be close.
    resid = add("contribution minus (standalone AUROC - 0.5), shuffled arm",
                t.kmer_gain - (t.kmer_alone - 0.5))
    log(f"  contribution - (standalone - 0.5) = {resid:+.5f}: with an uninformative baseline "
        f"the increment IS the model's own AUROC")

    # AND THE COMPARISON WITH THE THREE MATCHED ARMS, which is why this is in the paper.
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    key = "dataset" if "dataset" in pub.columns else pub.columns[0]
    j = t.merge(pub, on=key, how="inner", suffixes=("", "_pub"))
    log(f"\n  against the three matched arms, on the {len(j)} shared datasets:")
    for arm, lbl in (("dn", "dinucleotide-matched"), ("gc", "GC-matched"),
                     ("neg2", "bias-aware")):
        cc = f"comp_{arm}"
        gg = f"gain_{arm}"
        if cc not in j.columns or gg not in j.columns:
            continue
        add(f"composition AUROC, {arm} arm, shared datasets", j[cc])
        add(f"4-mer contribution, {arm} arm, shared datasets", j[gg])
        log(f"    {lbl:22s} baseline {j[cc].mean():.4f}   contribution "
            f"{j[gg].mean():+.4f}")
    log(f"    {'dinucleotide-SHUFFLED':22s} baseline {j.comp_auroc.mean():.4f}   contribution "
        f"{j.kmer_gain.mean():+.4f}")

    # THE DIRECTION, WHICH IS THE PAPER'S CLAIM AND ITS BOUNDARY. Across the three matched
    # arms a lower baseline goes with a LARGER contribution, which is the inversion the title
    # is about. The shuffled arm has the lowest baseline of all AND the largest contribution,
    # so it lies on the SAME side of the relation, not the opposite one. State it that way:
    # shuffling is not a fourth point on the matching axis, it is a different operation.
    if "comp_dn" in j.columns:
        lowest = j[["comp_gc", "comp_dn", "comp_neg2"]].mean().min()
        out.append({"check": "shuffled baseline below the lowest matched baseline",
                    "value": float(lowest - j.comp_auroc.mean()), "n": len(j)})
        biggest = j[["gain_gc", "gain_dn", "gain_neg2"]].mean().max()
        out.append({"check": "shuffled contribution above the largest matched contribution",
                    "value": float(j.kmer_gain.mean() - biggest), "n": len(j)})
        log(f"    -> the shuffled arm's baseline is {lowest - j.comp_auroc.mean():.4f} below "
            f"the lowest matched one and its contribution {j.kmer_gain.mean() - biggest:+.4f} "
            f"above the largest")

        # THE FOUR-PROTOCOL SPAN. The paper's headline is a 5.42-fold span across three
        # matching protocols. Shuffling is in wider use than any of them, and adding it as a
        # fourth point widens the span again. This is what the recommendation rests on: not
        # that one protocol is wrong, but that the quantity is not comparable across the four
        # constructions this literature actually uses.
        arms = [j[f"gain_{a}"].to_numpy(float) for a in ("gc", "dn", "neg2")]
        allarms = arms + [j.kmer_gain.to_numpy(float)]
        means = {a: float(j[f"gain_{a}"].mean()) for a in ("gc", "dn", "neg2")}
        means["shuffled"] = float(j.kmer_gain.mean())
        lo_arm, hi_arm = min(means, key=means.get), max(means, key=means.get)
        b = np.array([max(v[i].mean() for v in allarms) / min(v[i].mean() for v in allarms)
                      for i in draws])
        span = means[hi_arm] / means[lo_arm]
        out.append({"check": "four-protocol span of the 4-mer contribution",
                    "value": float(span), "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(j),
                    "note": f"{hi_arm} over {lo_arm}"})
        three = max(v.mean() for v in arms) / min(v.mean() for v in arms)
        out.append({"check": "three-protocol span, same datasets, for comparison",
                    "value": float(three), "n": len(j)})
        log(f"\n  span across the three MATCHED protocols {three:.2f}x; adding the shuffled "
            f"arm makes it {span:.2f}x ({hi_arm} over {lo_arm})")

    # WHAT SHUFFLING COSTS IN DATA, which is a practical objection to the protocol and is
    # rarely reported: a dinucleotide shuffle of a low-complexity window can come back a near
    # copy of its source, and those pairs have to be discarded.
    d1, d2 = int(t.dropped_failed.sum()), int(t.dropped_similar.sum())
    out.append({"check": "positives discarded, shuffle failed", "value": d1, "n": len(t)})
    out.append({"check": "positives discarded, shuffle too similar to its source",
                "value": d2, "n": len(t)})
    log(f"\n  discarded {d1} positives whose shuffle failed and {d2} whose shuffle came back "
        f"more than 90% identical to its source")

    pd.DataFrame(out).to_csv(TABLES / "shuffled_arm.csv", index=False)
    log("\nwrote shuffled_arm.csv and shuffled_arm_per_dataset.csv")


if __name__ == "__main__":
    main()
