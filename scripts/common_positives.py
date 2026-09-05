"""Does the protocol effect survive restricting both arms to the positives they share?

    python scripts/common_positives.py --store ../rbp-store
    python scripts/common_positives.py --from-cache

THE OBJECTION. The design claim is that only the negatives change. It is very nearly true and
not exactly true: a positive is retained only when its matcher finds an acceptable negative, the
GC and dinucleotide matchers fail on different windows, and the two arms end up with positive
sets whose median Jaccard is 0.9972 and whose minimum is 0.9164
(scripts/positive_set_overlap.py). An external review put it plainly: "almost the only
difference" is more accurate than "only", and a qualitative reason for declining a
common-positive analysis is not the analysis.

The Discussion's reason for declining was real and is worth keeping: the windows the matchers
disagree on are not a random sample. A matcher fails where a positive's composition is extreme,
so intersecting discards the hardest windows preferentially, and the intersection answers a
slightly different question than the panel does. That is an argument for reporting BOTH, not for
reporting neither.

WHAT THIS DOES. For each dataset, take the positives present in both composition-matched arms,
keyed by genomic coordinate, keep each arm's own matched negative for those positives, and
recompute the nested contribution. Both arms then run on identical positive sets, so the
negatives are exactly the only difference. The comparison of interest is whether the
dinucleotide-over-GC contrast holds, and by how much it moves.

The 4-mer only, because the neural arms' per-window scores are committed for the published
positive sets and rescoring an intersection means retraining. That limitation is the same one
scripts/cross_fitting.py has and is stated for the same reason.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "gc", "dn": "dinuc"}
KEY = ["chrom", "start", "end"]


def restrict(d, keep):
    """Rows of one arm whose POSITIVE is in `keep`, with that positive's paired negative.

    Negatives are matched 1:1 to positives, so dropping a positive must drop its partner or the
    classes stop being balanced and the AUROC is no longer comparable to the published one. The
    pairing is by row order within the file, which is how the matcher wrote it: positive i and
    negative i are partners.
    """
    pos = d[d.label == 1].reset_index(drop=True)
    neg = d[d.label == 0].reset_index(drop=True)
    if len(pos) != len(neg):
        return None
    sel = pd.MultiIndex.from_frame(pos[KEY]).isin(keep)
    return pd.concat([pos[sel], neg[sel]], ignore_index=True)


def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    rows = []
    for n, ds in enumerate(list(pub.dataset)[:limit or None], 1):
        protein, cell = ds.split(":")
        d = {}
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                d = None
                break
            d[arm] = pd.read_csv(f, sep="\t")
        if not d:
            continue
        keys = [pd.MultiIndex.from_frame(d[a][d[a].label == 1][KEY]) for a in ARMS]
        shared = keys[0].intersection(keys[1])
        rec = {"dataset": ds, "protein": protein, "cell": cell,
               "n_shared_positives": len(shared)}
        ok = True
        for arm in ARMS:
            # BOTH SIDES GO THROUGH restrict(), including the one that drops nothing.
            # restrict() returns all positives then all negatives, which is not the file's own
            # row order, and reordering alone moves the gain by about 4e-6 through
            # floating-point summation order in the logistic fit. Differencing the restricted
            # run against the PUBLISHED value would fold that into the answer -- it showed up as
            # a nonzero shift on 83 datasets when only 23 had lost a positive. Differencing two
            # runs that are identical but for the intersection leaves only the intersection.
            full = restrict(d[arm], pd.MultiIndex.from_frame(d[arm][d[arm].label == 1][KEY]))
            r = restrict(d[arm], shared)
            if r is None or r.label.nunique() < 2 or len(r) < 200:
                ok = False
                break
            for tag, x in (("", r), ("_full", full)):
                sc, _, _ = kmer_oof(x.seq_rna.values, x.label.values, x.fold.values, k=4)
                g = gain_over_composition(x.seq_rna.values, sc, x.label.values, x.fold.values)
                rec[f"gain{tag}_{arm}"] = float(g.delta)
                rec[f"comp{tag}_{arm}"] = float(g.auroc_composition)
            rec[f"n_{arm}"] = int(len(r))
            rec[f"n_full_{arm}"] = int(len(full))
        if not ok:
            continue
        rows.append(rec)
        log(f"[{n:3d}/94] {ds:18s} shared {len(shared):6d}  "
            f"gc {rec['gain_gc']:+.4f}  dn {rec['gain_dn']:+.4f}")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "common_positives_per_dataset.csv"
    t = pd.read_csv(per) if a.from_cache else build(a.store, a.n)
    if not a.from_cache:
        t.to_csv(per, index=False)

    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv").set_index("dataset")
    t = t.set_index("dataset")
    common = t.index.intersection(pub.index)
    t, pub = t.loc[common], pub.loc[common]

    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(4000)]
    out = []

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})

    add("datasets", np.full(len(t), len(t)))
    add("positives dropped by intersecting, fraction of the GC arm's",
        1 - t.n_shared_positives.to_numpy() / (pub.n_gc.to_numpy() / 2)
        if "n_gc" in pub else np.zeros(len(t)))
    for arm in ARMS:
        add(f"contribution on common positives, {arm} arm", t[f"gain_{arm}"])
        add(f"contribution on all retained positives, {arm} arm", t[f"gain_full_{arm}"])
        add(f"shift from intersecting, {arm} arm",
            t[f"gain_{arm}"].to_numpy() - t[f"gain_full_{arm}"].to_numpy())
        # The published value, carried alongside so the reconstruction can be checked against
        # it, but NOT the thing the shift is measured from.
        add(f"contribution as published, {arm} arm", pub[f"gain_{arm}"])
    add("dinucleotide-over-GC contrast on common positives",
        t.gain_dn.to_numpy() - t.gain_gc.to_numpy())
    add("dinucleotide-over-GC contrast on all retained positives",
        t.gain_full_dn.to_numpy() - t.gain_full_gc.to_numpy())

    ratio_c = t.gain_dn.mean() / t.gain_gc.mean()
    ratio_f = t.gain_full_dn.mean() / t.gain_full_gc.mean()
    out.append({"check": "dn/gc ratio on common positives", "value": float(ratio_c),
                "ci_low": "", "ci_high": "", "n": len(t), "note": ""})
    out.append({"check": "dn/gc ratio on all retained positives", "value": float(ratio_f),
                "ci_low": "", "ci_high": "", "n": len(t),
                "note": "the comparison of record; both sides processed identically"})
    out.append({"check": "datasets where the contrast stays positive on common positives",
                "value": float((t.gain_dn > t.gain_gc).sum()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "the published count is 88 of 94"})

    r = pd.DataFrame(out)
    r.to_csv(TABLES / "common_positives.csv", index=False)
    log("")
    for _, x in r.iterrows():
        ci = f"  [{x.ci_low:+.4f}, {x.ci_high:+.4f}]" if isinstance(x.ci_low, float) else ""
        log(f"  {x['check']:60s} {x['value']:+.4f}{ci}")
    log("\nwrote common_positives.csv and common_positives_per_dataset.csv")


if __name__ == "__main__":
    main()
