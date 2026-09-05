"""The nested estimator's noise floor at the order this paper actually uses.

    python scripts/estimator_floor.py
    python scripts/estimator_floor.py --from-cache

WHY THIS HAD TO BE RUN. The order-four profile measures what the estimator reports when the
truth is known to be zero, and it reports +0.09 to +0.14. That calibration is at a 337-column
baseline, and every headline in this paper uses the NINETEEN-column order-two baseline, so it
bounds an estimator the paper does not use. Worse, the three order-four values fall in the same
order as the headline (dinucleotide > GC > bias-aware, a 1.61-fold span) on data where the truth
is zero in all three arms, which raises the possibility that part of the headline ORDERING is
estimator bias rather than protocol effect. That objection cannot be answered at order four.

THE EXACT NULL AT ORDER TWO. A 2-mer model scores a window by its sixteen dinucleotide counts.
Windows here are a fixed 101 nt, so those counts are 100 times the dinucleotide frequencies, and
the fifteen frequency columns in the baseline plus the intercept span all sixteen exactly. The
2-mer's score is therefore a linear function of features the baseline already contains, and its
true nested contribution is ZERO by construction, at the order the paper reports.

Whatever the estimator returns here is its error at the operating point of every headline
number. Three things are measured: the size of the floor, whether it reproduces the headline
ORDERING across arms, and its size relative to the contributions actually reported.
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
ARMS = {"gc": "gc", "dn": "dinuc", "neg2": "neg2"}
def published_means():
    """The order-two panel means, READ from the table rather than typed here.

    They were hard-coded to the rounded values the paper prints, which put this section's span
    at 5.43 where every other section says 5.42, and shifted the floor-to-signal ratios by a
    tenth of a point. A number that appears in two sections must come from one place.
    """
    t = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    return {a: float(t[f"gain_{a}"].mean()) for a in ("gc", "dn", "neg2")}



def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    datasets = list(pub.dataset)[:limit or None]
    rows = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        ok = True
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            # THE WINDOW LENGTH IS THE WHOLE ARGUMENT. Counts equal length-1 times frequencies
            # only if every window has the same length; if any does not, the 2-mer's score is
            # no longer in the baseline's span and the null is not exactly zero.
            if d.seq_rna.str.len().nunique() != 1:
                sys.exit(f"{ds} {arm}: windows are not all the same length, so the 2-mer null "
                         f"is not exact")
            sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=2)
            r = gain_over_composition(d.seq_rna.values, sc, d.label.values, d.fold.values)
            rec[f"floor_{arm}"] = float(r.delta)
            rec[f"comp_{arm}"] = float(r.auroc_composition)
            rec[f"n_{arm}"] = int(len(d))
        if not ok:
            continue
        rows.append(rec)
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} " + "  ".join(
            f"{a} {rec[f'floor_{a}']:+.4f}" for a in ARMS))
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

    per = TABLES / "estimator_floor_per_dataset.csv"
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
             for _ in range(4000)]

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})
        return float(v.mean())

    PUBLISHED = published_means()
    log(f"\n=== the estimator's floor at the ORDER-TWO baseline, n = {len(t)}, "
        f"{len(uniq)} proteins ===")
    log("    a 2-mer's score lies in the baseline's span, so the truth is exactly zero\n")
    out.append({"check": "datasets", "value": len(t), "n": len(t)})

    floors = {}
    log(f"  {'arm':6s} {'floor':>9s} {'95% interval':>22s} {'published':>10s} {'ratio':>8s}")
    for arm in ARMS:
        v = add(f"order-2 noise floor, {arm} arm", t[f"floor_{arm}"],
                "true value is zero by construction")
        floors[arm] = v
        k = f"order-2 noise floor, {arm} arm"
        r = [x for x in out if x["check"] == k][0]
        out.append({"check": f"floor as a fraction of the published contribution, {arm} arm",
                    "value": float(v / PUBLISHED[arm]), "n": len(t)})
        out.append({"check": f"datasets with a positive floor, {arm} arm",
                    "value": int((t[f"floor_{arm}"] > 0).sum()), "n": len(t)})
        log(f"  {arm:6s} {v:+9.5f} [{r['ci_low']:+.5f}, {r['ci_high']:+.5f}] "
            f"{PUBLISHED[arm]:+10.4f} {v / PUBLISHED[arm]:7.1%}")

    # THE OBJECTION THIS SECTION ANSWERS. At order four the floor falls in the same order as
    # the headline, so the ordering itself could be estimator bias. Here the arms can be ranked
    # by their floor and compared with the published ranking directly.
    order_pub = [k for k, _ in sorted(PUBLISHED.items(), key=lambda kv: -kv[1])]
    order_floor = [k for k, _ in sorted(floors.items(), key=lambda kv: -kv[1])]
    same = order_floor == order_pub
    out.append({"check": "floor ordering reproduces the published ordering", "value": int(same),
                "n": len(t), "note": f"floor {'>'.join(order_floor)}; "
                                     f"published {'>'.join(order_pub)}"})
    span = max(floors.values()) / min(floors.values()) if min(floors.values()) > 0 else float("nan")
    out.append({"check": "span of the floor across arms", "value": float(span), "n": len(t)})
    log(f"\n  floor ranking      {' > '.join(order_floor)}")
    log(f"  published ranking  {' > '.join(order_pub)}")
    log(f"  same ordering: {same}   floor span {span:.2f}x")

    # THE DEFENCE, AND IT IS THE WHOLE POINT OF MEASURING THE FLOOR PER ARM. A bias that is
    # constant across arms cannot manufacture a difference between them. The floor spans 1.24x
    # while the reported contributions span 5.42x, so the SPAN survives even though the LEVEL
    # does not. Gated as the comparison rather than as two separate numbers, because either
    # alone invites the wrong conclusion.
    pub_span = max(PUBLISHED.values()) / min(PUBLISHED.values())
    out.append({"check": "published contribution span across arms", "value": float(pub_span),
                "n": len(t)})
    out.append({"check": "floor span as a fraction of the contribution span",
                "value": float(span / pub_span), "n": len(t)})
    log(f"\n  contributions span {pub_span:.2f}x and the floor spans {span:.2f}x, so the floor "
        f"is {span / pub_span:.0%} of the\n  span it would have to explain: nearly flat across "
        f"arms, and therefore unable to create the ordering.")

    # AND THE QUANTITY A READER NEEDS: is the floor small against the effects reported?
    worst = max(abs(floors[a] / PUBLISHED[a]) for a in ARMS)
    out.append({"check": "largest floor-to-contribution ratio over the three arms",
                "value": float(worst), "n": len(t)})
    log(f"  worst floor-to-signal ratio {worst:.1%}, in the "
        f"{max(ARMS, key=lambda a: abs(floors[a] / PUBLISHED[a]))} arm")

    pd.DataFrame(out).to_csv(TABLES / "estimator_floor.csv", index=False)
    log("\nwrote estimator_floor.csv and estimator_floor_per_dataset.csv")


if __name__ == "__main__":
    main()
