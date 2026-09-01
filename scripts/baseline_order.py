"""R1o: how much of "what the model contributes" is just one order of composition up?

    python scripts/baseline_order.py --n 30

THE OBJECTION, AND IT IS THE SHARPEST ONE ANYONE MADE ABOUT WHAT THE QUANTITY IS. A bag of
4-mer counts carries NO positional information. The composition baseline stops at order 2. So
the headline "nested contribution of a sequence model over composition" is, definitionally,
order-3-and-4 composition beyond order-1-and-2 composition -- and where the baseline stops is
an arbitrary choice doing all the work. `nested.py` defends stopping at order 2 as "exactly
what a dinucleotide shuffle would have preserved", which is a good reason to choose it and not
a reason to believe the choice does not matter.

Nobody had measured it. This does: refit the baseline at order 3 (mono + di + tri, 84
features) and recompute the same nested contribution.

WHAT IT COSTS THE PAPER. Roughly half of every arm's contribution and most of the R1 contrast.
The word "sequence model" has to be qualified: what is measured is short-range compositional
signal one order above wherever the matcher stopped, not motif recognition and not anything
positional.

WHAT IT DOES NOT COST. The FOLD RANGE, which is the paper's actual claim. Everything shrinks
roughly proportionally, so the protocol dependence is not an artefact of where the baseline
stops. That is the result worth having: the magnitude is a function of an analyst's choice,
the protocol dependence is not.

Subsampled by default because the order-3 fit is 82 features x 5 folds x 3 arms per dataset.
The subsample is size-stratified and its own R1 contrast is reported next to the published one
so the reader can see whether the subsample is representative.

THE BASELINE HAS TO BE THE PAPER'S BASELINE. An earlier version of this script built its own
composition block -- all 4 mono + all 16 dinucleotide frequencies, unstandardised, fed to
sklearn's default LogisticRegression. That is three departures from `nested.composition_features`
at once: a singular design matrix, no entropy column, and an L2 penalty applied to features on a
0.06 scale, which crushes the fit. It underfits the baseline and therefore OVERSTATES the gain,
by 1.22x on the GC arm and **2.31x on neg2**, where the published baseline is highest and has
the most to lose. The order-2/order-3 comparison was internally consistent, but the arm-level
numbers were not the paper's quantities and the fold range was not the paper's fold range.

So the composition block here is `nested.composition_features` generalised to arbitrary order
-- one column dropped per frequency family, entropy appended, every column standardised -- fit
with the same `_oof_scores` and compared with the same DeLong test. At order 2 it is that
function, so `gain2_*` MUST reproduce `three_arm_per_dataset.csv`. That is asserted per dataset
rather than hoped for, because "I re-implemented the baseline" is exactly the class of error it
was written to catch.
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
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import _counts, _oof_scores, entropy  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "processed/gc", "dn": "processed/dinuc", "neg2": "processed/neg2"}
REPRO_TOL = 5e-3


def log(m):
    print(m, flush=True)


def composition(seqs, max_k):
    """`nested.composition_features`, generalised past order 2. Identical at max_k=2."""
    seqs = list(seqs)
    cols = []
    for k in range(1, max_k + 1):
        C = _counts(seqs, k).astype(float)
        C /= np.maximum(C.sum(axis=1, keepdims=True), 1)
        cols.append(C[:, :-1])
    cols.append(entropy(seqs)[:, None])
    X = np.column_stack(cols)
    return np.column_stack([standardise(X[:, j]) for j in range(X.shape[1])])


def gain(X, score, y, folds):
    """DeLong gain of `score` over design `X`, out-of-fold on the study's own folds."""
    full = np.column_stack([X, standardise(score)])
    s_comp = _oof_scores(X, y, folds)
    s_full = _oof_scores(full, y, folds)
    ok = np.isfinite(s_comp) & np.isfinite(s_full)
    r = delong_test(s_full[ok], s_comp[ok], y[ok])
    return float(r["diff"]), float(r["auc_b"])  # auc_a is the FULL model, auc_b the baseline


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=30, help="size-stratified subsample")
    p.add_argument("--from-cache", action="store_true",
                   help="rebuild the summary from the committed per-dataset table. THE DEFAULT "
                        "PATH IN run.sh, because the window store is not published and running "
                        "without it on a clean clone silently produced an EMPTY table over the "
                        "committed evidence and then crashed.")
    a = p.parse_args()
    warnings.filterwarnings("ignore")
    store = Path(a.store)

    per = TABLES / "baseline_order_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it with --store")
    else:
        # FAIL LOUDLY RATHER THAN SILENTLY. The window store is ~3 GB and is not published, and
        # missing files used to be skipped one by one -- so on a clean clone this wrote an EMPTY
        # table over the committed evidence and then crashed on a zero-length bootstrap.
        if not (store / "processed" / "gc").exists():
            sys.exit(f"no window store at {store}. It is ~3 GB and is not published; use "
                     f"--from-cache to rebuild the summary from the committed per-dataset "
                     f"table instead.")
        panel = pd.read_csv(TABLES / "rehearsal_binding_gc.csv").sort_values("pairs")
        step = max(len(panel) // a.n, 1)
        sub = panel.iloc[::step].head(a.n)
        pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv").set_index("dataset")

        rows, worst = [], 0.0
        for i, r in enumerate(sub.itertuples(), 1):
            rec, ok = {"dataset": r.dataset, "protein": r.protein, "cell": r.cell}, True
            for arm, sd in ARMS.items():
                f = store / sd / r.cell / r.protein / "dataset.tsv"
                if not f.exists():
                    ok = False
                    break
                d = pd.read_csv(f, sep="\t")
                sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
                m = ~np.isnan(sc)
                seqs, y, fo = d.seq_rna.values[m], d.label.values[m], d.fold.values[m]
                for order in (2, 3):
                    g, comp = gain(composition(seqs, order), sc[m], y, fo)
                    rec[f"gain{order}_{arm}"] = g
                    rec[f"comp{order}_{arm}"] = comp
                # ANCHOR: order 2 IS the published baseline, so it must reproduce the
                # published gain. A re-implemented baseline that drifts is the bug this
                # script shipped with; it is now a per-dataset assertion.
                if r.dataset in pub.index:
                    ref = float(pub.loc[r.dataset, f"gain_{arm}"])
                    rec[f"published_{arm}"] = ref
                    worst = max(worst, abs(rec[f"gain2_{arm}"] - ref))
            if ok:
                rows.append(rec)
                log(f"[{i:3d}/{len(sub)}] {r.dataset:18s} "
                    + "  ".join(f"{k} {rec[f'gain2_{k}']:+.4f}->{rec[f'gain3_{k}']:+.4f}"
                                for k in ARMS))
        t = pd.DataFrame(rows)
        if t.empty:
            sys.exit("no dataset could be built; refusing to overwrite the committed table")
        log(f"\norder-2 vs published, max |difference| = {worst:.2e}")
        if worst > REPRO_TOL:
            sys.exit(f"order-2 gain does not reproduce three_arm_per_dataset.csv "
                     f"(max {worst:.2e} > {REPRO_TOL}). The composition block here is no "
                     f"longer the paper's, and every number below would be a different "
                     f"quantity. Refusing to write.")
        t.to_csv(per, index=False)

    out = []
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(t), size=(2000, len(t)))

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in idx])
        lo, hi = np.percentile(b, [2.5, 97.5])
        out.append({"check": check, "value": float(v.mean()), "ci_low": float(lo),
                    "ci_high": float(hi), "n": len(t), "note": note})
        return float(v.mean())

    log(f"\n=== R1o: baseline order, n = {len(t)} size-stratified datasets ===\n")

    # The anchor, re-asserted on every run including --from-cache: order 2 is the paper's
    # baseline, so it reproduces the paper's gain or nothing here is comparable to the paper.
    if all(f"published_{arm}" in t for arm in ARMS):
        w = max(float((t[f"gain2_{arm}"] - t[f"published_{arm}"]).abs().max()) for arm in ARMS)
        out.append({"check": "max |order-2 gain - published gain|", "value": w, "n": len(t)})
        log(f"  order-2 reproduces three_arm_per_dataset.csv to {w:.2e}")
        if w > REPRO_TOL:
            sys.exit(f"order-2 gain does not reproduce the published gain ({w:.2e})")
    else:
        sys.exit("no published_* anchor columns; regenerate with --store")

    log(f"\n  {'arm':6s} {'comp o2':>8s} {'comp o3':>8s} "
        f"{'over mono+di':>14s} {'over mono+di+tri':>18s} {'removed':>9s}")
    for arm in ARMS:
        add(f"composition AUROC, order-2 baseline, {arm} arm", t[f"comp2_{arm}"])
        add(f"composition AUROC, order-3 baseline, {arm} arm", t[f"comp3_{arm}"])
        g2 = add(f"gain over order-2 baseline, {arm} arm", t[f"gain2_{arm}"])
        g3 = add(f"gain over order-3 baseline, {arm} arm", t[f"gain3_{arm}"])
        frac = 1 - g3 / g2 if g2 else np.nan
        out.append({"check": f"fraction removed by order 3, {arm} arm", "value": float(frac),
                    "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
        log(f"  {arm:6s} {t[f'comp2_{arm}'].mean():8.4f} {t[f'comp3_{arm}'].mean():8.4f} "
            f"{g2:+14.4f} {g3:+18.4f} {100 * frac:8.0f}%")

    c2 = add("R1 contrast (dn-gc), order-2 baseline", t.gain2_dn - t.gain2_gc)
    c3 = add("R1 contrast (dn-gc), order-3 baseline", t.gain3_dn - t.gain3_gc)
    out.append({"check": "fraction of the R1 contrast removed by order 3",
                "value": float(1 - c3 / c2), "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
    log(f"\n  R1 contrast  order-2 {c2:+.4f}   order-3 {c3:+.4f}   "
        f"{100 * (1 - c3 / c2):.0f}% removed")

    # THE POINT: the range is what the paper claims, and it must survive.
    for order in (2, 3):
        means = [t[f"gain{order}_{arm}"].mean() for arm in ARMS]
        fr = max(means) / min(means) if min(means) > 0 else np.nan
        out.append({"check": f"fold range across protocols, order-{order} baseline",
                    "value": float(fr), "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
        log(f"  fold range across protocols, order-{order} baseline: {fr:.2f}x")

    pd.DataFrame(out).to_csv(TABLES / "baseline_order.csv", index=False)
    log("\nwrote baseline_order.csv and baseline_order_per_dataset.csv")


if __name__ == "__main__":
    main()
