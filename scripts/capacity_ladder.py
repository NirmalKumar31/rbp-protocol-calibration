"""How does the outer-fold channel scale with base-model capacity?

    python scripts/capacity_ladder.py --store ../rbp-store
    python scripts/capacity_ladder.py --from-cache

Specified in `docs/SENSITIVITY_SPEC.md` section C, committed before this was run.

What is being measured. The paper reports the outer-fold channel at two points, k = 4 and
k = 2, and asserts a DIRECTION for how it scales with capacity: a higher-capacity base model has
more scope to encode the withheld fold, so the channel should grow. That assertion carries real
weight, because it is the argument for why the unmeasured neural channel is unlikely to be
smaller than the k-mer one. Two points do not establish a direction. This measures five.

The channel at each rung is `two-stage minus cross-fitted`, both computed by
`scripts/cross_fitting.py` unmodified: ten complement fits per dataset per arm, one per
unordered fold pair.

Why a subsample of datasets, and the rule was fixed before the count. A k = 6 design has 4096
columns against a few thousand rows, so one rung costs far more than k = 4 and the full panel at
five rungs does not fit a sensible local budget. The rule in the specification is systematic
sampling by pair rank, every m-th dataset from the panel ordered by pair count, with m the
smallest integer whose projected runtime fits 90 minutes. That is the same rule
`scripts/select_panel.py` uses to pick the panel itself, and for the same reason: it spans the
size range rather than taking the biggest or the smallest, and dataset size correlates with
AUROC at r = +0.53 to +0.67 so a size-biased subsample would confound the ladder with size.

WHAT THIS CANNOT DO, and it is stated here as well as in the specification. It is a diagnostic
and NOT a substitute for neural cross-fitting. No extrapolation from this ladder to the CNN or
SpliceBERT is made or reported: a k-mer logistic and a fine-tuned transformer differ in more
than the number of parameters. If the channel does not grow with k, the paper's directional
assertion is withdrawn rather than defended.
"""

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rbp.eval.baseline import kmer_matrix  # noqa: E402
from rbp.eval.nested import composition_features  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "capacity_ladder.csv"
PER = TABLES / "capacity_ladder_per_dataset.csv"

ARMS = {"dn": "dinuc", "gc": "gc", "neg2": "neg2"}
LADDER = (2, 3, 4, 5, 6)
BUDGET_SECONDS = 90 * 60          # the specification's cap on the whole run


def rung(seqs, y, folds, k):
    """(two-stage, cross-fitted, composition) contributions for one k, one arm, one dataset."""
    import cross_fitting as cf
    comp, _ = composition_features(seqs, True, standardise_cols=True)
    a_comp = cf._pooled_comp(comp, y, folds)
    X, _ = kmer_matrix(list(seqs), k)
    pub, crossfit = cf.base_scores(X, y, folds)
    naive = {int(i): pub for i in np.unique(folds)}
    return (cf._pooled(comp, naive, y, folds) - a_comp,
            cf._pooled(comp, crossfit, y, folds) - a_comp,
            a_comp)


def sample_datasets(store, budget=BUDGET_SECONDS):
    """Systematic every-m-th by pair rank, with m set by a TIMED projection, not a guess.

    Timing is measured on the two datasets at the extremes of the ordered panel, so the estimate
    covers the cheap and the expensive end rather than one of them. Nothing about a
    CONTRIBUTION value enters the choice of m; only wall-clock seconds do.
    """
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    col = "pairs" if "pairs" in pub.columns else "n"
    order = pub.sort_values([col, "dataset"], kind="mergesort").reset_index(drop=True)
    have = [ds for ds in order.dataset
            if all((Path(store) / "processed" / s / ds.split(":")[1] / ds.split(":")[0]
                    / "dataset.tsv").exists() for s in ARMS.values())]
    if not have:
        sys.exit(f"no dataset.tsv found under {store}/processed; nothing to sample")

    probe = [have[0], have[-1]] if len(have) > 1 else have[:1]
    t0 = time.time()
    for ds in probe:
        protein, cell = ds.split(":")
        for s in ARMS.values():
            d = pd.read_csv(Path(store) / "processed" / s / cell / protein / "dataset.tsv",
                            sep="\t")
            for k in LADDER:
                rung(d.seq_rna.values, d.label.values, d.fold.values, k)
    per_dataset = (time.time() - t0) / len(probe)
    m = max(1, int(np.ceil(len(have) * per_dataset / budget)))
    picked = have[::m]
    log(f"  timing probe: {per_dataset:.1f}s per dataset for all five rungs and three arms")
    log(f"  {len(have)} datasets available, every {m}th fits the {budget // 60}-minute budget "
        f"-> {len(picked)} datasets")
    return picked, m, per_dataset


def build(store, limit=0):
    picked, m, per_dataset = sample_datasets(store)
    if limit:
        picked = picked[:limit]
    rows, skipped = [], 0
    for n, ds in enumerate(picked, 1):
        protein, cell = ds.split(":")
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        for arm, sub in ARMS.items():
            d = pd.read_csv(Path(store) / "processed" / sub / cell / protein / "dataset.tsv",
                            sep="\t")
            y, folds = d.label.values, d.fold.values
            if any(len(np.unique(y[folds == f])) < 2 for f in np.unique(folds)):
                skipped += 1
                continue
            rec[f"rows_{arm}"] = int(len(d))
            for k in LADDER:
                ts, cfv, _ = rung(d.seq_rna.values, y, folds, k)
                rec[f"k{k}_2s_{arm}"] = ts
                rec[f"k{k}_cf_{arm}"] = cfv
                rec[f"k{k}_channel_{arm}"] = ts - cfv
        rows.append(rec)
        log(f"  [{n:3d}/{len(picked)}] {ds:18s} " + "  ".join(
            f"k{k} ch {rec.get(f'k{k}_channel_dn', float('nan')):+.4f}" for k in LADDER))
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to write an empty table over committed evidence")
    return t, m, skipped, per_dataset


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    ap.add_argument("--n", type=int, default=0, help="smoke test; writes .partial.csv")
    ap.add_argument("--from-cache", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")

    meta = None
    if a.from_cache:
        if not PER.exists():
            sys.exit(f"no {PER.name}; run with --store first")
        t = pd.read_csv(PER)
    else:
        t, m, skipped, per_dataset = build(a.store, a.n)
        meta = (m, skipped, per_dataset)
        if a.n:
            p = PER.with_suffix(".partial.csv")
            t.to_csv(p, index=False)
            log(f"\n  --n {a.n}: wrote {p.name}, not the committed table")
            return
        t.to_csv(PER, index=False)

    out = []

    def add(check, value, n="", note=""):
        out.append({"check": check, "value": value, "ci_low": "", "ci_high": "",
                    "n": n, "note": note})

    add("datasets on the ladder", len(t), n=len(t),
        note="systematic every-m-th by pair rank; a subsample, and the rule is in "
             "docs/SENSITIVITY_SPEC.md")
    if meta:
        add("systematic sampling interval m", meta[0], n=len(t),
            note="chosen by a timed projection against a 90-minute budget, not by outcome")
        add("dataset-arms skipped, a fold lacked both classes", meta[1], n=len(t))

    log("")
    channels = {}
    for k in LADDER:
        for arm in ARMS:
            for kind, lab in (("2s", "two-stage"), ("cf", "cross-fitted"),
                              ("channel", "outer-fold channel")):
                col = f"k{k}_{kind}_{arm}"
                if col not in t.columns:
                    continue
                v = t[col].dropna()
                if not len(v):
                    continue
                mean = float(v.mean())
                add(f"panel-mean {lab}, k={k}, {arm} arm", mean, n=len(v))
                if kind == "channel":
                    channels[(k, arm)] = mean
        got = {arm: channels.get((k, arm)) for arm in ARMS}
        if all(v is not None for v in got.values()):
            log(f"  k={k}  channel  dn {got['dn']:+.5f}  gc {got['gc']:+.5f}  "
                f"neg2 {got['neg2']:+.5f}")

    # The direction the paper asserts, turned into a number per arm.
    for arm in ARMS:
        seq = [channels.get((k, arm)) for k in LADDER]
        if any(v is None for v in seq):
            continue
        mono = all(b >= a for a, b in zip(seq, seq[1:]))
        add(f"channel is monotone non-decreasing in k, {arm} arm", float(mono), n=len(t),
            note="the paper asserts the channel grows with capacity; reported either way")
        add(f"channel at k=6 over k=2, {arm} arm",
            seq[-1] / seq[0] if seq[0] > 0 else np.nan, n=len(t),
            note="how far the channel moves across the ladder")

    # And the span at each rung, under both estimators, since that is the paper's headline shape.
    for k in LADDER:
        for kind, lab in (("2s", "two-stage"), ("cf", "cross-fitted")):
            means = {}
            for arm in ARMS:
                col = f"k{k}_{kind}_{arm}"
                if col in t.columns and t[col].notna().any():
                    means[arm] = float(t[col].dropna().mean())
            if len(means) == 3 and min(means.values()) > 0:
                add(f"three-arm span, k={k}, {lab}",
                    max(means.values()) / min(means.values()), n=len(t))

    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {PER.name}")


if __name__ == "__main__":
    main()
