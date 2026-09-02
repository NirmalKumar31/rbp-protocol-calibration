"""R1r: is the order-3 collapse a property of the 4-mer, or of every model class?

    python scripts/baseline_order_models.py

R1o found that raising the composition baseline from order 2 to order 3 removes ~80% of the
4-mer's measured contribution. That is close to tautological for a bag-of-4-mers model: its
only information beyond trinucleotide frequency IS order 4. The question that matters for the
paper is whether a model with positional structure survives the same baseline. If SpliceBERT
also collapses, the honest reading is that everything this literature measures as "what the
model adds over composition" is one order of composition. If it does not, then the order-3
baseline is a sharper instrument than the order-2 one and the paper can say which models it
separates.

This runs R1o's design on R1g's evidence: three model classes, both arms, 94 datasets, the
per-window scores already committed under data/evidence. No GPU, nothing refitted on the deep
side -- the CNN and SpliceBERT scores are exactly the ones R1g used.

THE ROW SET AND THE ANCHOR. Rows are intersected across all three models before anything is
fitted, by importing deep_model_contrast's own loader rather than reimplementing it, so this
table's order-2 column IS that table's published column. That is asserted per cell: if the
order-2 gain here does not reproduce deep_contrast_per_dataset.csv to 5e-3, the composition
block has drifted and every order-3 number would be measuring something else.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from baseline_order import REPRO_TOL, composition  # noqa: E402
from deep_model_contrast import MIN_COVERAGE, MODELS, arm_roots, oof  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import _oof_scores  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"
ORDERS = (2, 3)


def log(m):
    print(m, flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0, help="limit datasets, for smoke tests")
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "baseline_order_models_per_dataset.csv"
    pub = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv").set_index("dataset")

    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it with --store")
    else:
        roots = arm_roots(a.store)
        datasets = list(pub.index)[: a.n or None]
        rows, worst = [], 0.0
        for i, ds in enumerate(datasets, 1):
            protein, cell = ds.split(":")
            rec, ok = {"dataset": ds, "protein": protein, "cell": cell}, True
            for arm, (dataroot, scoreroot) in roots.items():
                f = dataroot / cell / protein / "dataset.tsv"
                if not f.exists():
                    ok = False
                    break
                d = pd.read_csv(f, sep="\t")
                ids, got = set(d.id), {}
                for model in MODELS:
                    if model == "kmer":
                        continue
                    s = oof(scoreroot, cell, protein, model)
                    if s is None:
                        ok = False
                        break
                    got[model] = s
                    ids &= set(s.id)
                if not ok:
                    break
                if len(ids) / len(d) < MIN_COVERAGE:
                    ok = False
                    break
                dd = d[d.id.isin(ids)].reset_index(drop=True)
                sc, _, _ = kmer_oof(dd.seq_rna.values, dd.label.values, dd.fold.values, k=4)
                got["kmer"] = pd.DataFrame({"id": dd.id.values, "score": sc})

                y, fo = dd.label.values, dd.fold.values
                for order in ORDERS:
                    X = composition(dd.seq_rna.values, order)
                    s_comp = _oof_scores(X, y, fo)
                    for model in MODELS:
                        m = dd.merge(got[model], on="id", how="inner")
                        s_full = _oof_scores(
                            np.column_stack([X, standardise(m.score.values)]), y, fo)
                        good = np.isfinite(s_comp) & np.isfinite(s_full)
                        r = delong_test(s_full[good], s_comp[good], y[good])
                        rec[f"{model}_gain{order}_{arm}"] = float(r["diff"])
                        rec[f"comp{order}_{arm}"] = float(r["auc_b"])
                # ANCHOR: order 2 is R1g's own baseline, so it is R1g's own number.
                for model in MODELS:
                    ref = float(pub.loc[ds, f"{model}_gain_{arm}"])
                    rec[f"{model}_published_{arm}"] = ref
                    worst = max(worst, abs(rec[f"{model}_gain2_{arm}"] - ref))
            if ok:
                rows.append(rec)
                log(f"[{i:3d}/{len(datasets)}] {ds:18s} " + "  ".join(
                    f"{m[:4]} {rec[f'{m}_gain2_gc']:+.4f}->{rec[f'{m}_gain3_gc']:+.4f}"
                    for m in MODELS))
        t = pd.DataFrame(rows)
        if t.empty:
            sys.exit("no dataset could be built; refusing to overwrite the committed table")
        log(f"\norder-2 vs deep_contrast_per_dataset.csv, max |difference| = {worst:.2e}")
        if worst > REPRO_TOL:
            sys.exit(f"order-2 gain does not reproduce R1g ({worst:.2e} > {REPRO_TOL}); "
                     f"the composition block has drifted. Refusing to write.")
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
        lo, hi = np.percentile(b, [2.5, 97.5])
        out.append({"check": check, "value": float(v.mean()), "ci_low": float(lo),
                    "ci_high": float(hi), "n": len(t), "note": note})
        return float(v.mean())

    w = max(float((t[f"{m}_gain2_{arm}"] - t[f"{m}_published_{arm}"]).abs().max())
            for m in MODELS for arm in ("gc", "dn"))
    out.append({"check": "max |order-2 gain - R1g published gain|", "value": w, "n": len(t)})
    if w > REPRO_TOL:
        sys.exit(f"order-2 does not reproduce R1g ({w:.2e})")

    log(f"\n=== R1r: baseline order across model classes, n = {len(t)}, "
        f"{len(uniq)} proteins ===")
    log(f"  order-2 reproduces deep_contrast_per_dataset.csv to {w:.2e}\n")
    for arm in ("gc", "dn"):
        add(f"composition AUROC, order-2 baseline, {arm} arm", t[f"comp2_{arm}"])
        add(f"composition AUROC, order-3 baseline, {arm} arm", t[f"comp3_{arm}"])
        log(f"  {arm} arm, composition {t[f'comp2_{arm}'].mean():.4f} -> "
            f"{t[f'comp3_{arm}'].mean():.4f}")
        log(f"  {'model':12s} {'over order 2':>13s} {'over order 3':>13s} {'surviving':>10s}")
        for model in MODELS:
            g2 = add(f"{model} gain over order-2 baseline, {arm} arm", t[f"{model}_gain2_{arm}"])
            g3 = add(f"{model} gain over order-3 baseline, {arm} arm", t[f"{model}_gain3_{arm}"])
            # RATIO OF MEANS, bootstrapped as a ratio of means. Per-dataset ratios are not
            # usable here: a handful of datasets have an order-2 gain near zero, and the mean
            # of ratios then diverges -- the first version of this printed a lower bound of
            # -0.64 for a quantity whose point estimate is +0.22, and a NaN for the CNN.
            v2 = t[f"{model}_gain2_{arm}"].to_numpy(dtype=float)
            v3 = t[f"{model}_gain3_{arm}"].to_numpy(dtype=float)
            b = np.array([v3[i].mean() / v2[i].mean() for i in draws])
            lo, hi = np.percentile(b, [2.5, 97.5])
            out.append({"check": f"{model} fraction surviving order 3, {arm} arm",
                        "value": float(g3 / g2), "ci_low": float(lo), "ci_high": float(hi),
                        "n": len(t), "note": "ratio of means, protein-clustered"})
            log(f"  {model:12s} {g2:+13.4f} {g3:+13.4f} {100 * g3 / g2:9.0f}%")
        log("")

    # THE SHARE IS A RATIO, AND A RATIO HAS A DENOMINATOR. Reporting only "the k-mer keeps 22%
    # and SpliceBERT keeps 75%" invites the reading that the k-mer is uniquely FRAGILE. It is
    # not: the order-3 baseline absorbs almost exactly the same ABSOLUTE amount from every
    # model, and the shares differ because the totals differ. That has to be measured and
    # printed next to the shares, or the shares are misleading.
    log("  absolute amount absorbed by raising the baseline one order:")
    for arm in ("gc", "dn"):
        for model in MODELS:
            a = t[f"{model}_gain2_{arm}"] - t[f"{model}_gain3_{arm}"]
            v = add(f"{model} absolute absorbed by order 3, {arm} arm", a)
            log(f"    {arm} {model:11s} {v:+.4f}")
        d = ((t[f"splicebert_gain2_{arm}"] - t[f"splicebert_gain3_{arm}"])
             - (t[f"kmer_gain2_{arm}"] - t[f"kmer_gain3_{arm}"]))
        add(f"splicebert minus kmer, absolute absorbed, {arm} arm", d)

    # AND THE PER-DATASET SIGN, which is what the share cannot show. A mean of +0.0058 is
    # consistent with "small everywhere" and with "positive in two thirds, negative in a
    # third". For the k-mer over a trinucleotide baseline it is the second.
    for arm in ("gc", "dn"):
        for model in MODELS:
            n_pos = int((t[f"{model}_gain3_{arm}"] > 0).sum())
            out.append({"check": f"{model} order-3 gain positive in, {arm} arm",
                        "value": n_pos, "n": len(t)})
        log(f"  {arm} arm, order-3 gain positive in: " + "  ".join(
            f"{m} {int((t[f'{m}_gain3_{arm}'] > 0).sum())}/{len(t)}" for m in MODELS))

    # THE QUESTION: does a model with positional structure survive what a bag of 4-mers does not?
    for arm in ("gc", "dn"):
        base = t[f"kmer_gain3_{arm}"]
        for model in ("cnn", "splicebert"):
            d = t[f"{model}_gain3_{arm}"] - base
            add(f"{model} minus kmer, order-3 baseline, {arm} arm", d)
        r = (t["splicebert_gain3_" + arm].mean() / t["kmer_gain3_" + arm].mean())
        out.append({"check": f"splicebert / kmer at order 3, {arm} arm", "value": float(r),
                    "n": len(t)})
        r2 = (t["splicebert_gain2_" + arm].mean() / t["kmer_gain2_" + arm].mean())
        out.append({"check": f"splicebert / kmer at order 2, {arm} arm", "value": float(r2),
                    "n": len(t)})
        log(f"  {arm} arm: SpliceBERT/k-mer ratio {r2:.2f}x at order 2 -> {r:.2f}x at order 3")

    for order in ORDERS:
        for model in MODELS:
            c = t[f"{model}_gain{order}_dn"] - t[f"{model}_gain{order}_gc"]
            add(f"{model} R1 contrast (dn-gc), order-{order} baseline", c)
            # The ratio scale, because R1g's "the contrast grows with capacity" was withdrawn
            # for reversing here and must not be resurrected from the order-3 table either.
            mult = t[f"{model}_gain{order}_dn"].mean() / t[f"{model}_gain{order}_gc"].mean()
            out.append({"check": f"{model} dn/gc multiplier, order-{order} baseline",
                        "value": float(mult), "n": len(t)})
        log("  order-%d multipliers: " % order + "  ".join(
            f"{m} {t[f'{m}_gain{order}_dn'].mean() / t[f'{m}_gain{order}_gc'].mean():.2f}x"
            for m in MODELS))

    pd.DataFrame(out).to_csv(TABLES / "baseline_order_models.csv", index=False)
    log("\nwrote baseline_order_models.csv and baseline_order_models_per_dataset.csv")


if __name__ == "__main__":
    main()
