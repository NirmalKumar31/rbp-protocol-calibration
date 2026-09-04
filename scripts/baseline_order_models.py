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

This runs R1o's design on R1g's evidence: three model classes, all three arms, 94 datasets,
the
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

    # THE THIRD ARM WAS COMPUTED AND NEVER REPORTED. Every loop below used to read
    # ("gc", "dn") while the per-dataset table carried neg2 columns for all three models, so
    # the bias-aware arm's order-3 numbers existed on disk and appeared in nothing. Read the
    # arm set off the table instead, and require the two the paper's headline needs.
    ARMS_HERE = [x for x in ("gc", "dn", "neg2") if f"comp2_{x}" in t.columns]
    for need in ("gc", "dn"):
        if need not in ARMS_HERE:
            sys.exit(f"{per} has no {need} arm; refusing to report a partial table")

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
            for m in MODELS for arm in ARMS_HERE)
    out.append({"check": "max |order-2 gain - R1g published gain|", "value": w, "n": len(t)})
    if w > REPRO_TOL:
        sys.exit(f"order-2 does not reproduce R1g ({w:.2e})")

    log(f"\n=== R1r: baseline order across model classes, n = {len(t)}, "
        f"{len(uniq)} proteins ===")
    log(f"  order-2 reproduces deep_contrast_per_dataset.csv to {w:.2e}\n")
    for arm in ARMS_HERE:
        add(f"composition AUROC, order-2 baseline, {arm} arm", t[f"comp2_{arm}"])
        add(f"composition AUROC, order-3 baseline, {arm} arm", t[f"comp3_{arm}"])
        # The headroom the order-3 baseline takes away. This is what the compression
        # correction below is correcting for, so it has to be a committed number.
        add(f"composition rise from order 3, {arm} arm", t[f"comp3_{arm}"] - t[f"comp2_{arm}"])
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
    for arm in ARMS_HERE:
        for model in MODELS:
            a = t[f"{model}_gain2_{arm}"] - t[f"{model}_gain3_{arm}"]
            v = add(f"{model} absolute absorbed by order 3, {arm} arm", a)
            log(f"    {arm} {model:11s} {v:+.4f}")
        d = ((t[f"splicebert_gain2_{arm}"] - t[f"splicebert_gain3_{arm}"])
             - (t[f"kmer_gain2_{arm}"] - t[f"kmer_gain3_{arm}"]))
        add(f"splicebert minus kmer, absolute absorbed, {arm} arm", d)

    # THE COMPRESSION CORRECTION, WHICH THIS SCRIPT OMITTED. deep_model_contrast.py calls it
    # "not optional" because a nested gain is bounded above by 1 - baseline, and the order-3
    # baseline RAISES composition (+0.0196 gc, +0.0517 dn), cutting the headroom every model is
    # working in. So part of the "absorbed" amount is arithmetic, not absorption -- and it is a
    # LARGER part for a model with a bigger gain, which is exactly the asymmetry that made the
    # raw absorption look constant across model classes.
    #
    # Transplant each model's order-2 d' increment onto the order-3 baseline, the paper's own
    # estimator, and the residual is what the moving ceiling does NOT explain.
    from scipy.stats import norm
    r2 = np.sqrt(2.0)

    def dprime(a):
        return r2 * norm.ppf(np.clip(np.asarray(a, dtype=float), 1e-6, 1 - 1e-6))

    log("  compression-corrected: residual after transplanting the order-2 d' increment")
    for arm in ARMS_HERE:
        res = {}
        for model in MODELS:
            c2, c3 = t[f"comp2_{arm}"], t[f"comp3_{arm}"]
            inc = dprime(c2 + t[f"{model}_gain2_{arm}"]) - dprime(c2)
            pred3 = norm.cdf((dprime(c3) + inc) / r2) - c3
            r = t[f"{model}_gain3_{arm}"] - pred3
            add(f"{model} compression-predicted order-3 gain, {arm} arm", pred3)
            v = add(f"{model} compression-corrected order-3 residual, {arm} arm", r)
            res[model] = r
            log(f"    {arm} {model:11s} predicted {pred3.mean():+.4f}  observed "
                f"{t[f'{model}_gain3_{arm}'].mean():+.4f}  residual {v:+.4f}")
        # THE COMPARISON THAT REVERSES THE RAW FINDING. Corrected for the moving ceiling, the
        # k-mer loses MORE than SpliceBERT -- so "absorption is a constant" is itself partly an
        # artefact of omitting the correction, just as the shares were an artefact of their
        # denominators. The truth is between the two framings.
        d = res["kmer"] - res["splicebert"]
        add(f"kmer minus splicebert, compression-corrected residual, {arm} arm", d)
        ratio = float(res["kmer"].mean() / res["splicebert"].mean())
        out.append({"check": f"kmer/splicebert corrected-residual ratio, {arm} arm",
                    "value": ratio, "n": len(t)})
        log(f"    {arm} -> k-mer loses {ratio:.2f}x what SpliceBERT does once the moving "
            f"ceiling is removed")

    # AND THE ABSORBED SPREAD ON THE DINUCLEOTIDE ARM, which the first gate left ungated.
    for arm in ARMS_HERE:
        vals = [float((t[f"{m}_gain2_{arm}"] - t[f"{m}_gain3_{arm}"]).mean()) for m in MODELS]
        out.append({"check": f"absorbed spread across model classes, {arm} arm",
                    "value": float(max(vals) - min(vals)), "n": len(t)})
    for arm in ARMS_HERE:
        d = ((t[f"cnn_gain2_{arm}"] - t[f"cnn_gain3_{arm}"])
             - (t[f"kmer_gain2_{arm}"] - t[f"kmer_gain3_{arm}"]))
        add(f"cnn minus kmer, absolute absorbed, {arm} arm", d)

    # AND THE PER-DATASET SIGN, which is what the share cannot show. A mean of +0.0058 is
    # consistent with "small everywhere" and with "positive in two thirds, negative in a
    # third". For the k-mer over a trinucleotide baseline it is the second.
    for arm in ARMS_HERE:
        for model in MODELS:
            n_pos = int((t[f"{model}_gain3_{arm}"] > 0).sum())
            out.append({"check": f"{model} order-3 gain positive in, {arm} arm",
                        "value": n_pos, "n": len(t)})
        log(f"  {arm} arm, order-3 gain positive in: " + "  ".join(
            f"{m} {int((t[f'{m}_gain3_{arm}'] > 0).sum())}/{len(t)}" for m in MODELS))

    # THE QUESTION: does a model with positional structure survive what a bag of 4-mers does not?
    for arm in ARMS_HERE:
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

    # THE THREE-ARM SPAN AT EACH ORDER, on the full panel. The Discussion's 7.16 came from
    # baseline_order.py: 4-mer only, 30 size-stratified datasets, and TWO arms plus a
    # bias-aware column read off a different row set. Now that all three arms are here for all
    # three models, the span is the paper's own headline quantity recomputed one order up.
    # Print the denominator and its positive count beside every span, because a span whose
    # denominator is +0.0019 is a near-zero-denominator artefact and not a magnitude.
    if "neg2" in ARMS_HERE:
        log("\n  three-arm span of the contribution, per order:")
        for order in ORDERS:
            for model in MODELS:
                m = {arm: float(t[f"{model}_gain{order}_{arm}"].mean()) for arm in ARMS_HERE}
                lo_arm = min(m, key=m.get)
                span = m[max(m, key=m.get)] / m[lo_arm]
                b = np.array([
                    max(t[f"{model}_gain{order}_{a}"].to_numpy(float)[i].mean()
                        for a in ARMS_HERE)
                    / min(t[f"{model}_gain{order}_{a}"].to_numpy(float)[i].mean()
                          for a in ARMS_HERE) for i in draws])
                out.append({"check": f"{model} three-arm span, order-{order} baseline",
                            "value": float(span),
                            "ci_low": float(np.percentile(b, 2.5)),
                            "ci_high": float(np.percentile(b, 97.5)), "n": len(t),
                            "note": f"smallest arm {lo_arm} at {m[lo_arm]:+.4f}"})
                out.append({"check": f"{model} smallest-arm contribution, order-{order}",
                            "value": m[lo_arm], "n": len(t), "note": lo_arm})
                log(f"    order-{order} {model:11s} {span:6.2f}x  "
                    f"smallest arm {lo_arm} {m[lo_arm]:+.4f} "
                    f"({int((t[f'{model}_gain{order}_{lo_arm}'] > 0).sum())}/{len(t)} positive)")

    # CONCENTRATION. The Discussion says 51% of the surviving positive mass sits in 3 of 30
    # datasets. On the full panel that has to be recomputed, and it is worth having for every
    # model, because concentration is the difference between "a small effect everywhere" and
    # "nothing, plus a few outliers".
    log("\n  share of surviving positive mass in the top 3 datasets:")
    for arm in ARMS_HERE:
        for model in MODELS:
            v = t[f"{model}_gain3_{arm}"].to_numpy(float)
            pos = np.sort(v[v > 0])[::-1]
            share = float(pos[:3].sum() / pos.sum()) if pos.size else float("nan")
            out.append({"check": f"{model} top-3 share of order-3 positive mass, {arm} arm",
                        "value": share, "n": len(t),
                        "note": f"{pos.size} datasets positive"})
            # THE SHARE ALONE IS PANEL-SIZE DEPENDENT and would let a bigger panel look less
            # concentrated for free: 3 of 30 is a tenth of the panel, 3 of 94 a thirty-second.
            # Divide by the share three datasets would hold if the mass were spread evenly, so
            # the number means "times over-represented" and is comparable across panel sizes.
            if pos.size:
                out.append({"check": f"{model} top-3 over-representation, order-3, {arm} arm",
                            "value": float(share / (3 / pos.size)), "n": len(t),
                            "note": "top-3 share / uniform share, panel-size invariant"})
        log(f"    {arm:5s} " + "  ".join(
            f"{m} {100 * np.sort(t[f'{m}_gain3_{arm}'].to_numpy(float)[t[f'{m}_gain3_{arm}'].to_numpy(float) > 0])[::-1][:3].sum() / t[f'{m}_gain3_{arm}'].to_numpy(float)[t[f'{m}_gain3_{arm}'].to_numpy(float) > 0].sum():.0f}%"
            for m in MODELS))

    pd.DataFrame(out).to_csv(TABLES / "baseline_order_models.csv", index=False)
    log("\nwrote baseline_order_models.csv and baseline_order_models_per_dataset.csv")


if __name__ == "__main__":
    main()
