"""B3: the contribution as a function of where the composition baseline stops, orders 1 to 4.

    python scripts/order_profile.py                 # recompute, ~40 min
    python scripts/order_profile.py --from-cache    # summary only

THE OBJECTION THIS RETIRES. Every increment in this paper is measured over an order-two
composition baseline, and Section on order three shows the magnitudes are a strong function of
that stopping point. Reporting two orders invites "why not another one?", and answering it one
order at a time invites it again. A profile over orders 1 to 4 is immune to the question: it
reports the whole function instead of two of its values, and a reader who prefers a different
baseline can read their own number off it.

WHY IT STOPS AT FOUR, and this turned out to be the most informative part. The 4-mer model's
features ARE order-four counts, so at an order-four baseline the composition block spans the
model's entire feature space and the true contribution is exactly zero by construction. That is
a prediction with no free parameters, and it FAILS: the estimator reports a large positive
contribution there. The reason is that 337 standardised columns overfit at these sample sizes,
so the order-four baseline's own out-of-fold AUROC FALLS below the order-three baseline's on
most datasets, and the pre-fit model score then recovers what the baseline lost. The measured
"contribution" at order four is the gap between two estimators of the same information, not
information.

So order four is not reported as a baseline. It is reported as a calibration of the instrument:
it says how much apparent contribution this estimator manufactures when the truth is zero, and
that quantity is larger than most of the contributions this literature reports. Orders one to
three are the profile; order four is the noise floor, and the diagnostic that separates them is
whether the baseline's own AUROC still rises.

THE DOUBLE ANCHOR. Order two must reproduce deep_contrast_per_dataset.csv and order three must
reproduce baseline_order_models_per_dataset.csv, both per cell. Two independently produced tables
pin two of the four points, so the profile cannot be a self-consistent re-implementation of
something else. Same rows, same folds, same committed per-window model scores, nothing retrained.
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
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ORDERS = (1, 2, 3, 4)
ARMS = ("gc", "dn", "neg2")



def build(store, limit):
    roots = arm_roots(store)
    pub = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv").set_index("dataset")
    o3 = pd.read_csv(TABLES / "baseline_order_models_per_dataset.csv").set_index("dataset")
    datasets = list(pub.index)[:limit or None]
    rows = []
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
                rec[f"ncol{order}"] = X.shape[1]
                for model in MODELS:
                    m = dd.merge(got[model], on="id", how="inner")
                    s_full = _oof_scores(
                        np.column_stack([X, standardise(m.score.values)]), y, fo)
                    good = np.isfinite(s_comp) & np.isfinite(s_full)
                    r = delong_test(s_full[good], s_comp[good], y[good])
                    rec[f"{model}_gain{order}_{arm}"] = float(r["diff"])
                    rec[f"comp{order}_{arm}"] = float(r["auc_b"])
            rec[f"nrows_{arm}"] = int(len(dd))
        if not ok:
            continue
        # THE DOUBLE ANCHOR, per dataset, before the row is kept. Two tables produced by two
        # other scripts pin orders 2 and 3; a profile that agrees with neither would still look
        # perfectly smooth.
        for arm in ARMS:
            for model in MODELS:
                rec[f"{model}_anchor2_{arm}"] = float(pub.loc[ds, f"{model}_gain_{arm}"])
                rec[f"{model}_anchor3_{arm}"] = float(o3.loc[ds, f"{model}_gain3_{arm}"])
        rows.append(rec)
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} kmer gc " + " ".join(
            f"{rec[f'kmer_gain{o}_gc']:+.4f}" for o in ORDERS))
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no dataset could be built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0, help="limit datasets, for smoke tests")
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "order_profile_per_dataset.csv"
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

    # THE ANCHORS FIRST, and they are hard stops rather than reported diagnostics: every number
    # below is meaningless if the profile does not pass through the two published points.
    for order, tag in ((2, "anchor2"), (3, "anchor3")):
        w = max(float((t[f"{m}_gain{order}_{arm}"] - t[f"{m}_{tag}_{arm}"]).abs().max())
                for m in MODELS for arm in ARMS)
        src = ("deep_contrast_per_dataset.csv" if order == 2
               else "baseline_order_models_per_dataset.csv")
        out.append({"check": f"max |order-{order} gain - {src}|", "value": w, "n": len(t)})
        if w > REPRO_TOL:
            sys.exit(f"order {order} does not reproduce {src} ({w:.2e} > {REPRO_TOL})")
        log(f"  order {order} reproduces {src} to {w:.2e}")

    log(f"\n=== B3: composition-order profile, orders {ORDERS}, n = {len(t)}, "
        f"{len(uniq)} proteins ===\n")
    for order in ORDERS:
        out.append({"check": f"composition columns at order {order}",
                    "value": int(t[f"ncol{order}"].iloc[0]), "n": len(t)})
    log("  baseline width: " + "  ".join(
        f"order {o} {int(t[f'ncol{o}'].iloc[0])} cols" for o in ORDERS))

    for arm in ARMS:
        log(f"\n  {arm} arm")
        log(f"    {'order':>5s} {'composition':>12s} " + " ".join(f"{m:>12s}" for m in MODELS))
        for order in ORDERS:
            c = add(f"composition AUROC at order {order}, {arm} arm", t[f"comp{order}_{arm}"])
            vals = [add(f"{m} gain at order {order}, {arm} arm", t[f"{m}_gain{order}_{arm}"])
                    for m in MODELS]
            for m in MODELS:
                n_pos = int((t[f"{m}_gain{order}_{arm}"] > 0).sum())
                out.append({"check": f"{m} gain positive at order {order}, {arm} arm",
                            "value": n_pos, "n": len(t)})
            log(f"    {order:5d} {c:12.4f} " + " ".join(f"{v:+12.4f}" for v in vals))

    # DOES THE BASELINE STILL IMPROVE? This is the diagnostic that decides whether an order is
    # a baseline or a noise floor. A composition block whose own out-of-fold AUROC FALLS when
    # given more columns is overfitting, and a nested increment over an overfitted baseline
    # measures the gap between two estimators rather than a gap in information.
    log("\n  does the baseline's own out-of-fold AUROC still rise with order?")
    for arm in ARMS:
        for order in ORDERS[1:]:
            fell = int((t[f"comp{order}_{arm}"] < t[f"comp{order - 1}_{arm}"]).sum())
            out.append({"check": f"baseline AUROC fell from order {order - 1} to {order}, "
                                 f"{arm} arm", "value": fell, "n": len(t)})
        log(f"    {arm:5s} fell in: " + "  ".join(
            f"order {o - 1}->{o} {int((t[f'comp{o}_{arm}'] < t[f'comp{o - 1}_{arm}']).sum())}"
            f"/{len(t)}" for o in ORDERS[1:]))

    # THE ZERO OF THE SCALE, WHICH IS NOT AT ZERO. A bag of 4-mer counts has no information
    # beyond an order-four composition block, so a correct estimator must report nothing there.
    # This one reports a large positive number, and that number is the instrument's error when
    # the truth is known. It is the most useful thing in this section: it is the scale on which
    # every contribution in this literature should be read.
    log("\n  the 4-mer over an order-four baseline, which SPANS its own feature space, so the"
        "\n  true contribution is zero by construction and anything measured is estimator error:")
    for arm in ARMS:
        v = t[f"kmer_gain4_{arm}"]
        n_pos = int((v > 0).sum())
        out.append({"check": f"kmer gain at order 4 positive in, {arm} arm", "value": n_pos,
                    "n": len(t), "note": "the estimator's noise floor, not a contribution"})
        floor = add(f"kmer noise-floor gain at order 4, {arm} arm", v,
                    "true value is zero by construction")
        # AND HOW IT COMPARES WITH WHAT THE PAPER REPORTS. A noise floor is only alarming
        # relative to the signal, so state the ratio rather than leaving the reader to divide.
        ratio = floor / float(t[f"kmer_gain2_{arm}"].mean())
        out.append({"check": f"noise floor as a fraction of the order-2 gain, {arm} arm",
                    "value": float(ratio), "n": len(t)})
        log(f"    {arm:5s} {floor:+.4f}, positive in {n_pos}/{len(t)}, which is "
            f"{ratio:.2f}x the order-2 contribution")

    # AND THE MECHANISM, TESTED RATHER THAN ASSERTED. If the floor is overfitting of a
    # 337-column block, it must shrink as the sample grows. Spearman, because the relation
    # need not be linear and one dataset is 20x another.
    from scipy.stats import spearmanr
    log("\n  the floor against sample size (overfitting predicts a negative correlation):")
    for arm in ARMS:
        r, pv = spearmanr(t[f"nrows_{arm}"], t[f"kmer_gain4_{arm}"])
        out.append({"check": f"spearman(rows, order-4 noise floor), {arm} arm",
                    "value": float(r), "n": len(t), "note": f"p = {pv:.2e}"})
        log(f"    {arm:5s} rho {r:+.3f}  p = {pv:.2e}")

    # THE PROFILE'S SHAPE, as one number per model per arm: what fraction of the order-one
    # contribution survives each further order. This is the reporting object the section is
    # for, and it is a ratio of means for the reason spelled out in baseline_order_models.py.
    log("\n  fraction of the order-1 contribution surviving each order (ratio of means):")
    for arm in ARMS:
        for m in MODELS:
            base = t[f"{m}_gain1_{arm}"].to_numpy(float)
            for order in ORDERS[1:]:
                v = t[f"{m}_gain{order}_{arm}"].to_numpy(float)
                b = np.array([v[i].mean() / base[i].mean() for i in draws])
                out.append({"check": f"{m} share of order-1 gain surviving order {order}, "
                                     f"{arm} arm", "value": float(v.mean() / base.mean()),
                            "ci_low": float(np.percentile(b, 2.5)),
                            "ci_high": float(np.percentile(b, 97.5)), "n": len(t),
                            "note": "ratio of means, protein-clustered"})
        log(f"    {arm:5s} " + "  ".join(
            f"{m} " + "/".join(f"{t[f'{m}_gain{o}_{arm}'].mean() / t[f'{m}_gain1_{arm}'].mean():.2f}"
                               for o in ORDERS[1:]) for m in MODELS))

    # AND THE PAPER'S OWN TWO QUANTITIES AT EVERY ORDER: the two-arm contrast and the
    # three-arm span. This is the claim the profile exists to test -- that protocol dependence
    # is not an artefact of where the baseline stops -- stated over the whole function rather
    # than at two of its points.
    log("\n  the paper's quantities at every order:")
    for order in ORDERS:
        for m in MODELS:
            c = add(f"{m} two-arm contrast (dn-gc) at order {order}",
                    t[f"{m}_gain{order}_dn"] - t[f"{m}_gain{order}_gc"])
            means = {arm: float(t[f"{m}_gain{order}_{arm}"].mean()) for arm in ARMS}
            lo_arm = min(means, key=means.get)
            b = np.array([max(t[f"{m}_gain{order}_{x}"].to_numpy(float)[i].mean()
                              for x in ARMS)
                          / min(t[f"{m}_gain{order}_{x}"].to_numpy(float)[i].mean()
                                for x in ARMS) for i in draws])
            span = means[max(means, key=means.get)] / means[lo_arm]
            out.append({"check": f"{m} three-arm span at order {order}", "value": float(span),
                        "ci_low": float(np.percentile(b, 2.5)),
                        "ci_high": float(np.percentile(b, 97.5)), "n": len(t),
                        "note": f"smallest arm {lo_arm} at {means[lo_arm]:+.4f}"})
            log(f"    order {order} {m:11s} contrast {c:+.4f}   span {span:6.2f}x   "
                f"smallest arm {lo_arm} {means[lo_arm]:+.4f}")

    pd.DataFrame(out).to_csv(TABLES / "order_profile.csv", index=False)
    log("\nwrote order_profile.csv and order_profile_per_dataset.csv")


if __name__ == "__main__":
    main()
