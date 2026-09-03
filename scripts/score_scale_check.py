"""R1t: the rungs of the R1g ladder were not measured on the same scale. Does it matter?

    python scripts/score_scale_check.py

THE DEFECT, found by a methods referee and confirmed. The nested contribution adds the model's
score as a covariate to a logistic regression. The 4-mer enters as a LOG-ODDS
(`baseline.oof_scores` returns `decision_function`); the CNN and SpliceBERT enter as
PROBABILITIES, because that is what their sweeps wrote to `data/evidence/scores*`. A logistic
regression is not invariant to a monotone NONLINEAR transform of a covariate, so the deep
models' gains were measured with a range-compressed version of their own score while the k-mer's
were not.

That is an asymmetry in R1g, whose entire point is that the three rungs are measured on the same
rows, the same folds and the same estimator. AUROC itself is invariant to monotone transforms,
so the standalone AUROCs are untouched; it is only the nested fit that can move.

WHAT THIS SCRIPT DOES. Recomputes every deep gain on BOTH scales -- as published (probability)
and on logit(p) -- and the k-mer on both, so the direction of the effect is measured rather
than argued. If the deltas are negligible the ladder stands and the paper says so with a
number. If they are not, the deep arms have to be reported on the logit scale.

WHY NOT JUST REFIT EVERYTHING ON THE LOGIT. Because the published numbers are gated in
golden.yaml and quoted throughout the manuscript, and silently changing the estimand would
substitute one quantity for another under the same name. Measure the difference first.
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

from baseline_order import composition  # noqa: E402
from deep_model_contrast import MIN_COVERAGE, arm_roots, oof  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import _oof_scores  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"
DEEP = ("cnn", "splicebert")
CLIP = 1e-6
REPRO_TOL = 5e-3


def log(m):
    print(m, flush=True)


def logit(p):
    p = np.clip(np.asarray(p, dtype=float), CLIP, 1 - CLIP)
    return np.log(p / (1 - p))


def gain(X, score, y, folds, s_comp):
    full = np.column_stack([X, standardise(score)])
    s_full = _oof_scores(full, y, folds)
    ok = np.isfinite(s_comp) & np.isfinite(s_full)
    return float(delong_test(s_full[ok], s_comp[ok], y[ok])["diff"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "score_scale_per_dataset.csv"
    pub = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv").set_index("dataset")

    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it with --store")
    else:
        roots = arm_roots(a.store)
        rows, worst = [], 0.0
        datasets = list(pub.index)[: a.n or None]
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
                for m in DEEP:
                    s = oof(scoreroot, cell, protein, m)
                    if s is None:
                        ok = False
                        break
                    got[m] = s
                    ids &= set(s.id)
                if not ok:
                    break
                if len(ids) / len(d) < MIN_COVERAGE:
                    ok = False
                    break
                dd = d[d.id.isin(ids)].reset_index(drop=True)
                y, fo = dd.label.values, dd.fold.values
                X = composition(dd.seq_rna.values, 2)
                s_comp = _oof_scores(X, y, fo)

                # the k-mer, both ways: published is the logit, so the counterfactual is p
                sc, _, _ = kmer_oof(dd.seq_rna.values, y, fo, k=4)
                rec[f"kmer_logit_{arm}"] = gain(X, sc, y, fo, s_comp)
                rec[f"kmer_prob_{arm}"] = gain(X, 1 / (1 + np.exp(-sc)), y, fo, s_comp)

                # the deep models, both ways: published is the PROBABILITY
                for m in DEEP:
                    mm = dd.merge(got[m], on="id", how="inner")
                    rec[f"{m}_prob_{arm}"] = gain(X, mm.score.values, y, fo, s_comp)
                    rec[f"{m}_logit_{arm}"] = gain(X, logit(mm.score.values), y, fo, s_comp)
                    ref = float(pub.loc[ds, f"{m}_gain_{arm}"])
                    rec[f"{m}_published_{arm}"] = ref
                    worst = max(worst, abs(rec[f"{m}_prob_{arm}"] - ref))
            if ok:
                rows.append(rec)
                log(f"[{i:3d}/{len(datasets)}] {ds:18s} " + "  ".join(
                    f"{m[:4]} {rec[f'{m}_prob_gc']:+.4f}->{rec[f'{m}_logit_gc']:+.4f}"
                    for m in DEEP))
        t = pd.DataFrame(rows)
        if t.empty:
            sys.exit("no dataset could be built")
        log(f"\nprobability-scale gain vs published, max |difference| = {worst:.2e}")
        if worst > REPRO_TOL:
            sys.exit(f"the probability-scale refit does not reproduce the published deep gain "
                     f"({worst:.2e}); the published numbers are not what this script models.")
        t.to_csv(per, index=False)

    out = []
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    rng = np.random.default_rng(0)
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(2000)]

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        lo, hi = np.percentile(b, [2.5, 97.5])
        out.append({"check": check, "value": float(v.mean()), "ci_low": float(lo),
                    "ci_high": float(hi), "n": len(t), "note": note})
        return float(v.mean())

    w = max(float((t[f"{m}_prob_{arm}"] - t[f"{m}_published_{arm}"]).abs().max())
            for m in DEEP for arm in ("gc", "dn"))
    out.append({"check": "max |probability-scale gain - published gain|", "value": w,
                "n": len(t)})
    log(f"\n=== R1t: the score scale, n = {len(t)}, {len(uniq)} proteins ===")
    log(f"  published deep gains reproduce on the probability scale to {w:.2e}\n")

    log(f"  {'':22s} {'as published':>13s} {'on logit(p)':>13s} {'difference':>26s}")
    for arm in ("gc", "dn"):
        for m in DEEP:
            pv = add(f"{m} gain, probability scale, {arm} arm", t[f"{m}_prob_{arm}"])
            lv = add(f"{m} gain, logit scale, {arm} arm", t[f"{m}_logit_{arm}"])
            dv = add(f"{m} logit minus probability scale, {arm} arm",
                     t[f"{m}_logit_{arm}"] - t[f"{m}_prob_{arm}"])
            r = out[-1]
            log(f"  {arm} {m:18s} {pv:+13.4f} {lv:+13.4f}   {dv:+.4f} "
                f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]")
        # and the k-mer the other way, so the direction is measured on both rungs
        kl = add(f"kmer gain, logit scale (as published), {arm} arm", t[f"kmer_logit_{arm}"])
        kp = add(f"kmer gain, probability scale, {arm} arm", t[f"kmer_prob_{arm}"])
        add(f"kmer probability minus logit scale, {arm} arm",
            t[f"kmer_prob_{arm}"] - t[f"kmer_logit_{arm}"])
        log(f"  {arm} {'kmer':18s} {kl:+13.4f} {kp:+13.4f}   (published is the LOGIT)")

    # THE LADDER, both ways. This is what the asymmetry could have distorted.
    for arm in ("gc", "dn"):
        for m in DEEP:
            c_pub = t[f"{m}_prob_dn"] - t[f"{m}_prob_gc"]
            c_log = t[f"{m}_logit_dn"] - t[f"{m}_logit_gc"]
            if arm == "gc":
                add(f"{m} R1 contrast, probability scale", c_pub)
                add(f"{m} R1 contrast, logit scale", c_log)
    kc_pub = t["kmer_logit_dn"] - t["kmer_logit_gc"]
    add("kmer R1 contrast, logit scale (as published)", kc_pub)
    log("\n  R1 contrast (dn-gc), as published vs all-logit:")
    for m in DEEP:
        log(f"    {m:12s} {(t[f'{m}_prob_dn'] - t[f'{m}_prob_gc']).mean():+.4f}  ->  "
            f"{(t[f'{m}_logit_dn'] - t[f'{m}_logit_gc']).mean():+.4f}")
    log(f"    {'kmer':12s} {kc_pub.mean():+.4f}  (unchanged, already the logit)")

    pd.DataFrame(out).to_csv(TABLES / "score_scale_check.csv", index=False)
    log("\nwrote score_scale_check.csv and score_scale_per_dataset.csv")


if __name__ == "__main__":
    main()
