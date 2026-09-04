"""F3: does removing the bias-aware arm's region confound change the answer for a NEURAL model?

    python scripts/region_matched_neural.py
    python scripts/region_matched_neural.py --from-cache

WHAT WAS MISSING. The bias-aware arm matches fold only, so region alone separates its classes at
a median AUROC of 0.7484, and the region-matched rebuild removes that: the arm's composition
baseline falls from 0.8248 to 0.8052 and the 4-mer's contribution from +0.0122 to +0.0092. Until
now that repair existed for the 4-mer ONLY, because nothing could dispatch a neural sweep for the
arm. So the paper's most-criticised result had a correction available for one model class out of
three, and the obvious objection -- that a convolutional network or a language model might depend
on region in a way a bag of 4-mers does not -- had no answer.

This closes it. Both neural models were trained on the region-matched arm across all 94 datasets
and five folds on the same hardware and code path as the published sweeps, and their per-window
out-of-fold scores are read here exactly as `deep_model_contrast.py` reads the published ones.
Nothing is refitted on the deep side.

THE COMPARISON THAT MATTERS is not the level but the ORDERING. The bias-aware arm's job in this
paper is to be the protocol with the highest composition baseline and the smallest contribution.
If region matching moves it out of that position for any model class, the three-arm span is
partly a region artefact. If it does not, the span survives a correction applied to the one arm
whose position the design does not imply.
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

from deep_model_contrast import MIN_COVERAGE, MODELS, oof  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402

TABLES = ROOT / "results" / "tables"
EVIDENCE_RM = ROOT / "data" / "evidence" / "scores_neg2_rm"


def log(m):
    print(m, flush=True)


def build(store, limit):
    store = Path(store)
    dataroot = store / "processed" / "neg2_rm"
    scoreroot = EVIDENCE_RM if EVIDENCE_RM.exists() else store / "runs" / "neg2_rm"
    if not scoreroot.exists():
        sys.exit(f"no per-window scores at {scoreroot}; run the neg2_rm sweep first")
    pub = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv").set_index("dataset")
    datasets = list(pub.index)[:limit or None]

    rows, incomplete = [], []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        f = dataroot / cell / protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        ids, got = set(d.id), {}
        ok = True
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
            incomplete.append(ds)
            continue
        if len(ids) / len(d) < MIN_COVERAGE:
            incomplete.append(ds)
            continue
        dd = d[d.id.isin(ids)].reset_index(drop=True)
        sc, _, _ = kmer_oof(dd.seq_rna.values, dd.label.values, dd.fold.values, k=4)
        got["kmer"] = pd.DataFrame({"id": dd.id.values, "score": sc})

        rec = {"dataset": ds, "protein": protein, "cell": cell, "n": len(dd),
               "coverage": len(ids) / len(d)}
        for model in MODELS:
            m = dd.merge(got[model], on="id", how="inner")
            # gain_over_composition IS the published estimator, imported rather than
            # reimplemented, so this column is comparable with deep_contrast_per_dataset.csv
            # cell for cell rather than merely similar to it.
            r = gain_over_composition(m.seq_rna.values, m.score.values,
                                      m.label.values, m.fold.values)
            rec[f"comp_{model}"] = float(r.auroc_composition)
            rec[f"full_{model}"] = float(r.auroc_with_score)
            rec[f"gain_{model}"] = float(r.delta)
            rec[f"published_neg2_{model}"] = float(pub.loc[ds, f"{model}_gain_neg2"])
            rec[f"published_gc_{model}"] = float(pub.loc[ds, f"{model}_gain_gc"])
            rec[f"published_dn_{model}"] = float(pub.loc[ds, f"{model}_gain_dn"])
        rows.append(rec)
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} " + "  ".join(
            f"{m[:4]} {rec[f'gain_{m}']:+.4f} (neg2 {rec[f'published_neg2_{m}']:+.4f})"
            for m in MODELS))
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no dataset could be built; refusing to overwrite the committed table")
    if incomplete:
        log(f"\n{len(incomplete)} datasets incomplete: {incomplete[:8]}"
            f"{' ...' if len(incomplete) > 8 else ''}")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "region_matched_neural_per_dataset.csv"
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

    log(f"\n=== F3: the region-matched bias-aware arm, three model classes, n = {len(t)}, "
        f"{len(uniq)} proteins ===\n")
    out.append({"check": "datasets with complete region-matched neural scores",
                "value": len(t), "n": len(t)})

    # THE 4-MER IS THE CONTROL, and it is a real one: its published region-matched contribution
    # was computed by region_asymmetry.py from a different code path on the same windows. If
    # this script's 4-mer column does not land on that number, the neural columns beside it are
    # measuring something else.
    log(f"  {'model':12s} {'baseline':>10s} {'contribution':>13s} "
        f"{'published neg2':>15s} {'change':>9s}")
    for model in MODELS:
        c = add(f"composition AUROC, region-matched arm, {model} rows", t[f"comp_{model}"])
        g = add(f"{model} contribution, region-matched arm", t[f"gain_{model}"])
        p_ = add(f"{model} contribution, published bias-aware arm",
                 t[f"published_neg2_{model}"])
        d = add(f"{model} change from matching region", t[f"gain_{model}"]
                - t[f"published_neg2_{model}"])
        log(f"  {model:12s} {c:10.4f} {g:+13.4f} {p_:+15.4f} {d:+9.4f}")

    # THE ORDERING IS THE CLAIM. The region-matched arm must still give the SMALLEST
    # contribution of the three for every model class, or the span is partly a region artefact.
    log("\n  is the region-matched arm still the smallest of the three, per model?")
    n_ok = 0
    for model in MODELS:
        rm = float(t[f"gain_{model}"].mean())
        gc = float(t[f"published_gc_{model}"].mean())
        dn = float(t[f"published_dn_{model}"].mean())
        ok = rm < gc and rm < dn
        n_ok += int(ok)
        out.append({"check": f"region-matched arm is smallest of three, {model}",
                    "value": int(ok), "n": len(t),
                    "note": f"rm {rm:+.4f} vs gc {gc:+.4f} vs dn {dn:+.4f}"})
        span = dn / rm if rm else float("nan")
        b = np.array([t[f"published_dn_{model}"].to_numpy(float)[i].mean()
                      / t[f"gain_{model}"].to_numpy(float)[i].mean() for i in draws])
        out.append({"check": f"three-arm span with the region-matched arm, {model}",
                    "value": float(span), "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t)})
        log(f"    {model:12s} {'YES' if ok else 'NO':4s}  rm {rm:+.4f}  gc {gc:+.4f}  "
            f"dn {dn:+.4f}   span {span:.2f}x")
    out.append({"check": "model classes for which the region-matched arm stays smallest",
                "value": n_ok, "n": len(MODELS)})

    # AND THE CONTRAST AGAINST THE DINUCLEOTIDE ARM, which is the paper's headline pair, with
    # the bias-aware arm replaced by its region-matched rebuild.
    for model in MODELS:
        add(f"dinucleotide minus region-matched contrast, {model}",
            t[f"published_dn_{model}"] - t[f"gain_{model}"])

    pd.DataFrame(out).to_csv(TABLES / "region_matched_neural.csv", index=False)
    log("\nwrote region_matched_neural.csv and region_matched_neural_per_dataset.csv")


if __name__ == "__main__":
    main()
