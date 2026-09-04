"""B1: the nested contribution on four estimands that are not AUROC.

    python scripts/estimands.py --store ../rbp-store
    python scripts/estimands.py --store ../rbp-store --only KHSRP:K562

WHY THIS EXISTS. Two long reviews made the same point independently: the whole sweep never
leaves the ROC, so "no protocol-free measure of contribution exists" is a claim about AUROC and
not about measurement in general. The obvious reply is that an AUROC increment is bounded above
by 1 - baseline, so part of any comparison between arms with unequal baselines is arithmetic --
which is exactly the objection this paper spends a subsection answering for AUROC. If the
protocol ordering survives on an UNBOUNDED scale, that answer gets much stronger; if it does
not, the paper's central claim is a property of one summary statistic and has to say so.

See rbp.eval.estimands for what each one is. Row selection and the k-mer refit are imported
from deep_model_contrast.py, so the estimands are computed on exactly the rows the published
AUROC increments use.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from deep_model_contrast import MIN_COVERAGE, MODELS, arm_roots, oof  # noqa: E402
from nested_scale import panel  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.estimands import alternatives  # noqa: E402
from rbp.eval.nested import _oof_scores, composition_features  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARM_ORDER = ("dn", "gc", "neg2")


def log(m):
    print(m, flush=True)


def run(store, datasets):
    roots = arm_roots(store)
    rows = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        for arm, (dataroot, scoreroot) in roots.items():
            f = dataroot / cell / protein / "dataset.tsv"
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="\t")
            ids, got = set(d.id), {}
            bad = False
            for model in MODELS:
                if model == "kmer":
                    continue
                s = oof(scoreroot, cell, protein, model)
                if s is None:
                    bad = True
                    break
                got[model] = s
                ids &= set(s.id)
            if bad or len(ids) / len(d) < MIN_COVERAGE:
                continue
            dd = d[d.id.isin(ids)].reset_index(drop=True)
            y, folds = dd.label.values, dd.fold.values
            sc, _, _ = kmer_oof(dd.seq_rna.values, y, folds, k=4)
            scores = {"kmer": sc}
            for model, s in got.items():
                scores[model] = dd[["id"]].merge(s, on="id", how="left").score.to_numpy()

            comp, _ = composition_features(dd.seq_rna.values, True)
            s_comp = _oof_scores(comp, y, folds, "l2")
            okc = np.isfinite(s_comp)
            for model, raw in scores.items():
                col = standardise(raw)
                s_full = _oof_scores(np.column_stack([comp, col]), y, folds, "l2")
                ok = okc & np.isfinite(s_full) & np.isfinite(raw)
                r = alternatives(y[ok], s_comp[ok], s_full[ok], comp[ok], col[ok])
                rows.append({"dataset": ds, "arm": arm, "model": model,
                             "n": int(ok.sum()), **r})
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} {len([r for r in rows if r['dataset']==ds])}")
    return pd.DataFrame(rows)


KEYS = ("auroc_gain", "delta_deviance", "mcfadden_gain", "ap_gain", "idi", "residual_auroc")


def summarise(t):
    out = []
    for key in KEYS:
        for model in MODELS:
            v = {}
            for arm in ARM_ORDER:
                sub = t[(t.arm == arm) & (t.model == model)][key].dropna()
                if not len(sub):
                    continue
                v[arm] = float(sub.mean())
                out.append({"check": f"{key}, {arm} arm, {model}", "value": v[arm],
                            "n": len(sub)})
            # THE DEVIANCE IS QUOTED AS AN INTEGER, so the integer is emitted. The float
            # haystack in audit_manuscript.py indexes values at 3 to 6 decimals and the
            # integer haystack accepts only exact integers, so a table holding 966.7 sources
            # neither "967" nor "966.7". A paper that rounds has to say what it rounded to.
            if key == "delta_deviance":
                for arm, val in v.items():
                    out.append({"check": f"{key} rounded, {arm} arm, {model}",
                                "value": float(round(val)), "n": len(t)})
            if len(v) == 3 and key != "residual_auroc":
                # THE ORDERING IS THE CLAIM. dn largest, neg2 smallest, on every estimand.
                ordered = v["dn"] > v["gc"] > v["neg2"]
                out.append({"check": f"protocol ordering holds, {key}, {model}",
                            "value": int(ordered), "n": 3,
                            "note": f"dn {v['dn']:.4g} gc {v['gc']:.4g} neg2 {v['neg2']:.4g}"})
    # And the headline: on how many of the five increment estimands does the ordering hold for
    # every model class?
    inc = [k for k in KEYS if k != "residual_auroc"]
    held = 0
    for key in inc:
        rows = [r for r in out if r["check"].startswith(f"protocol ordering holds, {key},")]
        if len(rows) == len(MODELS) and all(r["value"] == 1 for r in rows):
            held += 1
    out.append({"check": "estimands on which the protocol ordering holds for every model",
                "value": held, "n": len(inc)})
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--only", default="")
    a = p.parse_args()
    datasets = [x.strip() for x in a.only.split(",") if x.strip()] or panel(Path(a.store))
    log(f"=== B1: five estimands, {len(datasets)} datasets ===\n")
    t = run(Path(a.store), datasets)
    if t.empty:
        sys.exit("no rows produced")
    t.to_csv(TABLES / "estimands_per_dataset.csv", index=False)
    s = summarise(t)
    s.to_csv(TABLES / "estimands.csv", index=False)
    log("\n" + s[s.check.str.startswith(("protocol ordering", "estimands on"))]
        .to_string(index=False))
    log("\nwrote estimands_per_dataset.csv and estimands.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
