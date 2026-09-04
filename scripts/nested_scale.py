"""1b and 1c: the covariate scale and the standardisation window, measured as a 2x2.

    python scripts/nested_scale.py --store ../rbp-store               # the whole panel
    python scripts/nested_scale.py --store ../rbp-store --only KHSRP:K562,AQR:HepG2

WHY ONE SCRIPT FOR BOTH. They are two coordinates of the same design matrix and share every
expensive step, so measuring them separately would fit the composition baseline four times per
(dataset, arm) instead of twice. See rbp.eval.scale_sensitivity for what each one is.

WHAT IS DONE WITH THE ANSWERS IS NOT SYMMETRIC. The covariate scale is an INCONSISTENCY
between arms of the model-class comparison -- a log odds in one, a probability in the other two
-- so it is corrected in the primary numbers. The standardisation window is IMPROPER but
symmetric and label-free, so it is reported as a bound rather than adopted; re-deriving every
published composition AUROC to remove an effect this size is not a trade worth making, and the
size is exactly what this measures.

Row selection is `deep_model_contrast.py`'s, imported rather than reimplemented: the same
intersection of the three models' covered rows, the same k-mer refit on that intersection, the
same minimum-coverage floor. If those two scripts ever disagree about which rows a dataset has,
this comparison is meaningless.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from deep_model_contrast import MIN_COVERAGE, MODELS, arm_roots, oof  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.scale_sensitivity import dataset_rows  # noqa: E402

TABLES = ROOT / "results" / "tables"


def log(m):
    print(m, flush=True)


def panel(store):
    """Every dataset with a window table in all three arms."""
    roots = arm_roots(store)
    sets = []
    for _arm, (dataroot, _s) in roots.items():
        found = set()
        for cell in ("K562", "HepG2"):
            d = dataroot / cell
            if d.is_dir():
                found |= {f"{p.name}:{cell}" for p in d.iterdir()
                          if (p / "dataset.tsv").exists()}
        sets.append(found)
    return sorted(set.intersection(*sets)) if sets else []


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
                log(f"  SKIP {ds} {arm}")
                continue
            dd = d[d.id.isin(ids)].reset_index(drop=True)
            sc, _, _ = kmer_oof(dd.seq_rna.values, dd.label.values, dd.fold.values, k=4)
            scores = {"kmer": sc}
            for model, s in got.items():
                scores[model] = dd[["id"]].merge(s, on="id", how="left").score.to_numpy()
            rows += dataset_rows(ds, arm, dd.seq_rna.values, dd.label.values,
                                 dd.fold.values, scores)
        done = [r for r in rows if r["dataset"] == ds]
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} {len(done)} cells")
    return pd.DataFrame(rows)


def summarise(t):
    """Panel means per (arm, model, scale, standardisation), and the two deltas."""
    g = t.groupby(["arm", "model", "scale", "standardisation"]).gain.mean().reset_index()
    out = []
    for _, r in g.iterrows():
        out.append({"check": f"gain, {r.arm} arm, {r.model}, {r.scale}, {r.standardisation}",
                    "value": float(r.gain), "n": int((t.arm == r.arm).sum())})

    # 1b: probability -> logit, at the published standardisation
    w = t[t.standardisation == "whole_dataset"]
    for arm in sorted(w.arm.unique()):
        for model in ("cnn", "splicebert"):
            a = w[(w.arm == arm) & (w.model == model) & (w.scale == "probability")]
            b = w[(w.arm == arm) & (w.model == model) & (w.scale == "logit")]
            if len(a) and len(b):
                j = a.merge(b, on="dataset", suffixes=("_p", "_l"))
                out.append({"check": f"logit minus probability, {arm} arm, {model}",
                            "value": float((j.gain_l - j.gain_p).mean()), "n": len(j),
                            "note": f"max |delta| {(j.gain_l - j.gain_p).abs().max():.4f}"})

    # 1c: whole-dataset -> within-fold, at each model's own primary scale
    prim = t[((t.model == "kmer") & (t.scale == "native"))
             | ((t.model != "kmer") & (t.scale == "logit"))]
    for arm in sorted(prim.arm.unique()):
        for model in ("kmer", "cnn", "splicebert"):
            a = prim[(prim.arm == arm) & (prim.model == model)
                     & (prim.standardisation == "whole_dataset")]
            b = prim[(prim.arm == arm) & (prim.model == model)
                     & (prim.standardisation == "within_fold")]
            if len(a) and len(b):
                j = a.merge(b, on="dataset", suffixes=("_w", "_f"))
                out.append({"check": f"within-fold minus whole-dataset, {arm} arm, {model}",
                            "value": float((j.gain_f - j.gain_w).mean()), "n": len(j),
                            "note": f"max |delta| {(j.gain_f - j.gain_w).abs().max():.4f}"})
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--only", default="", help="comma-separated PROTEIN:CELL subset")
    p.add_argument("--out", default="nested_scale")
    a = p.parse_args()

    datasets = ([x.strip() for x in a.only.split(",") if x.strip()]
                or panel(Path(a.store)))
    log(f"=== nested scale 2x2: {len(datasets)} datasets ===\n")
    t = run(Path(a.store), datasets)
    if t.empty:
        sys.exit("no rows produced")
    TABLES.mkdir(parents=True, exist_ok=True)
    t.to_csv(TABLES / f"{a.out}_per_dataset.csv", index=False)
    s = summarise(t)
    s.to_csv(TABLES / f"{a.out}.csv", index=False)
    log("\n" + s.to_string(index=False))
    log(f"\nwrote {a.out}_per_dataset.csv and {a.out}.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
