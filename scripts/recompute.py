"""Recompute published AUROCs from committed per-example scores. No cloud, no credentials.

WHY THIS EXISTS, and it is the sharpest lesson in the project.

`scripts/verify.py` asserts 86 named numbers and passes. Then someone tested it properly:
subtract 0.05 from `rehearsal_binding_gc.composition_auroc`, or zero every variant delta and
every conservation value, and it still reports **all checks passed**. It reads summary tables
written by the same analysis pass, so it detects whether the OUTPUT changed, not whether the
output was correctly DERIVED. It is a regression detector described as a correctness proof.

This file is the correctness half, and it is deliberately small enough to audit in one sitting.
It reads the per-example model scores committed under `data/evidence/`, pools the five
out-of-fold test sets, computes AUROC with sklearn, and compares against the published table.
Nothing in between: no project code, no config, no bucket.

    python scripts/recompute.py

Two model classes, not one. The published claim is a LADDER across architectures
(composition < k-mer < CNN < SpliceBERT), so recomputing only the last rung would prove the
arm every objection targets and leave the comparison itself unevidenced. Both deep arms have
committed scores, so both are checked; the two logistic arms are cheap to refit but their
per-example scores are the rehearsal files, checked separately.

WHAT A PASS MEANS. The number in the paper is the number the model produced on held-out data.
It does NOT mean the split was correct, the negatives were sensible, or the claim is
interesting -- those are arguments elsewhere. It means the arithmetic is not invented, which
`verify.py` alone cannot establish.

WHAT A FAIL MEANS. Either the tables drifted from the evidence, or the evidence was corrupted.
Both are worth stopping for. Zeroing the score column produces max|diff| ~0.45, so the check
has real power rather than passing on anything.
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
EVIDENCE = ROOT / "data" / "evidence" / "scores"

# The deep arms, and the column each one is published under.
ARMS = (("splicebert", "splicebert"), ("cnn", "cnn"))
N_FOLDS = 5

# THE REHEARSAL ARM, WHICH IS THE ONE THE PRIMARY RESULT RESTS ON.
#
# This was missing until a referee pointed a gun at it. `recompute.py` covered the two deep
# arms only, so R1 -- the paper's headline, built entirely on the composition and k-mer
# rehearsal -- had no evidence-level check at all. The referee zeroed the `score` column in all
# 74 rehearsal files and this script reported "Every recomputed AUROC matches", while
# verify.py returned 105/105. The most-checked project I have worked on could not detect its
# own primary result being deleted.
#
# The lesson generalises past this repo: coverage follows attention, and attention had gone to
# the interesting models rather than the load-bearing ones.
REHEARSAL = ROOT / "data" / "evidence" / "rehearsal"
REHEARSAL_TABLE = "rehearsal_binding_dinuc.csv"
TOL = 1.0e-9          # sklearn on the same floats; anything above this is drift, not noise
MIN_DATASETS = 90     # of 95; a few missing folds is a broken mirror, not a broken claim


def log(m):
    print(m, flush=True)


def fold_files(cell, protein, model):
    return sorted(glob.glob(str(EVIDENCE / cell / protein / model / "fold*" / "scores.tsv.gz")))


def recompute_arm(published, model, col, zero=False):
    """One AUROC per dataset, pooled over the five out-of-fold test sets."""
    rows = []
    for ds, pub in published[col].items():
        protein, cell = ds.split(":")
        fs = fold_files(cell, protein, model)
        if len(fs) != N_FOLDS:
            rows.append({"dataset": ds, "model": model, "published": pub,
                         "recomputed": np.nan, "abs_diff": np.nan,
                         "note": f"{len(fs)}/{N_FOLDS} folds present"})
            continue
        d = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
        score = np.zeros(len(d)) if zero else d.score.to_numpy()
        if d.label.nunique() < 2:
            rows.append({"dataset": ds, "model": model, "published": pub,
                         "recomputed": np.nan, "abs_diff": np.nan, "note": "one class"})
            continue
        got = roc_auc_score(d.label, score)
        rows.append({"dataset": ds, "model": model, "published": pub, "recomputed": got,
                     "abs_diff": abs(got - pub), "note": ""})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--zero-scores", action="store_true",
                   help="corrupt the evidence on purpose; the check MUST fail")
    p.add_argument("--tol", type=float, default=TOL)
    a = p.parse_args()

    four = TABLES / "matched_four_models.csv"
    if not four.exists():
        raise SystemExit(f"{four} is absent; nothing to check against")
    published = pd.read_csv(four).set_index("dataset")

    if not EVIDENCE.exists():
        raise SystemExit(
            f"{EVIDENCE} is absent. The per-example scores are what make this check possible;\n"
            f"without them every published AUROC rests on a table that the analysis wrote.")

    parts, failed = [], False
    for model, col in ARMS:
        if col not in published.columns:
            log(f"  {model:11} SKIP: no `{col}` column in matched_four_models.csv")
            continue
        r = recompute_arm(published, model, col, zero=a.zero_scores)
        parts.append(r)
        ok = r.dropna(subset=["abs_diff"])
        worst = ok.abs_diff.max() if len(ok) else np.nan
        good = len(ok) >= MIN_DATASETS and worst <= a.tol
        failed |= not good
        log(f"  [{'PASS' if good else 'FAIL'}] {model:11} {len(ok):3d} datasets recomputed, "
            f"max|diff| = {worst:.2e}  (want <= {a.tol:.0e} on >= {MIN_DATASETS})")
        for _, x in r[r.note != ""].iterrows():
            log(f"           {x.dataset}: {x.note}")

    # --- the rehearsal (k-mer) arm: the evidence behind the PRIMARY result ----------------
    reh = TABLES / REHEARSAL_TABLE
    if reh.exists() and REHEARSAL.exists():
        pub = pd.read_csv(reh)
        rows = []
        for _, r in pub.iterrows():
            f = REHEARSAL / r.cell / f"{r.protein}.scores.tsv.gz"
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="\t")
            score = np.zeros(len(d)) if a.zero_scores else d.score.to_numpy()
            if d.label.nunique() < 2:
                continue
            got = roc_auc_score(d.label, score)
            rows.append({"dataset": f"{r.protein}:{r.cell}", "model": "kmer (rehearsal)",
                         "published": float(r.auroc), "recomputed": got,
                         "abs_diff": abs(got - float(r.auroc)), "note": ""})
        if rows:
            rr = pd.DataFrame(rows)
            parts.append(rr)
            worst = rr.abs_diff.max()
            good = len(rr) >= MIN_DATASETS and worst <= a.tol
            failed |= not good
            log(f"  [{'PASS' if good else 'FAIL'}] {'kmer':11} {len(rr):3d} datasets recomputed, "
                f"max|diff| = {worst:.2e}  (want <= {a.tol:.0e} on >= {MIN_DATASETS})")
        else:
            log("  [FAIL] kmer        no rehearsal evidence found; the PRIMARY result is "
                "unchecked")
            failed = True

    if parts:
        out = pd.concat(parts, ignore_index=True)
        out.to_csv(TABLES / "recompute.csv", index=False)
        n = int(out.abs_diff.notna().sum())
        log(f"\n  {n} published AUROCs across 3 arms, "
            f"recomputed from {len(glob.glob(str(EVIDENCE / '*/*/*/fold*/scores.tsv.gz')))} "
            f"committed score files")
        log(f"  wrote {TABLES / 'recompute.csv'}")

    if a.zero_scores:
        log("\n  --zero-scores was set, so FAIL above is the correct outcome. "
            "A check that cannot fail is not a check.")
        return 0 if failed else 1        # inverted: passing on zeros is the bug
    if failed:
        log("\n  THE PUBLISHED NUMBERS DO NOT MATCH THE EVIDENCE.")
        return 1
    log("\n  Every recomputed AUROC matches the published table.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
