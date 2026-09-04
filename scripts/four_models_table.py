"""Rebuild matched_four_models.csv from committed tables, with no bucket.

WHY THIS EXISTS. `recompute.py` -- the end-to-end check that rebuilds published AUROCs from
per-example scores -- compares against `matched_four_models.csv`, and the only thing that wrote
that table was `cloud_analysis.py::four_models`, which reads two objects out of GCS. That
project's billing account is closed, so the one table the strongest reproducibility check in the
repo validates against could not be regenerated at all. After the dinucleotide retrain it was
stale, `recompute.py` failed by 0.025 and 0.021, and there was no offline way to fix it.

Both inputs are committed: `rehearsal_binding_dinuc.csv` for the composition and k-mer arms, and
`sweep_dinuc.csv` for the neural pooled AUROCs. The arithmetic is the same as `four_models`, and
it is reproduced here rather than refactored so the cloud path stays byte-for-byte what it was.

NOT CIRCULAR. This writes the pooled AUROC that `cloud_train.py aggregate` computed with
DeLong's estimator; `recompute.py` then rebuilds the same quantity from the per-window scores
with scikit-learn. Two implementations, two code paths, one number.

A ROW-SET DRIFT THE RETRAIN EXPOSED, recorded here because this is where it is visible. The
committed k-mer rehearsal scores cover marginally FEWER rows than the store's current window
tables: 46,380 against 46,384 for KHSRP:K562, 22,202 against 22,216 for AQR:HepG2. They were
computed before `dataset.tsv` was last regenerated, and the pre-retrain neural scores matched
that older set, so `aggregate`'s exact-join guard was satisfied and nothing showed. The
retrained scores cover the current tables in full, so for 22 of 190 rows that guard now
declines to compute `delta_vs_kmer` rather than compare two AUROCs over different rows. That is
the guard behaving correctly. The affected column is read by nothing: `verify.py` never loads
`sweep_dinuc.csv`, and the only consumer is the `pooled_auroc` pivot below. The model-class
comparison in the paper is unaffected, because `deep_model_contrast.py` refits the k-mer itself
on whatever rows it scores, so all three models there are on one row set by construction. The
drift is 0.009% of rows at worst and moves no published number; it is not repaired only because
regenerating the rehearsal would re-derive the paper's headline k-mer figures to fix a
four-row discrepancy.
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def main():
    reh = pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv")
    sweep = pd.read_csv(TABLES / "sweep_dinuc.csv")
    base = reh.set_index("dataset")
    wide = sweep.pivot_table(index="dataset", columns="model", values="pooled_auroc")
    d = pd.DataFrame({
        "composition_auroc": base.composition_auroc,
        "kmer_auroc": base.auroc,
        "pairs": base.pairs,
        "protein": base.protein,
        "cell": base.cell,
    })
    for m in ("cnn", "splicebert"):
        if m in wide:
            d[m] = wide[m]
    need = [c for c in ("composition_auroc", "kmer_auroc", "cnn", "splicebert") if c in d]
    d = d.dropna(subset=need).reset_index()
    d["delta_auroc"] = d.kmer_auroc - d.composition_auroc

    out = TABLES / "matched_four_models.csv"
    if out.exists():
        prev = pd.read_csv(out)
        if len(prev) != len(d):
            print(f"  NOTE: row count changes {len(prev)} -> {len(d)}")
        j = prev.merge(d, on="dataset", suffixes=("_prev", ""))
        for c in ("cnn", "splicebert", "kmer_auroc", "composition_auroc"):
            if f"{c}_prev" in j:
                delta = (j[c] - j[f"{c}_prev"]).abs()
                moved = int((delta > 1e-9).sum())
                print(f"  {c:18} changed in {moved:3d} of {len(j)} datasets, "
                      f"max |delta| {delta.max():.6f}")
    d.to_csv(out, index=False)
    print(f"\n{len(d)} datasets, {len(need)} models, {d.protein.nunique()} proteins; means "
          + ", ".join(f"{c}={d[c].mean():.4f}" for c in need))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
