"""Supplementary Table S1: every dataset in the study panel, joined to its ENCODE accessions.

    python scripts/table_s1.py

Without this the paper is not reproducible from the manuscript alone: the results name proteins
and cell lines, and ENCODE is indexed by accession. A reader who wants to re-derive one row has
to guess which of several experiments per protein was used.

One row per (protein, cell line) in the panel, carrying the ENCFF file and ENCSR experiment
accessions, the replicate count, the window counts each protocol produced and the composition
baseline each one left. The last two columns are the point: they let a reader see, per dataset,
how much room the protocol left before any model is involved.

FAILS LOUDLY ON AN INCOMPLETE JOIN. A supplementary table with silent gaps is worse than none,
because it looks complete.
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
CONFIG = ROOT / "config"
OUT = TABLES / "supplementary_table_s1.csv"


def main():
    acc = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(CONFIG.glob("panel_full*.tsv"))],
                    ignore_index=True).drop_duplicates(["protein", "cell_line"])
    panel = pd.read_csv(TABLES / "panel_summary.csv")
    arms = pd.read_csv(TABLES / "three_arm_per_dataset.csv")

    d = panel.merge(acc, left_on=["protein", "cell"], right_on=["protein", "cell_line"],
                    how="left")
    missing = d[d.accession.isna()]
    if len(missing):
        sys.exit(f"{len(missing)} panel datasets have no ENCODE accession: "
                 f"{', '.join(missing.dataset.head(10))}. A supplementary table with silent "
                 f"gaps looks complete and is not; refusing to write.")

    d = d.merge(arms[["dataset", "comp_gc", "comp_dn", "comp_neg2",
                      "gain_gc", "gain_dn", "gain_neg2"]], on="dataset", how="left")
    d["in_three_arm_panel"] = d.comp_gc.notna()

    cols = ["dataset", "protein", "cell", "accession", "experiment", "n_replicates", "pairs",
            "in_three_arm_panel", "comp_gc", "comp_dn", "comp_neg2",
            "gain_gc", "gain_dn", "gain_neg2"]
    d = d[cols].sort_values(["protein", "cell"])
    d.to_csv(OUT, index=False)

    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"  {len(d)} datasets, {d.protein.nunique()} proteins, "
          f"{d.experiment.nunique()} ENCODE experiments")
    print(f"  {int(d.in_three_arm_panel.sum())} carry all three protocols; "
          f"{int((~d.in_three_arm_panel).sum())} are ladder-only (R2)")
    print(f"  replicates per experiment: {d.n_replicates.min()}-{d.n_replicates.max()}")


if __name__ == "__main__":
    main()
