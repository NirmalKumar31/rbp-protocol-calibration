"""A column dictionary for every released table, generated rather than hand-maintained.

    python scripts/column_dictionary.py

WHY GENERATED. results/tables/SCHEMA.md described "124 CSVs, two shapes", of which 44 were
summary tables. An audit counted the tree: 135 files, 88 distinct column schemas, 57 with a
`check` column. Every one of those numbers was hand-typed and every one had drifted, in the
document whose job is to tell a reader what the columns mean.

Two shapes cover most files and the prose in SCHEMA.md explains them. This emits the rest: one
row per file per column, with dtype and an example value, so a reader of an unusual table --
per-fold, transport, matching, variant, task -- is not left inferring the columns from code.
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "COLUMNS.csv"


def main():
    rows = []
    files = sorted(list(TABLES.glob("*.csv")) + list(TABLES.glob("*.tsv"))
                   + list(TABLES.glob("*/*.csv")), key=lambda q: str(q))
    for f in files:
        if f.name in ("COLUMNS.csv", "PROVENANCE.csv"):
            continue
        delim = "\t" if f.suffix == ".tsv" else ","
        with f.open() as fh:
            r = csv.reader(fh, delimiter=delim)
            header = next(r, [])
            first = next(r, [])
        for i, col in enumerate(header):
            v = first[i] if i < len(first) else ""
            try:
                float(v)
                kind = "int" if v.strip().lstrip("-").isdigit() else "float"
            except ValueError:
                kind = "string" if v else "empty"
            rows.append({"table": str(f.relative_to(TABLES)), "column": col,
                         "dtype": kind, "example": v[:40]})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["table", "column", "dtype", "example"])
        w.writeheader()
        w.writerows(rows)
    schemas = len({tuple(x["column"] for x in rows if x["table"] == t)
                   for t in {r["table"] for r in rows}})
    log(f"  {len({r['table'] for r in rows})} tables, {len(rows)} columns, "
        f"{schemas} distinct schemas -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
