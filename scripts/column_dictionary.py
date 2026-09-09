"""A semantic data dictionary for every released table, generated rather than hand-maintained.

    python scripts/column_dictionary.py
    python scripts/column_dictionary.py --check      # fail if COLUMNS.csv is stale

WHY GENERATED. results/tables/SCHEMA.md described "124 CSVs, two shapes", of which 44 were
summary tables. An audit counted the tree: 135 files, 88 distinct column schemas, 57 with a
`check` column. Every one of those numbers was hand-typed and every one had drifted, in the
document whose job is to tell a reader what the columns mean.

Why it was rewritten. A second audit read the output and found it was not a data dictionary.
It carried four fields -- table, column, dtype, example -- where the dtype came from THE FIRST
DATA ROW ALONE and there was no definition, unit, key or missingness rule anywhere in it, while
It was advertised as "column definitions for every released table". The first-row
inference was not merely thin, it was wrong: 98 columns were typed `empty` and 87 of those 98
are populated further down the same file. A summary table's `ci_low` is blank on its first row
whenever that row is a count, so the single most common interval column in the release was
typed as having no data.

What a column entry now carries. dtype inferred from every non-empty value in the column;
n_rows, n_missing and the missing rule that applies; n_distinct; min and max for numerics;
whether the column is part of the table's key; a unit; and a definition. Definitions come from
DEFS below, matched on the column name, because the release uses a small controlled vocabulary
across 90 schemas: the same `check`/`value`/`ci_low` shape and the same `{quantity}_{arm}`
convention recur everywhere. Columns DEFS does not match are emitted with an empty definition
rather than a guessed one, and `--check` reports how many, so the gap is visible and countable
instead of hidden behind an example value.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "COLUMNS.csv"
SUMMARY = TABLES / "COLUMNS_SUMMARY.csv"
# The auditors' own output is not part of the release's data, and including it made the
# dictionary churn against itself: audit_manuscript.py rewrites manuscript_orphans.csv, whose
# `source` column is empty when there are no orphans and populated when there are, so a run
# that found none typed the column `empty` and the next run that found one failed the
# staleness gate. Same set provenance.py calls META, for the same reason.
SKIP = {"COLUMNS.csv", "COLUMNS_SUMMARY.csv", "PROVENANCE.csv", "manuscript_orphans.csv",
        "release_facts.csv", "verify_summary.csv"}

FIELDS = ["table", "column", "dtype", "unit", "key", "n_rows", "n_missing", "n_distinct",
          "min", "max", "example", "definition", "producing_script"]

# Every column points at code, even where prose does not cover it. 512 of the column names in
# the release are analysis-specific and appear once or twice; writing a sentence for each by
# hand is how the counts in SCHEMA.md went stale in the first place. The controlled vocabulary
# above covers the recurring shapes, and PROVENANCE.csv supplies the producing script for the
# rest, so a reader of an unusual column has somewhere definite to look rather than a guess.
PROVENANCE = TABLES / "PROVENANCE.csv"

# The arm suffixes, spelled once. SCHEMA.md's table is the prose version of this.
ARMS = {
    "gc": "GC-matched negatives",
    "dn": "dinucleotide-matched negatives",
    "dinuc": "dinucleotide-matched negatives",
    "neg2": "bias-aware negatives (other RBPs' binding sites)",
    "neg2_rm": "bias-aware negatives, donor draw stratified on transcript region",
    "shuf": "dinucleotide-shuffled negatives",
}

# (regex on the column name, unit, definition). First match wins, so order matters: put the
# exact names before the patterns that would also catch them.
DEFS = [
    (r"^check$", "", "The asserted quantity, in words. Primary key within the file, and the "
                     "string scripts/verify.py looks the row up by, so renaming one is a "
                     "breaking change"),
    (r"^value$", "", "The asserted quantity"),
    (r"^ci_low$", "", "Lower bound of a 95% interval. Protein-clustered percentile bootstrap, "
                      "4000 draws, unless the note says otherwise. Empty where the quantity "
                      "is a count or a ratio of means and has no resampling distribution"),
    (r"^ci_high$", "", "Upper bound of the same interval as ci_low"),
    (r"^n$", "count", "Datasets, or rows, behind the value. Usually 94"),
    (r"^note$", "", "Scope, caveat, or the definition of a non-obvious quantity. Empty where "
                    "there is nothing to say"),
    (r"^dataset$", "", "PROTEIN:CELL. The primary key of a per-dataset table"),
    (r"^protein$", "", "ENCODE target. NOT unique: 15 proteins appear in both cell lines, "
                       "which is why intervals are protein-clustered"),
    (r"^cell$", "", "K562 or HepG2"),
    (r"^accession$", "", "ENCODE file accession"),
    (r"^experiment$", "", "ENCODE experiment accession"),
    (r"^(pairs|n_pairs)$", "count", "Matched positive/negative pairs; rows are twice this"),
    (r"^fold$", "", "Chromosome-blocked fold index, 0 to 4, from config/folds.tsv"),
    (r"^chrom$", "", "Chromosome, GRCh38 primary assembly"),
    (r"^(start|end)$", "bp", "Window coordinate, 0-based half-open, GRCh38"),
    (r"^strand$", "", "+ or -"),
    (r"^label$", "", "1 for a positive window, 0 for its matched negative"),
    (r"^score$", "", "Out-of-fold model score for this window"),
    (r"^seed$", "", "Random seed for this draw or run"),
    (r"^arm$", "", "Negative-set protocol; see the arm suffixes in SCHEMA.md"),
    (r"^model$", "", "Model class: kmer, cnn or splicebert"),
    (r"^k$", "", "k-mer order"),
    (r"^jaccard$", "", "Intersection over union of two positive sets. Not a count ratio; an "
                       "earlier version of this column was min(n)/max(n) and is corrected"),
    (r"^n_pos_", "count", "Retained positives in that arm"),
    (r"^comp_", "AUROC", "Composition-only pooled out-of-fold AUROC in that arm: 19 features, "
                         "refit on that arm's own windows and folds"),
    (r"^(gain|delta)_", "AUROC", "Nested contribution in that arm: AUROC(19 composition "
                                 "features + model score) minus AUROC(composition alone), both "
                                 "pooled out of fold. Dimensionless, on [-1, 1]"),
    (r"^auroc", "AUROC", "Pooled out-of-fold AUROC, not a mean over per-fold AUROCs"),
    (r"^n_", "count", "Row count for that arm, which is twice the pair count"),
    (r"^p(_value)?$", "", "Two-sided p-value"),
    (r"^(se|stderr)", "", "Standard error"),
    (r"_(usd|cost)$", "USD", "Cost in US dollars; see docs/COST.md"),
    (r"^(gpu_)?(sec|seconds|secs)", "s", "Wall-clock seconds"),
    (r"^(md5|sha256|checksum)", "", "Content hash of the file named in this row"),
    (r"^(size|bytes|size_bytes)$", "bytes", "File size in bytes"),
    # The order-of-baseline family, which recurs across the baseline_order tables.
    (r"^comp([23])_", "AUROC", "Composition-only pooled out-of-fold AUROC in that arm at the "
                               "stated baseline order (comp2 = order two, 19 features; comp3 = "
                               "order three, 84 features)"),
    (r"^composition_auroc", "AUROC", "Composition-only pooled out-of-fold AUROC"),
    (r"^full_", "AUROC", "Pooled out-of-fold AUROC of the full nested model in that arm: the "
                         "composition features plus the model score"),
    (r"^(kmer|cnn|splicebert)_gain([23])_", "AUROC",
     "Nested contribution of that model class over a composition baseline of the stated order, "
     "in that arm, pooled out of fold"),
    (r"^(kmer|cnn|splicebert)_published_", "AUROC",
     "The published nested contribution for that model class in that arm, as it appears in the "
     "manuscript, carried alongside a recomputation for comparison"),
    (r"^(kmer|cnn|splicebert)_", "", "Quantity for that model class; see the producing script"),
    (r"^gain$", "AUROC", "Nested contribution: AUROC(composition features + model score) minus "
                         "AUROC(composition alone), both pooled out of fold"),
    (r"^coef", "", "Fitted coefficient"),
    (r"^git_sha$", "", "Short commit hash of the tree that produced the row"),
    (r"^cost", "USD", "Cost in US dollars; see docs/COST.md"),
    (r"^(shift|diff|delta)$", "", "Signed change in the quantity this table reports"),
    (r"^(frac|fraction|share|pct|proportion)", "", "A proportion on [0, 1] unless the note "
                                                   "says the column is a percentage"),
    (r"^(spearman|pearson|rho|corr)", "", "Correlation coefficient"),
]

# Columns that identify a row rather than measure one.
KEYS = {"check", "dataset", "protein", "cell", "fold", "arm", "model", "seed", "k", "table",
        "column", "accession", "experiment", "chrom", "start", "end", "strand", "path"}


def describe(col):
    for pat, unit, defn in DEFS:
        if re.search(pat, col):
            m = re.search(r"_(gc|dn|dinuc|neg2_rm|neg2|shuf)$", col)
            if m and unit == "AUROC":
                return unit, f"{defn.split(' in that arm')[0]} under {ARMS[m.group(1)]}" \
                    if " in that arm" in defn else defn
            return unit, defn
    return "", ""


def profile(values):
    """dtype, missing count, distinct count, min and max, over the WHOLE column."""
    present = [v for v in values if v.strip() != ""]
    n_missing = len(values) - len(present)
    if not present:
        return "empty", n_missing, 0, "", ""
    nums, ints = [], True
    for v in present:
        try:
            nums.append(float(v))
        except ValueError:
            return "string", n_missing, len(set(present)), "", ""
        if not v.strip().lstrip("-").isdigit():
            ints = False
    lo, hi = min(nums), max(nums)
    fmt = (lambda x: str(int(x))) if ints else (lambda x: f"{x:.6g}")
    return ("int" if ints else "float"), n_missing, len(set(present)), fmt(lo), fmt(hi)


def producers():
    """table -> producing script, read from PROVENANCE.csv.

    THIS MAKES THE TWO GENERATORS CIRCULAR, so they need two passes. provenance.py hashes
    COLUMNS.csv, and this reads provenance.py's output, so a newly committed table comes out
    with an empty producing_script on the first pass and provenance.py then records the hash of
    that incomplete file. Run column_dictionary, provenance, column_dictionary, provenance: the
    producer map is fixed after the first provenance run, so the second pass converges.

    Adding ten draw tables at once is what surfaced it. CI catches the half-done state, which is
    the right place for it, but it costs a red build to find out. Run
    `scripts/refresh_manifests.sh` rather than calling either generator by hand.
    """
    if not PROVENANCE.exists():
        return {}
    with PROVENANCE.open(newline="") as fh:
        return {r["table"]: r.get("producing_script", "") for r in csv.DictReader(fh)}


def build():
    rows = []
    prod = producers()
    files = sorted(list(TABLES.glob("*.csv")) + list(TABLES.glob("*.tsv"))
                   + list(TABLES.glob("*/*.csv")), key=str)
    for f in files:
        if f.name in SKIP:
            continue
        delim = "\t" if f.suffix == ".tsv" else ","
        with f.open(newline="") as fh:
            r = csv.reader(fh, delimiter=delim)
            header = next(r, [])
            data = list(r)
        for i, col in enumerate(header):
            values = [(d[i] if i < len(d) else "") for d in data]
            dtype, n_missing, n_distinct, lo, hi = profile(values)
            unit, defn = describe(col)
            example = next((v for v in values if v.strip()), "")
            rows.append({
                "table": str(f.relative_to(TABLES)), "column": col, "dtype": dtype,
                "unit": unit, "key": "yes" if col in KEYS else "",
                "n_rows": len(data), "n_missing": n_missing, "n_distinct": n_distinct,
                "min": lo, "max": hi, "example": example[:40], "definition": defn,
                "producing_script": prod.get(str(f.relative_to(TABLES)), ""),
            })
    return rows


def counts(rows):
    tables = {r["table"] for r in rows}
    schemas = {tuple(x["column"] for x in rows if x["table"] == t) for t in tables}
    defined = sum(1 for r in rows if r["definition"])
    return {"tables": len(tables), "columns": len(rows), "schemas": len(schemas),
            "columns_with_a_definition": defined,
            "columns_without_a_definition": len(rows) - defined,
            "columns_with_a_producing_script": sum(1 for r in rows if r["producing_script"]),
            "columns_typed_empty": sum(1 for r in rows if r["dtype"] == "empty")}


def write(rows):
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    # The counts are an artefact, not prose. SCHEMA.md hand-typed "135 tables, 1393 columns,
    # 90 distinct schemas" in the same paragraph that said counts there are "either generated
    # or absent", and all three were wrong. SCHEMA.md now points at this file.
    c = counts(rows)
    with SUMMARY.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["check", "value"])
        for k, v in c.items():
            w.writerow([k.replace("_", " "), v])
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="fail if the committed dictionary differs from a fresh build")
    a = ap.parse_args()
    rows = build()
    c = counts(rows)
    if a.check:
        if not OUT.exists():
            log("  COLUMNS.csv is missing")
            sys.exit(1)
        with OUT.open(newline="") as fh:
            have = list(csv.DictReader(fh))
        fresh = [{k: str(v) for k, v in r.items()} for r in rows]
        if [dict(r) for r in have] != fresh:
            log(f"  COLUMNS.csv is stale: {len(have)} committed rows against {len(fresh)} "
                "fresh. Rerun scripts/column_dictionary.py")
            sys.exit(1)
        log(f"  COLUMNS.csv matches: {c['tables']} tables, {c['columns']} columns")
        return
    c = write(rows)
    log(f"  {c['tables']} tables, {c['columns']} columns, {c['schemas']} distinct schemas "
        f"-> {OUT.relative_to(ROOT)}")
    log(f"  {c['columns_with_a_definition']} columns carry a definition, "
        f"{c['columns_without_a_definition']} do not; {c['columns_typed_empty']} are empty "
        "in every row")


if __name__ == "__main__":
    main()
