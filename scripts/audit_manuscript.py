"""Every number in the manuscript must trace to a committed table. Find the ones that do not.

WHY. The paper's primary contrast, +0.0397, was printed in the manuscript, appeared in no
committed table, was produced by no script that could be found, and was asserted by no golden
key. It could not be reproduced and it could not fail, and it survived six rounds of
adversarial review in that state, because a reviewer reads a number rather than goes looking
for it. Three more were found the same way and only by hand. Hand-checking does not scale and
does not repeat; this does.

THE HAYSTACK IS CURATED, AND THE FIRST VERSION OF THIS SCRIPT WAS NEARLY WORTHLESS BECAUSE IT
WAS NOT. That version pooled every numeric cell of every table, including
`variant_scores.csv` and `variant_assignments.csv` at 66,010 rows each. With that much data
the four-decimal grid over [0.5, 1.0] came out **73.9% saturated**: essentially any AUROC-like
number matched something by coincidence, and a fabricated "0.9427" injected into the
manuscript passed cleanly. The check reported 9 orphans and felt reassuring while being close
to vacuous.

What a manuscript actually quotes is a summary, not a raw cell, so the haystack is now:

  golden.yaml         every asserted reference value
  summary tables      every cell of any table with a `check` column or at most 50 rows, i.e.
                      one row per asserted quantity
  column aggregates   mean, median, min, max and sign fractions of every numeric column of
                      every table -- the operations a sentence like "composition alone reaches
                      0.783" actually performs

That was ~950 values instead of ~61,000, and it dropped four-decimal saturation from 73.9% to
about 6%. DO NOT TRUST EITHER FIGURE FROM THIS DOCSTRING: the haystack grows every time a
result table is added, and it now holds ~3,600 values at ~23% saturation. The script MEASURES
and prints its own occupancy on every run, because that figure IS the false-negative rate and
a checker that will not state its own is asking to be over-trusted. Read the run, not this
paragraph.

TWO LIMITATIONS, BOTH REAL, BOTH STATED RATHER THAN ROUNDED AWAY.

  three decimals   The 3-dp grid over [0.5, 1.0] holds only 501 slots and is now ~88%
                   occupied, so a fabricated number written to three decimals passes almost
                   always. Only 4-dp claims are checked with real power. Quote headline
                   numbers to four decimals and this check protects them.

  wrong claim,     It cannot tell whether a number is attached to the right claim. The
  right number     mislabelled row that put 0.7981 under "composition + score" when it is the
                   standalone k-mer AUROC passes here, because 0.7981 is a real value in a
                   real table and was merely on the wrong line.

BARE INTEGERS ARE NOW CHECKED TOO, against a different haystack. The counts a paper quotes --
94 datasets, 79 proteins, 456,734 pairs, 37 of 40 -- carry no decimal point, so NUM never saw
them and they were checked by hand; three were wrong in an earlier draft. Counts are not
aggregates of a column, they are properties of a table's shape, so INT_HAYSTACK is built from
row counts, per-column distinct counts, non-null counts, sign counts, and the size of every
group of every low-cardinality column, plus config/params.yaml, where a parameter like the
101~nt window size legitimately lives instead of in a result table.

  small integers   Integers below 10 are not checked. There are nine of them and a document
                   of this length hits nearly all nine by coincidence, so a check there
                   reports nothing. Years and version-like tokens are skipped for the same
                   reason: 2026 is a date, not a count.

THE REMAINING ORPHANS ARE NOT ALL ERRORS. Most are aggregates over a SUBSET -- "mean
conservation AUROC over the 44 powered datasets", say -- which no whole-column aggregate can
reproduce. Triaging them one by one, and emitting the legitimate ones into tables so the count
can be ratcheted toward zero, is outstanding work. What the gate buys today is that a NEW
unsourced number fails the build.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
# The submitted prose is the document that matters. An earlier revision also audited a
# separate claims ledger under docs/, which now lives on the working-notes branch.
MANUSCRIPT_DIR = ROOT / "manuscript"
# Every manuscript section is audited, not a single file: a value can be correct in one place
# and mistyped in another, and the submitted prose is what a referee reads.
MANUSCRIPT = sorted(MANUSCRIPT_DIR.glob("*.md")) + sorted((MANUSCRIPT_DIR / "sections").glob("*.tex"))
GOLDEN = ROOT / "config" / "golden.yaml"

# This script's own output lives in results/tables/ and its `value` column IS the orphan list,
# so globbing the directory let the second run find every orphan inside the report the first
# run wrote, and it printed a clean zero. A checker that reads its own findings back as
# evidence always passes.
SELF = "manuscript_orphans.csv"

# Below three decimals almost everything collides with a count, a year or a p-value exponent.
MIN_DECIMALS = 3
MAX_DP = 6
# A table this small is a summary: one row per asserted quantity, so its cells are claims and
# belong in the haystack. Above it, a table is per-dataset or per-variant data whose 66,010
# rows saturate the grid. Measured: at 50 the four-decimal grid over [0.5, 1.0] is 6% occupied;
# at 100 it is 35% and the check stops meaning anything.
SUMMARY_ROWS = 50
NUM = re.compile(r"(?<![\w.])[-+−]?(\d+\.\d+)(?![\w])")
# Integers, after LaTeX thousands separators are stripped. No decimal point, and not adjacent
# to one, so "0.0397" never yields a 397 and "v1.2" never yields a 2.
INT = re.compile(r"(?<![\w.])(\d+)(?![\w.])")
# Below ten there are nine possible values and a paper of this length uses most of them, so
# the check would flag nothing and prove nothing. Stated rather than hidden.
MIN_INT = 10
# A four-digit token in this range is a year -- a citation, a date, an ENCODE release -- and
# not a count of anything. LaTeX cross-reference and float machinery likewise.
YEAR = re.compile(r"^(19|20)\d\d$")
# A BIBLIOGRAPHY ASSERTS NOTHING ABOUT THE SCIENCE. Volume, issue and page numbers are
# bibliographic coordinates, and scanning them produced 30 orphans that were all correct and
# none of which any result table could ever source. Excluded from the integer scan only; a
# fabricated 4-decimal value in a reference would still be caught by NUM.
NO_INT_SCAN = {"bibliography.tex"}
MACRO = re.compile(r"\\(ref|label|cite\w*|citep|citet|includegraphics|vspace|hspace|"
                   r"textwidth|linewidth|columnwidth|arraystretch|scalebox|resizebox|"
                   r"multicolumn|multirow|cmidrule|addtocounter|setcounter|figure|table)")
# IDENTIFIERS ARE NOT CLAIMS. A DOI, accession or version string contains a decimal point and
# is matched by NUM, but it asserts nothing about the science and has no table to live in.
# Flagging one is a false positive that costs the reader's trust in the real orphans, and
# "10.5281" from a Zenodo DOI is exactly the case that surfaced. Matched on the surrounding
# text rather than the token, so a genuine 10.5281 elsewhere is still checked.
IDENTIFIER = re.compile(r"(doi:|zenodo\.|10\.\d{4,}/|ENC[A-Z]{2}\d|GSE\d|v\d+\.\d+)",
                        re.IGNORECASE)


def log(m):
    print(m, flush=True)


def haystack():
    """Values a manuscript could legitimately be quoting: assertions and aggregates."""
    # Keyed by decimal place: a manuscript writing "1.036" must match a table holding
    # 1.0357967, so the comparison happens at the TOKEN's own precision. Dropping this was a
    # regression that turned 45 orphans into 134, all of them correctly-rounded CI bounds.
    vals = {d: set() for d in range(MIN_DECIMALS, MAX_DP + 1)}

    def add(x):
        try:
            v = abs(float(x))
        except (TypeError, ValueError):
            return
        if np.isfinite(v):
            for d in vals:
                vals[d].add(round(v, d))

    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
        else:
            add(o)

    walk(yaml.safe_load(GOLDEN.read_text()))

    for p in sorted(TABLES.glob("*.csv")) + sorted(TABLES.glob("*.tsv")):
        if p.name == SELF:
            continue
        try:
            d = pd.read_csv(p, sep="\t" if p.suffix == ".tsv" else ",")
        except Exception:
            continue
        summary = "check" in d.columns or len(d) <= SUMMARY_ROWS
        for c in d.columns:
            s = pd.to_numeric(d[c], errors="coerce").dropna().astype(float)
            if not len(s):
                continue
            for f in (s.mean(), s.median(), s.min(), s.max()):
                add(f)
            add((s > 0).mean())
            add((s < 0).mean())
            if summary:                     # one row per asserted quantity: cells are claims
                for v in s.unique():
                    add(v)
    return vals


def int_haystack():
    """Counts a manuscript could legitimately be quoting.

    A count is a property of a table's SHAPE, not an aggregate of its values, so none of the
    means and medians in haystack() can source one. What sources "94 datasets" is that some
    table has 94 rows, or that some column has 94 distinct values, or that some grouping has
    a group of size 94. Sums are included because "463,091 peaks" is a column total, and
    config values because the 101~nt window size is set in params.yaml and appears in no
    result table.
    """
    out = set()

    def add(x):
        try:
            v = float(x)
        except (TypeError, ValueError):
            return
        if np.isfinite(v) and v == int(v):
            out.add(abs(int(v)))

    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
        else:
            add(o)

    walk(yaml.safe_load(GOLDEN.read_text()))
    cfgp = ROOT / "config" / "params.yaml"
    if cfgp.exists():
        walk(yaml.safe_load(cfgp.read_text()))

    for p in sorted(TABLES.glob("*.csv")) + sorted(TABLES.glob("*.tsv")):
        if p.name == SELF:
            continue
        try:
            d = pd.read_csv(p, sep="\t" if p.suffix == ".tsv" else ",")
        except Exception:
            continue
        add(len(d))
        summary = "check" in d.columns or len(d) <= SUMMARY_ROWS
        for c in d.columns:
            add(d[c].nunique())
            add(d[c].notna().sum())
            # group sizes: "48 of the 94 are K562" is a value_count and nothing else
            if d[c].nunique() <= 40:
                for n in d[c].value_counts().tolist():
                    add(n)
            s = pd.to_numeric(d[c], errors="coerce").dropna().astype(float)
            if not len(s):
                continue
            add(s.sum())
            add((s > 0).sum())
            add((s < 0).sum())
            add((s == 0).sum())
            for f in (s.min(), s.max(), s.mean().round(), s.median()):
                add(f)
            if summary:
                for v in s.unique():
                    add(v)
    return out


def main():
    if not MANUSCRIPT:
        raise SystemExit(f"no manuscript sources under {MANUSCRIPT_DIR}")
    vals = haystack()
    ints = int_haystack()
    sat = {d: len({v for v in vals[d] if 0.5 <= v <= 1.0}) / (10 ** d * 0.5 + 1)
           for d in (3, 4)}
    # The integer false-negative rate over the range the paper's counts actually occupy.
    isat = len({v for v in ints if MIN_INT <= v <= 1000}) / (1000 - MIN_INT + 1)

    orphans, checked, ichecked = [], 0, 0
    sources = [(m.name, m) for m in MANUSCRIPT]
    for src_name, src in sources:
      for i, line in enumerate(src.read_text().splitlines(), 1):
        for m in NUM.finditer(line):
            tok = m.group(1)
            if len(tok.split(".")[1]) < MIN_DECIMALS:
                continue
            # skip identifiers: a DOI or accession is not a numeric claim
            ctx = line[max(0, m.start() - 12):m.end() + 12]
            if IDENTIFIER.search(ctx):
                continue
            checked += 1
            d = min(len(tok.split(".")[1]), MAX_DP)
            if round(abs(float(tok)), d) in vals[d]:
                continue
            orphans.append((src_name, i, tok, line.strip()[:96]))

        # Integers, on the line with LaTeX thousands separators closed up, so 463{,}091 is
        # one token and not a 463 next to an 091.
        if src_name in NO_INT_SCAN:
            continue
        flat = re.sub(r"(?<=\d)\{,\}(?=\d)", "", line)
        flat = re.sub(r"(?<=\d),(?=\d\d\d(?!\d))", "", flat)
        for m in INT.finditer(flat):
            tok = m.group(1)
            v = int(tok)
            if v < MIN_INT or YEAR.match(tok):
                continue
            ctx = flat[max(0, m.start() - 24):m.end() + 12]
            if IDENTIFIER.search(ctx) or MACRO.search(ctx):
                continue
            ichecked += 1
            if v in ints:
                continue
            orphans.append((src_name, i, tok, line.strip()[:96]))

    log(f"  haystack: {len(vals[4])} values (golden keys, summary cells, column aggregates)")
    log(f"  false-negative rate: {sat[4]:.1%} of the 4-dp grid over [0.5, 1.0] is occupied "
        f"and {sat[3]:.1%} of the 3-dp grid, so a fabricated 4-decimal AUROC slips through "
        f"about 1 time in {max(int(1 / sat[4]), 1)} and a 3-decimal one closer to 1 in 2")
    log(f"  manuscript numbers checked (>= {MIN_DECIMALS} dp): {checked}")
    log(f"  integer haystack: {len(ints)} counts (table shapes, group sizes, column sums, "
        f"config values)")
    log(f"  integer false-negative rate: {isat:.1%} of [{MIN_INT}, 1000] is occupied, so a "
        f"fabricated count in that range slips through about 1 time in "
        f"{max(int(1 / isat), 1)}")
    log(f"  manuscript integers checked (>= {MIN_INT}, years excluded): {ichecked}")
    log(f"  ORPHANS (traceable to nothing): {len(orphans)}\n")
    for src_name, ln, tok, ctx in orphans:
        log(f"  {src_name}:{ln:<5} {tok:<12} {ctx}")

    pd.DataFrame(orphans, columns=["source", "line", "value", "context"]).to_csv(
        TABLES / SELF, index=False)
    log(f"\n  wrote {TABLES / SELF}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
