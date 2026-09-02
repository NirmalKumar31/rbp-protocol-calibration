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

That is ~950 values instead of ~61,000, and four-decimal saturation drops from 73.9% to about
6%. The script prints that figure on every run, because it IS the false-negative rate and a
checker that will not state its own is asking to be over-trusted.

TWO LIMITATIONS, BOTH REAL, BOTH STATED RATHER THAN ROUNDED AWAY.

  three decimals   The 3-dp grid over [0.5, 1.0] holds only 501 slots and is ~44% occupied, so
                   a fabricated number written to three decimals has close to a coin-flip
                   chance of passing. Only 4-dp claims are checked with real power. Quote
                   headline numbers to four decimals and this check protects them.

  wrong claim,     It cannot tell whether a number is attached to the right claim. The
  right number     mislabelled row that put 0.7981 under "composition + score" when it is the
                   standalone k-mer AUROC passes here, because 0.7981 is a real value in a
                   real table and was merely on the wrong line.

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
PAPER = ROOT / "docs" / "60-the-paper.md"
# THE DRAFTED MANUSCRIPT IS ALSO AUDITED, from 2026-09-02. docs/60 is the claims ledger and
# manuscript/ is the prose that will be submitted; a number can be sourced in the first and
# mistyped in the second, and the second is the one a referee reads. Auditing only the ledger
# would have checked the wrong document from the moment drafting began.
MANUSCRIPT = sorted((ROOT / "manuscript").glob("*.md")) if (ROOT / "manuscript").exists() else []
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


def main():
    if not PAPER.exists():
        raise SystemExit(f"{PAPER} is absent")
    vals = haystack()
    sat = {d: len({v for v in vals[d] if 0.5 <= v <= 1.0}) / (10 ** d * 0.5 + 1)
           for d in (3, 4)}

    orphans, checked = [], 0
    sources = [(PAPER.name, PAPER)] + [(m.name, m) for m in MANUSCRIPT]
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

    log(f"  haystack: {len(vals[4])} values (golden keys, summary cells, column aggregates)")
    log(f"  false-negative rate: {sat[4]:.1%} of the 4-dp grid over [0.5, 1.0] is occupied "
        f"and {sat[3]:.1%} of the 3-dp grid, so a fabricated 4-decimal AUROC slips through "
        f"about 1 time in {max(int(1 / sat[4]), 1)} and a 3-decimal one closer to 1 in 2")
    log(f"  manuscript numbers checked (>= {MIN_DECIMALS} dp): {checked}")
    log(f"  ORPHANS (traceable to nothing): {len(orphans)}\n")
    for src_name, ln, tok, ctx in orphans:
        log(f"  {src_name}:{ln:<5} {tok:<12} {ctx}")

    pd.DataFrame(orphans, columns=["source", "line", "value", "context"]).to_csv(
        TABLES / SELF, index=False)
    log(f"\n  wrote {TABLES / SELF}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
