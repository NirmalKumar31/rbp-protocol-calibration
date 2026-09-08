"""Manuscript numbers written to fewer than three decimals, traced to their evidence row.

Why this exists, and it closes the hole that let fourteen wrong numbers into the paper.

`scripts/audit_manuscript.py` checks every manuscript number with **three or more decimal
places** against a haystack of table cells, and reports its own false-negative rate honestly:
96.0% of the 3-decimal grid over [0.5, 1.0] is already occupied, so a 3-decimal number slips
through about one time in two. Below three decimals it does not look at all, because the grid is
saturated and a match would mean nothing.

So `6.58`, `58.5\\%`, `57 of 94` and `$3.58$` were never traced by anything. On 2026-09-06 and
again on 2026-09-07, two commits bumped a stated page count with a blind substitution: `57`
became `58`, then `58` became `59`. In a regex, `\\b58\\b` matches the `58` inside `6.58`,
because `.` is not a word character. Fourteen unrelated scientific numbers were rewritten across
`README.md`, `docs/COST.md`, `paper.tex`, `methods.tex` and `results.tex`. Four of them were
bumped twice, so they ended two away from the truth. Every gate stayed green: 1086/1086 verify
checks, zero manuscript orphans, every derived count consistent, and CI green on all four jobs.
An external audit found the difference and concluded that the committed PDF was **stale**,
because the fresh build matched the source. The source was wrong and the stale PDF was right.

The fix cannot be "check numbers below three decimals against the haystack", which is what
audit_manuscript declines to do for a good reason. It has to name the row. Each entry below
pins one low-precision manuscript claim to the exact table row it comes from, recomputes the
value from the committed evidence, rounds it the way the sentence rounds it, and requires the
sentence to say that. Nothing here hard-codes a number: change the evidence and this test
follows it, which is the whole point.

Adding a claim here costs four lines and buys a number that cannot silently drift again.
"""

import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"


def _row(table, check, column="value"):
    """The single committed row for `check`. Fails loudly rather than skipping on a rename."""
    p = TABLES / table
    if not p.exists():
        pytest.skip(f"{table} not in this checkout")
    d = pd.read_csv(p)
    hit = d[d["check"] == check]
    assert len(hit) == 1, (
        f"{table}: {len(hit)} rows match {check!r}. A renamed row must be renamed here too, "
        "or this test silently stops covering the claim it exists to cover")
    return float(hit.iloc[0][column])


# check = the sentence in the manuscript, as a regex with ONE group capturing the number.
# The group must capture exactly what a reader sees, so the rounding is what is asserted.
CLAIMS = [
    # --- the two-stage span's upper confidence bound. Bumped once, in four places. ---
    ("README.md", r"\*\*5\.42-fold\*\* \(4\.43 to ([\d.]+)\)",
     "cross_fitting.csv", "4-mer three-arm span, as published", "ci_high", 1.0, 2),
    ("manuscript/paper.tex", r"5\.42-fold \(4\.43--([\d.]+)\) with",
     "cross_fitting.csv", "4-mer three-arm span, as published", "ci_high", 1.0, 2),
    ("manuscript/sections/methods.tex", r"= 5\.42\$ \(\$4\.43\$ to \$([\d.]+)\$\)",
     "cross_fitting.csv", "4-mer three-arm span, as published", "ci_high", 1.0, 2),
    ("manuscript/sections/results.tex", r"& 5\.42 \(4\.43--([\d.]+)\) \\\\",
     "cross_fitting.csv", "4-mer three-arm span, as published", "ci_high", 1.0, 2),
    # results.tex also states it in prose, and this instance ESCAPED both bumps, which is how
    # the paper came to disagree with itself: three places said 6.59 and this one said 6.58.
    ("manuscript/sections/results.tex", r"4\.43 to ([\d.]+), protein-clustered",
     "cross_fitting.csv", "4-mer three-arm span, as published", "ci_high", 1.0, 2),

    # --- the per-dataset decomposition interval. Bumped TWICE: 57 -> 58 -> 59. ---
    ("README.md", r"gives 63% \(CI (\d+) to 68\)",
     "protocol_transport.csv",
     "share of variance from the EVALUATION protocol, per-dataset weighting", "ci_low", 100.0, 0),

    # --- how well the composition baseline recovers the arm. Bumped once. ---
    ("manuscript/sections/results.tex", r"only ([\d.]+)\\% accuracy against a chance",
     "protocol_or_baseline.csv", "arm recovered from the composition baseline alone",
     "value", 100.0, 1),

    # --- the block-preserving null's own share. Bumped TWICE: 57.4 -> 58.4 -> 59.4. ---
    ("manuscript/sections/results.tex", r"absorbs ([\d.]+)\\% by itself",
     "multiplier_variance.csv", "block-preserving null share, protein", "value", 100.0, 1),

    # --- residual AUROC for the 4-mer on the bias-aware arm. Bumped once. ---
    ("manuscript/sections/results.tex", r"in every arm: 0\.63, 0\.60 and ([\d.]+) for",
     "estimands.csv", "residual_auroc, neg2 arm, kmer", "value", 1.0, 2),

    # --- datasets where difficulty and contribution move oppositely. Bumped TWICE. ---
    ("manuscript/sections/results.tex", r"and (\d+) of 94 for GC-to-bias-aware",
     "standalone_auroc.csv",
     "datasets where difficulty and contribution move oppositely, GC to bias-aware", "value",
     1.0, 0),
]


@pytest.mark.parametrize("path,pattern,table,check,column,scale,dp", CLAIMS,
                         ids=[f"{c[0].split('/')[-1]}:{c[3][:38]}" for c in CLAIMS])
def test_a_low_precision_manuscript_number_matches_its_evidence_row(
        path, pattern, table, check, column, scale, dp):
    f = ROOT / path
    if not f.exists():
        pytest.skip(f"{path} not in this checkout")
    m = re.search(pattern, f.read_text())
    assert m, (
        f"{path}: the sentence stating {check!r} no longer matches this pattern. Either it was "
        "reworded, in which case update the pattern, or it was deleted, in which case delete "
        "this claim. A pattern that stops matching must FAIL, not quietly cover nothing")
    stated = float(m.group(1))
    truth = round(_row(table, check, column) * scale, dp)
    assert stated == truth, (
        f"{path} states {stated} for {check!r}; the committed evidence in {table} gives "
        f"{truth}. This is the defect class two blind page-count substitutions created: a "
        "number below three decimals is traced by nothing else in this repository")


def test_the_disagreement_reductions_match_the_committed_pairs():
    """Three percentages in one sentence, each a reduction between two committed rows.

    The first was bumped once, 58 to 59. Kept separate from CLAIMS because the stated value is
    derived from a PAIR of rows rather than read from one.
    """
    f = ROOT / "manuscript" / "sections" / "results.tex"
    if not f.exists():
        pytest.skip("results.tex not in this checkout")
    m = re.search(r"in 3 of 3, by (\d+)\\%, (\d+)\\% and (\d+)\\%", f.read_text())
    assert m, "the disagreement-reduction sentence no longer matches"
    stated = [int(g) for g in m.groups()]

    d = pd.read_csv(TABLES / "recommendation_works.csv").set_index("check")
    truth = []
    for pair in ("gc vs dn", "gc vs neg2", "dn vs neg2"):
        raw = float(d.loc[f"scale-free disagreement, raw, {pair}", "value"])
        headroom = float(d.loc[f"scale-free disagreement, headroom, {pair}", "value"])
        truth.append(round(100 * (1 - headroom / raw)))
    assert stated == truth, f"results.tex says {stated}, the evidence gives {truth}"


def test_the_minimum_stratified_span_matches_the_sensitivity_suite():
    """"no model's span falls below X in either cell line". Bumped once, 3.58 to 3.59."""
    f = ROOT / "manuscript" / "sections" / "results.tex"
    if not f.exists():
        pytest.skip("results.tex not in this checkout")
    m = re.search(r"falls below \$([\d.]+)\$ in either", f.read_text())
    assert m, "the stratified-span sentence no longer matches"
    d = pd.read_csv(TABLES / "sensitivity_suite.csv")
    rows = d[d.check.str.contains("HepG2 only|K562 only") & d.check.str.startswith("span,")]
    assert len(rows) == 6, f"expected 6 cell-line spans, found {len(rows)}"
    truth = round(rows.value.astype(float).min(), 2)
    assert float(m.group(1)) == truth, f"results.tex says {m.group(1)}, evidence gives {truth}"


def test_one_sweep_of_both_neural_models_over_three_arms():
    """The dollar figure a reader can check, recomputed from committed pair counts.

    Bumped TWICE, 57 to 58 to 59, which also broke the cost argument: at $59 the "ten times
    that" clause implies $590, against the $115 derived elsewhere in the same section.
    """
    f = ROOT / "manuscript" / "sections" / "methods.tex"
    if not f.exists():
        pytest.skip("methods.tex not in this checkout")
    text = f.read_text()
    rates = re.search(r"\\\$([\d.]+) per million pairs for the\s*\n?\s*convolutional network "
                      r"and \\\$([\d.]+) for SpliceBERT", text)
    assert rates, "the measured-rate sentence no longer matches"
    cnn, sbert = float(rates.group(1)), float(rates.group(2))
    stated = re.search(r"arms is about \\\$(\d+)", text)
    assert stated, "the sweep-cost sentence no longer matches"

    mpairs = 0.0
    for t in ("rehearsal_binding_dinuc.csv", "rehearsal_binding_gc.csv", "neg2_build.csv"):
        p = TABLES / t
        if not p.exists():
            pytest.skip(f"{t} not in this checkout")
        mpairs += pd.read_csv(p)["pairs"].sum() / 1e6
    truth = round((cnn + sbert) * mpairs)
    assert int(stated.group(1)) == truth, (
        f"methods.tex says ${stated.group(1)} for one sweep of both models over three arms; "
        f"{cnn} + {sbert} per Mpair over {mpairs:.4f} Mpairs gives ${truth}")


# --- the neural cross-fitting cost, which five files quoted four different ways -------------

# Files that may legitimately contain an obsolete figure BECAUSE they record it as an error.
CORRECTION_RECORDS = ("CHANGELOG.md", "docs/EXTERNAL_CORRECTION_1.md")

OBSOLETE = ("$570", "$573", "$590", "four times the GPU sweep", "ten times that")


def test_the_neural_cross_fitting_cost_is_one_figure_everywhere():
    """P0-2. The repository quoted $115, $570, $590, "ten times that" and "four times the GPU
    sweep" for the same unrun analysis, in five files, at the same time.

    "Four times the GPU sweep" names no sweep and was converted to $76 on one reading of it.
    "Ten times that" sat twenty-two lines above the $115 derivation in the SAME Methods section.
    A reader checking the paper's justification for omitting an analysis found the justification
    contradicting itself by a factor of five.

    The supported figure is $115: ten extra base fits per dataset where a sweep already does
    five, so twice the base cost, on each of three arms. It is checked here against the $57
    one-sweep figure, which test_one_sweep_of_both_neural_models_over_three_arms derives
    independently from the measured rates and the committed pair counts.
    """
    for rel in ("README.md", "manuscript/sections/methods.tex",
                "manuscript/sections/results.tex"):
        f = ROOT / rel
        if not f.exists():
            continue
        text = f.read_text()
        assert "115" in text, f"{rel} discusses the cost but no longer states $115"
        if rel in CORRECTION_RECORDS:
            continue
        for bad in OBSOLETE:
            assert bad not in text, (
                f"{rel} still carries the obsolete cost claim {bad!r}. One calculation, "
                "everywhere, or the paper argues against itself")


def test_the_cross_fitting_cost_is_twice_one_sweep():
    """$115 must be twice the $57 the same section derives, or the two figures are unrelated."""
    f = ROOT / "manuscript" / "sections" / "methods.tex"
    if not f.exists():
        pytest.skip("methods.tex not in this checkout")
    text = f.read_text()
    sweep = re.search(r"arms is about \\\$(\d+), and cross-fitting them is about twice that", text)
    assert sweep, "the sweep/cross-fitting sentence no longer matches"
    xfit = re.search(r"would need about \\\$(\d+):", text)
    assert xfit, "the cross-fitting derivation no longer matches"
    assert abs(int(xfit.group(1)) - 2 * int(sweep.group(1))) <= 2, (
        f"${xfit.group(1)} is not twice ${sweep.group(1)}; the two cost statements in one "
        "section do not agree")


# --- the external benchmark's four cells, on every surface that states them -------------------
#
# Every one of these is two decimals, so audit_manuscript.py never looked at them, and the
# withdrawn values are still legitimately present in docs/EXTERNAL_CORRECTION_1.md, which
# preserves them on purpose. So when the repaired analysis moved the chromosome-blocked
# cross-fitted cell from 1.68 (1.41 to 2.02) to 1.67 (1.40 to 2.00), the Discussion went on
# quoting the withdrawn figure as the current result for two commits and nothing objected: no
# orphan, no stale count, 1156/1156, CI green on all four jobs. An external audit found it by
# reading the paper. Six surfaces state these numbers and they have to agree.

EXTERNAL = {
    "supplied_2s":
        "directional R, supplied folds, two-stage (the published analysis, recomputed)",
    "supplied_cf":
        "directional R, supplied folds, CROSS-FITTED (D3: like-for-like with our primary)",
    "chrom_2s":
        "directional R, chromosome-blocked folds, two-stage (D2: the literal criterion)",
    "chrom_cf":
        "directional R, chromosome-blocked folds, cross-fitted (D2 and D3 together)",
}

# (file, cell, pattern, the column each capture group must equal). Patterns cross line breaks
# with \s+ because the TeX is hard-wrapped and the sentence is what matters, not the wrapping.
EXT_CLAIMS = [
    # Results, Table 14: the canonical statement, all four cells with their intervals.
    ("manuscript/sections/results.tex", "supplied_2s",
     r"supplied & two-stage & ([\d.]+) \(([\d.]+)--([\d.]+)\)", ("value", "ci_low", "ci_high")),
    ("manuscript/sections/results.tex", "supplied_cf",
     r"supplied & cross-fitted & ([\d.]+) \(([\d.]+)--([\d.]+)\)",
     ("value", "ci_low", "ci_high")),
    ("manuscript/sections/results.tex", "chrom_2s",
     r"chromosome-blocked & two-stage & ([\d.]+) \(([\d.]+)--([\d.]+)\)",
     ("value", "ci_low", "ci_high")),
    ("manuscript/sections/results.tex", "chrom_cf",
     r"chromosome-blocked & cross-fitted & ([\d.]+) \(([\d.]+)--([\d.]+)\)",
     ("value", "ci_low", "ci_high")),
    # Results prose, which restates two of them without their intervals.
    ("manuscript/sections/results.tex", "supplied_cf",
     r"It gives \$([\d.]+)\$ against the two-stage", ("value",)),
    ("manuscript/sections/results.tex", "supplied_2s",
     r"against the two-stage \$([\d.]+)\$", ("value",)),
    ("manuscript/sections/results.tex", "chrom_2s",
     r"The result is \$([\d.]+)\$ two-stage", ("value",)),
    ("manuscript/sections/results.tex", "chrom_cf",
     r"two-stage and \$([\d.]+)\$ cross-fitted, slightly larger", ("value",)),
    # The abstract states the two cross-fitted estimates.
    ("manuscript/paper.tex", "supplied_cf",
     r"produced a ([\d.]+)-fold\s+cross-fitted span \(95\\% CI ([\d.]+)--([\d.]+)\)",
     ("value", "ci_low", "ci_high")),
    ("manuscript/paper.tex", "chrom_cf",
     r"or ([\d.]+)-fold \(([\d.]+)--([\d.]+)\) after chromosome-blocked",
     ("value", "ci_low", "ci_high")),
    # The Discussion states all four estimator/partition combinations.
    ("manuscript/sections/discussion.tex", "supplied_2s",
     r"produced a ([\d.]+)-fold two-stage span \(95\\% CI ([\d.]+)--([\d.]+)\)",
     ("value", "ci_low", "ci_high")),
    ("manuscript/sections/discussion.tex", "supplied_cf",
     r"and a ([\d.]+)-fold\s+cross-fitted span \(([\d.]+)--([\d.]+)\)",
     ("value", "ci_low", "ci_high")),
    ("manuscript/sections/discussion.tex", "chrom_2s",
     r"estimates of\s+([\d.]+) \(([\d.]+)--([\d.]+)\) and",
     ("value", "ci_low", "ci_high")),
    ("manuscript/sections/discussion.tex", "chrom_cf",
     r"and ([\d.]+) \(([\d.]+)--([\d.]+)\)\. These datasets",
     ("value", "ci_low", "ci_high")),
    # README, which a reader meets before the paper.
    ("README.md", "supplied_2s",
     r"constructions is \*\*([\d.]+)\*\* \(95% CI ([\d.]+) to ([\d.]+)\) two-stage",
     ("value", "ci_low", "ci_high")),
    ("README.md", "supplied_cf",
     r"and \*\*([\d.]+)\*\* \(([\d.]+) to ([\d.]+)\) cross-fitted",
     ("value", "ci_low", "ci_high")),
    ("README.md", "chrom_2s", r"does not do, gives ([\d.]+) and [\d.]+ with", ("value",)),
    ("README.md", "chrom_cf", r"does not do, gives [\d.]+ and ([\d.]+) with", ("value",)),
]


@pytest.mark.parametrize("path,cell,pattern,columns", EXT_CLAIMS,
                         ids=[f"{c[0].split('/')[-1]}:{c[1]}:{'+'.join(c[3])}"
                              for c in EXT_CLAIMS])
def test_an_external_benchmark_cell_agrees_with_its_evidence_row(path, cell, pattern, columns):
    f = ROOT / path
    if not f.exists():
        pytest.skip(f"{path} not in this checkout")
    m = re.search(pattern, f.read_text())
    assert m, (
        f"{path}: the sentence stating the {cell} external cell no longer matches this "
        "pattern. Reworded, in which case fix the pattern; or deleted, in which case delete "
        "this claim. It must not quietly stop covering the number")
    assert len(m.groups()) == len(columns), "the pattern and the column list disagree"
    for got, col in zip(m.groups(), columns):
        truth = round(_row("external_sensitivity.csv", EXTERNAL[cell], col), 2)
        assert float(got) == truth, (
            f"{path} states {got} for the {cell} {col}; external_sensitivity.csv gives "
            f"{truth}. The withdrawn values live in docs/EXTERNAL_CORRECTION_1.md and "
            "nowhere else: prose describing the current analysis must state the current one")
