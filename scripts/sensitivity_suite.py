"""Does the headline survive the choices nobody argued for? One table, one row per choice.

    python scripts/sensitivity_suite.py

WHAT THIS IS FOR. The paper's headline is a ratio of panel MEANS: the 4-mer's nested
contribution under dinucleotide matching over its contribution under the bias-aware protocol,
5.42-fold. Every part of that sentence is a choice. The mean rather than the median. All 94
datasets rather than a subset. Both cell lines pooled. Datasets weighted equally rather than by
size. None of those choices was preregistered and none was argued for in the text, which means a
referee is entitled to ask what happens under the others. This runs them.

WHAT IT DELIBERATELY DOES NOT DO. It does not pick a winner. Every row is the same estimand
computed under a different aggregation or a different subset, and the point is the spread, not
the best member. Anyone reading this table to find the largest span has misread it.

PLANNED OR POST-HOC, STATED PER ROW. Four of these analyses existed before this script and are
referenced rather than rerun: the unbounded estimand (estimands.csv, delta_deviance), the
baseline-order sensitivity (baseline_order.csv), the dataset-size relation (robustness.csv), and
the cross-cell-line replication on the 15 proteins measured in both (r1_robustness.csv). Those
are marked `planned`. Everything computed here is marked `post-hoc, added after the results were
known`, because it was, and a sensitivity analysis chosen after seeing the result is weaker
evidence than one chosen before. Saying so is the difference between a robustness check and a
specification search.

WHAT IS NOT HERE AND WHY. Class-ratio sensitivity and the k-mer capacity ladder need the
composition features refit on resampled rows, which needs the window sequences. Only the
bias-aware arm's windows are on this machine; the composition-matched arms' are in cloud storage.
Both are in the costed plan in docs/REDRAW_PLAN.md rather than half-done here, because a
sensitivity run on one arm of three would answer a different question from the one asked.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
PER = TABLES / "three_arm_models_per_dataset.csv"
MATCH = TABLES / "match_quality_per_dataset.csv"
OUT = TABLES / "sensitivity_suite.csv"

MODELS = ("kmer", "cnn", "splicebert")
ARMS = ("dn", "gc", "neg2")
PLANNED = "planned before this script; the analysis it refers to is cited in the note"
POSTHOC = "post-hoc, added after the results were known"


def span(d, model):
    """The published estimand: ratio of panel means, largest arm over smallest."""
    m = {a: d[f"{model}_gain_{a}"].mean() for a in ARMS}
    lo = min(m.values())
    return (max(m.values()) / lo) if lo > 0 else np.nan


def span_with(d, model, agg):
    m = {a: agg(d[f"{model}_gain_{a}"].to_numpy()) for a in ARMS}
    lo = min(m.values())
    return (max(m.values()) / lo) if lo > 0 else np.nan


def trimmed(x, frac=0.10):
    x = np.sort(np.asarray(x, dtype=float))
    k = int(len(x) * frac)
    return x[k:len(x) - k].mean() if len(x) - 2 * k > 0 else x.mean()


def ordering_holds(d, model, agg=np.mean):
    m = [agg(d[f"{model}_gain_{a}"].to_numpy()) for a in ARMS]
    return float(m[0] > m[1] > m[2])          # dn > gc > neg2


def clustered_ci(d, model, agg, draws=4000, seed=0):
    """Protein-clustered percentile interval for a span under an arbitrary aggregator."""
    rng = np.random.default_rng(seed)
    prot = d.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    out = []
    for _ in range(draws):
        idx = np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
        out.append(span_with(d.iloc[idx], model, agg))
    out = np.asarray(out, dtype=float)
    out = out[np.isfinite(out)]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=4000)
    a = ap.parse_args()

    d = pd.read_csv(PER)
    rows = []

    def add(check, value, note, plan=POSTHOC, n=None, ci=(None, None)):
        rows.append({"check": check, "value": value,
                     "ci_low": "" if ci[0] is None else ci[0],
                     "ci_high": "" if ci[1] is None else ci[1],
                     "n": len(d) if n is None else n, "planning": plan, "note": note})

    # ---- the reference point, so every row below is comparable to something ----------------
    for m in MODELS:
        add(f"span, {m}, published (mean, all datasets, both cell lines)", span(d, m),
            "the value in the paper", PLANNED,
            ci=clustered_ci(d, m, np.mean, a.draws))

    # ---- aggregation: the mean is a choice ------------------------------------------------
    for label, agg in (("median", np.median), ("10% trimmed mean", trimmed)):
        for m in MODELS:
            add(f"span, {m}, {label} instead of mean", span_with(d, m, agg),
                "the ratio is of panel aggregates; the mean is not privileged, and a "
                "right-skewed per-dataset distribution would inflate it",
                ci=clustered_ci(d, m, agg, a.draws))
            add(f"ordering dn > gc > neg2 holds, {m}, {label}", ordering_holds(d, m, agg),
                "1 if the three arms keep their order under this aggregator")

    # ---- cell line: two cell lines, pooled, and nobody checked them apart ------------------
    for cell, sub in d.groupby("cell"):
        for m in MODELS:
            add(f"span, {m}, {cell} only", span(sub, m),
                "stratified rather than pooled. The panel is 2 cell lines and 1 assay, so a "
                "span that existed in only one of them would not support the general claim",
                n=len(sub))
            add(f"ordering dn > gc > neg2 holds, {m}, {cell} only", ordering_holds(sub, m),
                "1 if the three arms keep their order within this cell line", n=len(sub))

    # ---- dataset size: do a few large datasets carry the panel mean? -----------------------
    # `pairs` is per (dataset, arm); use the GC arm's, which is the panel of record.
    if MATCH.exists():
        mq = pd.read_csv(MATCH)
        gcq = mq[mq.arm == "gc"].set_index("dataset")
        sz = d.dataset.map(gcq["pairs"]).astype(float)
        big = sz >= sz.quantile(0.75)
        for m in MODELS:
            add(f"span, {m}, largest quartile by pairs dropped", span(d[~big], m),
                "if the panel mean is carried by a few large datasets, removing them moves it",
                n=int((~big).sum()))
            add(f"span, {m}, largest quartile only", span(d[big], m),
                "the complement of the row above", n=int(big.sum()))
        top = sz.sort_values(ascending=False)
        add("share of total pairs in the largest 10 datasets",
            float(top.head(10).sum() / top.sum()),
            "a panel mean over datasets weights each dataset equally regardless of this, "
            "which is a defensible choice but should be visible")

        # ---- match quality: drop the arms' worst-matched datasets ---------------------------
        # THE MATCHING IS THE PROTOCOL. A dataset whose GC matching left a large residual gap is
        # a dataset where "GC-matched" describes the intent rather than the achieved design, and
        # its contribution partly reflects the failure to match rather than the protocol.
        worst = mq[mq.arm == "gc"].set_index("dataset")["gc_gap_p90"]
        gap = d.dataset.map(worst).astype(float)
        for q, lab in ((0.90, "worst decile"), (0.75, "worst quartile")):
            keep = gap <= gap.quantile(q)
            for m in MODELS:
                add(f"span, {m}, {lab} by achieved GC gap dropped", span(d[keep], m),
                    "achieved match quality, not nominal tolerance; the p90 absolute GC gap "
                    "in the GC arm", n=int(keep.sum()))

    # ---- leave one protein out ------------------------------------------------------------
    # 79 proteins, 15 of which contribute two datasets. What is reported is the extreme, not a
    # distribution, because the question is whether ANY single protein carries the result.
    prot = d.protein.to_numpy()
    for m in MODELS:
        base = span(d, m)
        vals = []
        for q in np.unique(prot):
            vals.append((span(d[prot != q], m), q))
        lo, hi = min(vals), max(vals)
        add(f"span, {m}, leave-one-protein-out minimum", lo[0],
            f"dropping {lo[1]} and all its datasets; base {base:.4f}",
            n=len(np.unique(prot)))
        add(f"span, {m}, leave-one-protein-out maximum", hi[0],
            f"dropping {hi[1]} and all its datasets; base {base:.4f}",
            n=len(np.unique(prot)))
        add(f"span, {m}, largest leave-one-protein-out displacement",
            max(abs(lo[0] - base), abs(hi[0] - base)),
            "absolute, in fold units; the question is whether any single protein carries the "
            "span, not how it is distributed", n=len(np.unique(prot)))

    # ---- pointers to the sensitivities that already existed --------------------------------
    add("unbounded estimand: ordering holds on delta_deviance, all three models", 1.0,
        "estimands.csv rows 'protocol ordering holds, delta_deviance, {kmer,cnn,splicebert}'. "
        "AUROC is bounded on [0,1] and compresses near its ceiling; deviance is not bounded "
        "above, and the ordering survives it", PLANNED, n=3)
    add("baseline order: the span does not collapse at order three", 1.0,
        "baseline_order.csv. Raising the composition baseline from order two to order three "
        "shrinks every contribution but does not close the span; the fold range across "
        "protocols is reported there per model", PLANNED, n=3)
    add("dataset size: the contrast is not explained by size", 1.0,
        "robustness.csv row 'R1 effect vs log10(dataset size)'", PLANNED)
    add("cell line: the contrast replicates across the 15 proteins measured in both", 0.9094,
        "r1_robustness.csv 'REPLICATION of the contrast across cell lines'; pearson, "
        "p=2.63e-06", PLANNED, n=15)

    r = pd.DataFrame(rows)
    r.to_csv(OUT, index=False)

    log(f"\n  {len(r)} sensitivity rows -> {OUT.relative_to(ROOT)}")
    log("\n  spans by aggregator (published / median / 10% trimmed):")
    for m in MODELS:
        log(f"    {m:11s} {span(d, m):6.3f}  {span_with(d, m, np.median):6.3f}  "
            f"{span_with(d, m, trimmed):6.3f}")
    log("\n  spans by cell line:")
    for cell, sub in d.groupby("cell"):
        log(f"    {cell:6s} n={len(sub):3d}  " +
            "  ".join(f"{m} {span(sub, m):6.3f}" for m in MODELS))
    log("\n  largest leave-one-protein-out displacement, in fold units:")
    for m in MODELS:
        v = [x["value"] for x in rows
             if x["check"] == f"span, {m}, largest leave-one-protein-out displacement"][0]
        log(f"    {m:11s} {v:6.3f}")
    ordering = [x for x in rows if x["check"].startswith("ordering dn > gc > neg2 holds")]
    log(f"\n  the ordering dn > gc > neg2 holds in {int(sum(x['value'] for x in ordering))} of "
        f"{len(ordering)} aggregator-and-stratum combinations")


if __name__ == "__main__":
    main()
