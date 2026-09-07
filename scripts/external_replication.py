"""Claim A on datasets this study never touched, under another group's negative sets.

    python scripts/external_replication.py            # the disjoint datasets
    python scripts/external_replication.py --from-cache

The protocol for this analysis was committed first, in docs/EXTERNAL_BENCHMARK_PROTOCOL.md at
e76a80c, BEFORE THE 135-DATASET DISJOINT SUBSET BELOW WAS SCORED. The eligibility criteria, the
estimand, the uncertainty procedure and the falsification thresholds are quoted from it and
were not chosen after seeing this output.

THAT IS NARROWER THAN A PROSPECTIVE SEARCH, and an external audit was right to say the earlier
wording overstated it. This benchmark was not unknown: 3f96e8e and 2b5843a scored its
45-dataset intersection with our panel on 2026-08-31, six days before the protocol was written,
and 28a5a1f analysed it further. What was genuinely fixed in advance is the rule applied to the
previously unscored complement. Call it a pre-specified held-out-subset analysis of an
already-known external construction, which is what it is.

What the search found, including that the protocol's own premise was wrong. The protocol
asserted that the Horlacher et al. 2023 benchmark could not test Claim A "because testing Claim
A needs a dataset released under more than one negative-set construction for the same positives
and Horlacher releases one". It releases two, per fold, for the same positives: `negative-1`,
uniform positions from transcripts carrying a site of the target, and `negative-2`, other RBPs'
crosslink sites. The premise was wrong and the correction is recorded here rather than quietly
absorbed, because a prospective protocol whose errors are edited out afterwards is not one.

WHY A SECOND SCRIPT AND NOT AN OPTION ON horlacher_arm.py. That script deliberately restricts
to the intersection with our own panel, so that the measurement is on the same proteins; it is
the right design for the question it asks and its output is published. But it means all 45 of
its datasets are also in our 94, so it fails criterion 1 of the protocol and is NOT an
independent replication. It is an independent CONSTRUCTION: another group's windows, negatives,
peak calls and folds over largely the same experiments. This script takes the complement, the
Horlacher ENCODE datasets that our panel does not contain, which is the sample criterion 1
asks for.

What this still does not establish. The complement is still ENCODE eCLIP in K562 and HepG2, and
it is NOT an independent sample of proteins: 31 of its 108 proteins also appear among our 79,
in other cell lines or other experiments. It is disjoint in DATASETS and in processing, not in
biology. The no-shared-protein sensitivity below removes all 31 and the verdict survives, so
the replication is not carried by the overlap; but "independent sample of proteins", which an
earlier version of this docstring said, is false and the audit that caught it was right.
Claim A's scope after this remains eCLIP-derived panels. Nothing here bears on Claim B, the
directional relation, which does not replicate on this benchmark and stays labelled that way.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))

TABLES = ROOT / "results" / "tables"
DATA = ROOT.parent / "rbp-store" / "external" / "samples" / "processed" / "ENCODE"
GENOME = ROOT.parent / "rna-binding-proteins" / "data" / "raw" / "GRCh38.primary_assembly.genome.fa"
S1 = TABLES / "supplementary_table_s1.csv"
PER = TABLES / "external_replication_per_dataset.csv"
OUT = TABLES / "external_replication.csv"

# Quoted from docs/EXTERNAL_BENCHMARK_PROTOCOL.md. Do not edit these to match an outcome.
DEPOSIT_ROW = "datasets in the Horlacher ENCODE deposit"
MIN_DATASETS = 20
SUPPORT_SPAN = 1.5
SUPPORT_CI_LOW = 1.2


def our_panel():
    if not S1.exists():
        return set()
    d = pd.read_csv(S1)
    return set(d[d.in_three_arm_panel.astype(str) == "True"].dataset)


def disjoint_datasets():
    """Their ENCODE datasets that our study panel does not contain."""
    theirs = {p.name.replace("_", ":"): p.name for p in DATA.iterdir() if p.is_dir()}
    ours = our_panel()
    return {k: v for k, v in sorted(theirs.items()) if k not in ours}


def build(limit=0):
    from horlacher_arm import build as build_arm
    from pyfaidx import Fasta
    if not GENOME.exists():
        sys.exit(f"no genome at {GENOME}")
    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    sel = disjoint_datasets()
    keys = list(sel)
    if limit:
        keys = keys[:limit]
    log(f"  {len(disjoint_datasets())} of their datasets are outside our panel; "
        f"building {len(keys)}")
    rows = []
    for i, k in enumerate(keys, 1):
        protein, cell = k.split(":")
        rec, ok = {"dataset": k, "protein": protein, "cell": cell}, True
        for neg, tag in (("negative-1", "n1"), ("negative-2", "n2")):
            d = build_arm(fa, sel[k], neg)
            if d is None:
                ok = False
                break
            sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
            m = ~np.isnan(sc)
            g = gain_over_composition(d.seq_rna.values[m], sc[m], d.label.values[m],
                                      d.fold.values[m])
            rec[f"comp_{tag}"] = g.auroc_composition
            rec[f"gain_{tag}"] = g.delta
            rec[f"n_{tag}"] = g.n
        if ok:
            rows.append(rec)
            log(f"  [{i:3d}/{len(keys)}] {k:20s} n1 gain {rec['gain_n1']:+.4f}   "
                f"n2 gain {rec['gain_n2']:+.4f}")
    return pd.DataFrame(rows)


def span(d):
    """Descriptive magnitude: how far apart the two arms are, direction discarded."""
    a, b = d.gain_n1.mean(), d.gain_n2.mean()
    lo, hi = min(a, b), max(a, b)
    return float(hi / lo) if lo > 0 else float("nan")


def ratio(d):
    """PRE-LABELLED ratio, negative-1 over negative-2, which is the one the rule needs.

    max/min is bounded below by 1 by construction, so the protocol's "fails to replicate if the
    interval contains 1.0" clause is close to unsatisfiable on it: an external audit was right
    to call the rule ill-posed. Fixing the numerator and denominator in advance lets the
    interval fall below 1 if the arms swap, which is what "no replication" would look like.
    Both are reported; the verdict is computed on this one. In this sample the ordering never
    swaps in 4000 draws, so the two agree and nothing about the published number moves.
    """
    a, b = d.gain_n1.mean(), d.gain_n2.mean()
    return float(a / b) if b > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-cache", action="store_true")
    ap.add_argument("--n", type=int, default=0, help="smoke test; writes .partial.csv")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")

    if a.from_cache:
        t = pd.read_csv(PER)
    else:
        t = build(a.n)
        if t.empty:
            sys.exit("no dataset built")
        # --n WRITES A PARTIAL FILE. A smoke run once truncated a released 94-row table to 5.
        t.to_csv(PER.with_suffix(".partial.csv") if a.n else PER, index=False)
        if a.n:
            log(f"  --n {a.n}: wrote {PER.with_suffix('.partial.csv').name}, not the "
                "committed table")
            return

    # PROTEIN-CLUSTERED, as the protocol specifies, and it matters here: their release covers
    # some proteins in both cell lines, exactly as ours does.
    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(4000)]

    out = []

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})
        return float(v.mean())

    # The deposit's own size, so the manuscript's count of it traces to a table rather than to
    # nothing. Read from the unpacked deposit when it is present; omitted rather than guessed
    # when it is not, because a hardcoded 223 here is the hand-maintained count this repository
    # keeps getting wrong.
    # IDEMPOTENT IN CACHE MODE. This row used to be emitted only when the unpacked deposit was
    # present, so `--from-cache` on a clean clone silently DROPPED it, changed the committed
    # table and staled the provenance manifest. An audit found exactly that. Counted from the
    # deposit when it is there; otherwise carried forward from the committed table, which is
    # released evidence, rather than omitted or hardcoded.
    n_dep = None
    if DATA.exists():
        n_dep = sum(1 for q in DATA.iterdir() if q.is_dir())
        src = "counted from the unpacked deposit; the rest overlap our panel"
    elif OUT.exists():
        prev = pd.read_csv(OUT).set_index("check")
        if DEPOSIT_ROW in prev.index:
            n_dep = float(prev.loc[DEPOSIT_ROW, "value"])
            src = "carried forward from the committed table; the deposit is not unpacked here"
    if n_dep is not None:
        out.append({"check": DEPOSIT_ROW, "value": n_dep, "ci_low": "", "ci_high": "",
                    "n": n_dep, "note": src})
    out.append({"check": "datasets", "value": len(t), "ci_low": "", "ci_high": "", "n": len(t),
                "note": "Horlacher ENCODE datasets outside our 94-dataset study panel"})
    out.append({"check": "proteins", "value": len(uniq), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "the resampled unit"})
    out.append({"check": "dataset overlap with our study panel", "value": 0, "ci_low": "",
                "ci_high": "", "n": len(t),
                "note": "zero by construction against the 94-dataset three-arm panel; this is "
                        "what criterion 1 asks for and what the 45-dataset intersection lacks"})

    # PROTEIN OVERLAP IS NOT ZERO, and an earlier version of this file called the complement an
    # independent sample of proteins. Read from the released panel table rather than a literal.
    ours = pd.read_csv(S1)
    ours = set(ours[ours.in_three_arm_panel].protein) if "in_three_arm_panel" in ours else set()
    shared = sorted(ours & set(uniq))
    out.append({"check": "protein overlap with our study panel", "value": len(shared),
                "ci_low": "", "ci_high": "", "n": len(uniq),
                "note": "the two panels share biology even where they share no dataset, so "
                        "this is disjoint in datasets and processing, not in proteins"})

    g1 = add("nested contribution, negative-1 (bias-agnostic)", t.gain_n1)
    g2 = add("nested contribution, negative-2 (bias-aware)", t.gain_n2)
    add("composition baseline, negative-1", t.comp_n1)
    add("composition baseline, negative-2", t.comp_n2)
    add("CONTRAST, negative-1 minus negative-2", t.gain_n1 - t.gain_n2)

    pt = span(t)
    b = np.array([span(t.iloc[i]) for i in draws])
    b = b[np.isfinite(b)]
    lo, hi = (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))) if len(b) else ("", "")
    out.append({"check": "SPAN across their two negative-set constructions", "value": pt,
                "ci_low": lo, "ci_high": hi, "n": len(t),
                "note": "the protocol's primary estimand, transported unchanged: ratio of the "
                        "two arms' panel means, protein-clustered percentile interval"})

    # The pre-fixed decision, applied mechanically. The thresholds are quoted from the protocol
    # and are not adjustable here; writing the comparison as code rather than as prose is the
    # point, because prose can be softened after the fact and an assertion cannot.
    rpt = ratio(t)
    rb = np.array([ratio(t.iloc[i]) for i in draws])
    rb = rb[np.isfinite(rb)]
    rlo, rhi = ((float(np.percentile(rb, 2.5)), float(np.percentile(rb, 97.5)))
                if len(rb) else ("", ""))
    out.append({"check": "PRE-LABELLED RATIO, negative-1 over negative-2", "value": rpt,
                "ci_low": rlo, "ci_high": rhi, "n": len(t),
                "note": "the quantity the falsification rule is evaluated on, because max/min "
                        "is bounded below by 1 and so cannot straddle it"})
    out.append({"check": "the arm ordering never swaps under resampling", "value":
                float((rb > 1.0).mean()), "ci_low": "", "ci_high": "", "n": len(rb),
                "note": "fraction of draws with negative-1 above negative-2; at 1.0 the "
                        "pre-labelled ratio and the max/min span agree exactly"})

    supported = (rpt > SUPPORT_SPAN and isinstance(rlo, float) and rlo > SUPPORT_CI_LOW)
    failed = isinstance(rlo, float) and rlo < 1.0 < rhi
    powered = len(t) >= MIN_DATASETS
    verdict = ("supported" if (supported and powered)
               else "fails to replicate" if (failed and powered)
               else "indeterminate")
    out.append({"check": "protocol verdict for Claim A on an independent sample", "value": 1.0
                if verdict == "supported" else 0.0, "ci_low": "", "ci_high": "", "n": len(t),
                "note": f"{verdict}; criteria fixed in docs/EXTERNAL_BENCHMARK_PROTOCOL.md at "
                        f"commit e76a80c: supported if span > {SUPPORT_SPAN} and CI low > "
                        f"{SUPPORT_CI_LOW}, fails if the interval contains 1.0, indeterminate "
                        f"otherwise or if fewer than {MIN_DATASETS} datasets qualify"})

    # THE NO-SHARED-PROTEIN SENSITIVITY. Gated rather than quoted from a report, because a
    # number copied out of an audit by hand is the class of unsourced value this repository
    # spent a fortnight removing.
    k = t[~t.protein.isin(shared)]
    if len(k) >= MIN_DATASETS:
        kp = k.protein.to_numpy()
        ku = np.unique(kp)
        km = [np.flatnonzero(kp == q) for q in ku]
        krng = np.random.default_rng(0)
        kd = [np.concatenate([km[j] for j in krng.integers(0, len(ku), len(ku))])
              for _ in range(4000)]
        kpt = span(k)
        kb = np.array([span(k.iloc[i]) for i in kd])
        kb = kb[np.isfinite(kb)]
        klo, khi = float(np.percentile(kb, 2.5)), float(np.percentile(kb, 97.5))
        for lab, v in (("negative-1", k.gain_n1), ("negative-2", k.gain_n2)):
            out.append({"check": f"no-shared-protein contribution, {lab}",
                        "value": float(v.mean()), "ci_low": "", "ci_high": "", "n": len(k),
                        "note": f"the {len(shared)} proteins shared with our panel removed"})
        out.append({"check": "SPAN excluding every protein shared with our panel",
                    "value": kpt, "ci_low": klo, "ci_high": khi, "n": len(k),
                    "note": f"{len(k)} datasets over {len(ku)} proteins; the replication must "
                            "not depend on the shared biology"})
        out.append({"check": "no-shared-protein span still meets the pre-fixed criteria",
                    "value": 1.0 if (kpt > SUPPORT_SPAN and klo > SUPPORT_CI_LOW) else 0.0,
                    "ci_low": "", "ci_high": "", "n": len(k),
                    "note": f"same thresholds as the primary verdict: span > {SUPPORT_SPAN} "
                            f"and CI low > {SUPPORT_CI_LOW}"})

    r = pd.DataFrame(out)
    r.to_csv(OUT, index=False)

    log(f"\n=== Claim A on {len(t)} datasets our panel does not contain, {len(uniq)} proteins "
        "===\n")
    log(f"  negative-1 (bias-agnostic)  contribution {g1:+.4f}")
    log(f"  negative-2 (bias-aware)     contribution {g2:+.4f}")
    log(f"  SPAN                        {pt:.3f}x" +
        (f"  95% CI [{lo:.3f}, {hi:.3f}]" if isinstance(lo, float) else ""))
    log(f"\n  PRE-FIXED VERDICT: Claim A {verdict.upper()}")
    log(f"  wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
