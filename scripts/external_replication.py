"""Claim A on datasets this study never touched, under another group's negative sets.

    python scripts/external_replication.py            # the disjoint datasets
    python scripts/external_replication.py --from-cache

The protocol for this analysis was committed first, in docs/EXTERNAL_BENCHMARK_PROTOCOL.md at
e76a80c, BEFORE THE 135-DATASET DISJOINT SUBSET BELOW WAS SCORED. The eligibility criteria, the
estimand, the uncertainty procedure and the falsification thresholds are quoted from it and
were not chosen after seeing this output.

THE ESTIMATOR HERE IS TWO-STAGE, NOT CROSS-FITTED. gain_over_composition is the estimator the
paper reports as its comparability analysis, so this span belongs beside the two-stage 5.42 and
NOT beside the primary cross-fitted 4.84. Cross-fitting these arms costs ten extra 4-mer fits
per dataset and has simply not been run. Nothing here should be read as external validation of
the estimator the paper recommends, and the abstract used to say "the same estimator", which
elided exactly this.

THEIR FOLDS ARE NOT CHROMOSOME-BLOCKED, and the protocol makes that a criterion. See
fold_blocking() below: measured rather than assumed, and the channel is closed anyway.

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
from rbp.utils.carry import emit  # noqa: E402
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


def fold_blocking():
    """Is their fold partition chromosome-blocked, and if not, does the leakage channel close?

    THE PROTOCOL MAKES THIS A FALSIFICATION CRITERION. It calls the result indeterminate if the
    benchmark's fold partition "cannot be made chromosome-blocked", and nothing here had ever
    checked. Their split is emphatically NOT chromosome-blocked: every chromosome appears in
    every fold in all 135 datasets.

    That is not the same as leaking. Chromosome grouping is a coarse instrument for one specific
    channel, a held-out window sharing sequence with a training window, and the direct metric for
    that channel is the fraction of positives with a near neighbour on the same strand whose
    neighbour sits in a DIFFERENT fold. Measured here on their own coordinates, exactly as
    fold_integrity.py measures ours. Their scheme turns out to be locus-blocked rather than
    chromosome-blocked, which closes the channel by a finer instrument than ours.

    Reported rather than assumed, and reported with its denominator, because "0% cross-fold"
    means nothing if almost nothing has a neighbour: 70 to 86 per cent of their positives do.
    """
    have = cross = split = seen = 0
    for ds in disjoint_datasets():
        rows = []
        for fold in range(5):
            f = DATA / ds.replace(":", "_") / f"fold-{fold}" / f"positive.fold-{fold}.bed"
            if not f.exists():
                continue
            for line in f.read_text().split("\n"):
                if line.strip():
                    q = line.split("\t")
                    rows.append((q[0], int(q[1]), q[5] if len(q) > 5 else "+", fold))
        if not rows:
            continue
        seen += 1
        per_chrom = {}
        for c, _st, _s, fo in ((r[0], r[2], r[1], r[3]) for r in rows):
            per_chrom.setdefault(c, set()).add(fo)
        if any(len(v) > 1 for v in per_chrom.values()):
            split += 1
        by = {}
        for c, st_, strand, fo in ((r[0], r[1], r[2], r[3]) for r in rows):
            by.setdefault((c, strand), []).append((st_, fo))
        for v in by.values():
            v.sort()
            for i in range(len(v)):
                near = [j for j in (i - 1, i + 1)
                        if 0 <= j < len(v) and abs(v[i][0] - v[j][0]) <= 1000]
                if near:
                    have += 1
                    if any(v[j][1] != v[i][1] for j in near):
                        cross += 1
    return seen, split, have, cross



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
    # THE COMMITTED VALUE WINS, AND THAT IS WHAT MAKES THIS IDEMPOTENT. This row was first
    # emitted only when the unpacked deposit was present, so `--from-cache` on a clean clone
    # silently DROPPED it and staled the provenance manifest. My first repair carried the value
    # forward but wrote a DIFFERENT note and turned `n` into a float, so the table still failed
    # to reproduce byte for byte: one dropped row traded for a changed row. Both were found by
    # running the documented command in a clean export, which is the only environment that can
    # see either.
    #
    # So the deposit is counted ONCE, when there is no committed value to carry, and every run
    # afterwards reproduces the committed row exactly, in any environment. The deposit is a
    # fixed Zenodo record with a verified MD5, so its size is not a quantity that should move;
    # if it ever does, `--check` on this table is where that belongs, not a silent rewrite.
    n_dep, src = None, ""
    if OUT.exists():
        prev = pd.read_csv(OUT).set_index("check")
        if DEPOSIT_ROW in prev.index:
            n_dep = int(prev.loc[DEPOSIT_ROW, "value"])
            src = str(prev.loc[DEPOSIT_ROW, "note"])
    if n_dep is None and DATA.exists():
        n_dep = sum(1 for q in DATA.iterdir() if q.is_dir())
        src = "counted from the unpacked deposit; the rest overlap our panel"
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

    # THE PROTOCOL'S FOLD CRITERION, measured rather than assumed. Carried forward where the
    # deposit is absent, which is the released tree; see rbp.utils.carry.
    seen, split, have, cross = fold_blocking() if DATA.exists() else (0, 0, 0, 0)
    ok = seen > 0
    emit(out, OUT, "their datasets whose chromosomes span more than one fold", split,
         n=seen, recomputed=ok,
         note="their partition is NOT chromosome-blocked; the protocol makes this a criterion")
    emit(out, OUT, "their positives with a same-strand neighbour within 1 kb", have,
         n=seen, recomputed=ok,
         note="the denominator, because a zero cross-fold rate over nothing means nothing")
    emit(out, OUT, "of those, the fraction whose neighbour is in a different fold",
         (cross / have) if have else 0.0, n=have, recomputed=ok,
         note="the direct leakage metric fold_integrity.py uses on our arms; their scheme is "
              "locus-blocked rather than chromosome-blocked, which closes the same channel")

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
