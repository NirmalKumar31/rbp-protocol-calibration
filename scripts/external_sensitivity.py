"""The two external sensitivities Amendment 1 specifies: cross-fitting, and chromosome folds.

    python scripts/external_sensitivity.py                 # both, 135 datasets
    python scripts/external_sensitivity.py --n 4            # smoke test, writes .partial.csv
    python scripts/external_sensitivity.py --from-cache     # summary only

What this answers. `docs/EXTERNAL_BENCHMARK_AMENDMENT.md`, committed at 02a2bac BEFORE this
script was run, records three defects an external audit found in the external-benchmark claim
and fixes the first by specification:

  D1  the prewritten estimand max(mean arm)/min(mean arm) is bounded below by 1, so the stated
      failure rule "the interval contains 1.0" was close to unsatisfiable. The amendment fixes a
      DIRECTIONAL estimand, R = M(negative-1) / M(negative-2), whose numerator and denominator
      are set by the deposit's own labels before any value is computed, so R < 1 is attainable
      and means the ordering is the reverse of ours.
  D2  eligibility criterion 4 wanted a chromosome-blocked partition. Horlacher's supplied folds
      are not chromosome-blocked, and were accepted on a post-hoc near-neighbour diagnostic.
      This script builds the chromosome-blocked partition the criterion asked for and reports
      both.
  D3  the external analysis used the two-stage estimator while the internal PRIMARY estimator
      is cross-fitted, so the external figure was comparable with the internal 5.42 and not
      with the internal 4.84. This script transports the cross-fitted estimator unchanged.

So there are four cells: {two-stage, cross-fitted} x {supplied folds, chromosome folds}. The
first is the published external analysis and is recomputed here as a control: if it does not
reproduce the committed 1.69, nothing else in this table can be trusted.

THE DECISION RULE IS THE AMENDMENT'S AND IS NOT RENEGOTIABLE HERE. supports = point >= 1.5 and
lower bound >= 1.2; fails = the interval contains 1.0 or the point is below 1.0; weak = neither;
indeterminate = the denominator arm's panel mean is not positive, or under 20 datasets survive.
Every cell's verdict is emitted whatever it says.

Why it costs nothing. A 4-mer logistic fit on one of these datasets is milliseconds.
Cross-fitting adds ten complement fits per dataset per arm. The whole run is one laptop core for
well under an hour, with no GPU and no cloud.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rbp.eval.baseline import kmer_matrix  # noqa: E402
from rbp.eval.nested import composition_features  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "external_sensitivity.csv"
PER = TABLES / "external_sensitivity_per_dataset.csv"

N_FOLDS = 5
BOOT = 4000
SEED = 20260907              # fixed in the amendment, before any of this was run
K = 4
MIN_DATASETS = 20            # the amendment's indeterminacy floor
SUPPORT_POINT, SUPPORT_LOW = 1.5, 1.2


# --- the chromosome-blocked partition, exactly as the amendment specifies -------------------

def chrom_folds(chrom, n_folds=N_FOLDS):
    """Assign folds by chromosome, largest chromosome first into the emptiest fold.

    Deterministic by construction and by tie-break: chromosomes are ordered by descending
    window count and then by NAME ASCENDING, so the result is a function of the data and not of
    dictionary insertion order. No chromosome is ever split, which is the property criterion 4
    asked for and the supplied folds do not have.

    Returns None if there are fewer than `n_folds` distinct chromosomes, which cannot fill the
    folds; the caller counts that as an exclusion rather than dropping it silently.
    """
    chrom = np.asarray(chrom)
    names, counts = np.unique(chrom, return_counts=True)
    if len(names) < n_folds:
        return None
    order = sorted(range(len(names)), key=lambda i: (-counts[i], str(names[i])))
    load = [0] * n_folds
    where = {}
    for i in order:
        f = min(range(n_folds), key=lambda j: (load[j], j))
        where[names[i]] = f
        load[f] += int(counts[i])
    return np.array([where[c] for c in chrom])


def one_fold_per_chrom(chrom, folds):
    """True when every chromosome sits in exactly one fold. The criterion, measured."""
    seen = {}
    for c, f in zip(np.asarray(chrom), np.asarray(folds)):
        seen.setdefault(c, set()).add(int(f))
    return all(len(v) == 1 for v in seen.values())


def leakage(chrom, start, strand, folds):
    """(with a neighbour, of those crossing a fold) at 1 kb on the same strand.

    The same measure `external_replication.fold_blocking` applies to the supplied folds, so the
    two partitions are comparable on the channel chromosome blocking is a proxy for. Reported
    with its denominator, because "0% cross-fold" is vacuous if nothing has a neighbour.
    """
    by = {}
    for c, s, st, f in zip(chrom, start, strand, folds):
        by.setdefault((c, st), []).append((int(s), int(f)))
    have = cross = 0
    for v in by.values():
        v.sort()
        for i in range(len(v)):
            near = [j for j in (i - 1, i + 1)
                    if 0 <= j < len(v) and abs(v[i][0] - v[j][0]) <= 1000]
            if near:
                have += 1
                if any(v[j][1] != v[i][1] for j in near):
                    cross += 1
    return have, cross


# --- the two estimators, on whatever fold column it is given --------------------------------

def contributions(seqs, y, folds):
    """(two-stage, cross-fitted) nested contribution of a 4-mer over the composition block.

    Both come from cross_fitting.py unmodified. The two-stage value feeds every base score from
    the model that saw every fold but the row's own; the cross-fitted value replaces the
    covariate on each outer-TRAINING row with a score from the model that saw neither the outer
    test fold nor the row's own.
    """
    import cross_fitting as cf
    comp, _ = composition_features(seqs, True, standardise_cols=True)
    a_comp = cf._pooled_comp(comp, y, folds)
    X, _ = kmer_matrix(list(seqs), K)
    pub, crossfit = cf.base_scores(X, y, folds)
    naive = {int(i): pub for i in np.unique(folds)}
    return (cf._pooled(comp, naive, y, folds) - a_comp,
            cf._pooled(comp, crossfit, y, folds) - a_comp,
            a_comp)


def usable(y, folds):
    """Every fold must carry both label classes, or an out-of-fold AUROC does not exist."""
    y, folds = np.asarray(y), np.asarray(folds)
    return all(len(np.unique(y[folds == f])) == 2 for f in np.unique(folds)) \
        and len(np.unique(folds)) == N_FOLDS


def build(limit=0):
    import external_replication as er
    from horlacher_arm import build as build_arm
    from pyfaidx import Fasta
    if not er.GENOME.exists():
        sys.exit(f"no genome at {er.GENOME}")
    if not er.DATA.exists():
        sys.exit(f"the Horlacher deposit is not unpacked at {er.DATA}")
    fa = Fasta(str(er.GENOME), as_raw=True, sequence_always_upper=True)

    sel = er.disjoint_datasets()
    keys = list(sel)[:limit or None]
    log(f"  {len(sel)} of their datasets are outside our panel; analysing {len(keys)}")

    rows, excluded = [], {"no_chrom_folds": 0, "fold_class": 0, "no_build": 0}
    for i, key in enumerate(keys, 1):
        protein, cell = key.split(":")
        rec = {"dataset": key, "protein": protein, "cell": cell}
        arms, ok = {}, True
        for neg, tag in (("negative-1", "n1"), ("negative-2", "n2")):
            d = build_arm(fa, sel[key], neg)
            if d is None:
                ok = False
                break
            arms[tag] = d
        if not ok:
            excluded["no_build"] += 1
            continue

        # Supplied folds, both estimators. This cell reproduces the published analysis.
        for tag, d in arms.items():
            ts, cfv, comp = contributions(d.seq_rna.values, d.label.values, d.fold.values)
            rec[f"supplied_2s_{tag}"] = ts
            rec[f"supplied_cf_{tag}"] = cfv
            rec[f"comp_{tag}"] = comp
            rec[f"n_{tag}"] = int(len(d))

        # Chromosome folds, both estimators. Excluded and COUNTED if the partition cannot be
        # made, which is what criterion 4's "cannot be made chromosome-blocked" means per
        # dataset, rather than being dropped quietly.
        chrom_ok = True
        for tag, d in arms.items():
            cf_folds = chrom_folds(d.chrom.values)
            if cf_folds is None:
                excluded["no_chrom_folds"] += 1
                chrom_ok = False
                break
            if not usable(d.label.values, cf_folds):
                excluded["fold_class"] += 1
                chrom_ok = False
                break
            arms[tag] = d.assign(chrom_fold=cf_folds)
        if chrom_ok:
            for tag, d in arms.items():
                ts, cfv, comp = contributions(d.seq_rna.values, d.label.values,
                                              d.chrom_fold.values)
                rec[f"chrom_2s_{tag}"] = ts
                rec[f"chrom_cf_{tag}"] = cfv
                rec[f"chrom_comp_{tag}"] = comp
            d = arms["n1"]
            rec["chrom_blocked"] = int(one_fold_per_chrom(d.chrom.values, d.chrom_fold.values))
            strand = d.strand.values if "strand" in d.columns else np.full(len(d), "+")
            have, cross = leakage(d.chrom.values, d.start.values, strand, d.chrom_fold.values)
            rec["chrom_neighbours"], rec["chrom_cross_fold"] = have, cross
            have, cross = leakage(d.chrom.values, d.start.values, strand, d.fold.values)
            rec["supplied_neighbours"], rec["supplied_cross_fold"] = have, cross

        rows.append(rec)
        log(f"  [{i:3d}/{len(keys)}] {key:20s} supplied 2s {rec['supplied_2s_n1']:+.4f}/"
            f"{rec['supplied_2s_n2']:+.4f}  cf {rec['supplied_cf_n1']:+.4f}/"
            f"{rec['supplied_cf_n2']:+.4f}"
            + (f"  chrom cf {rec['chrom_cf_n1']:+.4f}/{rec['chrom_cf_n2']:+.4f}"
               if "chrom_cf_n1" in rec else "  chrom EXCLUDED"))

    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to write an empty table over committed evidence")
    for k, v in excluded.items():
        log(f"  excluded, {k}: {v}")
    return t, excluded


# --- the estimand, its interval and its verdict ---------------------------------------------

def directional(d, prefix):
    """R = M(negative-1) / M(negative-2) and its protein-clustered percentile interval.

    Protein is the cluster, not dataset: a protein can appear in both cell lines and the two
    are not independent. Resampling is over PROTEINS with the panel means recomputed inside
    each draw, which is the amendment's specification and the paper's convention.
    """
    n1, n2 = f"{prefix}_n1", f"{prefix}_n2"
    t = d.dropna(subset=[n1, n2])
    if len(t) < MIN_DATASETS:
        return dict(n=len(t), point=np.nan, lo=np.nan, hi=np.nan,
                    m1=np.nan, m2=np.nan, verdict="indeterminate",
                    why=f"only {len(t)} datasets survive, the floor is {MIN_DATASETS}")
    m1, m2 = float(t[n1].mean()), float(t[n2].mean())
    if m2 <= 0:
        return dict(n=len(t), point=np.nan, lo=np.nan, hi=np.nan, m1=m1, m2=m2,
                    verdict="indeterminate",
                    why="the negative-2 panel mean is not positive, so the ratio is undefined")
    point = m1 / m2

    rng = np.random.default_rng(SEED)
    groups = {p: g.index.to_numpy() for p, g in t.groupby("protein")}
    names = np.array(sorted(groups))
    draws = []
    for _ in range(BOOT):
        pick = rng.choice(len(names), len(names), replace=True)
        idx = np.concatenate([groups[names[j]] for j in pick])
        a, b = t.loc[idx, n1].mean(), t.loc[idx, n2].mean()
        if b > 0:
            draws.append(a / b)
    draws = np.array(draws)
    lo, hi = (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))) \
        if len(draws) > BOOT // 2 else (np.nan, np.nan)

    if not np.isfinite(lo):
        verdict, why = "indeterminate", "over half the draws had a non-positive denominator"
    elif lo <= 1.0 <= hi or point < 1.0:
        verdict, why = "fails", "the interval contains 1.0 or the point estimate is below it"
    elif point >= SUPPORT_POINT and lo >= SUPPORT_LOW:
        verdict, why = "supports", f"point >= {SUPPORT_POINT} and lower bound >= {SUPPORT_LOW}"
    else:
        verdict, why = "weak", ("directionally consistent, under the amendment's "
                                f"{SUPPORT_POINT}/{SUPPORT_LOW} thresholds")
    return dict(n=len(t), point=point, lo=lo, hi=hi, m1=m1, m2=m2, verdict=verdict, why=why)


CELLS = (("supplied_2s", "supplied folds, two-stage (the published analysis, recomputed)"),
         ("supplied_cf", "supplied folds, CROSS-FITTED (D3: like-for-like with our primary)"),
         ("chrom_2s", "chromosome-blocked folds, two-stage (D2: the literal criterion)"),
         ("chrom_cf", "chromosome-blocked folds, cross-fitted (D2 and D3 together)"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=0, help="smoke test; writes .partial.csv")
    ap.add_argument("--from-cache", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")

    excluded = None
    if a.from_cache:
        if not PER.exists():
            sys.exit(f"no {PER.name}; run without --from-cache first")
        t = pd.read_csv(PER)
    else:
        t, excluded = build(a.n)
        if a.n:
            p = PER.with_suffix(".partial.csv")
            t.to_csv(p, index=False)
            log(f"\n  --n {a.n}: wrote {p.name}, not the committed table")
        else:
            t.to_csv(PER, index=False)

    out = []

    def add(check, value, lo="", hi="", n="", note=""):
        out.append({"check": check, "value": value, "ci_low": lo, "ci_high": hi,
                    "n": n, "note": note})

    add("datasets analysed", len(t), n=len(t),
        note="Horlacher ENCODE datasets absent from our 95-dataset panel")

    log("")
    for prefix, label in CELLS:
        r = directional(t, prefix)
        add(f"directional R, {label}", r["point"], r["lo"], r["hi"], r["n"], note=r["why"])
        add(f"verdict, {prefix}", r["verdict"], n=r["n"],
            note="amendment decision rule: supports/weak/fails/indeterminate")
        if np.isfinite(r["m1"]):
            add(f"panel mean negative-1, {prefix}", r["m1"], n=r["n"])
            add(f"panel mean negative-2, {prefix}", r["m2"], n=r["n"])
        log(f"  {label}")
        log(f"      R = {r['point']:.4f}  [{r['lo']:.4f}, {r['hi']:.4f}]   n={r['n']}"
            f"   -> {r['verdict'].upper()}"
            if np.isfinite(r["point"]) else
            f"      INDETERMINATE: {r['why']}")

    # D2's own verification: did the refold actually deliver the property, and did it reduce the
    # channel chromosome blocking is a proxy for?
    if "chrom_blocked" in t.columns:
        c = t.dropna(subset=["chrom_blocked"])
        add("datasets whose every chromosome falls in exactly one fold",
            int(c.chrom_blocked.sum()), n=len(c),
            note="must equal the count with a chromosome partition, or the refold is broken")
        for tag, lab in (("chrom", "chromosome-blocked"), ("supplied", "their supplied")):
            have = float(c[f"{tag}_neighbours"].sum())
            cross = float(c[f"{tag}_cross_fold"].sum())
            add(f"cross-fold near-neighbour fraction, {lab} folds",
                cross / have if have else np.nan, n=int(have),
                note=f"{int(cross)} of {int(have)} positives with a same-strand neighbour "
                     "within 1 kb sit in a different fold from it")
        log(f"\n  chromosome-blocked on {int(c.chrom_blocked.sum())}/{len(c)} datasets")

    if excluded is not None:
        for k, v in excluded.items():
            add(f"datasets excluded, {k}", v, n=len(t),
                note="reported rather than dropped silently")

    if a.n:
        log("  --n: summary NOT written")
        return
    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {PER.name}")


if __name__ == "__main__":
    main()
