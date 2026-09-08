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


def shared_map(arms):
    """ONE chromosome-to-fold map per dataset, identical for both arms. Returns None if the
    partition cannot be built.

    The defect this replaces. `chrom_folds` was called once per arm, inside the loop over
    negative-1 and negative-2. It balances fold sizes by the counts it is given; the arms share
    their positives but differ in their negatives, so the counts differ and the greedy
    assignment diverges. Measured across all 135 datasets, the two maps were identical in ZERO
    of them. Each arm was still individually chromosome-blocked, which is why the verifier
    passed: it asserted one fold per chromosome, true of each arm separately, and never
    asserted the maps agreed. So a fold-design difference was folded into a comparison whose
    entire purpose is to hold the folds fixed and vary only the negative construction.

    The repair needs an ARM-INVARIANT source, and the positives are one: they are identical in
    both arms, which this function verifies rather than assumes. The map is built from the
    positive chromosome counts alone, so it cannot depend on either negative set.

    Chromosomes appearing only among negatives are rare, zero or one per dataset, but they
    exist. They are assigned after the positives, in NAME order, each to the currently emptiest
    fold, and the set is taken as the UNION over both arms so that the result does not depend
    on which arm is processed first.
    """
    tags = sorted(arms)
    pos = {t: arms[t][arms[t].label == 1] for t in tags}

    IDENTITY = ("chrom", "start", "end", "strand")

    def _key(d):
        """A sorted MULTISET over the full window identity.

        Comparing SETS of (chrom, start) would ignore strand, ignore the end coordinate and
        ignore multiplicity, so two arms differing in any of those would still pass the check
        this assertion exists to make.

        The four columns are REQUIRED, not taken if present. An earlier version composed the
        key from whichever of `end` and `strand` happened to exist, which means a caller that
        stops supplying one gets a silently weaker check instead of an error: the same failure
        shape as the "+"-for-every-strand fallback this whole correction exists to undo.
        """
        missing = [c for c in IDENTITY if c not in d.columns]
        if missing:
            sys.exit(f"the window table is missing {missing}, so a full window identity cannot "
                     "be built. This check is not allowed to degrade silently")
        return sorted(map(tuple, d[list(IDENTITY)].astype(str).to_numpy().tolist()))

    if len({tuple(_key(pos[t])) for t in tags}) != 1:
        sys.exit("the two arms do not share an identical positive multiset, so no "
                 "arm-invariant map exists")

    base = chrom_folds(pos[tags[0]].chrom.values)
    if base is None:
        return None
    who = pos[tags[0]].chrom.values
    cmap = {}
    load = [0] * N_FOLDS
    for c, f in zip(who, base):
        cmap.setdefault(c, int(f))
        load[int(f)] += 1

    # Negative-only chromosomes, from the union over arms, in name order, deterministically.
    extra = sorted({c for t in tags for c in arms[t].chrom.unique()} - set(cmap))
    for c in extra:
        f = min(range(N_FOLDS), key=lambda j: (load[j], j))
        cmap[c] = f
        load[f] += int(sum((arms[t].chrom == c).sum() for t in tags))
    return cmap


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
        cmap = shared_map(arms)
        if cmap is None:
            excluded["no_chrom_folds"] += 1
            chrom_ok = False
        else:
            for tag, d in arms.items():
                cf_folds = np.array([cmap[c] for c in d.chrom.values])
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
            # The property the old gate did not check: one map, both arms, byte-identical.
            per_arm = [{c: int(f) for c, f in zip(a.chrom.values, a.chrom_fold.values)}
                       for a in arms.values()]
            shared_keys = set(per_arm[0]) & set(per_arm[1])
            rec["map_identical"] = int(all(per_arm[0][c] == per_arm[1][c] for c in shared_keys))
            rec["map_shared_chroms"] = len(shared_keys)
            rec["chrom_blocked"] = int(all(
                one_fold_per_chrom(a.chrom.values, a.chrom_fold.values) for a in arms.values()))
            d = arms["n1"]
            # No fallback. A missing strand column is a hard failure, not a "+" for every
            # window: the silent fallback is what turned this into a strand-agnostic number
            # published under a same-strand label.
            for tag, dd in arms.items():
                if "strand" not in dd.columns:
                    sys.exit("horlacher_arm.windows() did not return strand; a same-strand "
                             "metric cannot be computed without it")
                if not set(dd.strand.unique()) <= {"+", "-"}:
                    sys.exit(f"unexpected strand values: {sorted(set(dd.strand.unique()))}")
                for lab, folds in (("chrom", dd.chrom_fold.values), ("supplied", dd.fold.values)):
                    have, cross = leakage(dd.chrom.values, dd.start.values, dd.strand.values,
                                          folds)
                    rec[f"{lab}_neighbours_{tag}"], rec[f"{lab}_cross_fold_{tag}"] = have, cross
            rec["minus_strand_fraction"] = float((arms["n1"].strand == "-").mean())

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

    # P0-4: the no-shared-protein subset, for EVERY cell rather than only the original one.
    # 31 of the 108 external proteins also appear in our panel, in different cell lines or
    # experiments, so the two samples share biology even where they share no dataset.
    import external_replication as er
    # our_panel() returns "PROTEIN:CELL" dataset keys; the overlap that matters is at the
    # PROTEIN level, because the same protein in another cell line is still shared biology.
    ours = {k.split(":")[0] for k in er.our_panel()}
    excl = t[~t.protein.isin(ours)] if ours else t.iloc[0:0]
    log(f"  no-shared-protein subset: {len(excl)} datasets over "
        f"{excl.protein.nunique() if len(excl) else 0} proteins "
        f"({len(t) - len(excl)} dropped, sharing {len(set(t.protein) & ours)} proteins)")

    log("")
    for prefix, label in CELLS:
        r = directional(t, prefix)
        add(f"directional R, {label}", r["point"], r["lo"], r["hi"], r["n"], note=r["why"])
        add(f"verdict, {prefix}", r["verdict"], n=r["n"],
            note="amendment decision rule: supports/weak/fails/indeterminate")
        if np.isfinite(r["m1"]):
            add(f"panel mean negative-1, {prefix}", r["m1"], n=r["n"])
            add(f"panel mean negative-2, {prefix}", r["m2"], n=r["n"])
        if len(excl) >= MIN_DATASETS:
            e = directional(excl, prefix)
            add(f"no-shared-protein R, {label}", e["point"], e["lo"], e["hi"], e["n"],
                note="every protein appearing in our panel excluded; " + e["why"])
            add(f"no-shared-protein verdict, {prefix}", e["verdict"], n=e["n"])
        log(f"  {label}")
        log(f"      R = {r['point']:.4f}  [{r['lo']:.4f}, {r['hi']:.4f}]   n={r['n']}"
            f"   -> {r['verdict'].upper()}"
            if np.isfinite(r["point"]) else
            f"      INDETERMINATE: {r['why']}")

    # D2's own verification: did the refold actually deliver the property, and did it reduce the
    # channel chromosome blocking is a proxy for?
    if "chrom_blocked" in t.columns:
        c = t.dropna(subset=["chrom_blocked"])
        add("datasets where BOTH arms share one identical chromosome map",
            int(c.map_identical.sum()) if "map_identical" in c.columns else 0, n=len(c),
            note="the property the first version did not have, and the gate did not check: "
                 "the fold partition must be held FIXED while the negative construction varies")
        add("datasets whose every chromosome falls in exactly one fold",
            int(c.chrom_blocked.sum()), n=len(c),
            note="must equal the count with a chromosome partition, or the refold is broken")
        # PER ARM, and named exactly. The first version computed this for negative-1 only,
        # with strand discarded, and called it a same-strand rate over "positives and
        # negatives together" without naming which negative construction.
        for tag, lab in (("chrom", "chromosome-blocked"), ("supplied", "their supplied")):
            for arm, armlab in (("n1", "negative-1"), ("n2", "negative-2")):
                hcol, ccol = f"{tag}_neighbours_{arm}", f"{tag}_cross_fold_{arm}"
                if hcol not in c.columns:
                    continue
                have, cross = float(c[hcol].sum()), float(c[ccol].sum())
                add(f"cross-fold same-strand neighbour fraction, {lab} folds, {armlab}",
                    cross / have if have else np.nan, n=int(have),
                    note=f"{int(cross)} of {int(have)}. DENOMINATOR: windows, positives and "
                         f"{armlab} negatives together, that HAVE an eligible same-strand "
                         "neighbour within 1 kb. Not all windows, and not positives only, "
                         "which is what external_replication.fold_blocking counts")
        if "minus_strand_fraction" in c.columns:
            add("fraction of windows on the minus strand", float(c.minus_strand_fraction.mean()),
                n=len(c), note="near a half, which is why discarding strand changed the answer")
        log(f"\n  chromosome-blocked on {int(c.chrom_blocked.sum())}/{len(c)} datasets")

    # --from-cache cannot recompute an exclusion count: the excluded datasets are, by
    # definition, absent from the per-dataset table it reads. Dropping the rows would take a
    # clean checkout from passing to failing by running a documented command, which is the
    # defect rbp.utils.carry exists to close. The committed row goes back byte-identical.
    from rbp.utils.carry import emit as carry_emit
    for k in ("no_chrom_folds", "fold_class", "no_build"):
        if excluded is not None:
            add(f"datasets excluded, {k}", excluded[k], n=len(t),
                note="reported rather than dropped silently")
        else:
            carry_emit(out, OUT, f"datasets excluded, {k}", None, recomputed=False)

    if a.n:
        log("  --n: summary NOT written")
        return
    pd.DataFrame(out).to_csv(OUT, index=False)
    log(f"\n  wrote {OUT.name} and {PER.name}")


if __name__ == "__main__":
    main()
