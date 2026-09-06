"""Does the residual sequence homology across chromosome folds carry the result?

    python scripts/homology_folds.py --store ../rbp-store
    python scripts/homology_folds.py --from-cache

What was already known and what was not. Methods report exact 32-mer sharing between held-out
and training windows: a median of 2.0% of held-out windows share any 32-mer, maximum 24.3%.
Two limits on that figure were named by reviewers and both are real. It was computed holding out
FOLD 0 ONLY, so four fifths of the partition was never audited. And measuring leakage is not
controlling it: chromosome grouping prevents two windows at the same locus landing on opposite
sides of a split, and does nothing about a paralogue or a repeat on another chromosome.

Three things are reported, in increasing strength.

1. The audit over ALL FIVE held-out folds, so the quoted median and maximum describe the
   partition rather than one fifth of it.

2. A FILTER SENSITIVITY, which is the conservative control. Keep the published chromosome folds,
   delete every held-out window that shares more than half its 32-mers with a training window,
   and recompute. This removes the leakage without changing the partition, so nothing else moves.

3. THE EXHAUSTIVE FILTER, which removes ALL cross-fold exact sharing. Delete every held-out
   window that shares ANY 32-mer with a training window, at stride one. The partition is
   untouched, so chromosome blocking survives by construction, and the code asserts afterwards
   that not one indexed 32-mer crosses a fold.

Two attempts at a fold-regrouping control failed and the second failure is the interesting one.

The first assigned homology components to folds from scratch, which silently discarded
chromosome grouping: two peaks 10 kb apart share no 32-mer and became free to split. The second
tried to TIGHTEN the chromosome partition by pulling straddling components into one fold. An
audit tested it and it failed both of its own claims. It indexed a 32-mer only every eight
bases, so it linked two windows only when their shared sequence sat at the same offset modulo
eight in both, catching roughly one homologous pair in eight; and moving a component's windows
left the rest of their chromosomes behind, splitting 9 to 15 chromosomes per dataset that had
been whole. Measured on three datasets it left 560 to 981 exact 32-mers still crossing folds.

A correct regrouping is not merely unimplemented, it is impossible here, and that is worth
reporting. Requiring a partition to respect chromosome blocking AND to place every
32-mer-sharing pair in one fold means taking connected components of the graph whose nodes are
chromosomes and whose edges are shared 32-mers. Repetitive sequence connects almost everything:
on AQR:HepG2 and BCLAF1:HepG2 that graph has a SINGLE component holding 100% of windows, and on
HNRNPC:K562 the largest holds 99.6%. There is no five-fold split satisfying both constraints.
Deleting the offending windows, which is what 3 does, is the only way to satisfy both.
"""

import argparse
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.audit.leakage import build_reference, overlap_profile  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "gc", "dn": "dinuc", "neg2": "neg2"}
K = 32
CUTOFF = 0.5        # "more than half its 32-mers seen in training" = effectively a duplicate
N_FOLDS = 5


def components(seqs):
    """Connected components of the graph linking windows that share a 32-mer."""
    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    idx = defaultdict(list)
    for i, s in enumerate(seqs):
        for j in range(len(s) - K + 1):
            idx[s[j:j + K]].append(i)
    for v in idx.values():
        r0 = find(v[0])
        for x in v[1:]:
            rx = find(x)
            if rx != r0:
                parent[rx] = r0
    out = defaultdict(list)
    for i in range(len(seqs)):
        out[find(i)].append(i)
    return list(out.values())


def collapse_check(seqs, chrom):
    """Would a partition respecting BOTH chromosome blocking and homology exist?

    Nodes are chromosomes, edges join two chromosomes sharing a 32-mer, and the answer is the
    size of the largest connected component as a fraction of the dataset. Near 1 means no such
    five-fold partition exists, because everything is in one group.
    """
    chroms = sorted(set(chrom))
    ci = {c: i for i, c in enumerate(chroms)}
    parent = list(range(len(chroms)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    idx = defaultdict(set)
    for s, c in zip(seqs, chrom):
        for j in range(len(s) - K + 1):
            idx[s[j:j + K]].add(ci[c])
    for v in idx.values():
        v = list(v)
        r0 = find(v[0])
        for x in v[1:]:
            rx = find(x)
            if rx != r0:
                parent[rx] = r0
    grp = defaultdict(set)
    for i in range(len(chroms)):
        grp[find(i)].add(chroms[i])
    biggest = max(sum(1 for c in chrom if c in g) for g in grp.values())
    return len(grp), biggest / len(seqs)


def audit_all_folds(d):
    """Fraction of each held-out fold's windows echoed in its training folds, per fold."""
    seqs = d.seq_rna.tolist()
    folds = np.asarray(d.fold.values)
    per = []
    for f in np.unique(folds):
        te = [s for s, ff in zip(seqs, folds) if ff == f]
        tr = [s for s, ff in zip(seqs, folds) if ff != f]
        if not te or not tr:
            continue
        frac = overlap_profile(te, build_reference(tr, K), K)
        per.append({"fold": int(f), "any": float((frac > 0).mean()),
                    "gt50": float((frac > CUTOFF).mean())})
    return per


def pair_key(d, arm):
    """One key per matched positive/negative pair, so a filter can remove both members.

    Why this is needed. Both filters below delete held-out WINDOWS. Negatives are matched 1:1 to
    positives, so deleting one member of a pair leaves the other unpartnered, and the deletions
    are not symmetric between classes -- measured on three datasets the GC arm dropped 268
    positives against 73 negatives, and on the full panel the imbalance runs the other way. The
    model is refitted after filtering, so the resulting number confounds three changes: removal
    of sharing, destruction of the matching, and a class-balance shift. Only the first is the
    thing being measured.

    gc and dinuc encode the pair in the id: PROTEIN_pos_i partners PROTEIN_neg_i. The
    bias-aware arm does not -- its positives keep their original indices while its negatives are
    renumbered -- so the pair is recoverable only from position within a fold, which is how
    build_neg2.py constructs them. That is inferred rather than declared, so it is asserted.
    """
    ids = d.id.astype(str).to_numpy()
    y = d.label.to_numpy()
    if arm in ("gc", "dn"):
        key = np.array([i.rsplit("_", 1)[-1] for i in ids])
        if len(set(key[y == 1])) != int((y == 1).sum()):
            return None
        return key
    order = np.empty(len(d), dtype=object)
    fa = d.fold.to_numpy()
    for fv in np.unique(fa):
        for lab in (1, 0):
            sel = np.flatnonzero((fa == fv) & (y == lab))
            order[sel] = [f"{fv}:{k}" for k in range(len(sel))]
    for fv in np.unique(fa):
        if ((fa == fv) & (y == 1)).sum() != ((fa == fv) & (y == 0)).sum():
            return None
    return order


def drop_pairs(keep, key):
    """Extend a window mask to whole pairs: if either member goes, both go."""
    doomed = set(key[~keep])
    return np.array([k not in doomed for k in key])


def gain(seqs, y, folds):
    sc, _, _ = kmer_oof(seqs, y, folds, k=4)
    return float(gain_over_composition(seqs, sc, y, folds).delta)


def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    rows = []
    for n, ds in enumerate(list(pub.dataset)[:limit or None], 1):
        protein, cell = ds.split(":")
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        ok = True
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            seqs, y, folds = list(d.seq_rna), d.label.values, d.fold.values

            if arm == "dn":                      # the audit is a property of the windows
                per = audit_all_folds(d)
                rec["leak_any_max"] = max(p["any"] for p in per)
                rec["leak_any_mean"] = float(np.mean([p["any"] for p in per]))
                rec["leak_gt50_max"] = max(p["gt50"] for p in per)
                rec["leak_any_fold0"] = [p for p in per if p["fold"] == 0][0]["any"]

            rec[f"gain_{arm}"] = gain(seqs, y, folds)

            # 2. Filter: drop held-out windows echoed in their own training folds.
            keep = np.ones(len(d), dtype=bool)
            fa = np.asarray(folds)
            for fv in np.unique(fa):
                te = fa == fv
                ref = build_reference([s for s, m in zip(seqs, te) if not m], K)
                frac = overlap_profile([s for s, m in zip(seqs, te) if m], ref, K)
                keep[np.flatnonzero(te)[frac > CUTOFF]] = False
            key = pair_key(d, arm)
            if key is None:
                sys.exit(f"{ds} {arm}: matched pairs are not recoverable, so a window filter "
                         f"cannot be made pair-aware. Refusing to report a confounded number.")
            keep = drop_pairs(keep, key)
            if keep.sum() > 200 and len(np.unique(y[keep])) == 2:
                if int((y[keep] == 1).sum()) != int((y[keep] == 0).sum()):
                    sys.exit(f"{ds} {arm}: the >50% filter left the classes unbalanced")
                rec[f"gain_filtered_{arm}"] = gain([s for s, m in zip(seqs, keep) if m],
                                                   y[keep], fa[keep])
                rec[f"dropped_{arm}"] = float(1 - keep.mean())
            else:
                ok = False
                break

            # 3. The exhaustive filter: delete any held-out window sharing ANY 32-mer with
            # training. Same partition, so chromosome blocking is untouched, and the
            # assertion below proves the sharing is gone rather than assuming it.
            keep0 = np.ones(len(d), dtype=bool)
            for fv in np.unique(fa):
                te = fa == fv
                ref = build_reference([s for s, m in zip(seqs, te) if not m], K)
                frac = overlap_profile([s for s, m in zip(seqs, te) if m], ref, K)
                keep0[np.flatnonzero(te)[frac > 0]] = False
            keep0 = drop_pairs(keep0, key)
            if keep0.sum() > 200 and len(np.unique(y[keep0])) == 2:
                if int((y[keep0] == 1).sum()) != int((y[keep0] == 0).sum()):
                    sys.exit(f"{ds} {arm}: the exhaustive filter left the classes unbalanced")
                ks = [s for s, m in zip(seqs, keep0) if m]
                kf = fa[keep0]
                idx = defaultdict(set)
                for s, fv in zip(ks, kf):
                    for j in range(len(s) - K + 1):
                        idx[s[j:j + K]].add(int(fv))
                crossing = sum(1 for v in idx.values() if len(v) > 1)
                if crossing:
                    sys.exit(f"{ds} {arm}: {crossing} 32-mers still cross folds after the "
                             f"exhaustive filter. The filter is not doing what it claims.")
                rec[f"gain_nosharing_{arm}"] = gain(ks, y[keep0], kf)
                rec[f"dropped0_{arm}"] = float(1 - keep0.mean())
            else:
                ok = False
                break
            if arm == "dn":
                ngrp, frac = collapse_check(seqs, list(d.chrom))
                rec["collapse_groups"] = ngrp
                rec["collapse_largest"] = frac
        if not ok:
            continue
        rows.append(rec)
        log(f"[{n:3d}/94] {ds:18s} dn {rec['gain_dn']:+.4f} "
            f"filt {rec['gain_filtered_dn']:+.4f} nosh {rec['gain_nosharing_dn']:+.4f} "
            f"(-{100*rec['dropped0_dn']:.1f}%)")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "homology_folds_per_dataset.csv"
    t = pd.read_csv(per) if a.from_cache else build(a.store, a.n)
    if not a.from_cache:
        if a.n:
            per = per.with_suffix(".partial.csv")
            log(f"  --n {a.n} given: writing {per.name}, not the committed table")
        t.to_csv(per, index=False)

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

    add("datasets", np.full(len(t), len(t)))
    add("held-out windows sharing any 32-mer with training, mean over all five folds",
        t.leak_any_mean, "the published figure audited fold 0 only")
    add("the same, worst fold of each dataset", t.leak_any_max)
    add("held-out windows sharing more than half their 32-mers, worst fold", t.leak_gt50_max)
    add("fraction of windows the filter removes, dinucleotide arm", t.dropped_dn)
    add("fraction of windows the exhaustive filter removes, dinucleotide arm", t.dropped0_dn)
    out.append({"check": "chromosome-plus-homology groups, median over datasets",
                "value": float(t.collapse_groups.median()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "1 means no partition can satisfy both constraints"})
    out.append({"check": "largest chromosome-plus-homology group as a fraction of windows",
                "value": float(t.collapse_largest.median()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "median over datasets"})
    # The two extremes the Methods quote. Reported as plain statistics over datasets rather
    # than bootstrap means, because a maximum has no useful resampling distribution.
    out.append({"check": "median across datasets of the per-dataset mean 32-mer sharing",
                "value": float(t.leak_any_mean.median()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "over all five held-out folds"})
    out.append({"check": "maximum 32-mer sharing over every dataset and fold",
                "value": float(t.leak_any_max.max()), "ci_low": "", "ci_high": "",
                "n": len(t),
                "note": "the fold-0-only audit reported 24.3%, which understated this"})

    for arm in ARMS:
        add(f"contribution as published, {arm} arm", t[f"gain_{arm}"])
        add(f"contribution with echoed held-out windows removed, {arm} arm",
            t[f"gain_filtered_{arm}"])
        add(f"contribution with ALL cross-fold sharing removed, {arm} arm",
            t[f"gain_nosharing_{arm}"])

    def span(cols):
        m = [t[c].mean() for c in cols]
        return float(max(m) / min(m)) if min(m) > 0 else float("nan")

    for tag, pat in (("as published", "gain_{}"), ("filtered", "gain_filtered_{}"),
                     ("all cross-fold sharing removed", "gain_nosharing_{}")):
        out.append({"check": f"three-arm span, {tag}",
                    "value": span([pat.format(a_) for a_ in ARMS]),
                    "ci_low": "", "ci_high": "", "n": len(t), "note": ""})

    r = pd.DataFrame(out)
    summary = TABLES / ("homology_folds.partial.csv" if a.n else "homology_folds.csv")
    r.to_csv(summary, index=False)
    log("")
    for _, x in r.iterrows():
        log(f"  {x['check']:66s} {x['value']:+.4f}")
    log(f"\nwrote {summary.name} and {per.name}")


if __name__ == "__main__":
    main()
