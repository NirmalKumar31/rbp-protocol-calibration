"""Does the residual sequence homology across chromosome folds carry the result?

    python scripts/homology_folds.py --store ../rbp-store
    python scripts/homology_folds.py --from-cache

WHAT WAS ALREADY KNOWN AND WHAT WAS NOT. Methods report exact 32-mer sharing between held-out
and training windows: a median of 2.0% of held-out windows share any 32-mer, maximum 24.3%.
Two limits on that figure were named by reviewers and both are real. It was computed holding out
FOLD 0 ONLY, so four fifths of the partition was never audited. And measuring leakage is not
controlling it: chromosome grouping prevents two windows at the same locus landing on opposite
sides of a split, and does nothing about a paralogue or a repeat on another chromosome.

THREE THINGS ARE REPORTED, IN INCREASING STRENGTH.

1. The audit over ALL FIVE held-out folds, so the quoted median and maximum describe the
   partition rather than one fifth of it.

2. A FILTER SENSITIVITY, which is the conservative control. Keep the published chromosome folds,
   delete every held-out window that shares more than half its 32-mers with a training window,
   and recompute. This removes the leakage without changing the partition, so nothing else moves.

3. HOMOLOGY-TIGHTENED FOLDS, which is the control reviewers asked for. Windows are linked when
   they share a 32-mer; any connected component of that graph straddling two chromosome folds is
   pulled wholly into one of them. No homologous pair can then straddle a split, and chromosome
   grouping is kept everywhere the two do not conflict. Components are small -- the largest holds
   0.5% of a dataset's windows -- so this moves a fraction of a percent of windows.

THE FIRST VERSION OF 3 WAS BACKWARDS AND IS WORTH RECORDING. It assigned components to folds
from scratch, balancing size, which discarded chromosome grouping entirely: two peaks 10 kb
apart share no 32-mer, so they became free to land on opposite sides of a split. Contributions
rose, and presenting that as a leakage control would have invited the obvious reply that the
partition was weaker rather than stronger. Tightening the published partition instead is
strictly stronger than it, so if the span survives it, it survives both objections at once.
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
STRIDE = 8          # every 8th 32-mer: any real overlap of 101-nt windows still shares one
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
        for j in range(0, len(s) - K + 1, STRIDE):
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


def homology_folds(seqs, y, chrom_folds):
    """Chromosome folds with every straddling homology component pulled into ONE fold.

    THIS ADDS TO THE PUBLISHED PARTITION RATHER THAN REPLACING IT, and the first version of this
    function did the opposite. It assigned components to folds from scratch, balancing size,
    which silently discarded chromosome grouping: two peaks 10 kb apart share no 32-mer, so they
    became free to land on opposite sides of a split. Contributions duly rose, and calling that
    a leakage control would have been backwards -- the partition was weaker, not stronger, and a
    referee would say so in one sentence.

    A component that already sits inside one fold is left alone, which is almost all of them.
    A component that straddles folds is moved wholly into the fold holding most of it, breaking
    ties toward the smaller fold. The result still respects chromosome grouping everywhere the
    two do not conflict, and where they do conflict homology wins. It is therefore strictly
    stronger than the published partition: no homologous pair straddles a split, and no locus
    pair does either except where a homology component forced the move.
    """
    fold = np.asarray(chrom_folds).copy()
    moved = 0
    for c in components(seqs):
        if len(c) < 2:
            continue
        here = fold[c]
        if len(np.unique(here)) == 1:
            continue
        counts = np.bincount(here, minlength=N_FOLDS)
        best = np.flatnonzero(counts == counts.max())
        # Tie-break toward the fold that is currently smaller, so the moves do not all pile up.
        sizes = np.bincount(fold, minlength=N_FOLDS)
        target = int(best[np.argmin(sizes[best])])
        fold[c] = target
        moved += len(c)
    y = np.asarray(y)
    if any(len(np.unique(y[fold == f])) < 2 for f in np.unique(fold)):
        return None, moved
    return fold, moved


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
            if keep.sum() > 200 and len(np.unique(y[keep])) == 2:
                rec[f"gain_filtered_{arm}"] = gain([s for s, m in zip(seqs, keep) if m],
                                                   y[keep], fa[keep])
                rec[f"dropped_{arm}"] = float(1 - keep.mean())
            else:
                ok = False
                break

            # 3. Homology-grouped folds, replacing the chromosome partition.
            hf, moved = homology_folds(seqs, y, folds)
            if hf is None:
                ok = False
                break
            rec[f"gain_homology_{arm}"] = gain(seqs, y, hf)
            rec[f"moved_{arm}"] = moved / len(d)
        if not ok:
            continue
        rows.append(rec)
        log(f"[{n:3d}/94] {ds:18s} dn {rec['gain_dn']:+.4f} "
            f"filt {rec['gain_filtered_dn']:+.4f} homol {rec['gain_homology_dn']:+.4f}")
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
    add("fraction of windows the tightening moves, dinucleotide arm", t.moved_dn)
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
        add(f"contribution under homology-grouped folds, {arm} arm", t[f"gain_homology_{arm}"])

    def span(cols):
        m = [t[c].mean() for c in cols]
        return float(max(m) / min(m)) if min(m) > 0 else float("nan")

    for tag, pat in (("as published", "gain_{}"), ("filtered", "gain_filtered_{}"),
                     ("homology-grouped folds", "gain_homology_{}")):
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
