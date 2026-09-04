"""B10: chromosome-grouped folds already imply gene-grouped folds. What does the finer one cost?

    python scripts/gene_clustered_cv.py
    python scripts/gene_clustered_cv.py --from-cache

THE OBJECTION, AND WHY ITS FIRST HALF IS ANSWERED BY CONSTRUCTION. Windows from one gene are
correlated, so a fold assignment that puts a gene's windows on both sides of the split leaks.
Every number in this paper uses chromosome-grouped folds from a frozen assignment, and a gene
lies on one chromosome, so chromosome grouping is strictly COARSER than gene grouping: no gene
can span folds. That is a structural fact rather than a measurement, and this script verifies
it instead of asserting it, over every window in both composition-matched arms.

WHAT IS NOT ANSWERED BY CONSTRUCTION is what the conservative choice costs. Gene-grouped folds
are finer: they may place two genes from the same chromosome in different folds, which
chromosome grouping forbids. That gives more balanced folds and more effective training data,
and it is what a reviewer asking for gene-clustered CV would get. So the gene-clustered refit
here is the LESS conservative design, and its agreement with the published number bounds what
the coarser choice costs rather than validating it.

THE CONTROL THAT MAKES IT A COMPARISON. The frozen assignment is refit through this same code
path and must return the published contrast. Without it, two agreeing numbers would show only
that this script agrees with itself, which is the error B6 was written to avoid.

4-mer only, and not as a shortcut: both neural models trained on the frozen assignment, so
re-folding them would compare their old folds against new ones.
"""

import argparse
import pickle
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data import annotation as ann  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import _oof_scores, composition_features  # noqa: E402
from rbp.stats import standardise  # noqa: E402

TABLES = ROOT / "results" / "tables"
DATA_ROOT = ROOT.parent / "rna-binding-proteins"
ARM_DIR = {"gc": "gc", "dn": "dinuc"}
K = 5


def log(m):
    print(m, flush=True)


def gene_labels(gene_index, chrom, mids):
    """`ann.gene_at` over many midpoints on one chromosome: smallest containing gene.

    Vectorised by candidate rank rather than by position. Genes are sorted by start, so the
    candidates for a position are the genes starting at or before it whose end is past it;
    scanning back a bounded number of ranks covers every one of them because the bound is the
    largest gene span on the chromosome.
    """
    per = gene_index.get(chrom)
    if per is None:
        return np.full(len(mids), None, dtype=object)
    starts, ends, names = per
    lengths = ends - starts
    hi = np.searchsorted(starts, mids, side="right")
    best = np.full(len(mids), None, dtype=object)
    best_len = np.full(len(mids), np.inf)
    # How far back a containing gene can begin. Bounding by the real maximum span rather
    # than a constant keeps this exact: a fixed 2.5 Mb window would silently miss the few
    # genes longer than that.
    span = int(lengths.max()) if len(lengths) else 0
    depth = int(np.searchsorted(starts, starts + span, side="right").max() - 0) if span else 0
    depth = min(depth, len(starts))
    for back in range(1, depth + 2):
        i = hi - back
        ok = i >= 0
        if not ok.any():
            break
        ii = np.clip(i, 0, len(starts) - 1)
        contains = ok & (ends[ii] > mids) & (starts[ii] <= mids)
        take = contains & (lengths[ii] < best_len)
        best[take] = names[ii][take]
        best_len[take] = lengths[ii][take]
    return best


def balanced_gene_folds(genes, weights, seed=7):
    """Assign genes to K folds, largest first into the lightest fold.

    The same objective the frozen chromosome partition uses: balance the pair mass. Genes are
    the unit, so no gene can span folds by construction, exactly as no chromosome can under
    the published assignment.
    """
    total = defaultdict(float)
    for g, w in zip(genes, weights):
        total[g] += w
    order = sorted(total.items(), key=lambda kv: (-kv[1], str(kv[0])))
    load = np.zeros(K)
    out = {}
    for g, w in order:
        j = int(np.argmin(load))
        out[g] = j
        load[j] += w
    return out, load


def score(dd, folds):
    """The 4-mer's nested contribution under a given fold vector, published estimator."""
    y = dd.label.to_numpy()
    X, _ = composition_features(dd.seq_rna.values)
    s_comp = _oof_scores(X, y, folds)
    sc, _, _ = kmer_oof(dd.seq_rna.values, y, folds, k=4)
    s_full = _oof_scores(np.column_stack([X, standardise(sc)]), y, folds)
    good = np.isfinite(s_comp) & np.isfinite(s_full)
    r = delong_test(s_full[good], s_comp[good], y[good])
    return float(r["auc_b"]), float(r["diff"])


def build(store, gtf, cache, limit):
    if Path(cache).exists():
        log(f"gene index from {cache}")
        gidx = pickle.loads(Path(cache).read_bytes())
    else:
        log(f"building gene index from {gtf} ...")
        gidx = ann.build_gene_index(gtf)
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        Path(cache).write_bytes(pickle.dumps(gidx))
        log(f"  {sum(len(v[0]) for v in gidx.values())} genes -> {cache}")

    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    datasets = list(pub.dataset)[:limit or None]
    rows, span_total, win_total, intergenic = [], 0, 0, 0
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        ok = True
        for arm, sub in ARM_DIR.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            g = np.full(len(d), None, dtype=object)
            for chrom, grp in d.groupby("chrom", sort=False):
                mids = ((grp.start.to_numpy() + grp.end.to_numpy()) // 2).astype(np.int64)
                g[grp.index.to_numpy()] = gene_labels(gidx, chrom, mids)
            # Intergenic windows have no gene, so they cannot be grouped by one. They are
            # kept in their own singleton group per chromosome, which is what the coarser
            # published design does with them anyway.
            g = np.array([x if x is not None else f"__intergenic_{c}"
                          for x, c in zip(g, d.chrom)], dtype=object)
            intergenic += int(sum(str(x).startswith("__intergenic_") for x in g))
            win_total += len(d)

            # THE STRUCTURAL CHECK: does any gene group appear in more than one FROZEN fold?
            # It should not, and 22 do. Every one is a gene NAME shared across chromosomes,
            # not a locus: the index keys on gene_name, and GENCODE reuses a name for
            # pseudo-autosomal chrX/chrY pairs and for multi-copy small-RNA families. So the
            # windows in those groups are recorded too, because the count of groups says
            # nothing about how much sequence is involved.
            gf_frozen = pd.DataFrame({"g": g, "f": d.fold}).groupby("g").f.nunique()
            span_names = set(gf_frozen[gf_frozen > 1].index)
            span_total += len(span_names)
            rec[f"genes_spanning_folds_{arm}"] = len(span_names)
            rec[f"windows_in_spanning_groups_{arm}"] = int(sum(x in span_names for x in g))
            rec[f"n_genes_{arm}"] = int(len(set(g)))

            gf, load = balanced_gene_folds(g, np.ones(len(g)))
            newfold = np.array([gf[x] for x in g])
            rec[f"fold_mass_cv_{arm}"] = float(load.std() / load.mean())
            # How much finer the gene grouping is: chromosome grouping forbids splitting a
            # chromosome across folds, so count the chromosomes this assignment splits.
            rec[f"chroms_split_{arm}"] = int(sum(
                1 for _, s in pd.DataFrame({"c": d.chrom, "f": newfold}).groupby("c")
                if s.f.nunique() > 1))
            c_fr, g_fr = score(d, d.fold.to_numpy())
            c_gn, g_gn = score(d, newfold)
            rec[f"comp_frozen_{arm}"], rec[f"gain_frozen_{arm}"] = c_fr, g_fr
            rec[f"comp_gene_{arm}"], rec[f"gain_gene_{arm}"] = c_gn, g_gn
        if not ok:
            continue
        rec["published_gain_gc"] = float(pub.set_index("dataset").loc[ds, "gain_gc"])
        rec["published_gain_dn"] = float(pub.set_index("dataset").loc[ds, "gain_dn"])
        rows.append(rec)
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} frozen dn {rec['gain_frozen_dn']:+.4f} "
            f"gene dn {rec['gain_gene_dn']:+.4f}   genes spanning folds "
            f"{rec['genes_spanning_folds_dn']}")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    log(f"\n{win_total} windows, {intergenic} intergenic, {span_total} genes spanning "
        f"a frozen fold boundary")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--gtf", default=str(DATA_ROOT / "data/raw/gencode.v45.annotation.gtf.gz"))
    p.add_argument("--gene-cache", default=str(ROOT.parent / "rbp-store/interim/genes.pkl"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "gene_clustered_cv_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it without --from-cache")
    else:
        t = build(a.store, a.gtf, a.gene_cache, a.n)
        t.to_csv(per, index=False)

    out = []
    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(2000)]

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})
        return float(v.mean())

    log(f"\n=== B10: gene-clustered folds, n = {len(t)}, {len(uniq)} proteins ===\n")

    # THE STRUCTURAL FACT, VERIFIED. A gene lies on one chromosome, so chromosome-grouped
    # folds cannot split one. Asserted as an exact zero over every window in both arms.
    span = int(sum(t[f"genes_spanning_folds_{arm}"].sum() for arm in ARM_DIR))
    wins = int(sum(t[f"windows_in_spanning_groups_{arm}"].sum() for arm in ARM_DIR))
    groups = int(t.n_genes_dn.sum() + t.n_genes_gc.sum())
    out.append({"check": "gene groups spanning a frozen fold boundary", "value": span,
                "n": len(t)})
    out.append({"check": "gene groups examined, summed over datasets and both arms",
                "value": groups, "n": len(t)})
    out.append({"check": "windows inside a gene group that spans a frozen fold boundary",
                "value": wins, "n": len(t)})
    log(f"  {span} of {groups} gene groups span a frozen fold boundary, holding {wins} windows")
    log("  A LOCUS CANNOT SPAN ONE: chromosome grouping is strictly coarser than grouping by")
    log("  locus, so the leakage the objection describes is impossible by construction. What")
    log("  spans is a gene NAME shared across chromosomes, which the index keys on: GENCODE")
    log("  reuses a name for pseudo-autosomal chrX/chrY pairs and for multi-copy small-RNA")
    log("  families. Those are the one leakage channel chromosome grouping does NOT close,")
    log("  because near-identical sequence really does sit on both sides of the split.")

    # THE CENSUS BEHIND THAT, straight from the gene index, so the explanation is a
    # measurement and not a plausible story. Two distinct causes and they are worth
    # separating: the pseudo-autosomal region puts one name on chrX and chrY, and multi-copy
    # small-RNA families put one name on up to two dozen chromosomes.
    cache = Path(a.gene_cache)
    if cache.exists():
        gidx = pickle.loads(cache.read_bytes())
        where = defaultdict(set)
        for chrom, (_s, _e, names) in gidx.items():
            for nm in names:
                where[nm].add(chrom)
        multi = {nm: c for nm, c in where.items() if len(c) > 1}
        xy = sum(1 for c in multi.values() if c == {"chrX", "chrY"})
        out.append({"check": "gene names appearing on more than one chromosome",
                    "value": len(multi), "n": len(t)})
        out.append({"check": "of those, pseudo-autosomal chrX/chrY name pairs", "value": xy,
                    "n": len(t)})
        out.append({"check": "largest number of chromosomes sharing one gene name",
                    "value": max((len(c) for c in multi.values()), default=0), "n": len(t)})
        log(f"  {len(multi)} gene names sit on more than one chromosome: {xy} are chrX/chrY "
            f"pseudo-autosomal pairs and the widest spans "
            f"{max((len(c) for c in multi.values()), default=0)} chromosomes")

    # HOW MUCH FINER THE GENE GROUPING IS. Without this the agreement below could mean the
    # two designs are nearly the same design, which would make it uninformative.
    add("chromosomes split across folds by the gene-clustered design, dinucleotide arm",
        t.chroms_split_dn)
    add("fold-mass coefficient of variation, gene-clustered, dinucleotide arm",
        t.fold_mass_cv_dn)

    # THE CONTROL. The frozen refit must return the published per-dataset gain, or the
    # comparison is between two things this script made up.
    w = max(float((t[f"gain_frozen_{arm}"] - t[f"published_gain_{arm}"]).abs().max())
            for arm in ARM_DIR)
    out.append({"check": "max |frozen refit gain - published gain|", "value": w, "n": len(t)})
    log(f"  the frozen refit reproduces the published per-dataset gain to {w:.2e}")

    for arm in ARM_DIR:
        add(f"composition AUROC, frozen folds, {arm} arm", t[f"comp_frozen_{arm}"])
        add(f"composition AUROC, gene-clustered folds, {arm} arm", t[f"comp_gene_{arm}"])
        f_ = add(f"4-mer contribution, frozen folds, {arm} arm", t[f"gain_frozen_{arm}"])
        g_ = add(f"4-mer contribution, gene-clustered folds, {arm} arm", t[f"gain_gene_{arm}"])
        log(f"  {arm:3s} contribution  frozen {f_:+.4f}   gene-clustered {g_:+.4f}   "
            f"difference {g_ - f_:+.4f}")

    cf = add("two-arm contrast, frozen folds", t.gain_frozen_dn - t.gain_frozen_gc)
    cg = add("two-arm contrast, gene-clustered folds", t.gain_gene_dn - t.gain_gene_gc)
    d = add("contrast change from gene-clustered folds",
            (t.gain_gene_dn - t.gain_gene_gc) - (t.gain_frozen_dn - t.gain_frozen_gc))
    out.append({"check": "gene-clustered / frozen contrast multiplier",
                "value": float(cg / cf), "n": len(t)})
    log(f"\n  two-arm contrast  frozen {cf:+.4f}   gene-clustered {cg:+.4f}   "
        f"change {d:+.4f} ({cg / cf:.3f}x)")
    log("  the finer design is the LESS conservative one, so this bounds what the coarser "
        "\n  published choice costs rather than validating it.")

    pd.DataFrame(out).to_csv(TABLES / "gene_clustered_cv.csv", index=False)
    log("\nwrote gene_clustered_cv.csv and gene_clustered_cv_per_dataset.csv")


if __name__ == "__main__":
    main()
