"""B11: does the matching ALGORITHM change what the protocol measures?

    python scripts/matching_robustness.py --n 15
    python scripts/matching_robustness.py --from-cache

TWO FREE PARAMETERS IN THE DINUCLEOTIDE ARM, both chosen once and never varied. Candidates are
drawn at 8 per positive (`pool_multiple=8`), and each positive takes the nearest unused
candidate in a GREEDY pass over positives in bucket order. Greedy is defended in
`negatives.py` on the grounds that an exact assignment "buys very little here and costs a great
deal", which was an assertion. This measures it, over three pool sizes and both algorithms.

WHAT THE EXACT VERSION IS AND IS NOT. Assignment decomposes exactly by (region, chromosome),
because a candidate can only serve a positive drawn from the same pool, so the global problem is
a set of small independent ones and `linear_sum_assignment` solves each to optimality. It is
optimal within the sampled candidate pool, not over the genome: a larger pool is a different
problem, which is why pool size is varied alongside it rather than held fixed.

THE PREDICTION, and it is the paper's own thesis applied to itself. Better matching raises the
composition baseline, and a higher baseline leaves less for a model to add. So a better
assignment should REDUCE the measured contribution rather than leave it alone, and the size of
that reduction is how much of the headline is an artefact of one implementation choice. A
protocol label like "dinucleotide-matched" does not pin the number down.

Positives come from the committed window tables, unchanged, so only the negatives differ between
settings. Folds come with the positives, so the chromosome grouping is identical throughout.
4-mer only: the neural models trained on the frozen negatives and cannot be re-scored here.
"""

import argparse
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data import negatives as neg  # noqa: E402
from rbp.data import windows as win  # noqa: E402
from rbp.data.encode import peak_path  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import composition_features  # noqa: E402
from rbp.eval.nested import _oof_scores  # noqa: E402
from rbp.stats import standardise  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

TABLES = ROOT / "results" / "tables"
DATA_ROOT = ROOT.parent / "rna-binding-proteins"
# (pool_multiple, pool_min). THE FLOOR HAD TO BE SCALED TOO, and finding that out is half of
# what this section reports. n_want is max(pool_min, pool_multiple * bucket size), and buckets
# are small: at the published pool_min of 1500 the multiple binds only for a bucket holding
# more than 188 positives. Varying 4x/8x/16x with the floor fixed returned byte-identical
# negative sets on the first dataset tested, because the floor bound in every bucket. So the
# floor is scaled with the multiple, and how often it binds at the published setting is
# measured and reported.
POOLS = ((4, 750), (8, 1500), (16, 3000))
ALGOS = ("greedy", "optimal")
MAX_COST_CELLS = 5e7          # per bucket; above this the candidate set is subsampled


def log(m):
    print(m, flush=True)


def buckets_of(positives):
    """The matcher's own grouping: one candidate pool per (region, chromosome)."""
    b = {}
    for i, p in enumerate(positives):
        b.setdefault((p["region"], p["chrom"]), []).append(i)
    return b


def assign_optimal(positives, peaks, fasta, region_index, size, min_peak_distance,
                   seed, pool_multiple, pool_min=1500, pools=None):
    """build_negatives_dinuc's construction with an exact assignment per bucket.

    Deliberately a near-copy of build_negatives_dinuc rather than a flag on it. The published
    function produced every committed dataset; adding a branch inside it would put the paper's
    primary artefact and an exploratory variant on one code path, and the greedy arm here has
    to be the published function itself for the comparison to mean anything.
    """
    from scipy.optimize import linear_sum_assignment

    rng = np.random.default_rng(seed)
    if pools is None:
        excl = neg.exclusion_zones(peaks, min_peak_distance)
        regions = {p["region"] for p in positives}
        pools = {r: neg.available(region_index, r, excl, size) for r in regions}
    buckets = {}
    for i, p in enumerate(positives):
        buckets.setdefault((p["region"], p["chrom"]), []).append(i)

    rows = [None] * len(positives)
    dists = np.full(len(positives), np.nan)
    dropped = {"no_pool": 0, "no_match": 0, "subsampled_buckets": 0}
    used_global = set()
    n_di = size - 1

    for (region, chrom), idxs in buckets.items():
        pool = pools.get(region, {})
        if chrom not in pool:
            dropped["no_pool"] += len(idxs)
            continue
        n_want = max(pool_min, pool_multiple * len(idxs))
        starts, seqs, cand = neg.candidate_pool(fasta, pool, chrom, size, n_want, rng,
                                                drop_n=True, normalise=False)
        if len(starts) == 0:
            dropped["no_pool"] += len(idxs)
            continue
        # Drop candidates already taken by an earlier bucket, for the same reason the greedy
        # path keeps a global `used` set: one interval can sit in two region pools.
        free = np.array([(chrom, int(s)) not in used_global for s in starts])
        if not free.any():
            dropped["no_match"] += len(idxs)
            continue
        starts, seqs = starts[free], [s for s, f in zip(seqs, free) if f]
        cand = cand[free]
        # MEMORY GUARD, reported rather than silent. The cost matrix is n x M dense; a bucket
        # with thousands of positives and a 16x pool would allocate gigabytes. Subsampling
        # candidates keeps the assignment exact for the problem it is given and is recorded so
        # the "optimal" label cannot be read as unconditional.
        if len(idxs) * len(starts) > MAX_COST_CELLS:
            keep = rng.choice(len(starts), int(MAX_COST_CELLS // len(idxs)), replace=False)
            keep.sort()
            starts, cand = starts[keep], cand[keep]
            seqs = [seqs[j] for j in keep]
            dropped["subsampled_buckets"] += 1
        target = neg.dinuc_matrix([positives[i]["seq_dna"] for i in idxs], normalise=False)
        cost = np.abs(target[:, None, :] - cand[None, :, :]).sum(axis=2)
        if cost.shape[1] < cost.shape[0]:
            dropped["no_match"] += cost.shape[0] - cost.shape[1]
        r, c = linear_sum_assignment(cost)
        for row, col in zip(r, c):
            i = idxs[row]
            j = int(col)
            used_global.add((chrom, int(starts[j])))
            p = positives[i]
            dna = seqs[j]
            dists[i] = float(cost[row, col]) / n_di
            rows[i] = {"chrom": chrom, "start": int(starts[j]),
                       "end": int(starts[j]) + size, "strand": p["strand"],
                       "region": region, "gc": round(win.gc_content(dna), 4),
                       "seq_dna": dna, "seq_rna": win.to_rna(dna, p["strand"])}
    return rows, dropped, dists


def score(pos_df, negs):
    """The 4-mer's nested contribution on one set of negatives, published estimator."""
    keep = [i for i, m in enumerate(negs) if m is not None]
    if len(keep) < 200:
        return None
    p = pos_df.iloc[keep].reset_index(drop=True)
    n = p.copy()
    n["seq_rna"] = [negs[i]["seq_rna"] for i in keep]
    n["label"] = 0
    dd = pd.concat([p, n], ignore_index=True)
    y, fo = dd.label.values, dd.fold.values
    X, _ = composition_features(dd.seq_rna.values)
    s_comp = _oof_scores(X, y, fo)
    sc, _, _ = kmer_oof(dd.seq_rna.values, y, fo, k=4)
    s_full = _oof_scores(np.column_stack([X, standardise(sc)]), y, fo)
    good = np.isfinite(s_comp) & np.isfinite(s_full)
    r = delong_test(s_full[good], s_comp[good], y[good])
    return {"pairs": len(p), "comp": float(r["auc_b"]), "full": float(r["auc_a"]),
            "gain": float(r["diff"])}


def build(store, index_path, fasta_path, want):
    from pyfaidx import Fasta

    cfg = cfgmod.load()
    size = cfg.windows["size"]
    mpd = cfg.negatives["min_peak_distance"]
    log("loading region index and genome ...")
    index = pickle.loads(Path(index_path).read_bytes())
    fasta = Fasta(str(fasta_path))

    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    pub = pub.dropna(subset=["gain_dn"]).sort_values("n_dn")
    # SIZE-STRATIFIED, because the pool multiple is a multiple of n and its effect cannot be
    # separated from dataset size on a subset chosen by convenience.
    pick = np.linspace(0, len(pub) - 1, want).round().astype(int)
    datasets = list(pub.iloc[pick].dataset)

    rows = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        f = Path(store) / "processed" / "dinuc" / cell / protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        pos_df = d[d.label == 1].reset_index(drop=True)
        positives = [{"chrom": r.chrom, "start": int(r.start), "end": int(r.end),
                      "strand": r.strand, "region": r.region, "gc": r.gc,
                      "seq_dna": r.seq_dna, "seq_rna": r.seq_rna}
                     for r in pos_df.itertuples()]
        try:
            peaks = list(win.read_peaks(peak_path(DATA_ROOT, protein, cell)))
        except (FileNotFoundError, RuntimeError) as e:
            log(f"  {ds}: {e}")
            continue
        # The region pools depend on the peaks and the annotation, not on the pool multiple
        # or the assignment rule, so they are computed once for all six settings. Without
        # this the same subtraction ran six times and was most of the runtime.
        excl = neg.exclusion_zones(peaks, mpd)
        pools = {r: neg.available(index, r, excl, size)
                 for r in {q["region"] for q in positives}}
        for algo in ALGOS:
            for mult, pmin in POOLS:
                t0 = time.time()
                if algo == "greedy":
                    negs, dr, dist = neg.build_negatives_dinuc(
                        positives, peaks, fasta, index, size, min_peak_distance=mpd,
                        seed=cfg.seed, drop_n=True, pool_multiple=mult, pool_min=pmin,
                        pools=pools)
                else:
                    negs, dr, dist = assign_optimal(
                        positives, peaks, fasta, index, size, mpd, cfg.seed, mult,
                        pool_min=pmin, pools=pools)
                s = score(pos_df, negs)
                if s is None:
                    continue
                rows.append({"dataset": ds, "protein": protein, "cell": cell,
                             "algo": algo, "pool_multiple": mult, "pool_min": pmin,
                             "floor_bound_buckets": int(sum(
                                 pmin > mult * len(v) for v in buckets_of(positives).values())),
                             "n_buckets": len(buckets_of(positives)),
                             "mean_l1": float(np.nanmean(dist)),
                             "median_l1": float(np.nanmedian(dist)),
                             "matched": int(np.isfinite(dist).sum()),
                             "subsampled_buckets": int(dr.get("subsampled_buckets", 0)),
                             "seconds": round(time.time() - t0, 1), **s})
                log(f"[{i:3d}/{len(datasets)}] {ds:18s} {algo:8s} pool {mult:2d}x  "
                    f"pairs {s['pairs']:6d}  L1 {np.nanmean(dist):.4f}  comp {s['comp']:.4f}  "
                    f"gain {s['gain']:+.4f}  {round(time.time() - t0, 1)}s")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--index", default=str(DATA_ROOT / "data/interim/regions.pkl"))
    p.add_argument("--fasta",
                   default=str(DATA_ROOT / "data/raw/GRCh38.primary_assembly.genome.fa"))
    p.add_argument("--n", type=int, default=15)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "matching_robustness_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it without --from-cache")
    else:
        t = build(a.store, a.index, a.fasta, a.n)
        t.to_csv(per, index=False)

    out = []
    ds = sorted(t.dataset.unique())
    log(f"\n=== B11: matching-algorithm robustness, {len(ds)} datasets, "
        f"{len(t)} rebuilds ===\n")
    out.append({"check": "datasets rebuilt", "value": len(ds), "n": len(t)})
    out.append({"check": "negative sets rebuilt", "value": len(t), "n": len(t)})
    out.append({"check": "settings per dataset", "value": len(ALGOS) * len(POOLS),
                "n": len(t)})

    # THE PUBLISHED SETTING MUST REPRODUCE. Greedy at 8x IS the committed dinucleotide arm, so
    # its gain has to come back equal to three_arm_per_dataset.csv on the same datasets. This
    # is the anchor: without it, six agreeing numbers would show only that this script is
    # self-consistent, exactly as in the partition and aggregation sections.
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv").set_index("dataset")
    base = t[(t.algo == "greedy") & (t.pool_multiple == 8)].set_index("dataset")
    err = (base.gain - pub.loc[base.index, "gain_dn"]).abs()
    out.append({"check": "max |greedy 8x gain - published dinucleotide gain|",
                "value": float(err.max()), "n": len(base),
                "note": "greedy at 8x IS the published construction"})
    log(f"  greedy at 8x reproduces the published gain to {err.max():.2e} "
        f"over {len(base)} datasets")

    log(f"\n  {'algo':9s} {'pool':>5s} {'mean L1':>9s} {'composition':>12s} "
        f"{'contribution':>13s} {'matched':>8s}")
    cells = {}
    for algo in ALGOS:
        for mult, _pmin in POOLS:
            s = t[(t.algo == algo) & (t.pool_multiple == mult)]
            if s.empty:
                continue
            cells[(algo, mult)] = s
            for name, col in (("mean L1 distance", "mean_l1"),
                              ("composition AUROC", "comp"),
                              ("4-mer contribution", "gain")):
                out.append({"check": f"{name}, {algo} at {mult}x", "value": float(s[col].mean()),
                            "n": len(s)})
            out.append({"check": f"pairs matched, {algo} at {mult}x",
                        "value": int(s.matched.sum()), "n": len(s)})
            # HOW OFTEN THE MULTIPLE IS INERT. At the published setting the floor binds in
            # almost every bucket, which is why 4x, 8x and 16x with a fixed floor return the
            # same negatives. This is the number that makes the parameter's scaling honest.
            out.append({"check": f"fraction of buckets where the pool floor binds, {algo} "
                                 f"at {mult}x",
                        "value": float((s.floor_bound_buckets / s.n_buckets).mean()),
                        "n": len(s)})
            log(f"  {algo:9s} {mult:4d}x {s.mean_l1.mean():9.4f} {s.comp.mean():12.4f} "
                f"{s.gain.mean():+13.4f} {int(s.matched.sum()):8d}")

    # DOES THE EXACT ASSIGNMENT MATCH BETTER? If it does not, the greedy defence is right for
    # a reason nobody had checked; if it does, the interesting question is whether the better
    # match moves the measurement.
    for mult, _pmin in POOLS:
        g, o = cells.get(("greedy", mult)), cells.get(("optimal", mult))
        if g is None or o is None:
            continue
        m = g.set_index("dataset").join(o.set_index("dataset"), rsuffix="_opt", how="inner")
        out.append({"check": f"L1 improvement from exact assignment at {mult}x",
                    "value": float((m.mean_l1 - m.mean_l1_opt).mean()), "n": len(m)})
        out.append({"check": f"contribution change from exact assignment at {mult}x",
                    "value": float((m.gain_opt - m.gain).mean()), "n": len(m)})
        out.append({"check": f"composition change from exact assignment at {mult}x",
                    "value": float((m.comp_opt - m.comp).mean()), "n": len(m)})
        log(f"  exact assignment at {mult}x: L1 "
            f"{(m.mean_l1_opt - m.mean_l1).mean():+.4f}, composition "
            f"{(m.comp_opt - m.comp).mean():+.4f}, contribution "
            f"{(m.gain_opt - m.gain).mean():+.4f}")

    # THE RANGE OF THE HEADLINE OVER THE SIX SETTINGS, which is what a reader needs: how much
    # of the dinucleotide arm's contribution is the protocol and how much is one afternoon's
    # implementation choice.
    means = {k: float(v.gain.mean()) for k, v in cells.items()}
    comps = {k: float(v.comp.mean()) for k, v in cells.items()}
    rng_gain = max(means.values()) - min(means.values())
    out.append({"check": "range of the 4-mer contribution over all matching settings",
                "value": rng_gain, "n": len(t)})
    out.append({"check": "range of the composition baseline over all matching settings",
                "value": max(comps.values()) - min(comps.values()), "n": len(t)})
    out.append({"check": "smallest 4-mer contribution over all matching settings",
                "value": min(means.values()), "n": len(t)})
    log(f"\n  over all {len(means)} settings the contribution ranges {rng_gain:.4f} "
        f"({min(means.values()):+.4f} to {max(means.values()):+.4f}) and the composition "
        f"baseline {max(comps.values()) - min(comps.values()):.4f}")

    # AND THE DIRECTION, which is the paper's thesis applied to its own free parameters: a
    # better match raises the baseline and a higher baseline leaves a smaller contribution.
    # Measured across all six settings as a correlation, so it is one number and not a story
    # told about six.
    from scipy.stats import pearsonr
    xs = np.array([comps[k] for k in cells])
    ys = np.array([means[k] for k in cells])
    if len(xs) > 2:
        r, pv = pearsonr(xs, ys)
        out.append({"check": "correlation across settings between baseline and contribution",
                    "value": float(r), "n": len(xs), "note": f"p = {pv:.3f}"})
        log(f"  across settings, baseline vs contribution r = {r:+.3f} (p = {pv:.3f}): the "
            f"free parameters move the measurement THROUGH the baseline")

    sub = int(t.subsampled_buckets.sum())
    out.append({"check": "buckets whose candidate set was subsampled for memory",
                "value": sub, "n": len(t)})
    if sub:
        log(f"  {sub} buckets had their candidate set subsampled, so 'optimal' there is exact "
            f"for a reduced pool")

    pd.DataFrame(out).to_csv(TABLES / "matching_robustness.csv", index=False)
    log("\nwrote matching_robustness.csv and matching_robustness_per_dataset.csv")


if __name__ == "__main__":
    main()
