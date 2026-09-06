"""R1k: a third negative-set arm, built the way the field's own benchmark builds one.

    python scripts/build_neg2.py --store ../rbp-store

WHY A THIRD ARM. Both existing arms draw negatives from unbound genomic windows matched on
composition, and that shared design carries every objection a referee raised: 40% of the
negatives sit in loci no transcript is produced from (R1j), roughly half are antisense to any
transcript (R1c), and how well the matcher does its job is set by undocumented sampler
parameters. A two-point contrast between two variants of one flawed design can only say so
much.

Horlacher et al. 2023 (Brief Bioinform 24(5):bbad307) define the bias-aware alternative, and
this arm implements their negative-2 verbatim: sample negatives from *the binding sites of
other RBPs assayed in the same cell line*, excluding any that overlap the target protein's own
sites, 1:1. Their window length is 101 nt, identical to ours.

WHAT IT BUYS, and it is four objections at once. Such negatives are

  * 100% transcribed, because they are somebody's crosslink site -- R1j cannot apply;
  * 100% CLIP-accessible, so "not sequenced" cannot masquerade as "not bound";
  * strand-correct by construction, each inheriting the strand of the peak it came from,
    which is its own true strand -- R1c cannot apply;
  * free of the composition matcher entirely, so the contrast cannot be a readout of
    `pool_multiple` or of greedy assignment.

WHAT IT COSTS. The negatives are no longer composition-matched at all, so this arm is not a
third point on the same axis -- it is a different axis. It answers "what happens to measured
model value under the protocol the field's own benchmark recommends", which is the question
Horlacher's -0.065 to -0.085 AUROC drop calibrates against our -0.1095.

TWO DESIGN CHOICES WORTH ARGUING WITH.

  1. FOLD IS MATCHED, region and composition are NOT. Our CV assigns folds by chromosome and
     every existing positive/negative pair lands in the same fold; the leakage audit checks
     that and finds 0 violations across 3.59M rows. Sampling a negative from a different fold
     would break that invariant for a reason unrelated to the science, so the pool is
     restricted to the positive's own fold. Region and composition are deliberately left
     free -- matching them is precisely what this arm exists to avoid.
  2. EXCLUSION IS BY DISTANCE, not by overlap. Horlacher says "do not overlap"; we use the
     project's own `negatives.min_peak_distance` (500 nt), which is stricter and consistent
     with how the other two arms define "unbound".
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils.log import log  # noqa: E402

CELLS = ("K562", "HepG2")



def load_positives(gc_root, cell, proteins):
    """Every panel protein's positive windows in one cell line, tagged with their source."""
    frames = []
    for p in proteins:
        f = gc_root / cell / p / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        d = d[d.label == 1].copy()
        d["source"] = p
        frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else None


def too_close(pool, target, margin):
    """Mask of pool rows within `margin` of any target window on the same chromosome."""
    bad = np.zeros(len(pool), dtype=bool)
    for chrom, tg in target.groupby("chrom"):
        idx = np.flatnonzero(pool.chrom.values == chrom)
        if not len(idx):
            continue
        starts = np.sort(tg.start.values)
        ends = np.sort(tg.end.values)
        ps, pe = pool.start.values[idx], pool.end.values[idx]
        # nearest target start at or before each pool end, and nearest end at or after start
        i = np.searchsorted(starts, pe + margin, side="right")
        j = np.searchsorted(ends, ps - margin, side="left")
        bad[idx] = i > j            # some target interval falls inside the padded window
    return bad


def draw(cand, grp, rng, match_region):
    """Draw one negative per positive from `cand`, optionally within region class.

    Without --match-region this is the arm as published: a uniform draw inside the fold.
    With it, the draw is stratified so the negatives reproduce the positives' region
    marginals, which is what the two composition-matched matchers do by construction.
    """
    if not match_region:
        n = min(len(grp), len(cand))
        return ([cand.iloc[rng.choice(len(cand), n, replace=False)]] if n else []), len(grp) - n
    out, short = [], 0
    for region, sub in grp.groupby("region"):
        pool_r = cand[cand.region.values == region]
        n = min(len(sub), len(pool_r))
        short += len(sub) - n
        if n:
            out.append(pool_r.iloc[rng.choice(len(pool_r), n, replace=False)])
    return out, short


def build(store, seed=7, match_region=False, arm=None):
    """Build the bias-aware arm. `arm` overrides the output directory name.

    THE OVERRIDE EXISTS SO A REDRAW CANNOT CLOBBER THE PUBLISHED ARM. This wrote
    unconditionally to processed/neg2, so running it with a different seed to measure
    draw-to-draw variability would have destroyed the negatives every published bias-aware
    number was computed from, in place, with no copy. scripts/negative_draws.py passes
    arm="neg2_seed11" and so on.
    """
    cfg = cfgmod.load()
    margin = int(cfg["negatives"]["min_peak_distance"])
    gc_root = Path(store) / "processed" / "gc"
    out_root = Path(store) / "processed" / (arm or ("neg2_rm" if match_region else "neg2"))
    panel = pd.read_csv(ROOT / "results" / "tables" / "rehearsal_binding_gc.csv")
    rng = np.random.default_rng(seed)

    made, report = 0, []
    for cell in CELLS:
        proteins = sorted(panel[panel.cell == cell].protein)
        allpos = load_positives(gc_root, cell, proteins)
        if allpos is None:
            continue
        log(f"{cell}: {len(proteins)} panel proteins, {len(allpos):,} pooled positive windows")
        for target in proteins:
            tgt = allpos[allpos.source == target]
            pool = allpos[allpos.source != target].copy()
            pool = pool[~too_close(pool, tgt, margin)]
            if pool.empty:
                continue

            picked = []
            short = 0
            for fold, grp in tgt.groupby("fold"):
                got, miss = draw(pool[pool.fold == fold], grp, rng, match_region)
                picked += got
                short += miss
            if not picked:
                continue
            neg = pd.concat(picked, ignore_index=True)
            neg["label"] = 0
            neg["id"] = [f"{target}_n2_{i}" for i in range(len(neg))]
            # Keep the positives that actually got a partner, 1:1 within the strata the draw
            # used -- fold alone, or fold and region when the draw was stratified.
            keys = ["fold", "region"] if match_region else ["fold"]
            # Tuple-normalise both sides: a one-element groupby list yields scalar keys on
            # some pandas versions and 1-tuples on others, and a silent mismatch here drops
            # every dataset while reporting success.
            def _k(x):
                return x if isinstance(x, tuple) else (x,)
            counts = {_k(k): v for k, v in neg.groupby(keys).size().items()}
            keep = []
            for key, grp in tgt.groupby(keys):
                k = int(counts.get(_k(key), 0))
                if k:
                    keep.append(grp.iloc[:k])
            if not keep:
                continue
            pos = pd.concat(keep, ignore_index=True)

            ds = pd.concat([pos, neg], ignore_index=True).drop(columns=["source"])
            d = out_root / cell / target
            d.mkdir(parents=True, exist_ok=True)
            ds.to_csv(d / "dataset.tsv", sep="\t", index=False)
            made += 1
            report.append({"dataset": f"{target}:{cell}", "protein": target, "cell": cell,
                           "pairs": len(pos), "pool_after_exclusion": len(pool),
                           "unmatched_positives": short,
                           "donor_proteins": int(neg.source.nunique())
                           if "source" in neg else np.nan})
            if made % 20 == 0:
                log(f"  {made} datasets built")
    r = pd.DataFrame(report)
    name = "neg2_rm_build.csv" if match_region else "neg2_build.csv"
    r.to_csv(ROOT / "results" / "tables" / name, index=False)
    log(f"\nbuilt {made} datasets -> {out_root}")
    log(f"  median pairs {int(r.pairs.median()):,}, total {int(r.pairs.sum()):,}")
    dropped = int(r.unmatched_positives.sum())
    pct = 100 * dropped / (r.pairs.sum() + dropped)
    log(f"  positives dropped for want of a same-fold donor: {dropped:,} ({pct:.2f}%)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--arm", default=None,
                   help="output directory name under processed/. Use this for a redraw so the "
                        "published arm is not overwritten.")
    p.add_argument("--match-region", action="store_true",
                   help="stratify the draw on transcript region, as a DIAGNOSTIC arm. This is "
                        "not Horlacher's protocol, which leaves region free; it exists to "
                        "measure how much of this arm's baseline is region mix.")
    a = p.parse_args()
    build(a.store, a.seed, a.match_region, a.arm)


if __name__ == "__main__":
    main()
