"""How much of the reported interval is one lucky draw of the negatives?

    python scripts/negative_draws.py --store ../rbp-store
    python scripts/negative_draws.py --from-cache

THE GAP THIS CLOSES, PARTLY. Every published protocol used ONE stochastic draw of its negatives,
and the protein-clustered bootstrap resamples proteins with that draw held fixed. So the stated
intervals do not contain negative-construction uncertainty, and three audits have said so. The
Limitations quantify a single redraw of the dinucleotide arm -- median absolute per-dataset
deviation 0.0031, maximum 0.0079 -- but a quantity that is measured and not propagated is still
missing from the interval a reader acts on.

WHY THE BIAS-AWARE ARM AND NOT ALL THREE. Redrawing a composition-matched arm means generating
new candidate windows from the genome, and the genome FASTA is a 3 GB download that is not part
of the released evidence. The bias-aware arm needs no such thing: its negatives ARE other
proteins' binding sites, every one of which is already in the window store. So this arm can be
redrawn exactly -- same procedure, same pool, different seed -- at no cost beyond CPU time, and
the other two cannot. That is a real limit on what follows and it is stated rather than papered
over: what is measured here is draw variability for one of three protocols.

The arm chosen is also the informative one for this question. It has the smallest contribution
of the three, so a draw effect of a given absolute size matters most there, and the Results
already say that arm's absolute level is not distinguishable from the estimator's floor.

WHAT IS REPORTED. Five independent draws including the published one. For each: the panel-mean
contribution, its protein-clustered interval, and then the two variance components separated --
between-protein, which the published bootstrap already carries, and between-draw, which it does
not. The combined interval adds them.
"""

import argparse
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from build_neg2 import build as build_neg2  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
PUBLISHED_SEED = 7
# Four extra draws. Prespecified here rather than chosen after looking: 4 gives a usable
# between-draw variance without a rebuild per seed dominating the run.
EXTRA_SEEDS = (11, 23, 42, 101)


def contribution(root, dataset):
    f = root / dataset.split(":")[1] / dataset.split(":")[0] / "dataset.tsv"
    if not f.exists():
        return None
    d = pd.read_csv(f, sep="\t")
    if d.label.nunique() < 2 or len(d) < 200:
        return None
    sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
    return float(gain_over_composition(d.seq_rna.values, sc, d.label.values,
                                       d.fold.values).delta)


def build(store, limit, keep):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    datasets = list(pub.dataset)[:limit or None]
    frames = {}
    for seed in (PUBLISHED_SEED,) + EXTRA_SEEDS:
        arm = "neg2" if seed == PUBLISHED_SEED else f"neg2_seed{seed}"
        root = Path(store) / "processed" / arm
        if seed != PUBLISHED_SEED and not root.exists():
            log(f"\ndrawing seed {seed} into processed/{arm}")
            build_neg2(store, seed=seed, match_region=False, arm=arm)
        log(f"\nscoring seed {seed} ({arm})")
        vals = {}
        for i, ds in enumerate(datasets, 1):
            v = contribution(root, ds)
            if v is not None:
                vals[ds] = v
            if i % 25 == 0:
                log(f"  [{i:3d}/{len(datasets)}]")
        frames[seed] = vals
        if seed != PUBLISHED_SEED and not keep:
            # 94 datasets of windows per seed is gigabytes. The contributions are what matter.
            shutil.rmtree(root, ignore_errors=True)

    common = set.intersection(*(set(v) for v in frames.values()))
    rows = []
    for ds in sorted(common):
        protein, cell = ds.split(":")
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        for seed in frames:
            rec[f"gain_seed{seed}"] = frames[seed][ds]
        rows.append(rec)
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--keep-draws", action="store_true",
                   help="do not delete the rebuilt window directories afterwards")
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "negative_draws_per_dataset.csv"
    t = pd.read_csv(per) if a.from_cache else build(a.store, a.n, a.keep_draws)
    if not a.from_cache:
        if a.n:
            per = per.with_suffix(".partial.csv")
            log(f"  --n {a.n} given: writing {per.name}, not the committed table")
        t.to_csv(per, index=False)

    seeds = [int(c.replace("gain_seed", "")) for c in t.columns if c.startswith("gain_seed")]
    M = t[[f"gain_seed{s}" for s in seeds]].to_numpy(float)     # datasets x draws

    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(4000)]
    out = []

    def add(check, value, lo="", hi="", note=""):
        out.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                    "n": len(t), "note": note})

    add("datasets", len(t))
    add("independent draws", len(seeds))
    per_draw = M.mean(axis=0)
    for s, m in zip(seeds, per_draw):
        add(f"panel-mean contribution, seed {s}", m,
            note="the published draw" if s == PUBLISHED_SEED else "")

    # THE TWO COMPONENTS. Between-protein is what the published bootstrap resamples; the
    # published interval is its half-width. Between-draw is what it holds fixed.
    boot = np.array([M[i, 0].mean() for i in draws])
    se_protein = float(boot.std(ddof=1))
    se_draw = float(per_draw.std(ddof=1))
    add("between-protein standard error, the published draw", se_protein,
        note="what the protein-clustered bootstrap already carries")
    add("between-draw standard error of the panel mean", se_draw,
        note="what it holds fixed, and therefore omits")
    combined = float(np.hypot(se_protein, se_draw))
    add("combined standard error", combined, note="the two added in quadrature")
    add("ratio of combined to published interval width", combined / se_protein if se_protein else 0,
        note="how much wider an interval including draw uncertainty would be")

    # Per-dataset spread, which is the number the Limitations already quote for one redraw.
    spread = M.max(axis=1) - M.min(axis=1)
    add("median per-dataset range across draws", float(np.median(spread)))
    add("maximum per-dataset range across draws", float(spread.max()))
    # The ordering claim across draws is NOT reported here, deliberately. It needs the other
    # two arms' contributions under the same redraw, and those arms cannot be redrawn without
    # the genome. Emitting len(seeds) as though it were that count, which an earlier version of
    # this file did, is reporting a constant as a measurement.
    add("draws whose panel mean stays below the GC arm's published 0.0265",
        float((per_draw < 0.0265).sum()),
        note="the weaker ordering statement the available redraws can support")

    r = pd.DataFrame(out)
    summary = TABLES / ("negative_draws.partial.csv" if a.n else "negative_draws.csv")
    r.to_csv(summary, index=False)
    log("")
    for _, x in r.iterrows():
        log(f"  {x['check']:62s} {x['value']:+.5f}")
    log(f"\nwrote {summary.name} and {per.name}")


if __name__ == "__main__":
    main()
