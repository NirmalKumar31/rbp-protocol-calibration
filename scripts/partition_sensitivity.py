"""B6: how much of the headline contrast depends on the one chromosome partition we froze?

    python scripts/partition_sensitivity.py --store ../rbp-store --partitions 3

WHY. Every number in this paper uses `config/folds.tsv`, one balanced chromosome-to-fold
assignment chosen once. The fold-to-fold spreads reported throughout therefore describe
variability WITHIN that partition, and say nothing about variability ACROSS the choice of
partition -- which the Limitations had to state as an untested design parameter. This measures
it.

FOUR-MER ONLY, and that is not a shortcut. The CNN and SpliceBERT scores were produced by
sweeps that trained on the frozen partition, so re-partitioning would compare their old folds
against new ones. The 4-mer is refit here on whatever partition is supplied, which makes it the
only model that can answer the question at all.

THE ALTERNATIVES MEET THE SAME CRITERIA, not merely different seeds of a shuffle: the same
balanced-mass objective over the same chromosomes, with the same minimum chromosomes per fold,
optimised from different starts. A partition that failed the criteria would be a different
design rather than a re-draw of this one.

MASS IS POSITIVE-WINDOW COUNTS, taken from the committed window tables, rather than the peak
counts `optimize_folds.py` used. Those windows are what the folds actually have to balance, the
peak BEDs are not in this repository, and the two differ only in that one peak yields one
window.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data import splits  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "gc", "dn": "dinuc"}


def log(m):
    print(m, flush=True)


def frozen():
    d = pd.read_csv(ROOT / "config" / "folds.tsv", sep="\t")
    return dict(zip(d.chrom, d.fold, strict=True))


def mass_matrix(store, datasets):
    """(dataset x chromosome) positive-window counts, from the GC arm's tables."""
    root = Path(store) / "processed" / "gc"
    chroms, rows, names = set(), [], []
    per = {}
    for ds in datasets:
        protein, cell = ds.split(":")
        f = root / cell / protein / "dataset.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t", usecols=["chrom", "label"])
        c = d[d.label == 1].chrom.value_counts()
        per[ds] = c
        chroms |= set(c.index)
        names.append(ds)
    chroms = sorted(chroms)
    for ds in names:
        rows.append([int(per[ds].get(c, 0)) for c in chroms])
    return names, chroms, np.asarray(rows, dtype=np.int64)


def contrast(store, datasets, assign, label):
    """Panel-mean 4-mer contrast (dinuc minus GC) under one chromosome assignment."""
    gains = {"gc": [], "dn": []}
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        got = {}
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                break
            d = pd.read_csv(f, sep="\t", usecols=["seq_rna", "label", "chrom"])
            fold = d.chrom.map(assign)
            m = fold.notna().to_numpy()
            if m.sum() < 100 or d.label.values[m].std() == 0:
                break
            fold = fold[m].astype(int).to_numpy()
            if len(np.unique(fold)) < 2:
                break
            seqs, y = d.seq_rna.values[m], d.label.values[m]
            sc, _, _ = kmer_oof(seqs, y, fold, k=4)
            ok = np.isfinite(sc)
            g = gain_over_composition(seqs[ok], sc[ok], y[ok], fold[ok])
            got[arm] = g.delta
        if len(got) == 2:
            gains["gc"].append(got["gc"])
            gains["dn"].append(got["dn"])
        if i % 25 == 0:
            log(f"    [{label}] {i}/{len(datasets)}")
    n = len(gains["gc"])
    if not n:
        return None
    gc, dn = float(np.mean(gains["gc"])), float(np.mean(gains["dn"]))
    return {"partition": label, "n": n, "gain_gc": gc, "gain_dn": dn,
            "contrast": dn - gc, "multiplier": dn / gc if gc else np.nan}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--partitions", type=int, default=3)
    p.add_argument("--k", type=int, default=5)
    p.add_argument("--min-per-fold", type=int, default=3)
    a = p.parse_args()

    three = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    datasets = sorted(three.dataset)
    names, chroms, counts = mass_matrix(a.store, datasets)
    log(f"=== B6: partition sensitivity, {len(names)} datasets x {len(chroms)} chromosomes ===\n")

    froz = frozen()
    assigns = [("frozen (config/folds.tsv)", froz)]
    for s in range(1, a.partitions + 1):
        _loss, _worst, best = splits.optimize_folds(
            counts, k=a.k, min_per_fold=a.min_per_fold, restarts=40, iters=4000,
            seed=100 + s)
        cand = dict(zip(chroms, [int(x) for x in best], strict=True))
        # A re-draw only counts if it is actually a different partition.
        if all(cand.get(c) == froz.get(c) for c in chroms):
            log(f"  seed {100 + s} reproduced the frozen partition; skipped")
            continue
        sizes = [sum(1 for c in chroms if cand[c] == f) for f in range(a.k)]
        assert min(sizes) >= a.min_per_fold, sizes
        assigns.append((f"alternative {s}", cand))
        log(f"  alternative {s}: chromosomes per fold {sizes}")

    rows = []
    for label, assign in assigns:
        log(f"\n  computing under {label} ...")
        r = contrast(a.store, datasets, assign, label)
        if r:
            rows.append(r)
            log(f"    gc {r['gain_gc']:+.4f}  dn {r['gain_dn']:+.4f}  "
                f"contrast {r['contrast']:+.4f}  multiplier {r['multiplier']:.3f}x")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing computed")
    t.to_csv(TABLES / "partition_sensitivity_per_partition.csv", index=False)

    out = []
    for _, r in t.iterrows():
        out.append({"check": f"4-mer contrast under {r.partition}", "value": float(r.contrast),
                    "n": int(r.n)})
    out += [{"check": "partitions compared", "value": len(t), "n": len(t)},
            {"check": "minimum 4-mer contrast over partitions",
             "value": float(t.contrast.min()), "n": len(t)},
            {"check": "maximum 4-mer contrast over partitions",
             "value": float(t.contrast.max()), "n": len(t)},
            {"check": "range of the 4-mer contrast over partitions",
             "value": float(t.contrast.max() - t.contrast.min()), "n": len(t)},
            {"check": "contrast is positive under every partition",
             "value": int((t.contrast > 0).all()), "n": len(t)}]
    pd.DataFrame(out).to_csv(TABLES / "partition_sensitivity.csv", index=False)
    log(f"\n  contrast over {len(t)} partitions: {t.contrast.min():+.4f} to "
        f"{t.contrast.max():+.4f}, range {t.contrast.max() - t.contrast.min():.4f}")
    log("\nwrote partition_sensitivity.csv and partition_sensitivity_per_partition.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
