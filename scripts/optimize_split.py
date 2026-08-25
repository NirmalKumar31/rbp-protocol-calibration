"""Choose the chromosome -> split assignment that fits every protein best.

Optimises on peak counts only (never labels or performance), so it is stratification
rather than leakage. Prints the assignment to paste into config/params.yaml.

    python scripts/optimize_split.py
"""
import argparse, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rbp.data import annotation as ann, splits  # noqa: E402
from rbp.data.cobinding import peak_paths_from  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def report(name, names, counts, assign, target):
    loss, worst = splits.assignment_loss(counts, assign, target)
    totals = counts.sum(axis=1, keepdims=True); totals[totals == 0] = 1
    props = np.hstack([counts[:, assign == k].sum(axis=1, keepdims=True) / totals
                       for k in range(3)])
    print(f"\n=== {name}  (loss {loss:.4f}, worst deviation {worst:.3f}) ===")
    print(f"{'protein':9} {'train':>7} {'val':>7} {'test':>7}")
    for i, p in enumerate(names):
        print(f"{p:9} {props[i,0]:7.2f} {props[i,1]:7.2f} {props[i,2]:7.2f}")
    print(f"{'range':9} " + " ".join(
        f"{props[:,k].min():.2f}-{props[:,k].max():.2f}".rjust(7) for k in range(3)))
    return props


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--restarts", type=int, default=60)
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    target = (0.64, 0.16, 0.20)
    drop = set(cfg.encode.get("exclude_chroms", []))
    chroms = [c for c in ann.MAIN_CHROMS if c not in drop]
    if drop:
        print(f"excluding {sorted(drop)}")

    paths = {r["protein"]: next((ROOT / "data/raw/peaks").glob(f"{r['protein']}.*.bed.gz"))
             for r in cfgmod.proteins()}
    names, counts = splits.peak_counts(paths, chroms)
    print(f"{len(names)} proteins x {len(chroms)} chromosomes")

    cur = np.zeros(len(chroms), dtype=int)
    for j, c in enumerate(chroms):
        cur[j] = 2 if c in cfg.split["test"] else 1 if c in cfg.split["val"] else 0
    report("CURRENT assignment", names, counts, cur, target)

    loss, worst, best = splits.optimize_assignment(
        counts, target=target, restarts=a.restarts, seed=cfg.seed)
    report("OPTIMISED assignment", names, counts, best, target)

    out = {k: [chroms[j] for j in range(len(chroms)) if best[j] == i]
           for i, k in enumerate(("train", "val", "test"))}
    print("\npaste into config/params.yaml under `split:`")
    for k in ("test", "val"):
        print(f"  {k}:   [{', '.join(out[k])}]")
    print(f"  # train (implicit): {', '.join(out['train'])}")


if __name__ == "__main__":
    main()
