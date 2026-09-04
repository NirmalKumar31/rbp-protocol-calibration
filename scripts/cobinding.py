"""Stage 2b: measure how much the panel's RBPs bind the same windows.

    python scripts/cobinding.py --config config/params.yaml
"""
import argparse, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rbp.data import cobinding  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--peaks", default=str(ROOT / "data/raw/peaks"))
    p.add_argument("--outdir", default=str(ROOT / "results/tables"))
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    size = cfg.windows["size"]

    paths = cobinding.peak_paths_from(a.peaks)
    print(f"{len(paths)} proteins, {size}-nt windows\n")
    m, stats = cobinding.matrix(paths, window=size)

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    m.to_csv(out / "cobinding_matrix.csv")
    stats.to_csv(out / "cobinding_stats.csv")

    print("=== per protein: how much of its binding is shared with any other ===")
    s = stats.sort_values("frac_shared_with_any_other", ascending=False)
    for prot, r in s.iterrows():
        bar = "#" * int(round(r.frac_shared_with_any_other * 40))
        print(f"  {prot:9} {int(r.n_windows):6d} windows  "
              f"{r.frac_shared_with_any_other:5.1%} shared  {bar}")

    vals = m.to_numpy(dtype=float).copy()
    np.fill_diagonal(vals, np.nan)
    off = vals[~np.isnan(vals)]
    print(f"\nmean pairwise overlap (off-diagonal): {off.mean():.1%}")
    print(f"median: {np.median(off):.1%}   max: {off.max():.1%}")
    pairs = [(m.index[i], m.columns[j], vals[i, j])
             for i in range(len(m)) for j in range(len(m)) if not np.isnan(vals[i, j])]
    pairs.sort(key=lambda t: -t[2])
    print("\ntop 8 directed pairs (fraction of A's windows hitting B):")
    for a, b, v in pairs[:8]:
        print(f"  {a:9} -> {b:9} {v:5.1%}")
    print(f"mean shared-with-any-other: {stats.frac_shared_with_any_other.mean():.1%}")
    print(f"\nwrote {(out/'cobinding_matrix.csv').relative_to(ROOT)} and cobinding_stats.csv")


if __name__ == "__main__":
    main()
