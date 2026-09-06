"""Step 3: run the entire analysis end to end with k-mer scores, before any GPU.

This is the dress rehearsal. Every downstream stage -- out-of-fold scoring, the
composition control, DeLong ranking, the variant delta, the conservation control, the
transfer correlation -- runs here on a baseline that costs nothing. If the pipeline is
broken it breaks now for $0 rather than halfway through a paid sweep, and the paper's
tables and figures exist in draft before a model is trained.

It also produces a real result rather than only a smoke test: the composition control
across the whole panel needs no GPU at all, and how much of a reported binding AUROC is
nucleotide composition is a finding in its own right.

    python scripts/rehearsal.py --what binding      # per-dataset AUROC + composition
    python scripts/rehearsal.py --what ranking      # DeLong across k
    python scripts/rehearsal.py --what all
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.eval import baseline, delong, nested  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
CELLS = ("K562", "HepG2")


def datasets(cells=CELLS, min_pairs=None, arm="gc"):
    """Every prepared dataset that clears the threshold, keyed by 'PROTEIN:CELL'.

    The arm selects BOTH the panel file and the directory. This used to read
    panel_final_<cell>.tsv and load from data/processed/ unconditionally, which meant it
    filtered GC datasets using whichever arm's pair counts had been written last.
    """
    out, missing = panelmod.datasets(cells, arm, min_pairs)
    for cell in cells:
        if not panelmod.panel_path(cell, arm).exists():
            print(f"  no {panelmod.panel_path(cell, arm).name}; "
                  f"run prepare.py --match {arm} for {cell} first")
    if missing:
        print(f"  {len(missing)} in the panel but missing on disk: {missing[:5]}")
    return out


def load(path):
    return pd.read_csv(path, sep="\t")


def binding(paths, k, n_boot):
    """Per-dataset out-of-fold AUROC, the composition baseline, and the gain."""
    rows = []
    t0 = time.time()
    for i, (name, path) in enumerate(sorted(paths.items()), 1):
        df = load(path)
        y, folds, seqs = df.label.to_numpy(), df.fold.to_numpy(), df.seq_rna.tolist()
        res = baseline.evaluate(df, k=k)
        g = nested.gain_over_composition(seqs, res["scores"], y, folds)
        nt = nested.test_score(seqs, res["scores"], y, n_boot=n_boot, seed=7)
        rows.append({
            "dataset": name, "protein": name.split(":")[0], "cell": name.split(":")[1],
            "pairs": int(len(df) // 2), "n": res["n"],
            "auroc": res["auroc"], "ci_low": res["ci_low"], "ci_high": res["ci_high"],
            "composition_auroc": g.auroc_composition,
            "with_score_auroc": g.auroc_with_score,
            "delta_auroc": g.delta, "delta_ci_low": g.delta_ci_low,
            "delta_ci_high": g.delta_ci_high, "delta_p": g.delta_p, "helps": g.helps,
            "coef": nt.coef, "coef_ci_low": nt.ci_low, "coef_ci_high": nt.ci_high,
            "lr_p": nt.lr_p})
        if i % 10 == 0 or i == len(paths):
            el = time.time() - t0
            print(f"  [{i:3d}/{len(paths)}] {el:5.0f}s elapsed, "
                  f"~{el/i*(len(paths)-i):4.0f}s left", flush=True)
    return pd.DataFrame(rows)


def ranking(paths, ks):
    """DeLong comparisons between k values, standing in for architecture comparison.

    The mechanics are identical to comparing SpliceBERT against RNA-FM: several scores on
    the same test rows, correlated because they share the data. Exercising it on k values
    proves the machinery before there are architectures to compare.
    """
    rows = []
    for name, path in sorted(paths.items()):
        df = load(path)
        y, folds = df.label.to_numpy(), df.fold.to_numpy()
        scores = {}
        for k in ks:
            s, _, _ = baseline.oof_scores(df.seq_rna.tolist(), y, folds, k=k)
            scores[f"k{k}"] = s
        ok = np.all(np.isfinite(np.vstack(list(scores.values()))), axis=0)
        cmp = delong.pairwise({m: v[ok] for m, v in scores.items()}, y[ok])
        cmp.insert(0, "dataset", name)
        rows.append(cmp)
    return pd.concat(rows, ignore_index=True)


def summarise(res):
    print(f"\n{'':=<74}")
    print(f"{len(res)} datasets, {res.protein.nunique()} distinct proteins, "
          f"{res.pairs.sum():,} pairs")
    print(f"{'':=<74}\n")
    print(f"{'':22} {'median':>8} {'mean':>8} {'min':>8} {'max':>8}")
    for col, lab in (("auroc", "baseline AUROC"), ("composition_auroc", "composition alone"),
                     ("delta_auroc", "gain over composition")):
        c = res[col]
        print(f"{lab:22} {c.median():8.3f} {c.mean():8.3f} {c.min():8.3f} {c.max():8.3f}")

    print("\ngain over composition, distribution across datasets:")
    for lo, hi, lab in ((-1, 0.005, "<0.005 (composition explains it)"),
                        (0.005, 0.02, "0.005-0.02"), (0.02, 0.05, "0.02-0.05"),
                        (0.05, 1, ">0.05 (substantial)")):
        n = int(((res.delta_auroc >= lo) & (res.delta_auroc < hi)).sum())
        print(f"  {lab:36} {n:4d} ({100*n/len(res):4.1f}%)")
    print(f"\n  interval excludes zero: {int(res.helps.sum())}/{len(res)}")
    print(f"  lr_p < 0.05:            {int((res.lr_p < 0.05).sum())}/{len(res)}"
          f"   <- significance is not the finding; the gain is")

    print("\nthe 8 datasets composition explains most completely:")
    print(res.nsmallest(8, "delta_auroc")[
        ["dataset", "pairs", "composition_auroc", "auroc", "delta_auroc"]
    ].to_string(index=False))
    print("\nthe 8 where the model adds most:")
    print(res.nlargest(8, "delta_auroc")[
        ["dataset", "pairs", "composition_auroc", "auroc", "delta_auroc"]
    ].to_string(index=False))

    # The confound flagged during step 1: dataset size predicts AUROC, so it could drive
    # both sides of the transfer correlation. Measured here so it is not a surprise later.
    lp = np.log(res.pairs)
    print("\nconfound check, correlation with log(pairs):")
    for c in ("auroc", "composition_auroc", "delta_auroc"):
        print(f"  {c:20} r = {np.corrcoef(lp, res[c])[0,1]:+.3f}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--what", default="binding",
                   choices=["binding", "ranking", "all"])
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--ks", default="3,4,5,6")
    p.add_argument("--n-boot", type=int, default=200,
                   help="bootstrap resamples; 200 for the rehearsal, 2000 for the paper")
    p.add_argument("--limit", type=int, default=None, help="first N datasets, for a trial")
    p.add_argument("--arm", default="dinuc", choices=sorted(panelmod.ARMS),
                   help="negative arm; dinuc is the primary one")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    TABLES.mkdir(parents=True, exist_ok=True)

    paths = datasets(min_pairs=cfg.cv["min_pairs"], arm=a.arm)
    if a.limit:
        paths = dict(sorted(paths.items())[:a.limit])
    if not paths:
        raise SystemExit("no prepared datasets found")
    print(f"{len(paths)} datasets, k={a.k}, {a.n_boot} bootstrap resamples")

    # Output is keyed by arm for the same reason the panel file is: two arms writing to
    # one filename means whichever ran last is the one you read, silently.
    if a.what in ("binding", "all"):
        res = binding(paths, a.k, a.n_boot)
        out = TABLES / f"rehearsal_binding_{a.arm}.csv"
        res.to_csv(out, index=False)
        summarise(res)
        print(f"\nwrote {out.relative_to(ROOT)}")

    if a.what in ("ranking", "all"):
        ks = [int(x) for x in a.ks.split(",")]
        cmp = ranking(paths, ks)
        out = TABLES / f"rehearsal_ranking_{a.arm}.csv"
        cmp.to_csv(out, index=False)
        print(f"\nDeLong: {len(cmp)} comparisons, "
              f"{int((cmp.q < 0.05).sum())} significant after FDR")
        print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
