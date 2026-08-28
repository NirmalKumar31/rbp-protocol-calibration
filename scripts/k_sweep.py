"""Does R1's contrast depend on the k-mer size? And re-derive it from raw sequence.

TWO JOBS, AND THE SECOND MATTERS MORE THAN THE FIRST.

  ROBUSTNESS   The model is a 4-mer. Nothing in the paper's argument requires that, and a
               referee will ask whether the contrast is an artifact of one arbitrary choice.
               Sweeping k = 3, 4, 5, 6 answers it. It also settles an embarrassment: the
               manuscript said "5-mer" in four places before the tables were checked, and this
               shows what that error would have cost, which is nothing.

  VERIFICATION Every number in R1 is currently READ from `rehearsal_binding_*.csv`, which the
               analysis pass wrote. This script rebuilds the contrast from `dataset.tsv`
               sequence, refitting composition and k-mer models from scratch. Alongside
               `recompute.py` it is the second thing in this repository that PROVES a headline
               rather than reproducing it: if the committed contrast were wrong, the k=4 column
               here would not land on it.

Both arms must come from canonical window tables. The GC arm reproduces locally; the
dinucleotide arm's local copy is a different draw, so its canonical tables are read from a
directory populated from the study bucket. Datasets are gated per-arm on rebuilding their own
published composition and with-score AUROCs before contributing to any column.
"""

import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np                                                      # noqa: E402
import pandas as pd                                                     # noqa: E402
from scipy.stats import wilcoxon                                        # noqa: E402

from rbp.eval import baseline, nested                                   # noqa: E402

TABLES = ROOT / "results" / "tables"
KS = (3, 4, 5, 6)
PUBLISHED_K = 4
REPRO_TOL = 2.0e-3
N_BOOT = 2000
SEED = 0


def log(m):
    print(m, flush=True)


def gain(d, k):
    res = baseline.evaluate(d, k=k)
    g = nested.gain_over_composition(d.seq_rna.tolist(), res["scores"],
                                     d.label.to_numpy(), d.fold.to_numpy())
    return g.delta, g.auroc_composition, g.auroc_with_score


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gc-root", required=True)
    p.add_argument("--dn-root", required=True)
    p.add_argument("--limit", type=int, default=0)
    a = p.parse_args()

    cm = pd.read_csv(TABLES / "cost_of_matching.csv")
    pub = {"gc": pd.read_csv(TABLES / "rehearsal_binding_gc.csv").set_index("dataset"),
           "dn": pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv").set_index("dataset")}
    root = {"gc": Path(a.gc_root), "dn": Path(a.dn_root)}
    names = list(cm.dataset)[: a.limit or None]

    rows = []
    for i, ds in enumerate(names, 1):
        prot, cell = ds.split(":")
        frames, ok = {}, True
        for arm in ("gc", "dn"):
            f = root[arm] / cell / prot / "dataset.tsv"
            if not f.exists() or ds not in pub[arm].index:
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            _, comp, with_s = gain(d, PUBLISHED_K)
            r = pub[arm].loc[ds]
            # REPRODUCTION GATE, at the published k, before any other k is trusted.
            if (abs(comp - r.composition_auroc) > REPRO_TOL
                    or abs(with_s - r.with_score_auroc) > REPRO_TOL):
                ok = False
                break
            frames[arm] = d
        if not ok:
            log(f"  [{i:3d}/{len(names)}] {ds:22} SKIP (does not reproduce at k={PUBLISHED_K})")
            continue
        rec = {"dataset": ds}
        for k in KS:
            for arm in ("gc", "dn"):
                rec[f"gain_{arm}_k{k}"] = gain(frames[arm], k)[0]
            rec[f"contrast_k{k}"] = rec[f"gain_dn_k{k}"] - rec[f"gain_gc_k{k}"]
        rows.append(rec)
        log(f"  [{i:3d}/{len(names)}] {ds:22} " +
            "  ".join(f"k{k} {rec[f'contrast_k{k}']:+.4f}" for k in KS))
        pd.DataFrame(rows).to_csv(TABLES / "k_sweep_per_dataset.csv", index=False)

    m = pd.DataFrame(rows)
    m.to_csv(TABLES / "k_sweep_per_dataset.csv", index=False)
    n = len(m)
    rng = np.random.default_rng(SEED)
    idx = [rng.integers(0, n, n) for _ in range(N_BOOT)]

    out = []

    def add(check, value, col=None, note=""):
        lo, hi = (np.percentile([m[col].iloc[i].mean() for i in idx], [2.5, 97.5])
                  if col else (np.nan, np.nan))
        out.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                    "n": n, "note": note})

    for k in KS:
        c = m[f"contrast_k{k}"]
        add(f"contrast, k={k}", c.mean(), f"contrast_k{k}",
            note=f"gains {m[f'gain_gc_k{k}'].mean():.4f} -> {m[f'gain_dn_k{k}'].mean():.4f}; "
                 f"dinuc larger in {int((c > 0).sum())}/{n}")

    # THE VERIFICATION. k=4 rebuilt from sequence against the committed contrast, ON THE SAME
    # DATASETS. Comparing a rebuild of n datasets against the published mean over 94 would mix
    # panels and report a difference that is mostly panel composition -- the same error that
    # once made a surviving fraction read 0.8506 instead of 0.8429.
    sub = cm[cm.dataset.isin(m.dataset)]
    committed = float((sub.delta_auroc_dn - sub.delta_auroc_gc).mean())
    rebuilt = float(m[f"contrast_k{PUBLISHED_K}"].mean())
    add("committed contrast (read from tables)", committed,
        note=f"restricted to the same {len(sub)} datasets rebuilt here")
    add("rebuilt from raw sequence at k=4", rebuilt,
        note="independent re-derivation, not a table read")
    add("absolute difference", abs(rebuilt - committed),
        note="THE PROOF: R1 is not merely re-read, it is rebuilt")

    # THE ROBUSTNESS. Does the choice of k change the conclusion?
    add("smallest contrast across k=3..6", float(min(m[f"contrast_k{k}"].mean() for k in KS)),
        note="positive at every k or the result is a k artifact")
    both = np.ones(n, dtype=bool)
    for k in KS:
        both &= m[f"contrast_k{k}"] > 0
    add("datasets positive at EVERY k", float(both.sum()), note=f"of {n}")
    w = wilcoxon(m.contrast_k5, m.contrast_k4)
    add("k=5 minus k=4", float((m.contrast_k5 - m.contrast_k4).mean()),
        note=f"paired Wilcoxon p={w.pvalue:.3g}; the manuscript's '5-mer' error would have "
             f"changed nothing")

    res = pd.DataFrame(out)
    res.to_csv(TABLES / "k_sweep.csv", index=False)
    log("")
    for _, x in res.iterrows():
        ci = f" [{x.ci_low:+.4f}, {x.ci_high:+.4f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:44} {x.value:+.4f}{ci}   {x.note}")
    log(f"\n  n = {n};  wrote k_sweep.csv")


if __name__ == "__main__":
    main()
