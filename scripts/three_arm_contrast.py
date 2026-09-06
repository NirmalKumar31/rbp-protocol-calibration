"""R1k: the same measurement under three negative-set protocols, not two.

    python scripts/three_arm_contrast.py --store ../rbp-store
    python scripts/three_arm_contrast.py --from-cache

The paper's whole argument is that a benchmark AUROC's meaning is set by how its negatives were
built. Demonstrating that with TWO protocols invites the reply that both are variants of one
flawed design -- composition-matched unbound genomic windows, with everything R1c and R1j say
about them. Three protocols, one of which is the field's own, is a different kind of evidence.

  gc      unbound windows, matched on region and GC content
  dinuc   unbound windows, matched on region and the full dinucleotide profile
  neg2    OTHER RBPs' binding sites in the same cell line, 1:1, target's own sites excluded
          -- Horlacher et al. 2023's bias-aware protocol, built by scripts/build_neg2.py

Only the 4-mer is run here. It is the model the headline is about, it costs nothing, and the
deep models were never trained on the third arm. Reporting a k-mer-only third arm is honest;
quietly comparing a k-mer on one arm against SpliceBERT on another would not be.

WHAT TO EXPECT, AND WHY THE PREDICTION IS RISKY. Under neg2 both classes are real crosslink
sites, so composition should separate them far less well than it does against untranscribed
genomic background -- the composition baseline should fall toward 0.5. If the nested
contribution then RISES the way it does under dinucleotide matching, the paper's claim holds
under a protocol that shares none of the two existing arms' defects. If it does not, that is a
real limit on the claim and it belongs in the paper.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = ("gc", "dn", "neg2")
DIRS = {"gc": "processed/gc", "dn": "processed/dinuc", "neg2": "processed/neg2"}



def per_dataset(store, datasets):
    store = Path(store)
    out = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        row = {"dataset": ds, "protein": protein, "cell": cell}
        ok = True
        for arm in ARMS:
            f = store / DIRS[arm] / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
            m = pd.DataFrame({"seq": d.seq_rna.values, "y": d.label.values,
                              "f": d.fold.values, "s": sc}).dropna()
            g = gain_over_composition(m.seq.values, m.s.values, m.y.values, m.f.values)
            row[f"comp_{arm}"] = g.auroc_composition
            row[f"full_{arm}"] = g.auroc_with_score
            row[f"gain_{arm}"] = g.delta
            row[f"n_{arm}"] = g.n
        if not ok:
            continue
        out.append(row)
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} comp "
            + "  ".join(f"{a} {row[f'comp_{a}']:.3f}" for a in ARMS)
            + "   gain " + "  ".join(f"{a} {row[f'gain_{a}']:+.4f}" for a in ARMS))
    return pd.DataFrame(out)


def summarise(d, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    rows = []

    def add(check, series, note=""):
        v = float(np.mean(series))
        b = np.array([np.asarray(series)[i].mean() for i in idx])
        lo, hi = np.percentile(b, [2.5, 97.5])
        rows.append({"check": check, "value": v, "ci_low": lo, "ci_high": hi,
                     "n": len(d), "note": note})
        return v

    for a in ARMS:
        add(f"composition alone, {a} arm", d[f"comp_{a}"])
        add(f"apparent AUROC (composition + score), {a} arm", d[f"full_{a}"])
        add(f"nested contribution, {a} arm", d[f"gain_{a}"],
            "the quantity the paper claims")
    for a, b in (("dn", "gc"), ("neg2", "gc"), ("neg2", "dn")):
        add(f"CONTRAST, {a} minus {b}", d[f"gain_{a}"] - d[f"gain_{b}"])
        pos = int(((d[f"gain_{a}"] - d[f"gain_{b}"]) > 0).sum())
        rows.append({"check": f"datasets with a positive contrast, {a} minus {b}",
                     "value": pos, "ci_low": np.nan, "ci_high": np.nan, "n": len(d),
                     "note": ""})
    # THE TRANSPLANT MATRIX. Is the whole three-arm pattern just compression? Carry each
    # arm's d' increment onto each other arm's baseline; the residual is the protocol effect
    # for that pair. Compression alone would make every residual zero.
    R2 = np.sqrt(2.0)
    dp = lambda a: R2 * norm.ppf(np.clip(a, 1e-6, 1 - 1e-6))   # noqa: E731
    au = lambda z: norm.cdf(z / R2)                            # noqa: E731
    for src in ARMS:
        inc = dp(d[f"full_{src}"]) - dp(d[f"comp_{src}"])
        for tgt in ARMS:
            if src == tgt:
                continue
            pred = au(dp(d[f"comp_{tgt}"]) + inc) - d[f"comp_{tgt}"]
            add(f"protocol effect, {src} increment on {tgt} baseline",
                d[f"gain_{tgt}"] - pred,
                "observed minus what compression alone predicts")
    rho, pv = spearmanr(np.r_[tuple(d[f"comp_{a}"] for a in ARMS)],
                        np.r_[tuple(d[f"gain_{a}"] for a in ARMS)])
    rows.append({"check": "spearman(composition baseline, nested gain), all arms pooled",
                 "value": float(rho), "ci_low": np.nan, "ci_high": np.nan,
                 "n": 3 * len(d), "note": f"p={pv:.2e}"})

    # The multiplier, on datasets where both arms have a positive gain.
    for a, b in (("dn", "gc"), ("neg2", "gc")):
        m = (d[f"gain_{a}"] > 0) & (d[f"gain_{b}"] > 0)
        lr = np.log(d.loc[m, f"gain_{a}"] / d.loc[m, f"gain_{b}"])
        rows.append({"check": f"protocol multiplier, {a} over {b}",
                     "value": float(np.exp(lr.mean())), "ci_low": np.nan,
                     "ci_high": np.nan, "n": int(m.sum()), "note": ""})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--from-cache", action="store_true")
    p.add_argument("--n-boot", type=int, default=2000)
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "three_arm_per_dataset.csv"
    if a.from_cache:
        d = pd.read_csv(per)
    else:
        panel = sorted(pd.read_csv(TABLES / "rehearsal_binding_gc.csv").dataset)
        d = per_dataset(a.store, panel)
        if d.empty:
            sys.exit("no dataset has all three arms")
        d.to_csv(per, index=False)

    s = summarise(d, n_boot=a.n_boot)
    s.to_csv(TABLES / "three_arm_contrast.csv", index=False)
    log(f"\n=== R1k: three negative-set protocols, 4-mer model, n = {len(d)} ===\n")
    q = s.set_index("check")
    log(f"  {'arm':6s} {'composition':>13s} {'apparent':>10s} {'nested contribution':>22s}")
    for arm, name in (("gc", "GC"), ("dn", "dinuc"), ("neg2", "neg2")):
        log(f"  {name:6s} {q.loc[f'composition alone, {arm} arm','value']:13.4f}"
            f" {q.loc[f'apparent AUROC (composition + score), {arm} arm','value']:10.4f}"
            f" {q.loc[f'nested contribution, {arm} arm','value']:+22.4f}")
    log("")
    for a_, b_ in (("dn", "gc"), ("neg2", "gc"), ("neg2", "dn")):
        r = q.loc[f"CONTRAST, {a_} minus {b_}"]
        n = int(q.loc[f"datasets with a positive contrast, {a_} minus {b_}", "value"])
        log(f"  contrast {a_:4s} - {b_:4s}  {r['value']:+.4f} "
            f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]   positive in {n}/{len(d)}")
    log("")
    for a_, b_ in (("dn", "gc"), ("neg2", "gc")):
        r = q.loc[f"protocol multiplier, {a_} over {b_}"]
        log(f"  multiplier {a_:4s} / {b_:4s}  {r['value']:.2f}x  (n={int(r['n'])})")
    log("\nwrote three_arm_per_dataset.csv and three_arm_contrast.csv")


if __name__ == "__main__":
    main()
