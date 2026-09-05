"""Did every committed score file actually use the study's chromosome folds? For 20 of them, no.

    python scripts/fold_integrity.py            # needs the window store
    python scripts/fold_integrity.py --from-cache

THE DEFECT, found by an adversarial statistician and confirmed three ways. `config/folds.tsv`
is a frozen chromosome-to-fold map, and chromosome grouping is the whole point: a model must
never be tested on a chromosome it trained on. For **20 of the 94 dinucleotide-arm datasets**,
the committed per-window CNN and SpliceBERT scores were produced under a DIFFERENT partition --
a stratified random split that preserved fold SIZES (so it is invisible to any size check) but
put **up to 23 chromosomes in every fold**. The GC and bias-aware arms are clean, 94 of 94.

The bias-aware arm was added here late. The manuscript already asserted it was clean, on the
grounds that its donors are sampled within fold, and that was the same kind of argument-from-
construction that the dinucleotide arm's docstring made and got wrong. It is the denominator of
every reported span, so it is the arm least safe to leave unmeasured.

WHY IT MATTERED AND WAS INVISIBLE. `deep_model_contrast.py` asserted in its own docstring that
"both arms use the same chromosome folds, the same seed ... the only difference is how the
negative windows were chosen". That was false for those 20, and none of the 549 verifier
assertions checked it -- the harness gated the VALUES the scores produced and never what
produced them, which is this project's oldest and most expensive lesson.

WHAT IS MEASURED HERE, per (dataset, arm, model):
  * agreement between the score file's own `fold` column and `dataset.tsv`'s;
  * the maximum number of distinct chromosomes inside any one fold, which is the structural
    tell -- under the study's design it is 4 or 5 and can never exceed it;
  * the DIRECT leakage metric: the fraction of scored rows having a same-strand genomic
    neighbour within 1 kb that was assigned to a different fold. Under chromosome grouping
    this is exactly 0 by construction.

AND THE SENSITIVITY, because the question a reader asks is not "is there a defect" but "does
any conclusion move". The k-mer is refit in-script on the study folds for every dataset, so it
is a clean internal control, and the R1g contrasts are recomputed with the 20 dropped.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rbp.utils.log import log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": ("processed/gc", "data/evidence/scores_gc"),
        "dn": ("processed/dinuc", "data/evidence/scores"),
        "neg2": ("processed/neg2", "data/evidence/scores_neg2")}
MODELS = ("cnn", "splicebert")
NEIGHBOUR_NT = 1000
MAX_CHROMS_PER_FOLD = 5      # the frozen map's largest fold holds 5 chromosomes
MODELS_ALL = ("kmer", "cnn", "splicebert")



def cross_fold_neighbours(chrom, strand, start, fold):
    """Fraction of rows with a same-strand neighbour within 1 kb in a DIFFERENT fold.

    Zero by construction under chromosome grouping, because a neighbour on the same
    chromosome is necessarily in the same fold. Non-zero is direct evidence that the
    partition is not chromosome-grouped, independent of any fold-label comparison.
    """
    df = pd.DataFrame({"chrom": chrom, "strand": strand, "start": start, "fold": fold})
    n_hit = 0
    for _, g in df.groupby(["chrom", "strand"], sort=False):
        s = g.start.to_numpy()
        f = g.fold.to_numpy()
        o = np.argsort(s, kind="mergesort")
        s, f = s[o], f[o]
        lo = np.searchsorted(s, s - NEIGHBOUR_NT, side="left")
        hi = np.searchsorted(s, s + NEIGHBOUR_NT, side="right")
        # Vectorised: a row has a cross-fold neighbour iff its window holds fewer rows of its
        # OWN fold than the window's total size. Prefix sums per fold label make that O(n) --
        # the obvious per-row loop is 14M Python iterations over the full panel.
        labels = np.unique(f)
        pre = {k: np.concatenate([[0], np.cumsum(f == k)]) for k in labels}
        own = np.empty(len(s), dtype=np.int64)
        for k in labels:
            m = f == k
            own[m] = pre[k][hi[m]] - pre[k][lo[m]]
        n_hit += int(np.count_nonzero(own < (hi - lo)))
    return float(n_hit / len(df)) if len(df) else 0.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--from-cache", action="store_true")
    p.add_argument("--n", type=int, default=0)
    a = p.parse_args()

    per = TABLES / "fold_integrity_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it with --store")
    else:
        store = Path(a.store)
        if not (store / "processed" / "gc").exists():
            sys.exit(f"no window store at {store}. It is ~3 GB and is not published; the audit "
                     f"table is committed, so use --from-cache to re-gate it.")
        panel = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv")
        if a.n:
            panel = panel.head(a.n)
        rows = []
        for i, r in enumerate(panel.itertuples(), 1):
            for arm, (sd, ev) in ARMS.items():
                f = store / sd / r.cell / r.protein / "dataset.tsv"
                if not f.exists():
                    continue
                d = pd.read_csv(f, sep="\t",
                                usecols=["id", "fold", "chrom", "strand", "start"])
                for model in MODELS:
                    parts, ok = [], True
                    for k in range(5):
                        sp = ROOT / ev / r.cell / r.protein / model / f"fold{k}" / "scores.tsv.gz"
                        if not sp.exists():
                            ok = False
                            break
                        parts.append(pd.read_csv(sp, sep="\t", usecols=["id", "fold"]))
                    if not ok:
                        continue
                    s = pd.concat(parts, ignore_index=True)
                    m = d.merge(s, on="id", suffixes=("_tab", "_score"))
                    if m.empty:
                        continue
                    rows.append({
                        "dataset": r.dataset, "protein": r.protein, "cell": r.cell,
                        "arm": arm, "model": model, "n": len(m),
                        "fold_agreement": float((m.fold_tab == m.fold_score).mean()),
                        "max_chroms_per_score_fold":
                            int(m.groupby("fold_score").chrom.nunique().max()),
                        "max_chroms_per_table_fold":
                            int(m.groupby("fold_tab").chrom.nunique().max()),
                        "cross_fold_neighbours_score": cross_fold_neighbours(
                            m.chrom, m.strand, m.start, m.fold_score),
                        "cross_fold_neighbours_table": cross_fold_neighbours(
                            m.chrom, m.strand, m.start, m.fold_tab),
                    })
            if i % 10 == 0:
                log(f"  [{i}/{len(panel)}]")
        t = pd.DataFrame(rows)
        if t.empty:
            sys.exit("nothing audited; refusing to overwrite the committed table")
        t.to_csv(per, index=False)

    t["grouped"] = t.max_chroms_per_score_fold <= MAX_CHROMS_PER_FOLD
    out = []
    log(f"\n=== fold integrity: {len(t)} (dataset, arm, model) score sets ===\n")
    log(f"  {'arm':5s} {'sets':>5s} {'chromosome-grouped':>19s} {'aligned to folds.tsv':>21s}"
        f" {'max chroms/fold':>16s}")
    for arm in ARMS:
        s = t[t.arm == arm]
        if s.empty:
            continue
        n_ds = s.dataset.nunique()
        bad = s[~s.grouped].dataset.nunique()
        aligned = s[s.fold_agreement > 0.999].dataset.nunique()
        out += [
            {"check": f"score sets audited, {arm} arm", "value": len(s), "n": n_ds},
            {"check": f"datasets NOT chromosome-grouped, {arm} arm", "value": bad, "n": n_ds},
            {"check": f"datasets aligned to folds.tsv, {arm} arm", "value": aligned, "n": n_ds},
            {"check": f"max chromosomes in any score fold, {arm} arm",
             "value": int(s.max_chroms_per_score_fold.max()), "n": n_ds},
            {"check": f"max cross-fold 1kb neighbour fraction, {arm} arm",
             "value": float(s.cross_fold_neighbours_score.max()), "n": n_ds},
            {"check": f"max cross-fold 1kb neighbour fraction under folds.tsv, {arm} arm",
             "value": float(s.cross_fold_neighbours_table.max()), "n": n_ds}]
        log(f"  {arm:5s} {len(s):5d} {n_ds - bad:11d}/{n_ds:<7d} {aligned:13d}/{n_ds:<7d}"
            f" {int(s.max_chroms_per_score_fold.max()):16d}")

    leaky = sorted(t.loc[~t.grouped, "dataset"].unique())
    out.append({"check": "leaky datasets", "value": len(leaky), "n": len(leaky),
                "note": ";".join(leaky)})
    if leaky:
        log(f"\n  NOT chromosome-grouped ({len(leaky)} datasets, dinucleotide arm):")
        w = t[~t.grouped].groupby("dataset").agg(
            agree=("fold_agreement", "min"),
            nbr=("cross_fold_neighbours_score", "max")).sort_values("nbr", ascending=False)
        for ds, r in w.head(8).iterrows():
            log(f"    {ds:16s} fold agreement {r.agree:.3f}   "
                f"cross-fold 1kb neighbours {100 * r.nbr:.1f}%")
        log(f"    ... and {max(len(w) - 8, 0)} more")

    # THE SENSITIVITY. Does any conclusion move? The k-mer is refit on the study folds for
    # every dataset, so it is a clean internal control against which the deep arms are read.
    d = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv")
    d["clean"] = ~d.dataset.isin(leaky)
    out.append({"check": "clean datasets for the R1g sensitivity", "value": int(d.clean.sum()),
                "n": len(d)})
    log(f"\n  R1g contrast (dn - gc), all {len(d)} vs the {int(d.clean.sum())} clean:")
    for model in MODELS_ALL:
        c = d[f"{model}_gain_dn"] - d[f"{model}_gain_gc"]
        full, clean = float(c.mean()), float(c[d.clean].mean())
        out += [{"check": f"{model} R1 contrast, all datasets", "value": full, "n": len(d)},
                {"check": f"{model} R1 contrast, chromosome-grouped only", "value": clean,
                 "n": int(d.clean.sum())},
                {"check": f"{model} R1 contrast shift when leaky sets are dropped",
                 "value": clean - full, "n": len(d)}]
        log(f"    {model:11s} {full:+.4f} -> {clean:+.4f}   shift {clean - full:+.4f}")

    # And the difference-in-differences, which isolates the leakage from dataset selection:
    # the deep-minus-kmer gap in the dn arm relative to the same datasets' gc arm.
    #
    # ONLY MEANINGFUL WHILE THERE IS SOMETHING LEAKY TO COMPARE AGAINST. After the 20 stale
    # datasets were retrained the leaky set is empty, and a mean over an empty selection is
    # NaN: the block printed "leaky +nan ... DiD +nan" and a median of nan pairs, which reads
    # as a broken pipeline rather than as a repaired one. Say what happened instead.
    if not leaky:
        log("\n  no datasets are off-partition, so the leaky-versus-clean sensitivity has "
            "nothing to compare and is not computed.")
        log("  the 20 datasets that were off-partition were retrained on the study folds; "
            "see cloud/modal/retrain_dinuc_20.txt for the frozen list.")
    else:
        log("\n  difference-in-differences (deep minus k-mer gap, dn vs gc), leaky vs clean:")
        for model in MODELS:
            gap = ((d[f"{model}_gain_dn"] - d["kmer_gain_dn"])
                   - (d[f"{model}_gain_gc"] - d["kmer_gain_gc"]))
            did = float(gap[~d.clean].mean() - gap[d.clean].mean())
            out.append({"check": f"{model} difference-in-differences, leaky minus clean",
                        "value": did, "n": len(d),
                        "note": "leaky datasets are also the largest; size-confounded"})
            log(f"    {model:11s} leaky {gap[~d.clean].mean():+.4f}  "
                f"clean {gap[d.clean].mean():+.4f}  DiD {did:+.4f}")
        # WINDOWS, not pairs: n_dn is len() of the scored frame, which holds one row per
        # window and so two per pair. Reported as pairs for two drafts, and quoted that way in
        # the manuscript, where every other size is in pairs.
        log("\n  NOTE: the leaky datasets are the panel's largest (median "
            f"{d.loc[~d.clean, 'n_dn'].median()/2:.0f} vs "
            f"{d.loc[d.clean, 'n_dn'].median()/2:.0f} "
            "pairs), so the DiD confounds leakage with size and is an UPPER bound.")

    pd.DataFrame(out).to_csv(TABLES / "fold_integrity.csv", index=False)
    log("\nwrote fold_integrity.csv and fold_integrity_per_dataset.csv")


if __name__ == "__main__":
    main()
