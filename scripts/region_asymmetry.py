"""The bias-aware arm does not match transcript region. How much of its result is that?

    python scripts/region_asymmetry.py              # needs the window store
    python scripts/region_asymmetry.py --from-cache

THE ASYMMETRY, and it was undisclosed until three referees found it independently. Both
composition-matched matchers draw candidates of the SAME region class -- `gc` requires it and
`dinuc` buckets on `(region, chrom)` -- so the region label carries nothing about the label
there. `build_neg2.py` matches fold only and says so in its own docstring, which the manuscript
did not echo. That matters because the bias-aware arm's high composition baseline is the one
result the paper calls not-implied-by-design, and it was explained as biology: "distinct RBPs
occupy compositionally distinct sites".

THREE MEASUREMENTS:

  1. How much the region label alone separates the classes, per arm. Scored by the likelihood
     ratio over the five region marginals, which is the optimal region-only score.
  2. A region-matched bias-aware arm built through the pipeline, `build_neg2.py
     --match-region`, which stratifies the draw on region instead of drawing uniformly inside
     the fold. It keeps every pair (456,734, the same as the published arm) and achieves the
     match exactly: region-only AUROC 0.5000, TV 0.0000.
  3. The same thing by reweighting the donors ALREADY drawn, as a second estimate from a
     different construction. It discards 30% of rows, so it is the weaker of the two and is
     reported only because agreement between two constructions is worth more than either.

THE ANSWER: region alone separates the bias-aware classes at median AUROC 0.748 against exactly
0.5000 in both other arms. Matching it lowers the arm's composition baseline from 0.8248 to
0.8052 and its contribution from +0.0122 to +0.0092, so 46% of the arm's baseline excess over
the GC arm is region mix. The arm STILL carries the highest baseline and the lowest contribution
of the three, and the span widens from 5.42 to 7.20. So the ordering is not a region artefact,
but the mechanism sentence was overstated and the magnitude is part annotation.

NOTE WHAT THIS ARM IS NOT. Horlacher's negative-2 does not match region either, so a
region-matched version is not the field's protocol and is not a fourth point on the same axis.
It is a diagnostic that measures how much of one arm's baseline is transcript annotation.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from rbp.utils.log import log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
DIRS = {"dn": "dinuc", "gc": "gc", "neg2": "neg2"}
SEED = 7
EPS = 1e-9



def region_lr_auroc(label, region):
    """AUROC of the optimal region-only score, and the marginals' total variation distance."""
    p = pd.Series(region[label == 1]).value_counts(normalize=True)
    q = pd.Series(region[label == 0]).value_counts(normalize=True)
    keys = set(p.index) | set(q.index)
    lr = {k: np.log((p.get(k, 0.0) + EPS) / (q.get(k, 0.0) + EPS)) for k in keys}
    tv = 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)
    return float(roc_auc_score(label, pd.Series(region).map(lr))), float(tv)


def region_matched_rows(d, rng):
    """All positives, plus the negatives that reproduce their region marginals, per fold."""
    keep = []
    for _, g in d.groupby("fold", sort=True):
        pos, neg = g[g.label == 1], g[g.label == 0]
        if pos.empty or neg.empty:
            continue
        keep.append(pos.index.to_numpy())
        share = pos.region.value_counts(normalize=True)
        avail = neg.region.value_counts()
        n_take = int(np.floor(min(avail.get(r, 0) / s for r, s in share.items())))
        for r, s in share.items():
            idx = neg.index[neg.region.values == r].to_numpy()
            k = min(int(round(n_take * s)), len(idx))
            if k > 0:
                keep.append(rng.choice(idx, size=k, replace=False))
    return np.concatenate(keep) if keep else np.array([], dtype=int)


def measure(store, panel, rng):
    from rbp.eval.baseline import oof_scores as kmer_oof
    from rbp.eval.nested import gain_over_composition

    rows = []
    for i, r in enumerate(panel.itertuples(), 1):
        row = {"dataset": r.dataset, "protein": r.protein, "cell": r.cell}
        ok = True
        for arm, sub in DIRS.items():
            f = store / sub / r.cell / r.protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t", usecols=["label", "region", "fold", "seq_rna"])
            a, tv = region_lr_auroc(d.label.to_numpy(), d.region.to_numpy())
            row[f"region_auroc_{arm}"], row[f"region_tv_{arm}"] = a, tv
            if arm != "neg2":
                continue

            # THE PRIMARY REGION-MATCHED ESTIMATE: the arm rebuilt through the pipeline with
            # the draw stratified on region. Every pair is kept, so this is not a subsample of
            # the published arm and carries no selection caveat.
            f2 = store / "neg2_rm" / r.cell / r.protein / "dataset.tsv"
            if f2.exists():
                d2 = pd.read_csv(f2, sep="\t", usecols=["label", "region", "fold", "seq_rna"])
                a3, tv3 = region_lr_auroc(d2.label.to_numpy(), d2.region.to_numpy())
                sc2, _, _ = kmer_oof(d2.seq_rna.values, d2.label.values, d2.fold.values, k=4)
                g2 = np.isfinite(sc2)
                r2 = gain_over_composition(d2.seq_rna.values[g2], sc2[g2],
                                           d2.label.values[g2], d2.fold.values[g2])
                row.update({"region_auroc_built": a3, "region_tv_built": tv3,
                            "retained_built": len(d2) / len(d),
                            "comp_built": r2.auroc_composition,
                            "full_built": r2.auroc_with_score, "gain_built": r2.delta})

            m = d.loc[region_matched_rows(d, rng)].reset_index(drop=True)
            if m.label.nunique() < 2 or len(m) < 200:
                ok = False
                break
            a2, tv2 = region_lr_auroc(m.label.to_numpy(), m.region.to_numpy())
            sc, _, _ = kmer_oof(m.seq_rna.values, m.label.values, m.fold.values, k=4)
            g = np.isfinite(sc)
            res = gain_over_composition(m.seq_rna.values[g], sc[g], m.label.values[g],
                                        m.fold.values[g])
            row.update({"region_auroc_matched": a2, "region_tv_matched": tv2,
                        "retained": len(m) / len(d), "comp_matched": res.auroc_composition,
                        "full_matched": res.auroc_with_score, "gain_matched": res.delta})
        if ok:
            rows.append(row)
        if i % 10 == 0:
            log(f"  [{i}/{len(panel)}]")
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    per = TABLES / "region_asymmetry_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
    else:
        store = Path(a.store) / "processed"
        if not (store / "neg2").exists():
            sys.exit(f"no window store at {store}; use --from-cache to re-gate the table")
        t = measure(store, pd.read_csv(TABLES / "three_arm_per_dataset.csv"),
                    np.random.default_rng(SEED))
        if t.empty:
            sys.exit("nothing measured; refusing to overwrite the committed table")
        t.to_csv(per, index=False)

    m = t.merge(pd.read_csv(TABLES / "three_arm_per_dataset.csv"),
                on=["dataset", "protein", "cell"])
    out = []
    log(f"\n=== region-only separability, n={len(m)} datasets ===\n")
    for arm in DIRS:
        c = m[f"region_auroc_{arm}"]
        out += [{"check": f"region-only AUROC, {arm} arm", "value": float(c.median()),
                 "n": len(m), "note": f"mean {c.mean():.4f}, max {c.max():.4f}, "
                                      f"above 0.70 in {int((c > 0.70).sum())}"},
                {"check": f"region marginal TV distance, {arm} arm",
                 "value": float(m[f'region_tv_{arm}'].median()), "n": len(m), "note": ""}]
        log(f"  {arm:6s} median {c.median():.4f}  mean {c.mean():.4f}  "
            f"IQR {c.quantile(.25):.3f}-{c.quantile(.75):.3f}  "
            f"above 0.70 in {int((c > 0.70).sum())}/{len(m)}   "
            f"TV {m[f'region_tv_{arm}'].median():.4f}")

    built = "comp_built" in m.columns and m.comp_built.notna().all()
    log("\n=== bias-aware arm with region matched ===\n")
    log(f"  {'construction':26s} {'kept':>6s} {'region AUROC':>13s} {'composition':>12s}"
        f" {'contribution':>13s}")
    log(f"  {'published (region free)':26s} {'100%':>6s} "
        f"{m.region_auroc_neg2.median():13.4f} {m.comp_neg2.mean():12.4f} "
        f"{m.gain_neg2.mean():+13.4f}")
    if built:
        log(f"  {'pipeline, stratified draw':26s} {m.retained_built.median():6.0%} "
            f"{m.region_auroc_built.median():13.4f} {m.comp_built.mean():12.4f} "
            f"{m.gain_built.mean():+13.4f}")
    log(f"  {'reweighting the drawn arm':26s} {m.retained.median():6.0%} "
        f"{m.region_auroc_matched.median():13.4f} {m.comp_matched.mean():12.4f} "
        f"{m.gain_matched.mean():+13.4f}")
    log(f"  {'GC arm':26s} {'':>6s} {m.region_auroc_gc.median():13.4f} "
        f"{m.comp_gc.mean():12.4f} {m.gain_gc.mean():+13.4f}")
    log(f"  {'dinucleotide arm':26s} {'':>6s} {m.region_auroc_dn.median():13.4f} "
        f"{m.comp_dn.mean():12.4f} {m.gain_dn.mean():+13.4f}")

    out += [{"check": "rows retained under region matching",
             "value": float(m.retained.median()), "n": len(m), "note": "reweighting"},
            {"check": "region-only AUROC after matching, neg2 arm",
             "value": float(m.region_auroc_matched.median()), "n": len(m), "note": ""},
            {"check": "composition alone, neg2 arm region-matched",
             "value": float(m.comp_matched.mean()), "n": len(m), "note": "reweighting"},
            {"check": "nested contribution, neg2 arm region-matched",
             "value": float(m.gain_matched.mean()), "n": len(m), "note": "reweighting"}]
    if built:
        excess = (m.comp_neg2.mean() - m.comp_gc.mean())
        share = (m.comp_neg2.mean() - m.comp_built.mean()) / excess if excess else float("nan")
        out += [{"check": "rows retained, pipeline region-matched arm",
                 "value": float(m.retained_built.median()), "n": len(m), "note": ""},
                {"check": "region-only AUROC, pipeline region-matched arm",
                 "value": float(m.region_auroc_built.median()), "n": len(m), "note": ""},
                {"check": "composition alone, pipeline region-matched arm",
                 "value": float(m.comp_built.mean()), "n": len(m), "note": "PRIMARY"},
                {"check": "nested contribution, pipeline region-matched arm",
                 "value": float(m.gain_built.mean()), "n": len(m), "note": "PRIMARY"},
                {"check": "share of the neg2 baseline excess over gc that is region mix",
                 "value": float(share), "n": len(m), "note": ""}]
        log(f"\n  region accounts for {100 * share:.0f}% of the bias-aware arm's baseline "
            f"excess over the GC arm")

    key_comp = m.comp_built if built else m.comp_matched
    key_gain = m.gain_built if built else m.gain_matched
    still_high = bool(key_comp.mean() > max(m.comp_gc.mean(), m.comp_dn.mean()))
    still_low = bool(key_gain.mean() < min(m.gain_gc.mean(), m.gain_dn.mean()))
    span = float(m.gain_dn.mean() / key_gain.mean())
    out += [{"check": "region-matched neg2 still has the highest baseline",
             "value": int(still_high), "n": len(m), "note": ""},
            {"check": "region-matched neg2 still has the lowest contribution",
             "value": int(still_low), "n": len(m), "note": ""},
            {"check": "span, dinuc over region-matched neg2", "value": span, "n": len(m),
             "note": f"published dinuc over neg2 {m.gain_dn.mean() / m.gain_neg2.mean():.2f}"}]
    log(f"\n  still the highest baseline of the three?     {still_high}")
    log(f"  still the lowest contribution of the three?  {still_low}")
    log(f"  span dinuc / region-matched neg2: {span:.2f}   "
        f"(published {m.gain_dn.mean() / m.gain_neg2.mean():.2f})")

    # The dose-response, which is what makes region a candidate explanation rather than a
    # coincidence: the more region separates, the more the baseline rises and the less the
    # model adds.
    for lab, y in (("baseline rise, neg2 over gc", m.comp_neg2 - m.comp_gc),
                   ("contribution deficit, neg2 minus gc", m.gain_neg2 - m.gain_gc)):
        rho, pv = spearmanr(m.region_auroc_neg2, y)
        out.append({"check": f"spearman(region-only AUROC, {lab})", "value": float(rho),
                    "n": len(m), "note": f"p={pv:.4f}"})
        log(f"\n  spearman(region-only AUROC, {lab}) = {rho:+.3f}  p={pv:.4f}")

    # And the half of the panel where region separates least, as a second check that does not
    # depend on the reweighting.
    lo = m.region_auroc_neg2 <= m.region_auroc_neg2.median()
    out += [{"check": "composition alone, neg2 arm, least region-separable half",
             "value": float(m.loc[lo, "comp_neg2"].mean()), "n": int(lo.sum()), "note": ""},
            {"check": "composition alone, gc arm, same least region-separable datasets",
             "value": float(m.loc[lo, "comp_gc"].mean()), "n": int(lo.sum()), "note": ""},
            {"check": "contribution ratio neg2 over gc, least region-separable half",
             "value": float(m.loc[lo, "gain_neg2"].mean() / m.loc[lo, "gain_gc"].mean()),
             "n": int(lo.sum()), "note": ""}]
    log(f"\n  least region-separable half (n={int(lo.sum())}): comp_neg2 "
        f"{m.loc[lo, 'comp_neg2'].mean():.4f} vs comp_gc {m.loc[lo, 'comp_gc'].mean():.4f}; "
        f"gain ratio {m.loc[lo, 'gain_neg2'].mean() / m.loc[lo, 'gain_gc'].mean():.3f}")

    pd.DataFrame(out).to_csv(TABLES / "region_asymmetry.csv", index=False)
    log("\nwrote region_asymmetry.csv and region_asymmetry_per_dataset.csv")


if __name__ == "__main__":
    main()
