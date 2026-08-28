"""Q1, properly: is R1's contrast a strand artifact? Restriction plus a matched placebo.

WHY THE OBVIOUS VERSION OF THIS TEST IS WRONG, AND WHY THE PLACEBO IS THE WHOLE EXPERIMENT.

`negatives.py:328` gives each negative the POSITIVE's strand, because `annotation.py:126`
deliberately drops region strand ("A window's strand comes from its peak, so the region's own
strand is never needed"). True for positives, false for negatives: only ~55% of negatives end
up on the strand their own gene is transcribed from. The obvious control is to keep only pairs
whose negative is unambiguously sense and recompute the contrast. Run alone, that control
LIES. Restricting to sense-only pairs discards roughly half the training data, and a 256-feature
k-mer model loses more from that than a 19-feature composition baseline does -- in both arms.
So the contrast shrinks whether or not strand matters, and the naive reading attributes the
whole shrinkage to strand.

The fix is a placebo: drop the SAME NUMBER of pairs at random, several times, and compare. The
strand-specific effect is the difference between the two, and everything else cancels.

An earlier attempt at this question regressed the per-dataset contrast on each dataset's sense
fraction and found rho = -0.24 [-0.54, +0.11]. That test was weak by construction: `frac_sense`
spans only 0.433-0.615 across datasets, so a BETWEEN-dataset regression has almost no power
against a bias that is present WITHIN every dataset, and `frac_sense` is not exogenous anyway
(it correlates +0.427 with GC-arm AUROC, so it proxies region mix). Restriction moves
`frac_sense` to 1.0 by construction and has full leverage. See `scripts/strand_contrast.py`
for the weak version, retained because the contrast between the two designs is the point.

PRE-REGISTERED in docs/61 before this was run: sign retained, CI excluding zero, and at least
60% of the point estimate surviving.

REPRODUCTION IS CHECKED PER DATASET. Local window tables are used only where recomputing the
published composition and with-score AUROCs from them lands on the committed rehearsal row.
The GC arm reproduces 40/40 locally; the dinucleotide arm reproduces 13/40, because the local
copy is a different draw, so its canonical tables are read from the study bucket instead.
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

from rbp.eval import baseline, nested                                   # noqa: E402
from strand_audit import gene_index, own_strands                        # noqa: E402

TABLES = ROOT / "results" / "tables"
REPRO_TOL = 2.0e-3
N_PLACEBO = 5
N_BOOT = 2000
SEED = 0


def log(m):
    print(m, flush=True)


def pair_key(ids):
    """`PROT_pos_17` / `PROT_neg_17` -> 17. The two members share a key."""
    return np.array([s.rsplit("_", 1)[-1] for s in ids])


def nested_gain(d):
    """The published quantity: out-of-fold AUROC(composition + score) - AUROC(composition)."""
    res = baseline.evaluate(d, k=4)
    g = nested.gain_over_composition(d.seq_rna.tolist(), res["scores"],
                                     d.label.to_numpy(), d.fold.to_numpy())
    return g.delta, g.auroc_composition, g.auroc_with_score


def sense_pairs(d, idx):
    """Keys of pairs whose NEGATIVE sits unambiguously on a gene of its assigned strand."""
    neg = d[d.label == 0]
    keep = []
    for k, c, s, e, a in zip(pair_key(neg.id), neg.chrom, neg.start, neg.end, neg.strand):
        ss = own_strands(idx, c, int(s), int(e))
        if len(ss) == 1 and next(iter(ss)) == a:        # unambiguous AND matching
            keep.append(k)
    return set(keep)


def subset(d, keys):
    k = pair_key(d.id)
    return d[np.isin(k, list(keys))]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--gc-root", required=True)
    p.add_argument("--dn-root", required=True)
    p.add_argument("--limit", type=int, default=0)
    # Rebuild the summary from the committed per-dataset table without redoing the ~30 minutes
    # of refits. The per-dataset table IS the evidence; the summary is arithmetic on it.
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    if a.from_cache:
        return summarise(pd.read_csv(TABLES / "strand_placebo_per_dataset.csv"))

    log("building gene index")
    idx = gene_index(a.gtf)

    audited = list(pd.read_csv(TABLES / "strand_audit.csv").dataset)
    if a.limit:
        audited = audited[:a.limit]
    pub = {"gc": pd.read_csv(TABLES / "rehearsal_binding_gc.csv").set_index("dataset"),
           "dn": pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv").set_index("dataset")}
    root = {"gc": Path(a.gc_root), "dn": Path(a.dn_root)}

    rows = []
    for i, ds in enumerate(audited, 1):
        prot, cell = ds.split(":")
        per = {}
        ok = True
        for arm in ("gc", "dn"):
            f = root[arm] / cell / prot / "dataset.tsv"
            if not f.exists() or ds not in pub[arm].index:
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            r = pub[arm].loc[ds]
            full, comp, with_s = nested_gain(d)
            # REPRODUCTION GATE. A local table that does not rebuild the published row is a
            # different draw, and differencing against it measures the draw, not the strand.
            if (abs(comp - r.composition_auroc) > REPRO_TOL
                    or abs(with_s - r.with_score_auroc) > REPRO_TOL):
                ok = False
                break
            per[arm] = {"d": d, "full": full}
        if not ok:
            log(f"  [{i:2d}/{len(audited)}] {ds:22} SKIP (does not reproduce)")
            continue

        rec = {"dataset": ds}
        for arm in ("gc", "dn"):
            d = per[arm]["d"]
            allk = set(pair_key(d.id))
            sk = sense_pairs(d, idx)
            n_keep = len(sk)
            rec[f"n_pairs_{arm}"] = len(allk)
            rec[f"n_sense_{arm}"] = n_keep
            rec[f"full_{arm}"] = per[arm]["full"]
            rec[f"sense_{arm}"] = nested_gain(subset(d, sk))[0] if 0 < n_keep < len(allk) \
                else np.nan
            # PLACEBO: same number of pairs, chosen at random, several seeds.
            pl = []
            for s in range(N_PLACEBO):
                rng = np.random.default_rng(1000 + s)
                pick = set(rng.choice(sorted(allk), size=n_keep, replace=False))
                pl.append(nested_gain(subset(d, pick))[0])
            rec[f"placebo_{arm}"] = float(np.mean(pl))
            rec[f"placebo_sd_{arm}"] = float(np.std(pl))
        rows.append(rec)
        log(f"  [{i:2d}/{len(audited)}] {ds:22} sense {rec['n_sense_gc']}/{rec['n_pairs_gc']} gc, "
            f"{rec['n_sense_dn']}/{rec['n_pairs_dn']} dn   "
            f"contrast full {rec['full_dn'] - rec['full_gc']:+.4f} "
            f"sense {rec['sense_dn'] - rec['sense_gc']:+.4f} "
            f"placebo {rec['placebo_dn'] - rec['placebo_gc']:+.4f}")
        pd.DataFrame(rows).to_csv(TABLES / "strand_placebo_per_dataset.csv", index=False)

    m = pd.DataFrame(rows).dropna()
    if not len(m):
        raise SystemExit("no datasets reproduced; nothing to report")
    return summarise(m)


def summarise(m):
    if "c_full" not in m.columns:
        m["c_full"] = m.full_dn - m.full_gc
        m["c_sense"] = m.sense_dn - m.sense_gc
        m["c_placebo"] = m.placebo_dn - m.placebo_gc
        m["excess"] = m.c_sense - m.c_placebo
    # THE STRAND-CORRECTED CONTRAST. Only the strand-specific part is removed; the shrinkage
    # the placebo also shows is an artifact of discarding pairs, not of strand, so subtracting
    # it would be double-counting and would understate the effect the paper claims.
    m["corrected"] = m.c_full + m.excess
    m.to_csv(TABLES / "strand_placebo_per_dataset.csv", index=False)

    rng = np.random.default_rng(SEED)
    n = len(m)
    boots = {k: [] for k in ("c_full", "c_sense", "c_placebo", "excess",
                             "d_sense", "d_placebo", "corrected")}
    for _ in range(N_BOOT):
        s = m.iloc[rng.integers(0, n, n)]
        boots["c_full"].append(s.c_full.mean())
        boots["c_sense"].append(s.c_sense.mean())
        boots["c_placebo"].append(s.c_placebo.mean())
        boots["excess"].append(s.excess.mean())
        boots["d_sense"].append((s.c_sense - s.c_full).mean())
        boots["d_placebo"].append((s.c_placebo - s.c_full).mean())
        boots["corrected"].append(s.corrected.mean())

    out = []

    def add(check, value, key=None, note=""):
        lo, hi = np.percentile(boots[key], [2.5, 97.5]) if key else (np.nan, np.nan)
        out.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                    "n": n, "note": note})

    add("contrast, full data", m.c_full.mean(), "c_full",
        note="published contrast on all 94 is +0.0397")
    add("contrast, sense-only pairs", m.c_sense.mean(), "c_sense",
        note=f"mean {m.n_sense_gc.sum() / m.n_pairs_gc.sum():.1%} of GC pairs retained")
    add("contrast, PLACEBO (same n, random)", m.c_placebo.mean(), "c_placebo",
        note=f"{N_PLACEBO} seeds per dataset per arm")
    add("change from restriction", (m.c_sense - m.c_full).mean(), "d_sense",
        note="what the naive test would have reported as strand")
    add("change from placebo", (m.c_placebo - m.c_full).mean(), "d_placebo",
        note="the same shrinkage, with no strand involved")
    add("STRAND-SPECIFIC EXCESS", m.excess.mean(), "excess",
        note="THE ANSWER: restriction minus placebo")
    add("strand-CORRECTED contrast", m.corrected.mean(), "corrected",
        note="full contrast with only the strand-specific part removed")
    # ON THIS PANEL'S OWN CONTRAST, not the n=94 published one. Inserting +0.0397 into an
    # n=40 computation mixes panels and reported 0.8506 where the honest figure is 0.8429.
    add("fraction of the contrast surviving",
        m.corrected.mean() / m.c_full.mean(),
        note=f"floor was 0.60, pre-registered; panel's own contrast is "
             f"{m.c_full.mean():+.4f}, not the n=94 +0.0397")

    res = pd.DataFrame(out)
    res.to_csv(TABLES / "strand_placebo.csv", index=False)
    log("")
    for _, x in res.iterrows():
        ci = f" [{x.ci_low:+.4f}, {x.ci_high:+.4f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:38} {x.value:+.4f}{ci}   {x.note}")
    log(f"\n  n = {n} datasets;  wrote strand_placebo.csv")


if __name__ == "__main__":
    main()
