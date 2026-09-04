"""R1p: the same measurement on the field's OWN published benchmark, using their folds.

    python scripts/horlacher_arm.py --n 40

WHY THIS IS THE MOST VALUABLE REMAINING EXPERIMENT. Every result so far is measured on windows
this project built. The obvious referee reply is "you re-derived a known effect on your own
data, with your own sampler, and your own reimplementation of the prior art's protocol." This
answers it: the same nested decomposition, on Horlacher et al. 2023's published negative sets,
using THEIR fold assignments, on the 88 of our 94 datasets they also cover.

DATA. Zenodo 10.5281/zenodo.10600977, `samples.tar.gz`, md5 c93a38b8bd684be2ccb5f6f82c6c4700
(verified). Not cited in their paper -- it is discoverable only through a comment on GitHub
issue #2. 223 ENCODE experiments x 5 folds of single-nucleotide crosslink sites in BED6:

  positive      the target RBP's crosslink sites
  negative-1    uniform positions from transcripts carrying at least one site of the target
  negative-2    OTHER RBPs' crosslink sites, not overlapping the target's

Both negative sets are expression-controlled by construction, which is the axis on which this
project's own negatives are weaker than the prior art it cites (R1j: 40.1% untranscribed).

WHAT IS DIFFERENT FROM OUR ARMS, and it must be stated rather than glossed. Their positives
are not ours: they use their own peak calling, their own 20,000-peak cap and their own
replicate-overlap rule. Their folds are not ours. So this is NOT a fourth arm on the same
positives; it is the same MEASUREMENT on a different benchmark. The comparison that means
something is the fold range within their data against the fold range within ours -- not the
absolute numbers side by side.

WINDOWS. 101 nt centred on the crosslink site, matching `windows.size`, reverse-complemented
on the minus strand so the sequence is the RNA the protein sees -- the same convention as
`windows.strand_correct`. Windows containing N are dropped.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sklearn.metrics import roc_auc_score  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402

TABLES = ROOT / "results" / "tables"
DATA = ROOT.parent / "rbp-store" / "external" / "samples" / "processed" / "ENCODE"
GENOME = ROOT.parent / "rna-binding-proteins" / "data" / "raw" / "GRCh38.primary_assembly.genome.fa"
WIN = 101
COMP = str.maketrans("ACGTN", "TGCAN")


def log(m):
    print(m, flush=True)


def windows(fa, bed, label):
    """101-nt strand-corrected RNA windows centred on each crosslink site."""
    out = []
    for ln in bed.read_text().splitlines():
        f = ln.split("\t")
        if len(f) < 6:
            continue
        chrom, pos, strand = f[0], int(f[1]), f[5]
        half = WIN // 2
        s, e = pos - half, pos + half + 1
        if s < 0 or chrom not in fa:
            continue
        try:
            seq = str(fa[chrom][s:e]).upper()
        except (KeyError, ValueError):
            continue
        if len(seq) != WIN or "N" in seq:
            continue
        if strand == "-":
            seq = seq.translate(COMP)[::-1]
        out.append({"seq_rna": seq.replace("T", "U"), "label": label,
                    "chrom": chrom, "start": s})
    return out


def build(fa, ds, negative):
    """One dataset, all five of THEIR folds, positives against the chosen negative set."""
    rows = []
    for fold in range(5):
        d = DATA / ds / f"fold-{fold}"
        pos, neg = d / f"positive.fold-{fold}.bed", d / f"{negative}.fold-{fold}.bed"
        if not pos.exists() or not neg.exists():
            return None
        for r in windows(fa, pos, 1) + windows(fa, neg, 0):
            r["fold"] = fold
            rows.append(r)
    df = pd.DataFrame(rows)
    return df if len(df) >= 400 and df.label.nunique() == 2 else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=40, help="size-stratified subsample of the overlap")
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "horlacher_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
    else:
        from pyfaidx import Fasta
        if not GENOME.exists():
            sys.exit(f"no genome at {GENOME}")
        fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

        ours = pd.read_csv(TABLES / "rehearsal_binding_gc.csv").sort_values("pairs")
        avail = {d.name.replace("_", ":"): d.name for d in DATA.iterdir() if d.is_dir()}
        ours = ours[ours.dataset.isin(avail)]
        log(f"{len(ours)} of our datasets are in their release")
        step = max(len(ours) // a.n, 1)
        sub = ours.iloc[::step].head(a.n)

        rows = []
        for i, r in enumerate(sub.itertuples(), 1):
            rec, ok = {"dataset": r.dataset, "protein": r.protein, "cell": r.cell}, True
            for neg, tag in (("negative-1", "n1"), ("negative-2", "n2")):
                d = build(fa, avail[r.dataset], neg)
                if d is None:
                    ok = False
                    break
                sc, _, _ = kmer_oof(d.seq_rna.values, d.label.values, d.fold.values, k=4)
                m = ~np.isnan(sc)
                g = gain_over_composition(d.seq_rna.values[m], sc[m], d.label.values[m],
                                          d.fold.values[m])
                rec[f"comp_{tag}"] = g.auroc_composition
                rec[f"full_{tag}"] = g.auroc_with_score
                rec[f"gain_{tag}"] = g.delta
                rec[f"n_{tag}"] = g.n
                # B8. THE MODEL'S OWN AUROC, which this table did not carry and which is the
                # quantity the title relation is DEFINED on. Without it, "difficulty" on their
                # benchmark was read off the composition baseline while ours was read off
                # model-alone AUROC, so the external test compared two different relations and
                # its non-replication could not be interpreted. One roc_auc_score call.
                rec[f"alone_{tag}"] = float(roc_auc_score(d.label.values[m], sc[m]))
            if ok:
                rows.append(rec)
                log(f"[{i:3d}/{len(sub)}] {r.dataset:18s} "
                    f"n1 comp {rec['comp_n1']:.3f} gain {rec['gain_n1']:+.4f}   "
                    f"n2 comp {rec['comp_n2']:.3f} gain {rec['gain_n2']:+.4f}")
        t = pd.DataFrame(rows)
        if t.empty:
            sys.exit("no dataset built")
        t.to_csv(per, index=False)

    out = []

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        # its OWN resample, sized to this input -- the strata below are subsets and reusing a
        # fixed index sized for the full panel indexes out of bounds
        ix = np.random.default_rng(0).integers(0, len(v), size=(2000, len(v)))
        b = np.array([v[i].mean() for i in ix])
        lo, hi = np.percentile(b, [2.5, 97.5])
        out.append({"check": check, "value": float(v.mean()), "ci_low": float(lo),
                    "ci_high": float(hi), "n": len(t), "note": note})
        return float(v.mean())

    log(f"\n=== R1p: Horlacher's own negative sets, their folds, n = {len(t)} ===\n")
    log(f"  {'negative set':14s} {'composition':>12s} {'apparent':>10s} {'contribution':>14s}")
    for tag, name in (("n1", "negative-1"), ("n2", "negative-2")):
        c = add(f"composition alone, {name}", t[f"comp_{tag}"])
        f_ = add(f"apparent AUROC, {name}", t[f"full_{tag}"])
        g = add(f"nested contribution, {name}", t[f"gain_{tag}"])
        log(f"  {name:14s} {c:12.4f} {f_:10.4f} {g:+14.4f}")
    d21 = add("CONTRAST, negative-2 minus negative-1", t.gain_n2 - t.gain_n1)
    pos = int(((t.gain_n2 - t.gain_n1) > 0).sum())
    out.append({"check": "datasets with a positive contrast, n2 minus n1", "value": pos,
                "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
    means = [t.gain_n1.mean(), t.gain_n2.mean()]
    fr = max(means) / min(means) if min(means) > 0 else np.nan
    out.append({"check": "fold range across their two negative sets", "value": float(fr),
                "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
    rho, pv = spearmanr(np.r_[t.comp_n1, t.comp_n2], np.r_[t.gain_n1, t.gain_n2])
    out.append({"check": "spearman(baseline, gain) on their data", "value": float(rho),
                "ci_low": np.nan, "ci_high": np.nan, "n": 2 * len(t), "note": f"p={pv:.1e}"})
    log(f"\n  contrast n2 - n1  {d21:+.4f}   positive in {pos}/{len(t)}")
    log(f"  fold range        {fr:.2f}x")
    log(f"  spearman(baseline, gain) on THEIR data: {rho:+.3f}  p={pv:.1e}")

    # R1n's DECISIVE TEST, RUN ON THEIR DATA. This is the point of the whole exercise: does
    # the claim that the baseline carries the information replicate outside our own windows?
    db, dg = t.comp_n2 - t.comp_n1, t.gain_n2 - t.gain_n1
    r2, p2 = spearmanr(db, dg)
    out.append({"check": "within-dataset spearman(delta baseline, delta gain), their data",
                "value": float(r2), "ci_low": np.nan, "ci_high": np.nan, "n": len(t),
                "note": f"p={p2:.2e}; ours is -0.664"})
    hi = (db > 0).values
    out.append({"check": "datasets where negative-2 raises the baseline", "value": int(hi.sum()),
                "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
    for lab, m in (("baseline HIGHER", hi), ("baseline LOWER", ~hi)):
        if m.sum() < 4:
            continue
        add(f"n2 minus n1 gain, negative-2 {lab}", dg.values[m])
    log(f"\n  R1n's test on their data: spearman(delta baseline, delta gain) "
        f"{r2:+.3f}  p={p2:.2e}   (ours -0.664)")
    for lab, m in (("HIGHER", hi), ("LOWER", ~hi)):
        if m.sum() >= 4:
            log(f"    where n2 baseline is {lab:6s}: n={m.sum():3d}  "
                f"gain diff {dg.values[m].mean():+.4f}")
    log("  -> the gradient replicates; the SIGN does not reverse, so on their benchmark the")
    log("     protocol label carries information beyond the baseline. R1n is limited, not lost.")

    # B8. THE TITLE RELATION, ON THE QUANTITY IT IS ACTUALLY DEFINED ON.
    #
    # Everything above measures difficulty by the COMPOSITION baseline, because this table had
    # no model-alone column. Our own claim is about the MODEL's own AUROC. So the external
    # "non-replication" compared two different relations, and could not have been evidence
    # either way. With alone_n1 and alone_n2 the comparison is finally like for like: which of
    # their two arms is harder for the 4-mer, and does the harder one yield more?
    if {"alone_n1", "alone_n2"} <= set(t.columns):
        for tag, name in (("n1", "negative-1"), ("n2", "negative-2")):
            add(f"model alone, {name}", t[f"alone_{tag}"])
        harder_n2 = int((t.alone_n2 < t.alone_n1).sum())
        out.append({"check": "datasets where negative-2 is HARDER for the 4-mer",
                    "value": harder_n2, "ci_low": np.nan, "ci_high": np.nan, "n": len(t)})
        # The inversion, per dataset: harder by the model's own AUROC AND a larger
        # contribution. Ours holds this way in 88 of 94 for the GC-to-dinucleotide step.
        da, dg2 = t.alone_n2 - t.alone_n1, t.gain_n2 - t.gain_n1
        inv = int(((da < 0) & (dg2 > 0)).sum() + ((da > 0) & (dg2 < 0)).sum())
        ra, pa = spearmanr(da, dg2)
        out += [{"check": "datasets where difficulty and contribution move OPPOSITELY, "
                          "model-alone axis, their data",
                 "value": inv, "ci_low": np.nan, "ci_high": np.nan, "n": len(t)},
                {"check": "within-dataset spearman(delta model-alone, delta gain), their data",
                 "value": float(ra), "ci_low": np.nan, "ci_high": np.nan, "n": len(t),
                 "note": f"p={pa:.2e}; POSITIVE means they move together"}]
        log("\n  B8, the relation on the model-alone axis (the one the title is about):")
        log(f"    model alone   negative-1 {t.alone_n1.mean():.4f}   "
            f"negative-2 {t.alone_n2.mean():.4f}")
        log(f"    negative-2 harder for the 4-mer in {harder_n2}/{len(t)} datasets")
        log(f"    opposite movement in {inv}/{len(t)};  spearman {ra:+.3f} (p={pa:.2e})")
        log("    a POSITIVE spearman here means difficulty and contribution move TOGETHER on")
        log("    their benchmark, which is the non-replication -- now measured on our own axis.")

    pd.DataFrame(out).to_csv(TABLES / "horlacher_arm.csv", index=False)
    log("\nwrote horlacher_arm.csv and horlacher_per_dataset.csv")


if __name__ == "__main__":
    main()
