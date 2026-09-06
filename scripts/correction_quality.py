"""Was the composition correction clean for this dataset, or did it destroy real signal?

*** THIS SCRIPT'S OUTPUT IS NOT INTERPRETABLE. KEPT AS A RECORD OF A FAILED ATTEMPT. ***

The stratification variable it computes (locality retention) rests on a measure that does
not do what it claims. A constructed test in tests/unit/test_locality.py shows the locality
probe returns a LARGE effect (Cohen's d about 1.8) for a signal that is pure global
composition with no local feature at all. The reason is structural: a bag-of-k-mers model IS
local by construction -- it represents composition as weights on k-mers -- and the
disruptive site is chosen as the most enriched k-mer, i.e. precisely where the model's
weight sits. The comparison is biased by construction and cannot separate motif from
composition.

What it produced, for the record: the AUROC drop was essentially flat across strata
(K562 0.083 / 0.107 / 0.092 for clean / partial / over-corrected; HepG2 0.076 / 0.095 /
0.079) with corr(retention, drop) = -0.077 and +0.135. That is CONSISTENT with the drop
being confound removal rather than signal destruction, but it is equally consistent with
the stratification variable being noise, and there is no way to tell from here.

The question can only be answered with a model that could in principle be NON-local, i.e.
the CNN or the transformers. Until then the honest statement is the one measured on the 9
proteins with literature motifs: the correction over-corrects for repeat motifs (median
Cohen's d 1.567 -> 0.807, TARDBP 1.72 -> 0.34), so 0.10 is an UPPER BOUND on the
composition confound.

Dinucleotide-matched negatives drop the median model AUROC by 0.10. Reading that as
"0.10 of a reported AUROC is composition" requires knowing that the correction removed the
confound and not the signal -- and the literature-motif positive control says it sometimes
removes both (median Cohen's d 1.567 -> 0.807, TARDBP 1.72 -> 0.34).

The literature control covers 9 of 131 proteins. This runs the DATA-DRIVEN locality control
(eval/locality.py) on every dataset in both arms, so the AUROC drop can be reported
stratified by whether the correction was clean.

Validated against the literature control on the 9 proteins where both exist:
corr(literature d, locality d) = +0.809 in the GC arm and +0.963 in the matched arm, and
corr of the RETENTION ratio = +0.964. n is 9, so the standard error of those correlations is
about 0.41 -- they are consistent rather than conclusive, and that is stated in the output.

    python scripts/correction_quality.py --cell K562
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.eval import locality as loc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cell", default="K562")
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--kmer-k", type=int, default=5)
    a = p.parse_args()

    cmp_path = TABLES / f"compare_negatives_{a.cell}.csv"
    if not cmp_path.exists():
        raise SystemExit(f"need {cmp_path}; run compare_negatives.py first")
    cmp = pd.read_csv(cmp_path)

    rows = []
    for i, r in enumerate(cmp.itertuples(), 1):
        out = {"protein": r.protein, "cell": a.cell, "pairs": r.gc_pairs,
               "gc_auroc": r.gc_auroc, "dn_auroc": r.dn_auroc,
               "gc_comp": r.gc_composition, "dn_comp": r.dn_composition,
               "gc_gain": r.gc_gain, "dn_gain": r.dn_gain}
        ok = True
        for tag, root in (("gc", "data/processed"), ("dn", "data/processed_dinucmatch")):
            f = ROOT / root / a.cell / r.protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            res = loc.locality(pd.read_csv(f, sep="\t"), k=a.k, kmer_k=a.kmer_k)
            if res is None:
                ok = False
                break
            out[f"{tag}_locality"] = res["cohens_d"]
            out[f"{tag}_kmer"] = res["kmers"][0] if res["kmers"] else None
        if ok:
            rows.append(out)
        if i % 20 == 0:
            print(f"  [{i}/{len(cmp)}]", flush=True)

    t = pd.DataFrame(rows)
    t["auroc_drop"] = t.gc_auroc - t.dn_auroc
    t["comp_drop"] = t.gc_comp - t.dn_comp
    # Retention only means anything when the GC-arm control detected something to retain.
    t["retention"] = np.where(t.gc_locality > 0.5, t.dn_locality / t.gc_locality, np.nan)
    t.to_csv(TABLES / f"correction_quality_{a.cell}.csv", index=False)

    print(f"\n{'':=<76}")
    print(f"WAS THE CORRECTION CLEAN? {a.cell}, {len(t)} datasets")
    print(f"{'':=<76}\n")
    print("locality (Cohen's d, disruptive vs neutral mutation at the same base):")
    print(f"  GC-matched arm    median {t.gc_locality.median():+.3f}")
    print(f"  dinuc-matched arm median {t.dn_locality.median():+.3f}")
    print(f"  datasets with a detectable local signal (d>0.5) in the GC arm: "
          f"{int((t.gc_locality > 0.5).sum())}/{len(t)}")

    v = t.dropna(subset=["retention"])
    print(f"\nretention of local signal after correction, {len(v)} datasets with d>0.5:")
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        print(f"  p{int(q*100):02d}  {v.retention.quantile(q):.3f}")

    print("\n--- THE STRATIFIED RESULT ---")
    print(f"{'stratum':34} {'n':>4} {'AUROC drop':>11} {'comp drop':>10} {'gain after':>11}")
    bands = [("retention >= 0.8  (clean)", v.retention >= 0.8),
             ("retention 0.5-0.8 (partial)", (v.retention >= 0.5) & (v.retention < 0.8)),
             ("retention < 0.5   (over-corrected)", v.retention < 0.5)]
    for lab, mask in bands:
        s = v[mask]
        if not len(s):
            continue
        print(f"{lab:34} {len(s):4d} {s.auroc_drop.median():11.4f} "
              f"{s.comp_drop.median():10.4f} {s.dn_gain.median():11.4f}")
    print(f"{'ALL (with detectable local signal)':34} {len(v):4d} "
          f"{v.auroc_drop.median():11.4f} {v.comp_drop.median():10.4f} "
          f"{v.dn_gain.median():11.4f}")
    nodet = t[t.gc_locality <= 0.5]
    if len(nodet):
        print(f"{'no local signal to begin with':34} {len(nodet):4d} "
              f"{nodet.auroc_drop.median():11.4f} {nodet.comp_drop.median():10.4f} "
              f"{nodet.dn_gain.median():11.4f}")

    if len(v) > 8:
        print(f"\ncorr(retention, AUROC drop) = "
              f"{np.corrcoef(v.retention, v.auroc_drop)[0,1]:+.3f}")
        print("  positive would mean: datasets that KEPT their signal lost MORE AUROC,")
        print("  i.e. the drop is confound removal rather than signal destruction")
    print(f"\nwrote results/tables/correction_quality_{a.cell}.csv")


if __name__ == "__main__":
    main()
