"""Step 3: does the delta machinery detect a disruption we planted ourselves?

Every null result in this project depends on this passing. If we mutate the base that
matters in a protein's known motif and the delta score does NOT move, then the pipeline
cannot detect binding disruption at all, and "no signal in ClinVar" says nothing about
ClinVar -- it says the instrument is broken.

The design is a matched pair per window. Take a window containing the motif. Make one
mutant that changes the critical motif base, and one that makes the SAME substitution at
a position at least 25 nt away from every motif occurrence. Both mutants change exactly
one base, from the same reference base to the same alternate base. The only difference is
whether the change lands on the motif. So any difference in delta is attributable to the
motif and not to the mutation itself.

Motifs are literature-derived and listed in params.yaml, which is why this covers 9
proteins rather than all 131: inventing a motif for a protein whose specificity is not
established would make the control circular.

    python scripts/positive_control.py --k 4
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.eval import baseline  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.variants import positive_control as pc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
CELLS = ("K562", "HepG2")


def run_one(protein, cell, motif, offset, replacement, cfg, k):
    f = ROOT / "data/processed" / cell / protein / "dataset.tsv"
    if not f.exists():
        return None
    df = pd.read_csv(f, sep="\t")
    pos = df[df.label == 1]
    if len(pos) < 50:
        return None

    pairs, dropped = pc.build_pairs(
        pos, motif, offset, replacement,
        min_distance=cfg.positive_control["min_motif_distance"], seed=cfg.seed)
    if pairs.empty:
        return {"protein": protein, "cell": cell, "motif": motif, "n_pairs": 0,
                "note": f"no usable windows ({dropped})"}
    cap = cfg.positive_control["max_windows"]
    if pairs.id.nunique() > cap:                 # cap for runtime, both kinds kept paired
        keep = pairs.id.drop_duplicates().sample(cap, random_state=cfg.seed)
        pairs = pairs[pairs.id.isin(keep)]

    # Score with the fold models, so a window is never scored by a model that trained on
    # its own chromosome. Using one model fit on everything would let it recognise the
    # window rather than the motif.
    models, vec = baseline.fit_fold_models(df.seq_rna.tolist(), df.label.to_numpy(),
                                           df.fold.to_numpy(), k=k)
    fold_of = dict(zip(df.chrom, df.fold))
    pairs = pairs.assign(fold=pairs.chrom.map(fold_of))
    pairs = pairs.dropna(subset=["fold"])
    if pairs.empty:
        return None

    d = baseline.variant_delta(models, vec, pairs.ref_seq, pairs.alt_seq,
                              pairs.fold.to_numpy())
    scored = pairs.assign(delta=np.abs(d)).dropna(subset=["delta"])
    if scored.empty or scored.kind.nunique() < 2:
        return None

    es = pc.effect_size(scored)
    dis = scored.loc[scored.kind == "disruptive", "delta"]
    neu = scored.loc[scored.kind == "neutral", "delta"]
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(dis, neu, alternative="greater")
    return {"protein": protein, "cell": cell, "motif": motif, **es,
            "auroc_disruptive_vs_neutral": round(
                float(u.statistic / (len(dis) * len(neu))), 4),
            "p": float(u.pvalue), "note": ""}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--k", type=int, default=4)
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    motifs = pc.motifs_from_config(cfg)

    rows = []
    for protein, (motif, off, rep) in sorted(motifs.items()):
        for cell in CELLS:
            r = run_one(protein, cell, motif, off, rep, cfg, a.k)
            if r:
                rows.append(r)
                print(f"  {protein:9} {cell:6} {motif:6} "
                      f"d={r.get('cohens_d', float('nan')):+.3f} "
                      f"auroc={r.get('auroc_disruptive_vs_neutral', float('nan')):.3f} "
                      f"n={r.get('n_pairs', 0)}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(TABLES / "positive_control.csv", index=False)
    ok = res[res.note == ""]
    print(f"\n{'':=<70}")
    print(f"POSITIVE CONTROL: {len(ok)} protein-cell-line pairs tested")
    print(f"{'':=<70}")
    passed = ok[(ok.cohens_d > 0) & (ok.p < 0.05)]
    print(f"  disruptive delta > neutral delta at p<0.05: {len(passed)}/{len(ok)}")
    print(f"  median Cohen's d: {ok.cohens_d.median():+.3f}")
    print(f"  median AUROC (disruptive vs neutral): "
          f"{ok.auroc_disruptive_vs_neutral.median():.3f}")
    if len(passed) < len(ok):
        print(f"\n  FAILED: "
              f"{', '.join(ok[~ok.index.isin(passed.index)].protein + ':' + ok[~ok.index.isin(passed.index)].cell)}")
        print("  a failure here means the delta score cannot detect a disruption we "
              "planted, so any null result for that protein is uninterpretable")
    print("\nwrote results/tables/positive_control.csv")


if __name__ == "__main__":
    main()
