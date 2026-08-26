"""Step 3, variant arm: the whole downstream analysis with k-mer scores, no GPU.

Runs the chain that produces the paper's central claim -- ClinVar variants attached to
each protein's binding sites, scored for predicted disruption, tested against
pathogenicity with conservation controlled, pooled across the panel, and finally
correlated across proteins to ask whether binding-prediction quality transfers to
variant-scoring usefulness.

Every number here is a rehearsal, because the scores come from a k-mer model rather than
a trained network. The point is that the pipeline, the sample sizes and the power are
real, so we learn whether the analysis is even askable before paying for the models.

Stages, each cached to disk so a rerun is cheap:

    --what assign   variants near peaks, per dataset          -> variant_assignments.csv
    --what score    k-mer delta per variant                   -> variant_scores.csv
    --what phylop   conservation, over HTTP range requests    -> variant_conservation.csv
    --what test     conservation control + pooled + transfer  -> variant_results.csv
"""

import argparse
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from pyfaidx import Fasta  # noqa: E402

from rbp.data import annotation as ann  # noqa: E402
from rbp.data import encode  # noqa: E402
from rbp.eval import baseline  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402
from rbp.variants import assign, clinvar, phylop  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
INTERIM = ROOT / "data" / "interim"
CELLS = ("K562", "HepG2")


def panel(cfg, arm="dinuc"):
    """(protein, cell, pairs) for every dataset in the final panel of one arm."""
    rows = []
    for cell in CELLS:
        rows += [(prot, cell, pairs) for prot, pairs in panelmod.read_panel(cell, arm)]
    return rows


def fold_map(cfg):
    lines = (ROOT / cfg.cv["folds"]).read_text().strip().splitlines()[1:]
    return {c: int(k) for c, k in (ln.split("\t") for ln in lines)}


# ---------------------------------------------------------------------------------------

def stage_assign(cfg):
    """Which ClinVar variants sit near which protein's peaks."""
    v = cfg["variants"]
    fmap = fold_map(cfg)
    keep = set(ann.MAIN_CHROMS) - set(cfg.encode.get("exclude_chroms", []))

    print("loading ClinVar (strict Pathogenic/Benign SNVs, noncoding) ...", flush=True)
    variants = list(clinvar.load(ROOT / "data/raw/clinvar.vcf.gz",
                                 v["pathogenic_labels"], v["benign_labels"],
                                 noncoding=v["noncoding_consequences"],
                                 chroms=keep, snv_only=True))
    npath = sum(x["label"] for x in variants)
    print(f"  {len(variants):,} variants  ({npath:,} pathogenic, "
          f"{len(variants)-npath:,} benign)")

    rows = []
    ds = panel(cfg)
    for i, (prot, cell, pairs) in enumerate(ds, 1):
        idx = assign.peak_index(encode.peak_path(ROOT, prot, cell), keep)
        got = assign.assign(variants, idx, v["peak_margin"], fmap)
        for g in got:
            rows.append({"protein": prot, "cell": cell, "pairs": pairs,
                         **{k: g[k] for k in ("vid", "chrom", "pos", "pos_vcf", "ref",
                                              "alt", "label", "fold", "strand",
                                              "peak_distance")}})
        if i % 20 == 0 or i == len(ds):
            print(f"  [{i:3d}/{len(ds)}] {len(rows):,} assignments", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "variant_assignments.csv", index=False)
    report_assign(df)
    return df


def report_assign(df):
    per = (df.groupby(["protein", "cell"])
             .agg(n=("vid", "size"), n_path=("label", "sum")).reset_index())
    per["n_benign"] = per.n - per.n_path
    print(f"\n{'':=<70}")
    print(f"pooled: {len(df):,} variant-dataset pairs, "
          f"{df.label.sum():,} pathogenic, {int((df.label == 0).sum()):,} benign")
    print(f"distinct variants: {df.vid.nunique():,}")
    print(f"{'':=<70}")
    print(f"\nper-dataset pathogenic count (this is what killed the per-protein analysis):")
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        print(f"  p{int(q*100):02d}  {per.n_path.quantile(q):6.0f}")
    for t in (10, 20, 30, 50):
        print(f"  datasets with >={t:3d} pathogenic: {int((per.n_path >= t).sum()):3d}"
              f"/{len(per)}")
    per.to_csv(TABLES / "variant_availability_panel.csv", index=False)
    print("\nwrote results/tables/variant_availability_panel.csv")


# ---------------------------------------------------------------------------------------

def stage_score(cfg, k):
    """k-mer delta per variant, using only the fold model that never saw its chromosome."""
    a = pd.read_csv(TABLES / "variant_assignments.csv")
    fasta = Fasta(str(ROOT / "data/raw/GRCh38.primary_assembly.genome.fa"))
    size, shifts = cfg.windows["size"], cfg["variants"]["shifts"]
    how = cfg["variants"]["delta"]

    out, dropped_tot = [], {"ref_mismatch": 0, "no_usable_window": 0}
    groups = list(a.groupby(["protein", "cell"]))
    t0 = time.time()
    for i, ((prot, cell), g) in enumerate(groups, 1):
        ds = pd.read_csv(ROOT / "data/processed" / cell / prot / "dataset.tsv", sep="\t")
        models, vec = baseline.fit_fold_models(ds.seq_rna.tolist(), ds.label.to_numpy(),
                                               ds.fold.to_numpy(), k=k)
        recs = g.to_dict("records")
        table, dropped = assign.build_scoring_table(recs, fasta, size, shifts)
        for key, n in dropped.items():
            dropped_tot[key] += n
        if not table:
            continue
        t = pd.DataFrame(table)
        d = baseline.variant_delta(models, vec, t.seq_ref, t.seq_alt, t.fold.to_numpy())
        vids, deltas = assign.collapse_delta(t.vid.to_numpy(), d, how=how)
        lab = dict(zip(t.vid, t.label))
        fld = dict(zip(t.vid, t.fold))
        for vid, dv in zip(vids, deltas):
            out.append({"protein": prot, "cell": cell, "vid": vid, "label": lab[vid],
                        "fold": fld[vid], "delta": dv})
        if i % 20 == 0 or i == len(groups):
            el = time.time() - t0
            print(f"  [{i:3d}/{len(groups)}] {len(out):,} scored, {el:5.0f}s, "
                  f"~{el/i*(len(groups)-i):4.0f}s left", flush=True)

    df = pd.DataFrame(out)
    df.to_csv(TABLES / "variant_scores.csv", index=False)
    print(f"\nscored {len(df):,} variant-dataset pairs")
    print(f"dropped: {dropped_tot}")
    if dropped_tot["ref_mismatch"]:
        print("  ref_mismatch means ClinVar's REF disagreed with the genome; those are "
              "excluded rather than scored against a wrong window")
    return df


# ---------------------------------------------------------------------------------------

def stage_phylop(cfg):
    """Conservation for every distinct variant position, cached."""
    s = pd.read_csv(TABLES / "variant_scores.csv")
    uniq = s.vid.drop_duplicates().to_frame()
    parts = uniq.vid.str.split(":", expand=True)
    uniq["chrom"], uniq["pos_vcf"] = parts[0], parts[1].astype(int)
    print(f"{len(uniq):,} distinct positions to annotate", flush=True)

    cache = INTERIM / "phylop_cache.tsv"
    have = phylop.load_cache(cache)
    print(f"  cache holds {len(have):,}")
    ann_df = phylop.annotate(uniq, cache=str(cache))
    print(f"  coverage: {phylop.assert_coverage(ann_df)}")
    ann_df[["vid", "chrom", "pos_vcf", "conservation"]].to_csv(
        TABLES / "variant_conservation.csv", index=False)
    return ann_df


# ---------------------------------------------------------------------------------------

def stage_test(cfg, n_boot):
    """The conservation control, pooled and per dataset, then the transfer correlation."""
    from rbp.variants import conservation as cons

    s = pd.read_csv(TABLES / "variant_scores.csv")
    c = pd.read_csv(TABLES / "variant_conservation.csv")[["vid", "conservation"]]
    df = s.merge(c, on="vid", how="left").dropna(subset=["delta", "conservation"])
    df["dataset"] = df.protein + ":" + df.cell
    print(f"{len(df):,} scored variant-dataset pairs with conservation")
    print(f"  {int(df.label.sum()):,} pathogenic, {int((df.label==0).sum()):,} benign")

    print("\n--- POOLED across the whole panel ---", flush=True)
    pooled = cons.run(df, ["delta"], group_col=None, n_boot=n_boot,
                      method=cfg.conservation["method"])
    cols = ["n", "n_pathogenic", "conservation_auroc", "delta_auroc",
            "corr_delta_conservation", "alone_coef", "alone_ci_low", "alone_ci_high",
            "controlled_coef", "controlled_ci_low", "controlled_ci_high",
            "attenuation", "controlled_survives"]
    print(pooled[[c for c in cols if c in pooled]].T.to_string())

    print("\n--- PER DATASET ---", flush=True)
    per = cons.run(df, ["delta"], group_col="dataset", n_boot=n_boot,
                   method=cfg.conservation["method"])
    per.to_csv(TABLES / "variant_results.csv", index=False)
    ok = per[per.note == ""]
    print(f"{len(ok)} datasets testable; "
          f"{int(ok.controlled_survives.sum())} survive the conservation control, "
          f"{int(ok.get('controlled_survives_fdr', pd.Series(dtype=bool)).sum())} after FDR")
    pooled.to_csv(TABLES / "variant_results_pooled.csv", index=False)

    transfer(ok)
    return per


def transfer(per):
    """The novel question: does binding-prediction quality predict variant usefulness?

    Correlated ACROSS datasets, with log(pairs) partialled out. Dataset size predicts
    binding AUROC and could predict the variant effect too, in which case a raw
    correlation would be size driving both rather than binding quality driving variant
    utility. The partial correlation is the pre-registered primary.
    """
    b = ROOT / "results/tables/rehearsal_binding.csv"
    if not b.exists():
        print("\n(transfer: needs rehearsal_binding.csv; run scripts/rehearsal.py first)")
        return
    bind = pd.read_csv(b)[["dataset", "auroc", "delta_auroc", "pairs"]]
    m = per.rename(columns={"group": "dataset"}).merge(bind, on="dataset")
    m = m.dropna(subset=["controlled_coef", "auroc"])
    if len(m) < 10:
        print(f"\n(transfer: only {len(m)} datasets matched; too few)")
        return

    print(f"\n--- TRANSFER, {len(m)} datasets ---")
    x, y, lp = m.auroc.values, m.controlled_coef.values, np.log(m.pairs.values)
    r = np.corrcoef(x, y)[0, 1]

    def resid(v):
        A = np.column_stack([np.ones(len(lp)), lp])
        return v - A @ np.linalg.lstsq(A, v, rcond=None)[0]

    pr = np.corrcoef(resid(x), resid(y))[0, 1]
    from scipy.stats import spearmanr
    print(f"  binding AUROC vs controlled variant coefficient:")
    print(f"    Pearson  r = {r:+.3f}")
    print(f"    Spearman   = {spearmanr(x, y).statistic:+.3f}")
    print(f"    PARTIAL (log pairs out) = {pr:+.3f}   <- pre-registered primary")
    print(f"  size confound: corr(log pairs, binding AUROC) = "
          f"{np.corrcoef(lp, x)[0,1]:+.3f}, "
          f"corr(log pairs, variant coef) = {np.corrcoef(lp, y)[0,1]:+.3f}")
    m.to_csv(ROOT / "results/tables/transfer.csv", index=False)
    print("  wrote results/tables/transfer.csv")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--what", required=True,
                   choices=["assign", "score", "phylop", "test", "all"])
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--n-boot", type=int, default=500)
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    TABLES.mkdir(parents=True, exist_ok=True)
    INTERIM.mkdir(parents=True, exist_ok=True)

    if a.what in ("assign", "all"):
        stage_assign(cfg)
    if a.what in ("score", "all"):
        stage_score(cfg, a.k)
    if a.what in ("phylop", "all"):
        stage_phylop(cfg)
    if a.what in ("test", "all"):
        stage_test(cfg, a.n_boot)


if __name__ == "__main__":
    main()
