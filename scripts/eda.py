"""Stage 4: dataset EDA, aimed at both the paper and model-design decisions.

    python scripts/eda.py --what fast      # summaries, integrity, motifs, widths, shift
    python scripts/eda.py --what kmer      # k-mer baselines (slower)
    python scripts/eda.py --what figures
    python scripts/eda.py --what all
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rbp.eval import eda  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results/tables"
FIG = ROOT / "results/figures"

# Literature motifs (core consensus). Used only as an external sanity check, never to
# build features, so it cannot leak into the models.
MOTIFS = {
    "RBFOX2": "GCAUG", "PUM2": "UGUA", "PUM1": "UGUA", "ELAVL1": "AUUUA",
    "QKI": "ACUAA", "TARDBP": "UGUGU", "PTBP1": "UCUU", "LIN28B": "GGAG",
    "TIA1": "UUUUU", "MBNL1": "UGCU", "MATR3": "CUCU", "FUS": "GUGGU",
}


def show(name, df, floatfmt="%.3f"):
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / f"{name}.csv", index=False)
    print(f"\n=== {name} ===")
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print(df.to_string(index=False))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--what", default="fast",
                   choices=["fast", "kmer", "descriptive", "figures", "all"])
    p.add_argument("--cell", default=None, help="default: encode.cell_line")
    p.add_argument("--arm", default=None, choices=sorted(panelmod.ARMS),
                   help="default: negatives.primary_arm")
    a = p.parse_args()
    cfg = cfgmod.load()
    # One cell line at a time: every table here is per-panel and mixing two cell lines
    # would average over a distinction the study is built on.
    cell = a.cell or cfg.encode["cell_line"]
    arm = panelmod.arm_of(cfg, a.arm)
    names = [n for n, pairs in panelmod.read_panel(cell, arm)
             if pairs >= cfg.cv["min_pairs"]]
    if not names:
        sys.exit(f"no panel for {cell} {arm}; run scripts/prepare.py first")
    print(f"loading {len(names)} datasets ({cell}, {arm} arm) ...", flush=True)
    ds = eda.load(panelmod.data_dir(cell, arm), names)
    print(f"loaded {sum(len(d) for d in ds.values()):,} rows")

    if a.what in ("fast", "all"):
        show("eda_panel_summary", eda.panel_summary(ds))
        show("eda_integrity", eda.integrity(ds))
        show("eda_motif_enrichment", eda.motif_enrichment(ds, MOTIFS))
        show("eda_nucleotide_composition", eda.nucleotide_composition(ds))
        paths = {n: next((ROOT / "data/raw/peaks" / cell).glob(f"{n}.*.bed.gz"))
                 for n in names}
        show("eda_peak_widths", eda.peak_width_stats(paths, cfg.windows["size"]))
        show("eda_split_shift", eda.split_shift(ds))
        show("eda_redundancy", eda.redundancy(ds))

    if a.what in ("descriptive", "all"):
        paths = {n: next((ROOT / "data/raw/peaks" / cell).glob(f"{n}.*.bed.gz"))
                 for n in names}
        show("eda_peak_quality", eda.peak_quality(paths))
        show("eda_complexity", eda.complexity(ds))
        show("eda_gc_by_region", eda.gc_by_region(ds))
        din = eda.dinucleotide_enrichment(ds)
        din.to_csv(TAB / "eda_dinucleotide_enrichment.csv", index=False)
        print("\n=== eda_dinucleotide_enrichment: strongest log2(pos/neg) per protein ===")
        cols = [c for c in din.columns if c.startswith("lr_")]
        for _, r in din.iterrows():
            v = r[cols].astype(float)
            top = v.abs().sort_values(ascending=False).index[:4]
            print(f"  {r.protein:9} " + "  ".join(f"{c[3:]} {v[c]:+.2f}" for c in top))
        chrom = eda.chromosome_distribution(
            paths, ROOT / "data/raw/GRCh38.primary_assembly.genome.fa.fai")
        chrom.to_csv(TAB / "eda_chromosome_distribution.csv", index=False)
        piv = chrom.pivot_table(index="chrom", values="per_mb", aggfunc="mean")
        print("\n=== peaks per Mb, averaged over proteins (top 8 / bottom 4) ===")
        piv = piv.sort_values("per_mb", ascending=False)
        print(piv.head(8).to_string())
        print("  ...")
        print(piv.tail(4).to_string())
        cb = TAB / "cobinding_matrix.csv"
        if cb.exists():
            cl, _ = eda.protein_clusters(pd.read_csv(cb, index_col=0))
            show("eda_protein_clusters", cl)

    if a.what in ("kmer", "all"):
        print("\nrunning k-mer baselines (this is the slow part) ...", flush=True)
        show("eda_kmer_sweep", eda.kmer_sweep(ds))
        rows = []
        for n in names:
            out = eda.kmer_baseline(ds[n], k=5)
            if out:
                rows.append({"protein": n, "auroc": out["auroc"],
                             "top_kmers": " ".join(out["top_positive"][:8])})
        show("eda_kmer_top_features", pd.DataFrame(rows))

    if a.what in ("figures", "all"):
        from rbp.eval import figures as fg
        FIG.mkdir(parents=True, exist_ok=True)
        summary = pd.read_csv(TAB / "eda_panel_summary.csv")
        fg.region_composition(summary, FIG / "eda_region_composition.png")
        fg.split_proportions(summary, FIG / "eda_split_proportions.png")
        fg.gc_matching(ds, FIG / "eda_gc_matching.png")
        fg.motif_enrichment(pd.read_csv(TAB / "eda_motif_enrichment.csv"),
                            FIG / "eda_motif_enrichment.png")
        fg.peak_widths(pd.read_csv(TAB / "eda_peak_widths.csv"),
                       cfg.windows["size"], FIG / "eda_peak_widths.png")
        for name, fn, src in [
            ("eda_dinucleotide_heatmap", fg.dinucleotide_heatmap, "eda_dinucleotide_enrichment"),
            ("eda_complexity", fg.complexity, "eda_complexity"),
            ("eda_gc_by_region", fg.gc_by_region, "eda_gc_by_region"),
            ("eda_peak_quality", fg.peak_quality, "eda_peak_quality"),
        ]:
            f = TAB / f"{src}.csv"
            if f.exists():
                fn(pd.read_csv(f), FIG / f"{name}.png")
        cd = TAB / "eda_chromosome_distribution.csv"
        if cd.exists():
            fg.chromosome_density(pd.read_csv(cd), FIG / "eda_chromosome_density.png")
        ps = eda.positional_signal(ds)
        fg.positional_signal(ps, FIG / "eda_positional_signal.png",
                             window=cfg.windows["size"])

        sw = TAB / "eda_kmer_sweep.csv"
        if sw.exists():
            fg.kmer_baseline(pd.read_csv(sw), FIG / "eda_kmer_baseline.png")
        cb = TAB / "cobinding_matrix.csv"
        if cb.exists():
            fg.cobinding_heatmap(pd.read_csv(cb, index_col=0),
                                 FIG / "eda_cobinding_heatmap.png")
        print(f"\nfigures -> {FIG.relative_to(ROOT)}")
        for f in sorted(FIG.glob("eda_*.png")):
            print(f"  {f.name:38} {f.stat().st_size/1e3:6.0f} KB")

    print(f"\ntables -> {TAB.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
