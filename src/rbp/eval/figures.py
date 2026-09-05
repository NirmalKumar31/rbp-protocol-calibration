"""EDA figures.

One palette, defined once. Categorical hues are assigned in a fixed order so a region
keeps its colour across every figure; magnitude uses a single hue ramp rather than a
rainbow. Grid and axes stay recessive so the data is what reads.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

INK = "#222222"
MUTED = "#666666"
GRID = "#E4E4E4"
POS = "#2E5C8A"
NEG = "#B4B4B4"
ACCENT = "#C1553B"
RAMP = plt.cm.Blues

# fixed order, so utr3 is always orange everywhere
REGION_COLOR = {
    "utr5": "#4C78A8", "utr3": "#F58518", "cds": "#54A24B",
    "exon_nc": "#EECA3B", "intron": "#B4B4B4",
}
REGION_ORDER = ("utr5", "utr3", "cds", "exon_nc", "intron")


def _style(ax, xlabel=None, ylabel=None, title=None):
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.grid(axis="y", visible=False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0, labelsize=9)
    if xlabel:
        ax.set_xlabel(xlabel, color=INK, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK, fontsize=10)
    if title:
        ax.set_title(title, color=INK, fontsize=12, fontweight="bold", loc="left", pad=12)
    return ax


def save(fig, path, dpi=150):
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def kmer_baseline(sweep, out, saturated=0.94):
    """Per-protein k-mer AUROC, with the headroom above each bar made explicit.

    The point of the figure is not the ranking but which proteins a deep model could
    still improve on, so the gap to 1.0 is drawn rather than left implicit.
    """
    d = sweep.sort_values("best_auroc")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(d) + 1.4))
    ax.barh(y, d.best_auroc, height=0.62,
            color=[NEG if v >= saturated else POS for v in d.best_auroc])
    ax.barh(y, 1 - d.best_auroc, left=d.best_auroc, height=0.62,
            color="#F2F2F2", zorder=0)
    for i, (v, k) in enumerate(zip(d.best_auroc, d.best_k)):
        ax.text(v + 0.006, i, f"{v:.3f}  (k={k})", va="center", fontsize=8.5, color=INK)
    ax.axvline(0.5, color=MUTED, lw=1, ls=":")
    ax.axvline(saturated, color=ACCENT, lw=1, ls="--")
    ax.text(saturated, len(d) - 0.2, "  effectively saturated", color=ACCENT, fontsize=8.5)
    ax.set_yticks(y, d.protein, fontsize=9)
    ax.set_xlim(0.5, 1.10)
    _style(ax, xlabel="test AUROC, logistic regression on k-mer counts",
           title="How far a bag of k-mers already gets")
    fig.text(0.0, -0.02, "grey = little headroom left   ·   blue = room for a better model",
             fontsize=8.5, color=MUTED)
    save(fig, out)


def motif_enrichment(df, out):
    d = df.sort_values("enrichment")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(6.8, 0.34 * len(d) + 1.4))
    ax.barh(y, d.enrichment, height=0.6,
            color=[ACCENT if v < 1.15 else POS for v in d.enrichment])
    ax.axvline(1.0, color=MUTED, lw=1)
    for i, v in enumerate(d.enrichment):
        ax.text(v + 0.08, i, f"{v:.2f}x", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(y, [f"{p}  {m}" for p, m in zip(d.protein, d.motif)], fontsize=9)
    _style(ax, xlabel="motif frequency in positives / matched negatives",
           title="Known motifs, against their own matched negatives")
    fig.text(0.0, -0.02, "red = no enrichment: the matched negatives absorbed the motif",
             fontsize=8.5, color=MUTED)
    save(fig, out)


def region_composition(summary, out):
    d = summary.set_index("protein")[[f"frac_{r}" for r in REGION_ORDER]]
    d = d.loc[d["frac_intron"].sort_values().index]
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.6, 0.34 * len(d) + 1.6))
    left = np.zeros(len(d))
    for r in REGION_ORDER:
        v = d[f"frac_{r}"].values
        ax.barh(y, v, left=left, height=0.66, color=REGION_COLOR[r], label=r,
                edgecolor="white", linewidth=0.8)
        left += v
    ax.set_yticks(y, d.index, fontsize=9)
    ax.set_xlim(0, 1)
    _style(ax, xlabel="fraction of positive windows",
           title="Where each protein binds")
    ax.legend(frameon=False, fontsize=8.5, ncol=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.09), labelcolor=INK)
    save(fig, out)


def gc_matching(datasets, out, examples=("TARDBP", "PUM2", "EWSR1")):
    have = [p for p in examples if p in datasets] or list(datasets)[:3]
    fig, axes = plt.subplots(1, len(have) + 1, figsize=(3.1 * (len(have) + 1), 3.0))
    for ax, p in zip(axes[:-1], have):
        df = datasets[p]
        bins = np.linspace(0.15, 0.85, 34)
        ax.hist(df[df.label == 1].gc, bins=bins, color=POS, alpha=0.85, label="bound")
        ax.hist(df[df.label == 0].gc, bins=bins, color=NEG, alpha=0.75, label="matched")
        _style(ax, xlabel="GC content", title=p)
        ax.grid(axis="x", visible=False)
        ax.grid(axis="y", color=GRID, lw=0.8)
        if ax is axes[0]:
            ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)
    gaps = np.concatenate([
        np.abs(d[d.label == 1].gc.values - d[d.label == 0].gc.values)
        for d in datasets.values()])
    ax = axes[-1]
    ax.hist(gaps, bins=40, color=POS)
    ax.axvline(0.05, color=ACCENT, lw=1, ls="--")
    ax.text(0.052, ax.get_ylim()[1] * 0.9, "tolerance", color=ACCENT, fontsize=8.5)
    _style(ax, xlabel="|GC gap| within a pair", title="all proteins pooled")
    ax.grid(axis="x", visible=False)
    ax.grid(axis="y", color=GRID, lw=0.8)
    fig.suptitle("Negatives are GC-matched to their positive", color=INK,
                 fontsize=12, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, out)


def split_proportions(summary, out, target=(0.64, 0.16, 0.20)):
    d = summary.set_index("protein")[["train_frac", "val_frac", "test_frac"]]
    y = np.arange(len(d))
    colors = (POS, "#7FA8CC", NEG)
    fig, ax = plt.subplots(figsize=(7.4, 0.34 * len(d) + 1.6))
    left = np.zeros(len(d))
    for (col, c, _t, name) in zip(d.columns, colors, target, ("train", "val", "test")):
        ax.barh(y, d[col], left=left, height=0.66, color=c, label=name,
                edgecolor="white", linewidth=0.8)
        left += d[col].values
    for t in np.cumsum(target)[:-1]:
        ax.axvline(t, color=ACCENT, lw=1, ls="--")
    ax.set_yticks(y, d.index, fontsize=9)
    ax.set_xlim(0, 1)
    _style(ax, xlabel="fraction of rows",
           title="Split proportions after optimising the chromosome assignment")
    ax.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.09), labelcolor=INK)
    fig.text(0.0, -0.04, "dashed = 64 / 16 / 20 target", fontsize=8.5, color=MUTED)
    save(fig, out)


def peak_widths(widths, window, out):
    d = widths.sort_values("median")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.0, 0.34 * len(d) + 1.4))
    ax.hlines(y, d.p10, d.p90, color=NEG, lw=4)
    ax.plot(d["median"], y, "o", color=POS, ms=6)
    ax.axvline(window, color=ACCENT, lw=1, ls="--")
    ax.text(window + 3, len(d) - 0.4, f"{window}-nt window", color=ACCENT, fontsize=8.5)
    for i, v in enumerate(d.frac_wider_than_window):
        ax.text(d.p90.iloc[i] + 6, i, f"{v:.0%} wider", va="center",
                fontsize=8, color=MUTED)
    ax.set_yticks(y, d.protein, fontsize=9)
    _style(ax, xlabel="peak width (nt): median with 10th-90th percentile",
           title="Peak widths against the fixed window")
    save(fig, out)


def cobinding_heatmap(matrix, out):
    m = matrix.copy()
    vals = np.array(m.to_numpy(dtype=float), copy=True)
    np.fill_diagonal(vals, np.nan)
    order = np.argsort(-np.nan_to_num(vals).sum(axis=1))
    names = [m.index[i] for i in order]
    vals = vals[np.ix_(order, order)]
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    im = ax.imshow(vals, cmap=RAMP, vmin=0, vmax=np.nanmax(vals))
    ax.set_xticks(range(len(names)), names, rotation=90, fontsize=8.5)
    ax.set_yticks(range(len(names)), names, fontsize=8.5)
    ax.tick_params(colors=MUTED, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.02)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, length=0, labelsize=8.5)
    cb.set_label("fraction of row's windows overlapping column", color=INK, fontsize=9)
    ax.set_title("Co-binding: which proteins share sites", color=INK, fontsize=12,
                 fontweight="bold", loc="left", pad=12)
    save(fig, out)


DIVERGE = plt.cm.RdBu_r   # log-ratios need a zero midpoint, not a single ramp


def positional_signal(profiles, out, window=101, highlight=("TARDBP", "ELAVL1", "FMR1")):
    """Per-position composition gap between positives and negatives.

    Flat means the discriminative signal is spread across the window, so the model has
    to be position invariant. A central spike would mean position itself is informative.
    """
    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    x = np.arange(window) - window // 2
    for p, v in profiles.items():
        hi = p in highlight
        ax.plot(x, v, lw=1.8 if hi else 0.8,
                color=POS if hi else "#D8D8D8", zorder=3 if hi else 1,
                label=p if hi else None)
    ax.axvline(0, color=ACCENT, lw=1, ls="--")
    _style(ax, xlabel="position relative to peak centre (nt)",
           ylabel="|positive - negative| base frequency",
           title="Where in the window the signal lives")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig.text(0.0, -0.02, "grey = the other proteins", fontsize=8.5, color=MUTED)
    save(fig, out)


def dinucleotide_heatmap(din, out):
    """log2(positive/negative) dinucleotide frequency: the residual composition shortcut."""
    cols = [c for c in din.columns if c.startswith("lr_")]
    m = din.set_index("protein")[cols]
    m.columns = [c[3:] for c in cols]
    m = m.loc[m.abs().max(axis=1).sort_values(ascending=False).index]
    v = np.nanmax(np.abs(m.to_numpy(dtype=float)))
    fig, ax = plt.subplots(figsize=(9.0, 0.34 * len(m) + 1.8))
    im = ax.imshow(m.to_numpy(dtype=float), cmap=DIVERGE, vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(m.columns)), m.columns, fontsize=8.5)
    ax.set_yticks(range(len(m)), m.index, fontsize=9)
    ax.tick_params(colors=MUTED, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, length=0, labelsize=8.5)
    cb.set_label("log2 positive / negative", color=INK, fontsize=9)
    ax.set_title("Dinucleotide composition still differs after GC matching",
                 color=INK, fontsize=12, fontweight="bold", loc="left", pad=12)
    save(fig, out)


def complexity(cx, out):
    """Homopolymer content and entropy, positives against matched negatives."""
    piv = cx.pivot(index="protein", columns="class",
                   values=["frac_with_run_ge", "entropy_mean"])
    piv.columns = [f"{a}_{b}" for a, b in piv.columns]
    piv["delta"] = piv.frac_with_run_ge_positive - piv.frac_with_run_ge_negative
    piv = piv.sort_values("delta")
    y = np.arange(len(piv))
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 0.34 * len(piv) + 1.6))

    ax = axes[0]
    ax.hlines(y, piv.frac_with_run_ge_negative, piv.frac_with_run_ge_positive,
              color=GRID, lw=2, zorder=1)
    ax.plot(piv.frac_with_run_ge_negative, y, "o", color=NEG, ms=6, label="matched")
    ax.plot(piv.frac_with_run_ge_positive, y, "o", color=POS, ms=6, label="bound")
    ax.set_yticks(y, piv.index, fontsize=9)
    _style(ax, xlabel="fraction of windows with a run of >=5",
           title="Homopolymer content")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="lower right")

    ax = axes[1]
    d = piv.entropy_mean_positive - piv.entropy_mean_negative
    ax.barh(y, d, height=0.6, color=[ACCENT if v < 0 else POS for v in d])
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_yticks(y, [""] * len(piv))
    _style(ax, xlabel="entropy(bound) - entropy(matched)",
           title="Sequence complexity")
    fig.text(0.0, -0.03,
             "red = bound windows are lower complexity than their matched negatives",
             fontsize=8.5, color=MUTED)
    fig.tight_layout()
    save(fig, out)


def gc_by_region(gcr, out):
    order = ("utr5", "cds", "intron", "exon_nc", "utr3")
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    for i, reg in enumerate(order):
        v = gcr[gcr.region == reg].gc_mean.values
        if not len(v):
            continue
        ax.scatter(np.full(len(v), i) + np.random.default_rng(0).normal(0, 0.05, len(v)),
                   v, s=26, color=REGION_COLOR[reg], alpha=0.85, zorder=3)
        ax.hlines(np.mean(v), i - 0.28, i + 0.28, color=INK, lw=2, zorder=4)
    ax.set_xticks(range(len(order)), order, fontsize=9)
    _style(ax, ylabel="mean GC of positive windows",
           title="GC differs systematically by region")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.grid(axis="x", visible=False)
    fig.text(0.0, -0.03, "each dot is one protein; bar is the cross-protein mean",
             fontsize=8.5, color=MUTED)
    save(fig, out)


def peak_quality(pq, out):
    d = pq.sort_values("signal_median")
    y = np.arange(len(d))
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 0.34 * len(d) + 1.5), sharey=True)
    ax = axes[0]
    ax.hlines(y, d.signal_median, d.signal_p90, color=GRID, lw=3)
    ax.plot(d.signal_median, y, "o", color=POS, ms=6)
    ax.plot(d.signal_p90, y, "o", color=NEG, ms=5)
    ax.set_yticks(y, d.protein, fontsize=9)
    _style(ax, xlabel="fold enrichment (median -> 90th pct)", title="Peak signal")
    ax = axes[1]
    ax.hlines(y, d.pval_median, d.pval_p90, color=GRID, lw=3)
    ax.plot(d.pval_median, y, "o", color=POS, ms=6)
    ax.plot(d.pval_p90, y, "o", color=NEG, ms=5)
    ax.set_xscale("log")
    _style(ax, xlabel="-log10 p (median -> 90th pct), log scale",
           title="Peak significance")
    fig.tight_layout()
    save(fig, out)


def chromosome_density(chrom, out):
    piv = chrom.pivot_table(index="chrom", values="per_mb", aggfunc="mean")
    order = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    piv = piv.reindex([c for c in order if c in piv.index])
    fig, ax = plt.subplots(figsize=(8.2, 3.4))
    ax.bar(range(len(piv)), piv.per_mb, color=POS, width=0.66)
    ax.set_xticks(range(len(piv)), [c[3:] for c in piv.index], fontsize=8.5)
    _style(ax, xlabel="chromosome", ylabel="peaks per Mb (mean over proteins)",
           title="Peak density tracks gene density")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.grid(axis="x", visible=False)
    save(fig, out)
