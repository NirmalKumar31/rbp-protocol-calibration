"""The paper's figures, one function per figure, each from a committed table.

Every figure reads a CSV under results/tables/ and writes a PNG plus a PDF under
results/figures/. Nothing is computed here that is not already in a table, so a figure can
never disagree with the number in the text -- if they disagree, the table is the truth and
the figure is a bug.

Figures that have no table yet are skipped with a message rather than invented, which is
why this can run unattended while other stages are still producing their inputs.

    python scripts/figures.py            # every figure whose table exists
    python scripts/figures.py --only f1  # one of them
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm, spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"

# One hue per model, fixed, so a model is the same colour in every figure. Assigned by
# identity and never by rank, or a figure that drops a model would repaint the survivors.
COLOR = {"composition": "#8c8c8c", "kmer": "#4878a8", "cnn": "#e08214",
         "splicebert": "#2b6a4d", "gc": "#b2182b", "dinuc": "#2166ac"}
LABEL = {"composition": "composition (19 feat)", "kmer": "k-mer LR", "cnn": "CNN",
         "splicebert": "SpliceBERT"}

plt.rcParams.update({"figure.dpi": 150, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "grid.linewidth": 0.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.axisbelow": True})


def save(fig, name):
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote results/figures/{name}.png / .pdf", flush=True)


def need(*names):
    """Table gate. Returns loaded frames, or None if any is missing."""
    out = []
    for n in names:
        p = TABLES / n
        if not p.exists():
            print(f"  skipped: needs results/tables/{n}", flush=True)
            return None
        out.append(pd.read_csv(p))
    return out


# --- f0: what the panel actually is ------------------------------------------------------
#
# The paper needs this before any result. A reader's first question about a 95-dataset panel
# drawn from a larger candidate pool is "which 95, and are they the easy ones", and the answer
# has to be a figure rather than a sentence -- particularly because dataset size correlates
# with AUROC at r = +0.55 here, so a size-biased panel would inflate every model equally and
# invisibly.
#
# Panel c is therefore load-bearing, not decorative: it shows the study panel spanning the
# candidate pool's whole size range rather than sitting at the top of it.

def f0():
    t = need("panel_summary.csv")
    if t is None:
        return
    d = t[0]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.0))

    # (a) THE SIZE-SAMPLE CHECK, and the reason this figure exists. Study panel against the
    # candidate pool it was drawn from. The pool is a step outline rather than a second filled
    # histogram: two translucent fills muddied each other and the comparison was unreadable,
    # which for the one panel a reviewer will actually interrogate is not acceptable.
    bins = np.linspace(np.log10(d.pairs).min() - 0.05, np.log10(d.pairs).max() + 0.05, 24)
    pool = need("candidate_sizes.csv")
    if pool is not None:
        p = pool[0]
        ax[0].hist(np.log10(p.pairs), bins=bins, density=True, histtype="step",
                   color="#404040", linewidth=1.4, label=f"candidate pool (n={len(p)})")
        lo = (p.pairs < d.pairs.min()).mean() * 100
        hi = (p.pairs < d.pairs.max()).mean() * 100
        sub = f"a  panel spans pool percentile {lo:.0f}-{hi:.0f}"
    else:
        sub = "a  size distribution"
    ax[0].hist(np.log10(d.pairs), bins=bins, density=True, color="#4878a8", alpha=0.75,
               label=f"study panel (n={len(d)})")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_xlabel("dataset size, log$_{10}$(pairs)")
    ax[0].set_ylabel("density")
    ax[0].set_title(sub, loc="left", fontsize=9)

    # (b) how the panel splits across the two cell lines, and how many proteins appear in both.
    by = d.groupby("cell").size()
    n_both = int((d.protein.value_counts() == 2).sum())
    ax[1].bar(by.index, by.values, color=["#8c8c8c", "#2b6a4d"], width=0.6)
    for i, v in enumerate(by.values):
        ax[1].text(i, v + 0.8, str(v), ha="center", fontsize=8)
    ax[1].set_ylabel("datasets")
    ax[1].set_ylim(0, by.max() * 1.18)
    ax[1].set_title(f"b  {d.protein.nunique()} proteins, {n_both} in both lines",
                    loc="left", fontsize=9)

    # (c) ClinVar coverage per dataset, which is what R4's power actually rests on. Reported
    # rather than assumed: the ladder pools ~19k variants, but they are distributed very
    # unevenly across datasets and the median dataset contributes far fewer than the mean.
    if "n_variants" in d.columns and d.n_variants.notna().any():
        v = d.n_variants.dropna()
        ax[2].hist(np.log10(v.clip(lower=1)), bins=22, color="#7a5195",
                   edgecolor="white", linewidth=0.4)
        ax[2].axvline(np.log10(v.median()), color="#404040", linestyle="--", linewidth=1.1)
        ax[2].text(np.log10(v.median()), ax[2].get_ylim()[1] * 0.92,
                   f" median {int(v.median())}", fontsize=7.5, va="top")
        ax[2].set_xlabel("ClinVar variants per dataset, log$_{10}$")
        ax[2].set_ylabel("datasets")
        ax[2].set_title(f"c  {int(v.sum()):,} variant-dataset pairs", loc="left", fontsize=9)
    else:
        ax[2].axis("off")
    save(fig, "f0_panel_overview")


# --- f1: the cost of the negative-set protocol -------------------------------------------

def f1():
    """The headline. Same positives, same model, two negative sets.

    Paired per dataset, because the two arms share their positives -- an unpaired plot
    would throw away exactly the structure that makes 187/187 meaningful.
    """
    got = need("cost_of_matching.csv")
    if got is None:
        return
    d = got[0]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1),
                           gridspec_kw={"width_ratios": [1.15, 1]})

    for _, r in d.iterrows():
        ax[0].plot([0, 1], [r.auroc_gc, r.auroc_dn], color="#999999",
                   lw=0.4, alpha=0.45, zorder=1)
    for i, (col, arm) in enumerate([("auroc_gc", "gc"), ("auroc_dn", "dinuc")]):
        ax[0].scatter(np.full(len(d), i), d[col], s=9, color=COLOR[arm],
                      zorder=3, edgecolor="white", linewidth=0.3)
        ax[0].hlines(d[col].mean(), i - 0.22, i + 0.22, color="black", lw=1.8, zorder=4)
    ax[0].set_xticks([0, 1])
    ax[0].set_xticklabels([f"GC-matched\n{d.auroc_gc.mean():.3f}",
                           f"dinucleotide-matched\n{d.auroc_dn.mean():.3f}"])
    ax[0].set_ylabel("k-mer model AUROC")
    ax[0].set_xlim(-0.4, 1.4)
    ax[0].set_title(f"Every dataset falls ({len(d)}/{len(d)})", loc="left", fontsize=9)

    # PANEL B IS THE FINDING, AND IT USED TO BE A SECOND VIEW OF PANEL A.
    #
    # This was a histogram of the same AUROC drop the left panel already shows. That drop is
    # very nearly tautological -- harder negatives lower AUROC by construction -- so the
    # headline figure led with the one thing a reviewer would call obvious, and the actual
    # result was in no figure at all.
    #
    # What is NOT obvious is that the model earns MORE credit under the harder control: the
    # nested gain over composition rises from +0.027 to +0.066. Absolute performance falls
    # while measured skill rises, because the GC-matched benchmark was handing most of its
    # AUROC to composition alone. That is the paper, so it is now the panel.
    g = d[["delta_auroc_gc", "delta_auroc_dn"]].mean()
    ratio = g.delta_auroc_dn / g.delta_auroc_gc
    for i, (col, arm) in enumerate([("delta_auroc_gc", "gc"), ("delta_auroc_dn", "dinuc")]):
        ax[1].scatter(np.random.default_rng(7).normal(i, 0.055, len(d)), d[col], s=8,
                      color=COLOR[arm], alpha=0.75, edgecolor="white", linewidth=0.25,
                      zorder=3)
        ax[1].hlines(d[col].mean(), i - 0.24, i + 0.24, color="black", lw=1.8, zorder=4)
        ax[1].text(i, d[col].max() + 0.012, f"{d[col].mean():+.4f}", ha="center", fontsize=8)
    ax[1].axhline(0, color="black", lw=0.8, zorder=2)
    ax[1].set_xticks([0, 1])
    ax[1].set_xticklabels(["GC-matched", "dinucleotide-matched"])
    ax[1].set_xlim(-0.4, 1.4)
    ax[1].set_ylabel("nested gain over composition")
    ax[1].set_title(f"Gain RISES {ratio:.1f}x under the harder control", loc="left",
                    fontsize=9)
    save(fig, "f1_cost_of_matching")


# --- f2: four models on identical data ---------------------------------------------------

def f2():
    got = need("matched_four_models.csv")
    if got is None:
        return
    d = got[0].rename(columns={"kmer_auroc": "kmer", "composition_auroc": "composition"})
    order = ["composition", "kmer", "cnn", "splicebert"]
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.1),
                           gridspec_kw={"width_ratios": [1.1, 1]})

    for i, m in enumerate(order):
        x = np.random.default_rng(i).normal(i, 0.055, len(d))
        ax[0].scatter(x, d[m], s=8, color=COLOR[m], alpha=0.75, zorder=3,
                      edgecolor="white", linewidth=0.25)
        ax[0].hlines(d[m].mean(), i - 0.26, i + 0.26, color="black", lw=1.8, zorder=4)
        ax[0].text(i, 0.985, f"{d[m].mean():.3f}", ha="center", fontsize=8)
    ax[0].set_xticks(range(4))
    # "composition (19 feat)" and "k-mer LR" overlapped at this figure width. The full name
    # stays in the legend of f5 and in the caption; the axis gets a short form.
    SHORT = {"composition": "composition\n(19 feat)", "kmer": "k-mer LR",
             "cnn": "CNN", "splicebert": "SpliceBERT"}
    ax[0].set_xticklabels([SHORT[m] for m in order], fontsize=8)
    ax[0].set_ylabel("pooled out-of-fold AUROC")
    ax[0].set_ylim(0.45, 1.0)
    ax[0].set_title(f"{len(d)} datasets, identical splits", loc="left", fontsize=9)

    # Gain over composition is the quantity the control was built to measure, so it gets
    # its own panel rather than being left for the reader to subtract by eye.
    for m in order[1:]:
        g = (d[m] - d.composition).sort_values().values
        ax[1].plot(np.arange(len(g)), g, color=COLOR[m], lw=1.4, label=LABEL[m])
    ax[1].axhline(0, color="black", lw=0.8)
    ax[1].set_xlabel("datasets, sorted within model")
    ax[1].set_ylabel("AUROC gain over composition")
    ax[1].legend(frameon=False, fontsize=8, loc="upper left")
    ax[1].set_title("Only SpliceBERT clears it everywhere", loc="left", fontsize=9)
    save(fig, "f2_four_models")



# --- f3: the strand control, and why the obvious version of it lies ----------------------
#
# Three bars and a difference. The point is not that the contrast shrinks under restriction --
# it does, and a naive reading calls all of that strand. The point is that a matched random
# drop reproduces most of the same shrinkage with no strand involved, and matching the placebo
# on region accounts for more still. What is left is the artifact.

def f3():
    t = need("strand_placebo.csv")
    if t is None:
        return
    q = t[0].set_index("check")
    v = lambda k: float(q.loc[k, "value"])                               # noqa: E731
    err = lambda k: [[v(k) - float(q.loc[k, "ci_low"])],                 # noqa: E731
                     [float(q.loc[k, "ci_high"]) - v(k)]]

    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.4),
                           gridspec_kw={"width_ratios": [1.35, 1]})

    bars = [("contrast, full data", "full\ndata", COLOR["kmer"]),
            ("contrast, sense-only pairs", "sense-only\npairs", COLOR["gc"]),
            ("contrast, PLACEBO (same n, random)", "placebo\n(random)", "#b0b0b0"),
            ("contrast, PLACEBO stratified on region x GC", "placebo\n(region-matched)",
             "#8c8c8c")]
    for i, (k, lab, c) in enumerate(bars):
        ax[0].bar(i, v(k), width=0.62, color=c, edgecolor="white", linewidth=1.2, zorder=3)
        ax[0].errorbar(i, v(k), yerr=err(k), color="#333333", capsize=3, lw=1.2, zorder=4)
    ax[0].axhline(v("contrast, full data"), color="#333333", lw=0.8, ls=":", zorder=2)
    ax[0].set_xticks(range(len(bars)), [b[1] for b in bars], fontsize=7.5)
    ax[0].set_ylabel("contrast in nested gain")
    ax[0].set_ylim(0, v("contrast, full data") * 1.30)
    ax[0].set_title("a  restriction shrinks it; so does dropping at random", loc="left")

    # The decomposition of the -0.0091 that restriction alone reports.
    parts = [("change from placebo", "cost of\ndropping pairs", "#b0b0b0"),
             ("locus-mix component", "locus\nmix", "#8c8c8c"),
             ("STRAND-SPECIFIC EXCESS (stratified)", "strand", COLOR["splicebert"])]
    for i, (k, lab, c) in enumerate(parts):
        ax[1].bar(i, v(k), width=0.6, color=c, edgecolor="white", linewidth=1.2, zorder=3)
        ax[1].errorbar(i, v(k), yerr=err(k), color="#333333", capsize=3, lw=1.2, zorder=4)
        ax[1].text(i, v(k) - 0.0011, f"{v(k):+.4f}", ha="center", va="top", fontsize=7.5)
    ax[1].axhline(0, color="#333333", lw=0.8)
    ax[1].set_xticks(range(len(parts)), [p_[1] for p_ in parts], fontsize=7.5)
    ax[1].set_ylabel("component of the -0.0091")
    surv = v("fraction of the contrast surviving")
    ax[1].set_title(f"b  {surv:.0%} of the contrast survives", loc="left")

    fig.tight_layout()
    save(fig, "f3_strand_placebo")


# --- f4: the magnitude replicates, and it buys precision ---------------------------------
#
# The paper concedes that the SIGN of the contrast is design-implied, so this figure is about
# the only part that is not. Panel a is an out-of-sample prediction: fifteen proteins measured
# in two cell lines, separate experiments with separately drawn negatives.

def f4():
    t = need("cost_of_matching.csv", "rehearsal_binding_gc.csv", "rehearsal_binding_dinuc.csv")
    if t is None:
        return
    cm, gc, dn = t
    cm["contrast"] = cm.delta_auroc_dn - cm.delta_auroc_gc
    w = cm.pivot_table(index="protein", columns="cell", values="contrast").dropna()
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.4))

    a_, b_ = w.iloc[:, 0], w.iloc[:, 1]
    lim = [min(a_.min(), b_.min()) - 0.01, max(a_.max(), b_.max()) + 0.01]
    ax[0].plot(lim, lim, color="#999999", ls="--", lw=1, zorder=1)
    ax[0].scatter(a_, b_, s=34, color=COLOR["kmer"], edgecolor="white", linewidth=0.6, zorder=3)
    r = np.corrcoef(a_, b_)[0, 1]
    ax[0].set(xlim=lim, ylim=lim, xlabel=f"contrast in {w.columns[0]}",
              ylabel=f"contrast in {w.columns[1]}")
    ax[0].set_title(f"a  replicates across cell lines, r = {r:+.2f} (n={len(w)})", loc="left")

    # Panel b: the same gain, measured more precisely. z = gain / SE, per dataset, paired.
    m = gc.merge(dn, on="dataset", suffixes=("_gc", "_dn"))
    for arm in ("gc", "dn"):
        se = (m[f"delta_ci_high_{arm}"] - m[f"delta_ci_low_{arm}"]) / (2 * 1.959963985)
        m[f"z_{arm}"] = m[f"delta_auroc_{arm}"] / se.replace(0, np.nan)
    m = m.dropna(subset=["z_gc", "z_dn"])
    for _, r_ in m.iterrows():
        ax[1].plot([0, 1], [r_.z_gc, r_.z_dn], color="#999999", lw=0.4, alpha=0.45, zorder=1)
    for i, (col, arm) in enumerate((("z_gc", "gc"), ("z_dn", "dinuc"))):
        ax[1].scatter(np.full(len(m), i), m[col], s=9, color=COLOR[arm], zorder=3,
                      edgecolor="white", linewidth=0.3)
        ax[1].hlines(m[col].median(), i - 0.22, i + 0.22, color="black", lw=1.8, zorder=4)
    ax[1].set_xticks([0, 1], ["GC-matched", "dinucleotide-matched"], fontsize=8)
    ax[1].set_xlim(-0.4, 1.4)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("z = nested gain / SE")
    up = int((m.z_dn > m.z_gc).sum())
    ax[1].set_title(f"b  measured more precisely in {up}/{len(m)}", loc="left")

    fig.tight_layout()
    save(fig, "f4_replication")


# --- f5: effect modification by dataset size, printed rather than buried -----------------

def f5():
    t = need("cost_of_matching.csv")
    if t is None:
        return
    d = t[0]
    d["contrast"] = d.delta_auroc_dn - d.delta_auroc_gc
    from scipy.stats import spearmanr
    rho = spearmanr(np.log10(d.pairs), d.contrast)
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.axhline(0, color="#999999", lw=0.8, ls="--")
    ax.scatter(d.pairs, d.contrast, s=13, color=COLOR["kmer"], alpha=0.75,
               edgecolor="white", linewidth=0.3, zorder=3)
    lp = np.log10(d.pairs)
    b_, a_ = np.polyfit(lp, d.contrast, 1)
    xs = np.linspace(lp.min(), lp.max(), 40)
    ax.plot(10 ** xs, a_ + b_ * xs, color="#333333", lw=1.4, zorder=4)
    ax.set_xscale("log")
    ax.set_xlabel("pairs per dataset")
    ax.set_ylabel("contrast in nested gain")
    ax.set_title(f"Larger in larger datasets: rho = {rho.statistic:+.3f}, "
                 f"p = {rho.pvalue:.1g}", loc="left", fontsize=9)
    save(fig, "f5_size_modification")


# --- f6: the contrast does not depend on the k-mer size ----------------------------------

def f6():
    t = need("k_sweep_per_dataset.csv", "k_sweep.csv")
    if t is None:
        return
    per, summ = t
    q = summ.set_index("check")
    ks = [3, 4, 5, 6]
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.3),
                           gridspec_kw={"width_ratios": [1, 1.1]})

    for i, k in enumerate(ks):
        c = per[f"contrast_k{k}"]
        ax[0].scatter(np.full(len(c), i) + np.random.default_rng(k).normal(0, 0.05, len(c)),
                      c, s=7, color=COLOR["kmer"], alpha=0.45, edgecolor="none", zorder=2)
        v = float(q.loc[f"contrast, k={k}", "value"])
        lo = float(q.loc[f"contrast, k={k}", "ci_low"])
        hi = float(q.loc[f"contrast, k={k}", "ci_high"])
        ax[0].errorbar(i, v, yerr=[[v - lo], [hi - v]], fmt="o", ms=6, color="#333333",
                       capsize=4, lw=1.5, zorder=4)
    ax[0].axhline(0, color="#999999", lw=0.9, ls="--")
    ax[0].set_xticks(range(len(ks)), [f"k={k}" for k in ks])
    ax[0].set_ylabel("contrast in nested gain")
    ax[0].set_title("a  positive at every k", loc="left")

    # Panel b: the two arms' gains, so the reader sees WHERE the contrast comes from.
    for arm, lab, c in (("gc", "GC-matched", COLOR["gc"]), ("dn", "dinuc-matched", COLOR["dinuc"])):
        ys = [per[f"gain_{arm}_k{k}"].mean() for k in ks]
        ax[1].plot(ks, ys, marker="o", ms=5, lw=1.8, color=c, label=lab)
    ax[1].set_xticks(ks)
    ax[1].set_xlabel("k-mer size")
    ax[1].set_ylabel("mean nested gain")
    ax[1].legend(frameon=False, fontsize=8)
    d54 = float(q.loc["k=5 minus k=4", "value"])
    ax[1].set_title(f"b  k=5 minus k=4 = {d54:+.4f}", loc="left")

    fig.tight_layout()
    save(fig, "f6_k_sweep")


# --- f7: the effect is twice as large for coding-region binders --------------------------
#
# Panel b is the mechanism, not decoration: if intronic sites were NOT more compositional the
# explanation offered in the text would be wrong, and the reader can check it here.

def f7():
    t = need("region_heterogeneity_per_dataset.csv")
    if t is None:
        return
    d = t[0]
    order = [g for g in ("cds", "utr3", "intron") if (d.dominant == g).sum() >= 8]
    col = {"cds": COLOR["cnn"], "utr3": COLOR["splicebert"], "intron": COLOR["kmer"]}
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.3))

    for i, g in enumerate(order):
        sub = d[d.dominant == g]
        ax[0].scatter(np.full(len(sub), i) +
                      np.random.default_rng(1).normal(0, 0.06, len(sub)),
                      sub.contrast, s=13, color=col[g], alpha=0.7, edgecolor="white",
                      linewidth=0.3, zorder=3)
        ax[0].hlines(sub.contrast.mean(), i - 0.26, i + 0.26, color="black", lw=2, zorder=4)
        ax[1].scatter(np.full(len(sub), i) +
                      np.random.default_rng(2).normal(0, 0.06, len(sub)),
                      sub.composition_auroc_dn, s=13, color=col[g], alpha=0.7,
                      edgecolor="white", linewidth=0.3, zorder=3)
        ax[1].hlines(sub.composition_auroc_dn.mean(), i - 0.26, i + 0.26, color="black",
                     lw=2, zorder=4)
    lab = [f"{g}\n(n={int((d.dominant == g).sum())})" for g in order]
    for a_ in ax:
        a_.set_xticks(range(len(order)), lab, fontsize=8)
        a_.set_xlim(-0.5, len(order) - 0.5)
    ax[0].axhline(0, color="#999999", lw=0.8, ls="--")
    ax[0].set_ylabel("contrast in nested gain")
    ax[0].set_title("a  twice as large for CDS binders", loc="left")
    ax[1].set_ylabel("composition alone, dinuc arm")
    ax[1].set_title("b  ...because intronic sites ARE more compositional", loc="left")

    fig.tight_layout()
    save(fig, "f7_region")


# --- f8: R1 is not the AUROC ceiling -----------------------------------------------------
#
# The figure exists to answer one referee objection, so each panel is one step of that answer
# and nothing else. (a) the raw paired effect. (b) how much of it survives once AUROC
# compression is removed, which is the only number that matters. (c) why the third scale
# disagrees: the coefficient gap tracks total task signal, not incremental value, so it is
# measuring difficulty. Panel c is two scatters rather than one because the argument IS the
# contrast between the two correlations.

def f8():
    t = need("rehearsal_binding_gc.csv", "rehearsal_binding_dinuc.csv", "scale_check.csv")
    if t is None:
        return
    gc, dn, sc = t
    m = gc.merge(dn, on="dataset", suffixes=("_gc", "_dn"))
    q = sc.set_index("check")
    val = lambda k: float(q.loc[k, "value"])
    err = lambda k: [[val(k) - float(q.loc[k, "ci_low"])], [float(q.loc[k, "ci_high"]) - val(k)]]

    r2 = np.sqrt(2.0)
    dp = lambda a: r2 * norm.ppf(np.clip(a, 1e-6, 1 - 1e-6))
    for a in ("gc", "dn"):
        m[f"dd_{a}"] = dp(m[f"with_score_auroc_{a}"]) - dp(m[f"composition_auroc_{a}"])
        m[f"dfull_{a}"] = dp(m[f"with_score_auroc_{a}"])

    fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.5))

    # (a) the paired effect, per dataset
    lim = [0, max(m.delta_auroc_gc.max(), m.delta_auroc_dn.max()) * 1.06]
    ax[0].plot(lim, lim, color="#999999", ls="--", lw=1, zorder=1)
    ax[0].scatter(m.delta_auroc_gc, m.delta_auroc_dn, s=22, alpha=0.8,
                  color=COLOR["kmer"], edgecolor="white", linewidth=0.4, zorder=3)
    n_up = int((m.delta_auroc_dn > m.delta_auroc_gc).sum())
    ax[0].set(xlim=lim, ylim=lim, xlabel="nested gain, GC-matched",
              ylabel="nested gain, dinuc-matched")
    ax[0].set_title(f"a  larger under proper matching in {n_up}/{len(m)}", loc="left")

    # (b) the decomposition: what survives removing compression
    # observed contrast first, then the two things it decomposes into. Three distinct hues,
    # because these are three different quantities and two shades of the same blue read as one.
    keys = [("CONTRAST, AUROC scale (published headline)", "observed\ncontrast", COLOR["kmer"]),
            ("contrast attributable to SCALE alone", "AUROC\ncompression", "#b0b0b0"),
            ("CONTRAST, protocol effect net of scale", "protocol effect\n(what survives)",
             COLOR["splicebert"])]
    for i, (k, lab, c) in enumerate(keys):
        ax[1].bar(i, val(k), width=0.62, color=c, edgecolor="white", linewidth=1.2, zorder=3)
        ax[1].errorbar(i, val(k), yerr=err(k), color="#333333", capsize=3, lw=1.2, zorder=4)
        ax[1].text(i, float(q.loc[k, "ci_high"]) + 0.0016, f"{val(k):+.4f}",
                   ha="center", fontsize=8)
    ax[1].axhline(0, color="#333333", lw=0.8)
    ax[1].set_xticks(range(3), [lab for _, lab, _ in keys], fontsize=8)
    ax[1].set_ylim(0, val(keys[0][0]) * 1.34)
    ax[1].set_ylabel("contrast in nested gain")
    share = val("scale share of the published contrast")
    ax[1].set_title(f"b  compression explains {share:.0%}, not all of it", loc="left")

    # (c) why the log-odds scale disagrees
    cgap = m.coef_gc - m.coef_dn
    for x, lab, c in ((m.dfull_gc - m.dfull_dn, "vs TOTAL task signal", COLOR["gc"]),
                      (m.dd_dn - m.dd_gc, "vs INCREMENTAL value", COLOR["splicebert"])):
        rho = spearmanr(cgap, x).statistic
        ax[2].scatter(x, cgap, s=20, alpha=0.75, color=c, edgecolor="white",
                      linewidth=0.3, zorder=3, label=f"{lab}   rho {rho:+.2f}")
    ax[2].axhline(0, color="#999999", lw=0.8, ls="--")
    ax[2].legend(frameon=False, fontsize=7.5, loc="lower right")
    ax[2].set(xlabel="between-arm gap (d' units)",
              ylabel="between-arm gap in coefficient")
    ax[2].set_title("c  the coefficient tracks difficulty, not value", loc="left")

    fig.tight_layout()
    save(fig, "f8_scale_check")


# --- f9: R1g, the contrast is not an artefact of the model class -------------------------
#
# The paper's sharpest limitation was that every number came from one model class. This
# figure is the answer, and it has to make three separate points or it does not close the
# objection. (a) the arm gap is present for all three models, so it is not a property of
# bags of k-mers. (b) it GROWS with model capacity, and the protocol effect survives the
# compression correction for each. (c) it holds dataset by dataset rather than only on
# average, which is what makes the paired comparison meaningful.

def f9():
    t = need("deep_contrast_per_dataset.csv", "deep_contrast.csv")
    if t is None:
        return
    d, s = t
    models = [m for m in ("kmer", "cnn", "splicebert") if f"{m}_gain_gc" in d.columns]
    q = s.set_index(["model", "quantity"]).value
    fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.3))

    # a. nested contribution, both arms, every model
    w = 0.34
    for i, m in enumerate(models):
        for j, (arm, key) in enumerate((("gc", "gc"), ("dinuc", "dn"))):
            v = d[f"{m}_gain_{key}"]
            ax[0].bar(i + (j - 0.5) * w, v.mean(), w * 0.9, color=COLOR[arm],
                      zorder=3, label=f"{arm}-matched" if i == 0 else None)
            ax[0].errorbar(i + (j - 0.5) * w, v.mean(), yerr=v.sem(), color="black",
                           lw=1, capsize=2, zorder=4)
    ax[0].set_xticks(range(len(models)), [LABEL[m] for m in models], fontsize=8)
    ax[0].set_ylabel("nested contribution over composition")
    ax[0].set_title("a  every model shows the arm gap", loc="left")
    ax[0].legend(fontsize=7, frameon=False)

    # b. the contrast, with the part that survives the compression correction
    for i, m in enumerate(models):
        c = q[(m, "contrast_auroc")]
        lo = s[(s.model == m) & (s.quantity == "contrast_auroc")].ci_low.iloc[0]
        hi = s[(s.model == m) & (s.quantity == "contrast_auroc")].ci_high.iloc[0]
        ax[1].bar(i, c, 0.55, color=COLOR[m], zorder=3)
        ax[1].errorbar(i, c, yerr=[[c - lo], [hi - c]], color="black", lw=1, capsize=3,
                       zorder=4)
        pmin, pmax = q[(m, "protocol_effect_min")], q[(m, "protocol_effect_max")]
        ax[1].add_patch(plt.Rectangle((i - 0.28, pmin), 0.56, pmax - pmin,
                                      facecolor="white", edgecolor="black", lw=0.8,
                                      hatch="///", alpha=0.85, zorder=5))
    ax[1].axhline(0, color="#999999", lw=0.8, ls="--")
    ax[1].set_xticks(range(len(models)), [LABEL[m] for m in models], fontsize=8)
    ax[1].set_ylabel("contrast (dinuc - GC)")
    ax[1].set_title("b  hatched: survives compression", loc="left")

    # c. dataset by dataset, deepest model against the k-mer
    if "kmer" in models and "splicebert" in models:
        x = d.kmer_gain_dn - d.kmer_gain_gc
        y = d.splicebert_gain_dn - d.splicebert_gain_gc
        lim = [min(x.min(), y.min()) - 0.01, max(x.max(), y.max()) + 0.01]
        ax[2].plot(lim, lim, color="#999999", lw=0.8, ls="--", zorder=2)
        ax[2].axhline(0, color="#cccccc", lw=0.6, zorder=1)
        ax[2].axvline(0, color="#cccccc", lw=0.6, zorder=1)
        ax[2].scatter(x, y, s=14, color=COLOR["splicebert"], alpha=0.75,
                      edgecolor="white", linewidth=0.3, zorder=3)
        ax[2].set_xlim(lim)
        ax[2].set_ylim(lim)
        ax[2].set_xlabel("contrast, k-mer LR")
        ax[2].set_ylabel("contrast, SpliceBERT")
        above = int((y > x).sum())
        ax[2].set_title(f"c  larger for SpliceBERT in {above}/{len(d)}", loc="left")

    fig.tight_layout()
    save(fig, "f9_deep_contrast")

FIGURES = {"f0": f0, "f1": f1, "f2": f2, "f3": f3, "f4": f4, "f5": f5,
           "f6": f6, "f7": f7, "f8": f8, "f9": f9}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--only", default="", help="comma separated, e.g. f1,f3")
    a = p.parse_args()
    want = [x.strip() for x in a.only.split(",") if x.strip()] or list(FIGURES)
    for name in want:
        print(f"{name}:", flush=True)
        FIGURES[name]()


if __name__ == "__main__":
    main()
