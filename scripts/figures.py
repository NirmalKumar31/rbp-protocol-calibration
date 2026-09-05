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
# COLOUR IS CHECKED, NOT CHOSEN. tests/unit/test_palette_cvd.py measures OKLab separation for
# every pair that shares a panel, under normal vision and under simulated protanopia and
# deuteranopia, and fails below the thresholds. Three pairs failed the first time it ran, all
# of them at NORMAL vision, meaning full-colour readers could not tell them apart either:
# k-mer against SpliceBERT at 14.8, the dinucleotide arm against Horlacher's first arm at 9.5,
# and the two greys of Figure 5b -- which are ADJACENT bars -- at 8.7. Do not edit these values
# without rerunning that test.
#
# "neg2" also used to be COLOR["cnn"], so one orange meant the bias-aware ARM in Figure 1 and
# the convolutional MODEL in Figure 4. Colour follows the entity, so the arm has its own now.
COLOR = {"composition": "#8c8c8c", "kmer": "#4878a8", "cnn": "#e08214",
         "splicebert": "#276419", "gc": "#b2182b", "dinuc": "#2166ac",
         "neg2": "#e7298a", "theirs": "#762a83",
         "grey_mid": "#969696", "grey_light": "#d9d9d9"}
LABEL = {"composition": "composition (19 feat)", "kmer": "4-mer LR", "cnn": "CNN",
         "splicebert": "SpliceBERT"}

plt.rcParams.update({"figure.dpi": 150, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "grid.linewidth": 0.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.axisbelow": True,
                     # TYPE 42 (TrueType), NOT matplotlib's default Type 3. Type 3 embeds
                     # glyphs as PostScript drawing programs: the text is not selectable, not
                     # searchable and not reliably extractable, and NAR, Bioinformatics and
                     # every IEEE venue reject it outright at submission. It is a one-line fix
                     # that otherwise surfaces as a desk rejection weeks later.
                     "pdf.fonttype": 42, "ps.fonttype": 42})
SAVE_DPI = 400          # >= 300 for print; the panel figures are dense


def save(fig, name):
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}", bbox_inches="tight", dpi=SAVE_DPI)
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


def clustered_mean_err(frame, col, n_boot=4000, seed=7):
    """Mean of `col` with an asymmetric error bar clustered on protein.

    WHY THIS EXISTS. Two panels drew `.sem()` over the 94 datasets, which treats them as 94
    independent observations. They are 79 proteins, fifteen of which contribute two datasets
    each at a within-protein correlation of 0.92 for the primary contrast, and EVERY headline
    interval in the paper resamples proteins. An error bar narrower than the inference it
    illustrates makes a figure look stronger than the text it belongs to, which is the one
    direction a figure must never err in.

    Returns (mean, [[lower], [upper]]) shaped for matplotlib's yerr.
    """
    rng = np.random.default_rng(seed)
    g = frame[["protein", col]].dropna()
    groups = [v.to_numpy() for _, v in g.groupby("protein")[col]]
    m = float(np.concatenate(groups).mean())
    draws = np.empty(n_boot)
    idx = np.arange(len(groups))
    for b in range(n_boot):
        pick = rng.choice(idx, len(idx), replace=True)
        draws[b] = np.concatenate([groups[i] for i in pick]).mean()
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return m, [[max(m - lo, 0.0)], [max(hi - m, 0.0)]]


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
    # From the palette rather than a literal: #2b6a4d was SpliceBERT's old green and became an
    # orphan when the CVD test moved it, leaving one figure on a hex nothing else referenced.
    ax[1].bar(by.index, by.values,
              color=[COLOR["composition"], COLOR["splicebert"]], width=0.6)
    for i, v in enumerate(by.values):
        ax[1].text(i, v + 0.8, str(v), ha="center", fontsize=8)
    ax[1].set_ylabel("datasets")
    ax[1].set_ylim(0, by.max() * 1.18)
    ax[1].set_title(f"b  {d.protein.nunique()} proteins, {n_both} in both lines",
                    loc="left", fontsize=9)

    # (c) THE THING THE PAPER IS ABOUT: what each protocol leaves the model to do. This panel
    # used to plot ClinVar coverage per dataset, for an analysis that has since been RETRACTED
    # -- so the overview figure was advertising a result the paper no longer makes. The three
    # composition-baseline distributions belong here instead: they are the paper's independent
    # variable, and their near-disjointness is R1l in one picture.
    arms = need("three_arm_per_dataset.csv")
    if arms is not None:
        a3 = arms[0]
        for key, col, lab in (("dn", COLOR["dinuc"], "dinucleotide-matched"),
                              ("gc", COLOR["gc"], "GC-matched"),
                              ("neg2", COLOR["neg2"], "bias-aware")):
            v = a3[f"comp_{key}"]
            ax[2].hist(v, bins=np.linspace(0.5, 1.0, 26), histtype="step", linewidth=1.5,
                       color=col, label=f"{lab}  {v.mean():.3f}")
        ax[2].legend(frameon=False, fontsize=6.8, loc="upper left")
        ax[2].set_ylim(0, ax[2].get_ylim()[1] * 1.35)   # room for the legend over the peak
        ax[2].set_xlabel("composition-only AUROC (the baseline the protocol leaves)")
        ax[2].set_ylabel("datasets")
        ax[2].set_title("c  same positives, three protocols", loc="left", fontsize=9)
    else:
        ax[2].axis("off")
    save(fig, "f0_panel_overview")


# --- f1: the cost of the negative-set protocol -------------------------------------------

def f1():
    """The headline. Same positives, same model, two negative sets.

    Paired per dataset, because the two arms share their positives -- an unpaired plot
    would throw away exactly the structure that makes 94/94 meaningful. The panel is 94 paired
    datasets; earlier counts refer to the candidate pool and do not apply here.
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
    for i, (k, _lab, c) in enumerate(bars):
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
    for i, (k, _lab, c) in enumerate(parts):
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
    def val(k):
        return float(q.loc[k, "value"])
    def err(k):
        return [[val(k) - float(q.loc[k, "ci_low"])], [float(q.loc[k, "ci_high"]) - val(k)]]

    r2 = np.sqrt(2.0)
    def dp(a):
        return r2 * norm.ppf(np.clip(a, 1e-6, 1 - 1e-6))
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
    for i, (k, _lab, c) in enumerate(keys):
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
            mu, err = clustered_mean_err(d, f"{m}_gain_{key}")
            ax[0].bar(i + (j - 0.5) * w, mu, w * 0.9, color=COLOR[arm],
                      zorder=3, label=f"{arm}-matched" if i == 0 else None)
            ax[0].errorbar(i + (j - 0.5) * w, mu, yerr=err, color="black",
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

# --- f10: R1k, the paper's central figure ------------------------------------------------
#
# Three protocols, one model, one set of positives. The figure has to carry the whole thesis,
# so it makes exactly three points and no others. (a) the measured contribution spans 5.4x and
# the ordering is NOT by negative-set hardness -- the field's own bias-aware protocol is the
# LOWEST. (b) why: the contribution tracks how much room the composition baseline leaves, which
# is the circularity critique adopted rather than deflected. (c) it is not an averaging
# artefact -- the reordering happens dataset by dataset.

def f10():
    # three_arm_contrast.csv was needed only by the panel that moved to Figure 3, but it stays
    # in the gate: if it is missing the run is incomplete and this figure should not be drawn
    # from a half-built results directory.
    t = need("three_arm_per_dataset.csv", "three_arm_contrast.csv")
    if t is None:
        return
    d = t[0]
    arms = [("dn", "dinucleotide\nmatched"), ("gc", "GC\nmatched"),
            ("neg2", "bias-aware\n(other RBPs' sites)")]
    col = {"dn": COLOR["dinuc"], "gc": COLOR["gc"], "neg2": COLOR["neg2"]}
    # TWO PANELS, NOT THREE. The old panel b -- contribution against the composition baseline
    # over all 282 cells -- is the same plot as Figure 3a, which belongs to the subsection that
    # argues from the gradient. Duplicating it here spent a third of this figure restating a
    # later one. Kept in Figure 3; the pointer below sends the reader there.
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4))

    # a. composition baseline and nested contribution, per protocol
    x = np.arange(len(arms))
    top = 0.0
    for i, (a, _) in enumerate(arms):
        m, err = clustered_mean_err(d, f"gain_{a}")
        ax[0].bar(i, m, 0.6, color=col[a], zorder=3)
        ax[0].errorbar(i, m, yerr=err, color="black", lw=1, capsize=3, zorder=4)
        ax[0].text(i, m + err[1][0] + 0.004, f"{m:+.4f}", ha="center", fontsize=7.5)
        top = max(top, m + err[1][0])
    ax[0].set_xticks(x, [lbl for _, lbl in arms], fontsize=7.5)
    ax[0].set_ylabel("nested contribution over composition")
    # From the whiskers, not a literal: the clustered interval is wider than the .sem() this
    # panel used to draw, and a hardcoded 0.082 put the tallest value label into the title.
    ax[0].set_ylim(0, top * 1.16)
    ax[0].set_title("a  same model, same positives: 5.4x", loc="left")

    # b. per dataset, so it is not an averaging artefact
    lim = [-0.02, max(d.gain_dn.max(), d.gain_gc.max()) + 0.01]
    ax[1].plot(lim, lim, color="#999999", lw=0.8, ls="--", zorder=2)
    ax[1].scatter(d.gain_gc, d.gain_dn, s=12, color=COLOR["dinuc"], alpha=0.7,
                  edgecolor="white", linewidth=0.25, zorder=3, label="dinuc vs GC")
    ax[1].scatter(d.gain_gc, d.gain_neg2, s=12, color=COLOR["neg2"], alpha=0.7,
                  edgecolor="white", linewidth=0.25, zorder=3, label="bias-aware vs GC")
    ax[1].set_xlim(lim)
    ax[1].set_ylim(lim)
    ax[1].set_xlabel("nested contribution, GC-matched")
    ax[1].set_ylabel("nested contribution, other protocol")
    ax[1].legend(fontsize=7, frameon=False, loc="upper left")
    up = int((d.gain_dn > d.gain_gc).sum())
    dn_ = int((d.gain_neg2 < d.gain_gc).sum())
    ax[1].set_title(f"b  above in {up}/94, below in {dn_}/94", loc="left")

    fig.tight_layout()
    save(fig, "f10_three_protocols")

# --- f11: no rescaling reaches protocol independence (R1m) -------------------------------

def f11():
    """The result that earns the title. Eight monotone transforms, none reaching 1.0x.

    A dot chart rather than bars: the quantity is a ratio with an interval, the coordinates are
    not commensurable with each other, and a bar from zero would invite reading the areas.
    """
    t = need("scale_sweep.csv")
    if t is None:
        return
    d = t[0].set_index("check")
    rows = [("raw AUROC gain", "fold range, raw AUROC gain"),
            ("Somers' D  (affine control)", "fold range, Somers' D gain"),
            ("arcsine increment", "fold range, arcsine increment"),
            ("cloglog increment", "fold range, cloglog increment"),
            ("d' increment (binormal)", "fold range, d' increment (binormal)"),
            ("logit increment", "fold range, logit increment"),
            ("headroom-normalised, g/(1-comp)", "fold range, headroom-normalised, g/(1-comp)")]
    rows = [(lab, k) for lab, k in rows if k in d.index]
    rows.sort(key=lambda r: -float(d.loc[r[1], "value"]))

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    for i, (_lab, k) in enumerate(rows):
        v = float(d.loc[k, "value"])
        lo, hi = float(d.loc[k, "ci_low"]), float(d.loc[k, "ci_high"])
        best = "headroom" in k
        col = COLOR["splicebert"] if best else COLOR["kmer"]
        ax.plot([lo, hi], [i, i], color=col, linewidth=2.0, solid_capstyle="round", zorder=2)
        ax.scatter([v], [i], s=34, color=col, zorder=3, edgecolor="white", linewidth=0.6)
        ax.text(hi + 0.15, i, f"{v:.2f}x", va="center", fontsize=7.5, color=col)
    ax.axvline(1.0, color="#404040", linestyle="--", linewidth=1.1)
    # Rotated onto the rule itself: horizontally it collided with the title, and nudging it
    # down would have run it through the shortest interval, which is the one that matters.
    ax.text(1.0, (len(rows) - 1) / 2, "  protocol independence", rotation=90, va="center",
            ha="left", fontsize=7, color="#404040")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for lab, _ in rows], fontsize=8)
    ax.set_xlabel("fold range in measured contribution across the three protocols")
    ax.set_xlim(0.6, None)
    ax.grid(axis="y", visible=False)
    ax.set_title("no monotone rescaling reaches a protocol-free quantity", loc="left",
                 fontsize=9)
    fig.tight_layout()
    save(fig, "f11_scale_sweep")


# --- f12: it is the baseline, not the protocol label (R1n) --------------------------------

def f12():
    """The thesis in two panels: the gradient, and the natural experiment that isolates it."""
    t = need("three_arm_per_dataset.csv", "protocol_or_baseline.csv")
    if t is None:
        return
    d, s = t[0], t[1].set_index("check")
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))

    # (a) every arm-dataset cell: the gain falls as the baseline the protocol leaves rises.
    for key, col, lab in (("dn", COLOR["dinuc"], "dinucleotide-matched"),
                          ("gc", COLOR["gc"], "GC-matched"),
                          # "neg2" is the internal arm key and had leaked into the legend as
                          # display text, while the panel beside it said "bias-aware". A reader
                          # has no way to know they are the same arm.
                          ("neg2", COLOR["neg2"], "bias-aware")):
        ax[0].scatter(d[f"comp_{key}"], d[f"gain_{key}"], s=13, color=col, alpha=0.7,
                      edgecolor="white", linewidth=0.25, label=lab, zorder=3)
    ax[0].axhline(0, color="#404040", linewidth=0.8)
    ax[0].set_xlabel("composition-only AUROC")
    ax[0].set_ylabel("nested contribution")
    ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_title("a  282 cells, one gradient", loc="left", fontsize=9)

    # (b) THE NATURAL EXPERIMENT. Where neg2 lowers the baseline, its deficit reverses. This is
    # the panel that separates "the baseline does it" from "the protocol label does it".
    hi = (d.comp_neg2 > d.comp_gc).values
    labels, vals, los, his = [], [], [], []
    for lab, key in (("bias-aware raises\nthe baseline", "concordant"),
                     ("bias-aware lowers\nthe baseline", "discordant")):
        k = f"neg2 minus gc gain, {key} datasets"
        if k not in s.index:
            ax[1].axis("off")
            break
        labels.append(f"{lab}\nn={int(hi.sum() if key == 'concordant' else (~hi).sum())}")
        vals.append(float(s.loc[k, "value"]))
        los.append(float(s.loc[k, "ci_low"]))
        his.append(float(s.loc[k, "ci_high"]))
    if vals:
        x = np.arange(len(vals))
        cols = [COLOR["gc"], COLOR["splicebert"]]
        for i, v in enumerate(vals):
            ax[1].plot([x[i], x[i]], [los[i], his[i]], color=cols[i], linewidth=2.2,
                       solid_capstyle="round", zorder=2)
            ax[1].scatter([x[i]], [v], s=44, color=cols[i], zorder=3, edgecolor="white",
                          linewidth=0.6)
            ax[1].text(x[i] + 0.12, v, f"{v:+.4f}", fontsize=8, va="center", color=cols[i])
        ax[1].axhline(0, color="#404040", linestyle="--", linewidth=1.1)
        ax[1].set_xticks(x)
        ax[1].set_xticklabels(labels, fontsize=8)
        ax[1].set_xlim(-0.5, len(vals) - 0.1)
        ax[1].set_ylabel("bias-aware minus GC, nested contribution")
        ax[1].grid(axis="x", visible=False)
        ax[1].set_title("b  the deficit follows the baseline, not the label", loc="left",
                        fontsize=9)
    fig.tight_layout()
    save(fig, "f12_protocol_or_baseline")


# --- f13: the order-3 collapse is the 4-mer's alone (R1r) ---------------------------------

def f13():
    """Raise the baseline one order and the k-mer loses four fifths; SpliceBERT loses a quarter."""
    t = need("baseline_order_models.csv")
    if t is None:
        return
    d = t[0].set_index("check")
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.3), sharey=True)
    models = ("kmer", "cnn", "splicebert")
    for j, (arm, title) in enumerate((("gc", "a  GC-matched"),
                                      ("dn", "b  dinucleotide-matched"))):
        for i, m in enumerate(models):
            k = f"{m} fraction surviving order 3, {arm} arm"
            if k not in d.index:
                continue
            v = float(d.loc[k, "value"])
            lo, hi = float(d.loc[k, "ci_low"]), float(d.loc[k, "ci_high"])
            ax[j].plot([i, i], [lo, hi], color=COLOR[m], linewidth=2.4,
                       solid_capstyle="round", zorder=2)
            ax[j].scatter([i], [v], s=46, color=COLOR[m], zorder=3, edgecolor="white",
                          linewidth=0.6)
            ax[j].text(i + 0.14, v, f"{100 * v:.0f}%", fontsize=8.5, va="center",
                       color=COLOR[m])
        ax[j].set_xticks(range(len(models)))
        ax[j].set_xticklabels([LABEL[m] for m in models], fontsize=8)
        ax[j].set_xlim(-0.5, len(models) - 0.35)
        ax[j].set_ylim(0, 1.0)
        ax[j].grid(axis="x", visible=False)
        ax[j].set_title(title, loc="left", fontsize=9)
    ax[0].set_ylabel("share of the contribution surviving\nan order-3 composition baseline")
    fig.tight_layout()
    save(fig, "f13_baseline_order_models")


# --- f14: external validation on an independent benchmark (R1p) --------------------------

def f14():
    """Their positives, their negatives, their folds. Only the measurement is ours.

    Two panels because the section makes two claims: the range travels, and the two-family
    mechanism travels. The second is the stronger one and was previously unplotted.
    """
    t = need("horlacher_per_dataset.csv", "three_arm_per_dataset.csv")
    if t is None:
        return
    h, d = t[0], t[1]
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.5))

    # (a) the two arms' contributions, paired per dataset.
    lim = (min(h.gain_n1.min(), h.gain_n2.min()) - 0.01,
           max(h.gain_n1.max(), h.gain_n2.max()) + 0.01)
    ax[0].plot(lim, lim, color="#404040", linewidth=0.9, linestyle="--", zorder=1)
    ax[0].scatter(h.gain_n1, h.gain_n2, s=16, color=COLOR["kmer"], alpha=0.8,
                  edgecolor="white", linewidth=0.3, zorder=3)
    ax[0].set_xlim(lim)
    ax[0].set_ylim(lim)
    ax[0].set_xlabel("contribution, negative-1 (transcript background)")
    ax[0].set_ylabel("contribution, negative-2 (other RBPs' sites)")
    below = int((h.gain_n2 < h.gain_n1).sum())
    ax[0].set_title(f"a  {below}/{len(h)} below the diagonal, "
                    f"{h.gain_n1.mean() / h.gain_n2.mean():.2f}x",
                    loc="left", fontsize=9)

    # (b) THE WITHIN-ARM GRADIENT, in their data and ours, WITH PROTEIN-CLUSTERED INTERVALS.
    #
    # This panel used to be five bare bars under a title asserting that the gradient is a
    # property of composition-matched negatives. That is the grouping the Results section
    # declines to endorse, because the statistic is not invariant to which term goes on the
    # horizontal axis and the dinucleotide arm changes family under two of the three choices.
    # Bars alone cannot show what the claim rests on. The intervals can: the three
    # composition-matched arms exclude zero and the two other-RBPs'-sites arms do not, in both
    # benchmarks independently, which is a statement about detectability and not a taxonomy.
    from scipy.stats import spearmanr
    rng = np.random.default_rng(7)

    def clustered(frame, cx, cy, n_boot=4000):
        """Spearman with a percentile interval over resampled proteins."""
        by = {q: g for q, g in frame.groupby("protein")}
        names = frame.protein.unique()
        draws = []
        for _ in range(n_boot):
            b = pd.concat([by[q] for q in rng.choice(names, len(names), replace=True)],
                          ignore_index=True)
            draws.append(spearmanr(b[cx], b[cy])[0])
        lo, hi = np.percentile(draws, [2.5, 97.5])
        return spearmanr(frame[cx], frame[cy])[0], lo, hi

    bars = [
        ("ours\nGC", *clustered(d, "comp_gc", "gain_gc"), COLOR["gc"]),
        ("ours\ndinuc", *clustered(d, "comp_dn", "gain_dn"), COLOR["dinuc"]),
        ("theirs\nneg-1", *clustered(h, "comp_n1", "gain_n1"), COLOR["theirs"]),
        ("ours\nbias-aware", *clustered(d, "comp_neg2", "gain_neg2"), COLOR["grey_mid"]),
        ("theirs\nneg-2", *clustered(h, "comp_n2", "gain_n2"), COLOR["grey_light"]),
    ]
    x = np.arange(len(bars))
    vals = [b[1] for b in bars]
    err = np.array([[b[1] - b[2] for b in bars], [b[3] - b[1] for b in bars]])
    ax[1].bar(x, vals, color=[b[4] for b in bars], width=0.68,
              edgecolor="white", linewidth=0.6)
    ax[1].errorbar(x, vals, yerr=np.abs(err), fmt="none", ecolor="#404040",
                   elinewidth=1.0, capsize=3.0, zorder=4)
    for i, b in enumerate(bars):
        # OUTSIDE the bar, past the far whisker: a label placed at the near whisker lands on
        # the coloured fill, where dark grey text on dark red or purple is unreadable.
        ax[1].text(i, b[2] - 0.055, f"{b[1]:+.2f}", ha="center", va="top", fontsize=7.5,
                   color="#404040")
    ax[1].axhline(0, color="#404040", linewidth=0.9)
    ax[1].set_ylim(min(b[2] for b in bars) - 0.15, max(0.12, max(b[3] for b in bars) + 0.05))
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([b[0] for b in bars], fontsize=7.5)
    ax[1].set_ylabel("within-arm Spearman(baseline, contribution)")
    ax[1].grid(axis="x", visible=False)
    ax[1].set_title("b  gradient detectable in composition-matched arms, not in the others",
                    loc="left", fontsize=9)
    fig.tight_layout()
    save(fig, "f14_external_validation")


# --- f15: the recommendation, and where it stops working (R1q + R1u) ----------------------

def f15():
    """The deliverable, tested in sample and out. It does not survive the external test.

    Plotted together deliberately: a figure that showed only the in-sample panel would be
    the strongest available misrepresentation of this paper's own evidence.
    """
    t = need("recommendation_works.csv", "transport_check.csv")
    if t is None:
        return
    r, tr = t[0].set_index("check"), t[1].set_index("check")
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))

    # (a) in sample: rank agreement raw -> headroom, three protocol pairs.
    pairs = [("gc", "dn"), ("gc", "neg2"), ("dn", "neg2")]
    labs = ["GC vs\ndinuc", "GC vs\nbias-aware", "dinuc vs\nbias-aware"]
    for i, (a, b) in enumerate(pairs):
        k1 = f"rank agreement, raw, {a} vs {b}"
        k2 = f"rank agreement, headroom, {a} vs {b}"
        if k1 not in r.index or k2 not in r.index:
            continue
        v1, v2 = float(r.loc[k1, "value"]), float(r.loc[k2, "value"])
        sig = i == 0          # only the first pair's improvement clears zero
        ax[0].plot([i - 0.16, i + 0.16], [v1, v2], color="#404040", linewidth=1.0, zorder=2)
        ax[0].scatter([i - 0.16], [v1], s=40, color=COLOR["grey_light"], zorder=3,
                      edgecolor="white", linewidth=0.5, label="raw" if i == 0 else None)
        ax[0].scatter([i + 0.16], [v2], s=40, color=COLOR["splicebert"], zorder=3,
                      edgecolor="white", linewidth=0.5,
                      label="headroom-normalised" if i == 0 else None)
        ax[0].text(i + 0.24, v2, "*" if sig else "n.s.", fontsize=8, va="center",
                   color="#404040")
    ax[0].set_xticks(range(len(labs)))
    ax[0].set_xticklabels(labs, fontsize=8)
    ax[0].set_xlim(-0.5, len(labs) - 0.3)
    ax[0].set_ylabel("cross-protocol rank agreement")
    ax[0].legend(frameon=False, fontsize=7.5, loc="lower left")
    ax[0].grid(axis="x", visible=False)
    ax[0].set_title("a  our data: improves in 3/3, one interval clear of zero", loc="left",
                    fontsize=9)

    # (b) OUT of sample, on their benchmark: both pre-registered criteria fire the wrong way.
    keys = [("rank agreement", "external rank agreement, raw",
             "external rank agreement, headroom"),
            ("disagreement", "external scale-free disagreement, raw",
             "external scale-free disagreement, headroom")]
    for i, (_lab, kr, kh) in enumerate(keys):
        if kr not in tr.index:
            continue
        v1, v2 = float(tr.loc[kr, "value"]), float(tr.loc[kh, "value"])
        worse = v2 < v1 if i == 0 else v2 > v1
        col = COLOR["gc"] if worse else COLOR["splicebert"]
        ax[1].plot([i - 0.16, i + 0.16], [v1, v2], color="#404040", linewidth=1.0, zorder=2)
        ax[1].scatter([i - 0.16], [v1], s=40, color=COLOR["grey_light"], zorder=3,
                      edgecolor="white", linewidth=0.5)
        ax[1].scatter([i + 0.16], [v2], s=40, color=col, zorder=3, edgecolor="white",
                      linewidth=0.5)
        ax[1].annotate("", xy=(i + 0.16, v2), xytext=(i - 0.16, v1),
                       arrowprops=dict(arrowstyle="->", color=col, lw=1.3))
        # Beside the arrow, not above or below it: above ran into the panel title and below
        # ran into the tick labels.
        ax[1].text(i + 0.26, (v1 + v2) / 2, "wrong\ndirection" if worse else "", ha="left",
                   va="center", fontsize=7.5, color=col)
    ax[1].set_xticks(range(len(keys)))
    ax[1].set_xticklabels([k[0] for k in keys], fontsize=8)
    ax[1].set_xlim(-0.5, len(keys) - 0.05)
    ax[1].set_ylabel("value on the independent benchmark")
    ax[1].grid(axis="x", visible=False)
    ax[1].set_title("b  independent benchmark: both criteria fail to replicate", loc="left",
                    fontsize=9)
    fig.tight_layout()
    save(fig, "f15_recommendation")


def f16():
    """B3: the contribution as a function of where the baseline stops, and where it breaks.

    THE ORDER-4 COLUMN IS DRAWN, NOT DROPPED. At order four the baseline spans the 4-mer's
    own feature space, so its true contribution is zero and the +0.09 to +0.14 the estimator
    reports is the instrument's error. Hiding that column would turn the figure into three
    tidy declining curves and lose the section's main result; drawing it without marking it
    would invite reading a noise floor as a contribution. So it is drawn beyond a rule, in
    a shaded panel region, with the baseline's own AUROC below it as the reason.
    """
    t = need("order_profile.csv")
    if t is None:
        return
    q = t[0].set_index("check")
    orders = [1, 2, 3, 4]
    arms = [("dn", "dinucleotide-matched"), ("gc", "GC-matched"), ("neg2", "bias-aware")]
    fig, ax = plt.subplots(2, 3, figsize=(9.4, 5.4), sharex=True,
                           gridspec_kw={"height_ratios": [2.0, 1.0]})

    for col, (arm, title) in enumerate(arms):
        a0, a1 = ax[0][col], ax[1][col]
        # The region where the baseline has stopped fitting. Shaded rather than cut.
        for a in (a0, a1):
            a.axvspan(3.5, 4.35, color=COLOR["grey_light"], alpha=0.55, zorder=0, lw=0)
        for model in ("kmer", "cnn", "splicebert"):
            ys, los, his = [], [], []
            for o in orders:
                k = f"{model} gain at order {o}, {arm} arm"
                ys.append(float(q.loc[k, "value"]))
                los.append(float(q.loc[k, "ci_low"]))
                his.append(float(q.loc[k, "ci_high"]))
            # Orders 1-3 joined; the step into 4 dashed, because it crosses from a baseline
            # into a noise floor and a solid line would assert a continuous quantity.
            a0.plot(orders[:3], ys[:3], color=COLOR[model], linewidth=2.0, zorder=3,
                    marker="o", markersize=5, markeredgecolor="white", markeredgewidth=0.6,
                    label=LABEL[model] if col == 0 else None)
            a0.plot(orders[2:], ys[2:], color=COLOR[model], linewidth=2.0, zorder=3,
                    linestyle=(0, (3, 2)), marker="o", markersize=5,
                    markeredgecolor="white", markeredgewidth=0.6)
            a0.fill_between(orders, los, his, color=COLOR[model], alpha=0.15, lw=0, zorder=2)
        a0.axhline(0, color="#404040", linewidth=0.8, zorder=1)
        a0.set_title(f"{'abc'[col]}  {title}", loc="left", fontsize=9)
        a0.grid(axis="x", visible=False)
        # THE BAND'S LABEL CARRIES ITS OWN REASON. Annotating the count beside the curve in
        # the panel below collided with the curve in two of three arms, and separating cause
        # from effect made the reader join them up. One label, both facts.
        fell = int(q.loc[f"baseline AUROC fell from order 3 to 4, {arm} arm", "value"])
        a0.text(3.93, a0.get_ylim()[0] + 0.97 * (a0.get_ylim()[1] - a0.get_ylim()[0]),
                f"baseline\nno longer fits\n(its own AUROC\nfalls on {fell}/94)",
                fontsize=6.5, ha="right", va="top", color="#404040")
        if col == 0:
            a0.set_ylabel("nested contribution")
            a0.legend(frameon=False, fontsize=7.5, loc="upper left")

        # THE DIAGNOSTIC, directly under the curves it explains: the baseline's own
        # out-of-fold AUROC. Where it stops rising, the panel above stops being a baseline.
        cs = [float(q.loc[f"composition AUROC at order {o}, {arm} arm", "value"])
              for o in orders]
        a1.plot(orders, cs, color=COLOR["composition"], linewidth=2.0, marker="s",
                markersize=4.5, markeredgecolor="white", markeredgewidth=0.6, zorder=3,
                label="baseline's own AUROC" if col == 0 else None)
        a1.set_xticks(orders)
        a1.set_xlabel("order of the composition baseline")
        a1.grid(axis="x", visible=False)
        if col == 0:
            a1.set_ylabel("composition AUROC")
            a1.legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout()
    save(fig, "f16_order_profile")

FIGURES = {"f0": f0, "f1": f1, "f2": f2, "f3": f3, "f4": f4, "f5": f5,
           "f6": f6, "f7": f7, "f8": f8, "f9": f9, "f10": f10,
           "f11": f11, "f12": f12, "f13": f13, "f14": f14, "f15": f15,
           "f16": f16}


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
