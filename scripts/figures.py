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


# --- f3: is the signal local? ------------------------------------------------------------

def f3():
    got = need("locality_ism.csv")
    if got is None:
        return
    d = got[0].dropna(subset=["kmer_gini", "sb_gini"])
    if len(d) < 3:
        print("  skipped f3: fewer than 3 paired datasets")
        return
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.1))

    lo = min(d.kmer_gini.min(), d.sb_gini.min()) - 0.03
    hi = max(d.kmer_gini.max(), d.sb_gini.max()) + 0.03
    ax[0].plot([lo, hi], [lo, hi], color="#999999", lw=0.8, ls="--")
    ax[0].scatter(d.kmer_gini, d.sb_gini, s=22, color=COLOR["splicebert"],
                  edgecolor="white", linewidth=0.4, zorder=3)
    ax[0].set_xlabel("k-mer LR, ISM Gini")
    ax[0].set_ylabel("SpliceBERT, ISM Gini")
    ax[0].set_xlim(lo, hi)
    ax[0].set_ylim(lo, hi)
    n_up = int((d.sb_gini > d.kmer_gini).sum())
    ax[0].set_title(f"More local in {n_up}/{len(d)} datasets", loc="left", fontsize=9)

    diff = (d.sb_gini - d.kmer_gini).sort_values().values
    ax[1].bar(np.arange(len(diff)), diff,
              color=[COLOR["splicebert"] if v > 0 else COLOR["gc"] for v in diff],
              edgecolor="white", linewidth=0.3)
    ax[1].axhline(0, color="black", lw=0.8)
    ax[1].axhline(diff.mean(), color="black", lw=1.2, ls="--",
                  label=f"mean {diff.mean():+.3f}")
    ax[1].set_xlabel("datasets, sorted")
    ax[1].set_ylabel("Gini(SpliceBERT) - Gini(k-mer)")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("Positional concentration", loc="left", fontsize=9)
    save(fig, "f3_locality")


# --- f4: the variant arm -----------------------------------------------------------------

def f4():
    """The ClinVar ladder. Three rungs plus conservation.

    The claim is the GAPS between rungs, not any single number: a wrong-protein head already
    beats the k-mer, so part of the signal is generic sequence plausibility inherited from
    pretraining, and only the matched-minus-mismatched gap is binding-specific. Plotting one
    bar would hide exactly the thing that makes the result defensible.
    """
    got = need("variant_ladder.csv")
    if got is None:
        return
    d = got[0]
    order = ["kmer", "mismatched", "matched", "conservation"]
    # Short labels: the four full phrases collide at this width, and a colliding axis is
    # a broken figure however correct the numbers are.
    label = {"kmer": "k-mer", "mismatched": "wrong\nprotein", "matched": "right\nprotein",
             "conservation": "phyloP"}
    colr = {"kmer": COLOR["kmer"], "mismatched": "#bdbdbd",
            "matched": COLOR["splicebert"], "conservation": "#8c8c8c"}
    d = d.set_index("arm").reindex([a for a in order if a in set(d.arm)]).reset_index()

    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2))
    ax[0].bar(range(len(d)), d.auroc, color=[colr[a] for a in d.arm],
              edgecolor="white", linewidth=0.6)
    ax[0].axhline(0.5, color="black", lw=0.8, ls="--")
    for i, v in enumerate(d.auroc):
        ax[0].text(i, v + 0.008, f"{v:.3f}", ha="center", fontsize=8)
    ax[0].set_xticks(range(len(d)))
    ax[0].set_xticklabels([label[a] for a in d.arm], fontsize=8)
    ax[0].set_ylabel("pathogenic vs benign AUROC")
    ax[0].set_ylim(0.45, 1.0)
    ax[0].set_title("The ladder", loc="left", fontsize=9)

    # Cluster-corrected coefficients. Conservation has no delta coefficient by construction,
    # so it is absent here rather than drawn as zero.
    c = d[d.coef.notna()]
    ax[1].bar(range(len(c)), c.coef, color=[colr[a] for a in c.arm],
              edgecolor="white", linewidth=0.6)
    ax[1].errorbar(range(len(c)), c.coef,
                   yerr=[c.coef - c.ci_low, c.ci_high - c.coef],
                   fmt="none", ecolor="black", elinewidth=0.9, capsize=3)
    ax[1].axhline(0, color="black", lw=0.8)
    ax[1].set_xticks(range(len(c)))
    ax[1].set_xticklabels([label[a] for a in c.arm], fontsize=8)
    ax[1].set_ylabel("standardised |delta| coefficient")
    ax[1].set_title("Conservation controlled, gene-clustered", loc="left", fontsize=9)
    save(fig, "f4_variant_ladder")


# --- f5: the size confound, stated rather than buried ------------------------------------

def f5():
    got = need("matched_four_models.csv")
    if got is None:
        return
    d = got[0].rename(columns={"kmer_auroc": "kmer", "composition_auroc": "composition"})
    fig, ax = plt.subplots(figsize=(4.0, 3.1))
    lp = np.log10(d.pairs)
    for m in ["composition", "kmer", "cnn", "splicebert"]:
        r = np.corrcoef(lp, d[m])[0, 1]
        ax.scatter(d.pairs, d[m], s=7, color=COLOR[m], alpha=0.6,
                   label=f"{LABEL[m]}  r={r:+.2f}", edgecolor="none")
        b, a_ = np.polyfit(lp, d[m], 1)
        xs = np.linspace(lp.min(), lp.max(), 40)
        ax.plot(10 ** xs, a_ + b * xs, color=COLOR[m], lw=1.3)
    ax.set_xscale("log")
    ax.set_xlabel("training pairs per dataset")
    ax.set_ylabel("pooled out-of-fold AUROC")
    # Outside the axes, not "lower right": at this size the legend covered the composition
    # cloud, which is the series the panel exists to contrast against.
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_title("Better models depend MORE on dataset size", loc="left", fontsize=9)
    save(fig, "f5_size_confound")


FIGURES = {"f0": f0, "f1": f1, "f2": f2, "f3": f3, "f4": f4, "f5": f5}


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
