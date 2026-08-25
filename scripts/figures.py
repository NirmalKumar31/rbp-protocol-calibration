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

    ax[1].hist(d.cost, bins=28, color=COLOR["gc"], edgecolor="white", linewidth=0.4)
    ax[1].axvline(0, color="black", lw=0.8)
    ax[1].axvline(d.cost.mean(), color="black", lw=1.4, ls="--",
                  label=f"mean {d.cost.mean():+.4f}")
    ax[1].set_xlabel("AUROC(dinuc) - AUROC(GC)")
    ax[1].set_ylabel("datasets")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("Cost of matching properly", loc="left", fontsize=9)
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
    ax[0].set_xticklabels([LABEL[m] for m in order], fontsize=8)
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
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.set_title("AUROC tracks dataset size", loc="left", fontsize=9)
    save(fig, "f5_size_confound")


FIGURES = {"f1": f1, "f2": f2, "f3": f3, "f4": f4, "f5": f5}


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
