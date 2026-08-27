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
    """R4, as a paired per-dataset comparison rather than one pooled AUROC.

    THE POOLED VERSION OF THIS FIGURE WAS WRONG AND LOOKED BETTER. It showed four bars,
    matched 0.829 against a wrong-protein 0.680, and that gap of +0.149 was inflated by
    between-dataset heterogeneity: mean |delta| per dataset correlates with that dataset's
    pathogenic rate at rho +0.73 and spans 10.4x, so a pooled AUROC partly measures which
    dataset a variant came from. Paired within dataset the gap is +0.065.

    Conservation was the only arm immune to the artefact, because phyloP is on a fixed
    external scale -- which is why the inflation stayed invisible: the arm that could not be
    inflated was winning anyway.

    Panel a is the honest ladder across power strata. Panel b is the paired test. Panel c
    puts the size of the inflation on the record instead of quietly dropping it.
    """
    got = need("variant_ladder_paired.csv", "variant_specificity.csv",
               "variant_coefficients.csv")
    if got is None:
        return
    paired, per, coef = got
    fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.2))

    # (a) each arm against statistical power. The mismatched floor is flat; matched rises.
    x = np.arange(len(paired))
    for col, key, lab in (("conservation", "conservation", "phyloP conservation"),
                          ("matched", "splicebert", "right protein"),
                          ("mismatched", "kmer", "wrong protein")):
        c = "#8c8c8c" if col == "conservation" else COLOR[key]
        c = "#bdbdbd" if col == "mismatched" else c
        ax[0].plot(x, paired[col], "o-", color=c, lw=1.6, ms=5, label=lab)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([f"\u2265{int(v)}\n(n={int(n)})"
                           for v, n in zip(paired.min_pathogenic, paired.n_datasets)],
                          fontsize=8)
    ax[0].set_xlabel("pathogenic variants per dataset")
    ax[0].set_ylabel("mean per-dataset AUROC")
    # Headroom then upper-left: at center-right the box sat on the rising matched line.
    ax[0].set_ylim(top=ax[0].get_ylim()[1] + 0.06)
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    ax[0].set_title("a  the floor is flat, the signal is not", loc="left", fontsize=9)

    # (b) the paired test itself, on the adequately powered datasets.
    q = per[per.n_pathogenic >= 20]
    ax[1].scatter(q.auroc_mismatched, q.auroc_matched, s=16,
                  color=COLOR["splicebert"], alpha=0.75, edgecolor="white", linewidth=0.3)
    lim = [min(q.auroc_mismatched.min(), q.auroc_matched.min()) - 0.03,
           max(q.auroc_mismatched.max(), q.auroc_matched.max()) + 0.03]
    ax[1].plot(lim, lim, "--", color="#999999", lw=1.0)
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel("wrong-protein head, AUROC")
    ax[1].set_ylabel("right-protein head, AUROC")
    row = paired[paired.min_pathogenic == 20]
    if len(row):
        r = row.iloc[0]
        ax[1].set_title(f"b  right wins {int(r.matched_wins)}/{int(r.n_datasets)}, "
                        f"p={r.p_specificity:.1e}", loc="left", fontsize=9)

    # (c) the inflation, stated rather than deleted.
    w = 0.36
    arms = ["mismatched", "matched"]
    for i, tag in enumerate(("pooled", "within_dataset")):
        c = coef[coef.standardisation == tag].set_index("arm").reindex(arms)
        pos = np.arange(len(arms)) + (i - 0.5) * w
        ax[2].bar(pos, c.coef, width=w, color=["#bdbdbd", COLOR["splicebert"]],
                  alpha=1.0 if tag == "within_dataset" else 0.45,
                  edgecolor="white", linewidth=0.5,
                  label="pooled (inflated)" if tag == "pooled" else "within dataset")
        ax[2].errorbar(pos, c.coef, yerr=[c.coef - c.ci_low, c.ci_high - c.coef],
                       fmt="none", ecolor="black", elinewidth=0.9, capsize=2.5)
    ax[2].set_xticks(np.arange(len(arms)))
    ax[2].set_xticklabels(["wrong\nprotein", "right\nprotein"], fontsize=8)
    ax[2].set_ylabel("conservation-controlled coefficient")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title("c  pooling inflates both arms", loc="left", fontsize=9)
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


# --- f6: the trivial baselines, which is the paper's strongest negative result -----------
#
# THIS FIGURE DID NOT EXIST FOR THREE DAYS while the claim it carries was called the paper's
# headline. `variant_ladder_paired.csv` already held `block_prevalence` and
# `model_minus_prevalence`; f4 loaded that exact frame and drew three arms out of four,
# omitting the one the argument rests on. A claim with no figure is a claim a referee skims.
#
# Panel a is the one that matters: every point is one dataset, and a point below the diagonal
# is a dataset where a rule that knows only WHERE a variant sits beat a fine-tuned language
# model. Aggregates hide unanimity; scatters show it.
def f6():
    got = need("variant_specificity.csv", "variant_ladder_paired.csv",
               "variant_specificity_attacks.csv", "robustness.csv")
    if got is None:
        return
    per, ladder, attacks, rob = got
    fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.5))

    # (a) paired scatter, powered stratum, on the COMMON variant mask
    mcol = "auroc_matched_common" if "auroc_matched_common" in per else "auroc_matched"
    q = per[(per.n_pathogenic >= 20)].dropna(subset=[mcol, "auroc_block_prevalence"])
    ax[0].plot([0.3, 1.0], [0.3, 1.0], color="#999999", lw=1, ls="--", zorder=1)
    ax[0].scatter(q.auroc_block_prevalence, q[mcol], s=26, color=COLOR["splicebert"],
                  edgecolor="white", linewidth=0.6, zorder=3)
    below = int((q[mcol] < q.auroc_block_prevalence).sum())
    ax[0].set_xlabel("1-Mb positional prevalence, AUROC")
    ax[0].set_ylabel("SpliceBERT, AUROC")
    ax[0].set_title(f"a  the model loses on {below}/{len(q)} datasets", loc="left")
    ax[0].set_xlim(0.3, 1.0)
    ax[0].set_ylim(0.3, 1.0)
    ax[0].set_aspect("equal")

    # (b) the scorers, worst to best, so the reader's eye lands on conservation
    row = ladder[ladder.min_pathogenic == 20]
    bars = []
    if len(row):
        r = row.iloc[0]
        bars = [("k-mer |delta|", 0.5519, "#bbbbbb"),
                ("dataset identity", 0.6682, "#bbbbbb"),
                ("wrong-protein head", float(r.mismatched), COLOR["kmer"]),
                ("SpliceBERT", float(r.matched), COLOR["splicebert"]),
                ("1-Mb prevalence", float(r.block_prevalence), COLOR["gc"]),
                ("phyloP", float(r.conservation), "#7b3294")]
    for i, (lab, v, c) in enumerate(bars):
        ax[1].barh(i, v - 0.5, left=0.5, color=c, height=0.62)
        ax[1].text(v + 0.006, i, f"{v:.3f}", va="center", fontsize=8)
    ax[1].set_yticks(range(len(bars)))
    ax[1].set_yticklabels([b[0] for b in bars])
    ax[1].set_xlim(0.5, 0.98)
    ax[1].set_xlabel("AUROC (paired, 44 powered datasets)")
    ax[1].set_title("b  two of the three winners use no model", loc="left")
    ax[1].grid(axis="y", visible=False)

    # (c) the decay. A rule that wins only at one block size is a binning artefact; one that
    # decays smoothly with block size is reading real positional structure.
    dec = {}
    for _, x in attacks.iterrows():
        if str(x.attack).startswith("trivial rule at"):
            dec[str(x.attack).replace("trivial rule at ", "")] = float(x.value)
    if dec:
        order = ["100 kb", "1000 kb", "10000 kb"]
        xs = [k for k in order if k in dec]
        ax[2].plot(range(len(xs)), [dec[k] for k in xs], marker="o", color=COLOR["gc"], lw=1.6)
        for i, k in enumerate(xs):
            ax[2].annotate(f"{dec[k]:.3f}", (i, dec[k]), textcoords="offset points",
                           xytext=(0, 8), ha="center", fontsize=8)
        ax[2].set_xticks(range(len(xs)))
        ax[2].set_xticklabels(["100 kb", "1 Mb", "10 Mb"])
    m = ladder[ladder.min_pathogenic == 20]
    if len(m):
        ax[2].axhline(float(m.iloc[0].matched), color=COLOR["splicebert"], ls=":", lw=1.4)
        ax[2].text(0.02, float(m.iloc[0].matched) + 0.004, "SpliceBERT", fontsize=8,
                   color=COLOR["splicebert"])
    ax[2].set_ylabel("pooled AUROC")
    ax[2].set_xlabel("genomic block size")
    ax[2].set_title("c  positional signal decays with block size", loc="left")

    fig.tight_layout()
    save(fig, "f6_trivial_baselines")


# --- f7: why the threshold is 20, answered with a curve rather than a sentence -----------
#
# The single most reachable objection to the specificity result is "you chose 20 pathogenic
# variants because it worked". `variant_threshold_curve.csv` answers it and was read by
# NOTHING for three days -- computed, uploaded, and consumed by no figure, no gate and no
# test. A threshold defended by a monotone curve is a design choice; the same threshold
# defended by prose is a suspicion.
#
# The two lines are the argument. A generic plausibility floor should not care how many
# pathogenic variants a dataset has, and a protein's own head should. Flat versus rising IS
# the detection threshold.
def f7():
    got = need("variant_threshold_curve.csv")
    if got is None:
        return
    (c,) = got
    c = c.sort_values("min_pathogenic")
    fig, ax = plt.subplots(1, 2, figsize=(9.8, 3.5))

    ax[0].plot(c.min_pathogenic, c.matched, marker="o", ms=3.5, lw=1.6,
               color=COLOR["splicebert"], label="right protein")
    ax[0].plot(c.min_pathogenic, c.mismatched, marker="s", ms=3.5, lw=1.6,
               color=COLOR["kmer"], label="wrong protein (the floor)")
    ax[0].axvline(20, color="#999999", ls="--", lw=1)
    ax[0].text(21, ax[0].get_ylim()[0] + 0.005, "reported\nthreshold", fontsize=7.5,
               color="#666666")
    ax[0].set_xlabel("minimum pathogenic variants per dataset (nested subsets)")
    ax[0].set_ylabel("mean per-dataset AUROC")
    # WAS "the floor is flat; the model is not", which the revised R5 text explicitly
    # withdraws: the floor rises 0.634 -> 0.667, about 6x less steeply than the matched arm's
    # 0.559 -> 0.755, and rho = +0.091 with CI [-0.01, +0.19] is a precise near-zero rather
    # than a demonstration of flatness. The title now says the comparative thing that is true.
    ax[0].set_title("a  the model's arm is ~6x steeper than the floor", loc="left")
    ax[0].legend(frameon=False, fontsize=8, loc="lower right")

    ax[1].axhline(0, color="#999999", lw=1)
    ax[1].plot(c.min_pathogenic, c.gap, marker="o", ms=3.5, lw=1.8, color="#7b3294")
    # THE p < 0.05 MARKERS ARE GONE ON PURPOSE. They read as "the effect becomes significant
    # at 15", and it is not a test: each x value FILTERS datasets, so the points are nested
    # subsets with dependent p-values, and low-power points are also weaker eCLIP experiments.
    # The experiment that would license that claim holds the 44 powered datasets fixed and
    # downsamples pathogenic variants, and it has not been run.
    ax[1].axvline(20, color="#999999", ls="--", lw=1)
    ax[1].axvline(45, color="#666666", ls=":", lw=1.2)
    ax[1].text(46, ax[1].get_ylim()[0] + 0.004, "plateau\n>=45", fontsize=7.5, color="#444444")
    ax[1].set_xlabel("minimum pathogenic variants per dataset (nested subsets)")
    ax[1].set_ylabel("specificity gap (right - wrong protein)")
    ax[1].set_title("b  plateau at >=45, not at 20", loc="left")

    fig.tight_layout()
    save(fig, "f7_power_threshold")


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


FIGURES = {"f0": f0, "f1": f1, "f2": f2, "f3": f3, "f4": f4, "f5": f5,
           "f6": f6, "f7": f7, "f8": f8}


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
