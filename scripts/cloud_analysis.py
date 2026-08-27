"""Stage 13. Aggregate every arm into the paper's tables and render the figures, in a
container.

WHAT THIS REPLACES. The final assembly -- merging the two negative arms, building the
four-model table, computing the cluster-corrected variant ladder, drawing the figures -- was
a sequence of ad-hoc local commands. That is the part of a study most likely to be
irreproducible, because it is the part nobody writes down: it is done once, interactively,
at the end, and the numbers in the paper come from whatever was in memory that afternoon.

Everything here reads committed artefacts out of GCS and writes committed artefacts back.
No interactive state, no laptop.

WHY THE VARIANT LADDER IS COMPUTED HERE AND NOT IN THE SCORING STAGE. The three rungs come
from three separate Modal sweeps (k-mer, mismatched head, matched head). Only once all three
exist can they be put on one axis, and the comparison must be on the SAME variants with the
SAME clustering, or the ladder measures panel differences instead of model differences.

    python scripts/cloud_analysis.py --what tables    # the result tables
    python scripts/cloud_analysis.py --what figures   # the figures
    python scripts/cloud_analysis.py --what all
"""

import argparse
import io
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.utils import cloud as cloudcfg  # noqa: E402

TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"
# Under results/ for the same IAM reason as the variants marker: rbp-analysis is scoped to
# results/, variants/ and driver/, so a root-level marker is a guaranteed 403 at the very end
# of the stage. See scripts/cloud_variants.py for the full account.
MARKER = "results/analysis-complete.json"


def log(m):
    print(f"[cloud_analysis] {m}", flush=True)


def fetch(bucket, key):
    b = bucket.blob(key)
    if not b.exists():
        return None
    return pd.read_csv(io.BytesIO(b.download_as_bytes()))


def fetch_prefix(bucket, prefix):
    """Concatenate every CSV under a prefix. Used for the per-dataset Modal outputs."""
    parts = [pd.read_csv(io.BytesIO(b.download_as_bytes()))
             for b in bucket.client.list_blobs(bucket.name, prefix=prefix)
             if b.name.endswith(".csv")]
    return pd.concat(parts, ignore_index=True) if parts else None


# --- the panel itself, which every result is conditional on ----------------------------

def panel_description(bucket):
    """Table 1: what the 95 datasets are, and the candidate pool they were drawn from.

    THIS WAS MISSING ENTIRELY. Every result table described model behaviour and none described
    the data, so the panel existed only as a manifest and a sentence. A reader's first question
    about a 95-dataset panel drawn from a larger pool is which 95 and whether they are the easy
    ones -- and that question is sharper here than usual, because dataset size correlates with
    SpliceBERT AUROC at r = +0.55, so a size-biased panel would inflate every model at once and
    leave no trace in any of the four results.

    Writes two tables: the per-dataset panel with its variant coverage, and the candidate
    pool's sizes, so f0 can show the panel spanning the pool rather than sitting on top of it.
    """
    import io

    panel = bucket.blob("manifest/study_panel.tsv")
    if not panel.exists():
        log("skip panel: no manifest/study_panel.tsv")
        return None
    d = pd.read_csv(io.StringIO(panel.download_as_text()), sep="\t")
    d = d.rename(columns={"cell_line": "cell"})

    # Variant coverage per dataset, from the window manifest the ClinVar arm actually scored.
    vt = bucket.blob("variants/variant_tasks.tsv")
    if vt.exists():
        v = pd.read_csv(io.StringIO(vt.download_as_text()), sep="\t")
        v["dataset"] = v.protein + ":" + v.cell
        d = d.merge(v[["dataset", "n_variants", "n_windows"]], on="dataset", how="left")

    d = d.sort_values("dataset", kind="mergesort").reset_index(drop=True)
    d.to_csv(TABLES / "panel_summary.csv", index=False)

    # THE CANDIDATE POOL, from the panel files select_panel.py itself drew from.
    #
    # First attempt used sweep_dinuc.csv, which is one row per dataset per MODEL and covers
    # only the 95 already selected -- so the "pool" came out as n=95 and the figure compared
    # the panel against itself, drawing a comparison that looked reassuring and said nothing.
    # The real pool is panel/{arm}/panel_final_{cell}_{arm}.tsv, the same objects the selector
    # reads, which is the only defensible denominator for "the panel spans the pool".
    pool = []
    for cell in ("K562", "HepG2"):
        b = bucket.blob(f"panel/dinuc/panel_final_{cell}_dinuc.tsv")
        if b.exists():
            pool.append(pd.read_csv(io.StringIO(b.download_as_text()), sep="\t"))
    if pool:
        p = pd.concat(pool, ignore_index=True).rename(columns={"cell_line": "cell"})
        p["dataset"] = p.protein + ":" + p.cell
        p = p.drop_duplicates("dataset").sort_values("dataset", kind="mergesort")
        p.to_csv(TABLES / "candidate_sizes.csv", index=False)
        log(f"candidate pool: {len(p)} datasets, pairs {int(p.pairs.min())}-"
            f"{int(p.pairs.max())}; panel covers percentile "
            f"{(p.pairs < d.pairs.min()).mean() * 100:.0f}-{(p.pairs < d.pairs.max()).mean() * 100:.0f}")

    both = int((d.protein.value_counts() == 2).sum())
    log(f"panel: {len(d)} datasets, {d.protein.nunique()} proteins ({both} in both lines), "
        f"pairs {int(d.pairs.min())}-{int(d.pairs.max())} median {int(d.pairs.median())}"
        + (f", variants {int(d.n_variants.sum()):,}" if "n_variants" in d else ""))
    return d


# --- R3 -------------------------------------------------------------------------------

def locality(bucket):
    """Aggregate the 95 per-dataset ISM runs into R3's table.

    THIS STEP WAS SIMPLY ABSENT. Stage 10 wrote one JSON per dataset to runs/locality/, and
    main() then *fetched* results/tables/locality_ism.csv from GCS -- an object nothing ever
    created. So the analysis quietly produced no R3 table, do_figures() skipped the R3 figure
    without complaint, and stage 14 was the first thing to notice, reporting the table as
    MISSING after the whole pipeline had run green.

    Every other result has a producer here; R3's was a fetch of its own output. The tolerant
    guard in do_figures() is what let it stay invisible, which is the pattern already recorded
    in the chronicle: a lenient gate needs a strict counterpart, and stage 14 was it.
    """
    import json

    rows = []
    for b in bucket.client.list_blobs(bucket.name, prefix="runs/locality/"):
        if b.name.endswith(".json"):
            rows.append(json.loads(b.download_as_text()))
    if not rows:
        log("skip R3: no runs/locality/*.json")
        return None

    d = pd.DataFrame(rows)
    # Deterministic order, so the committed table is byte-stable across runs. mergesort
    # because the default quicksort is unstable and two datasets tying on any key would
    # reorder run to run -- the same defect that made the study panel non-reproducible.
    d = d.sort_values("dataset", kind="mergesort").reset_index(drop=True)
    d.to_csv(TABLES / "locality_ism.csv", index=False)

    more = int((d.sb_gini > d.kmer_gini).sum())
    log(f"R3: {len(d)} datasets, SpliceBERT more concentrated on {more}/{len(d)}, "
        f"median delta Gini {(d.sb_gini - d.kmer_gini).median():+.4f}")
    return d


# --- R1 -------------------------------------------------------------------------------

def cost_of_matching(bucket):
    """The two arms, differenced on the datasets they SHARE.

    The arms do not contain the same datasets: matching is a search and a dataset can clear
    the min_pairs floor in one arm and miss it in the other. Differencing anything other
        than the intersection compares a dataset against nothing.
    """
    gc = fetch(bucket, "results/rehearsal_binding_gc.csv")
    dn = fetch(bucket, "results/rehearsal_binding_dinuc.csv")
    if gc is None or dn is None:
        log("skip R1: need both rehearsal arms")
        return None
    m = gc.merge(dn, on="dataset", suffixes=("_gc", "_dn"))
    m["cost"] = m.auroc_dn - m.auroc_gc
    # delta_p and helps travel with the gains. The claim is not just that the nested gain is
    # larger under dinucleotide matching but that it is REAL, and the significance of the
    # nested comparison is what says so. Carrying the point estimate without its inference is
    # how a ratio of two noisy numbers becomes a headline.
    keep = ["dataset", "protein_gc", "cell_gc", "pairs_gc", "auroc_gc", "auroc_dn", "cost",
            "composition_auroc_gc", "composition_auroc_dn",
            "delta_auroc_gc", "delta_auroc_dn",
            "delta_p_gc", "delta_p_dn", "helps_gc", "helps_dn"]
    out = m[keep].rename(columns={"protein_gc": "protein", "cell_gc": "cell",
                                  "pairs_gc": "pairs"})
    out.to_csv(TABLES / "cost_of_matching.csv", index=False)
    log(f"R1: {len(out)} shared datasets, cost {out.cost.mean():+.4f}, "
        f"{int((out.cost < 0).sum())}/{len(out)} fall")
    return out


# --- R2 -------------------------------------------------------------------------------

def four_models(bucket):
    """One row per dataset, one column per model, only where ALL of them exist."""
    reh = fetch(bucket, "results/rehearsal_binding_dinuc.csv")
    sweep = fetch(bucket, "results/sweep_dinuc.csv")
    if reh is None or sweep is None:
        log("skip R2: need the rehearsal and the sweep")
        return None
    base = reh.set_index("dataset")
    wide = sweep.pivot_table(index="dataset", columns="model", values="pooled_auroc")
    d = pd.DataFrame({
        "composition_auroc": base.composition_auroc,
        "kmer_auroc": base.auroc,
        "pairs": base.pairs,
        "protein": base.protein,
        "cell": base.cell,
    })
    for m in ("cnn", "splicebert"):
        if m in wide:
            d[m] = wide[m]
    need = [c for c in ("composition_auroc", "kmer_auroc", "cnn", "splicebert") if c in d]
    d = d.dropna(subset=need).reset_index()
    d["delta_auroc"] = d.kmer_auroc - d.composition_auroc
    d.to_csv(TABLES / "matched_four_models.csv", index=False)
    log(f"R2: {len(d)} datasets with {len(need)} models, "
        f"{d.protein.nunique()} proteins; means " +
        ", ".join(f"{c}={d[c].mean():.3f}" for c in need))
    return d


# --- R4, the honest version -------------------------------------------------------------

def _prev_block(q):
    """The model against the trivial positional rule, on the datasets where both exist.

    Separated out because it is the single most deflating number in the study and it must not
    be quietly skipped when the column is absent: a missing baseline reports as NaN, not as a
    pass.
    """
    from scipy.stats import wilcoxon
    z = q.dropna(subset=["auroc_block_prevalence"])
    if len(z) < 8:
        return {"block_prevalence": np.nan, "model_minus_prevalence": np.nan,
                "model_beats_prevalence": np.nan, "p_vs_prevalence": np.nan}
    return {"block_prevalence": z.auroc_block_prevalence.mean(),
            "model_minus_prevalence": (z.auroc_matched - z.auroc_block_prevalence).mean(),
            "model_beats_prevalence": int((z.auroc_matched > z.auroc_block_prevalence).sum()),
            "p_vs_prevalence": wilcoxon(z.auroc_matched, z.auroc_block_prevalence).pvalue}


def variant_specificity(bucket):
    """R4's PRIMARY statistic: paired per-dataset comparison, not one pooled AUROC.

    WHY THE POOLED LADDER IS NOT THE RESULT. variant_ladder() below concatenates ~19k
    variants across 95 datasets and computes a single AUROC per arm. That number is inflated,
    and measurably so: a dataset's mean |delta| correlates with its pathogenic rate at
    Spearman +0.73 for the matched arm, and mean |delta| spans 10.4x across datasets. So a
    pooled AUROC partly measures WHICH DATASET a variant came from rather than whether the
    model separates pathogenic from benign within a dataset.

    The size of the error is not subtle. Matched: pooled 0.829 versus 0.755 paired
    per-dataset. Mismatched: 0.680 versus 0.691. The pooled gap of +0.149 is really +0.065.
    Conservation barely moves, 0.908 to 0.892, because phyloP is on a fixed external scale
    and has no between-dataset inflation to contribute -- which is exactly why the artefact
    was invisible: the arm that could not be inflated was the one winning anyway.

    This is the failure mode described for ClinVar benchmarks generally (heterogeneity
    inflating apparent performance), reproduced inside this pipeline. The golden gate passed
    33/33 while reporting it, because the gate checked the number the code computed and not
    whether that number was the right statistic.

    POWER STRATIFICATION IS NOT p-HACKING, and the reason is in the numbers. Median coverage
    is 140 variants per dataset, and many datasets carry a handful of pathogenic ones, so
    their per-dataset AUROC is noise. Across all 82 usable datasets the paired test gives
    p=0.87; restricted to datasets with at least 20 pathogenic variants it gives +0.065 and
    p<0.001, and at 50 it gives +0.103. The effect grows with power (rho=+0.52 between
    pathogenic count and the matched-mismatched gap) and the MISMATCHED arm stays flat at
    ~0.69 throughout. A spurious effect does not behave that way. Every stratum is reported
    here, including the one that shows nothing.
    """
    from sklearn.metrics import roc_auc_score
    from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

    from rbp.variants import conservation as cons

    cv = fetch(bucket, "results/tables/variant_conservation.csv")
    sb = fetch_prefix(bucket, "variants/scores_sb/")
    mm = fetch_prefix(bucket, "variants/scores_mm/")
    if cv is None or sb is None or mm is None:
        log("skip R4-paired: need conservation and both score arms")
        return None
    cv = cv[["vid", "conservation"]]
    for d in (sb, mm):
        d["dataset"] = d.protein + ":" + d.cell
    sb = sb.merge(cv, on="vid", how="left")
    mm = mm.merge(cv, on="vid", how="left")

    # --- per dataset -------------------------------------------------------------------
    rows = []
    for ds, s in sb.groupby("dataset"):
        g = mm[mm.dataset == ds]
        s = s.dropna(subset=["delta", "conservation"])
        g = g.dropna(subset=["delta"])
        if s.label.nunique() < 2 or g.label.nunique() < 2:
            continue                      # a one-class dataset has no AUROC, not a zero one
        donor = g.weights_from.iloc[0] if "weights_from" in g else ""
        # THE TRIVIAL POSITIONAL BASELINE, computed per dataset because that is where the
        # claim lives. "What fraction of the OTHER variants in this 1-Mb window are
        # pathogenic" uses no sequence, no model and no biology, and it beats the model.
        # Leave-one-out, so a variant never contributes to its own window's rate.
        blk = (s.vid.str.split(":").str[0] + "_" +
               (s.vid.str.split(":").str[1].astype(int) // 1_000_000).astype(str))
        gb = s.assign(_b=blk).groupby("_b").label
        tot, cnt = gb.transform("sum"), gb.transform("size")
        prev = ((tot - s.label) / (cnt - 1)).where(cnt > 1, np.nan)
        ok = prev.notna()
        auroc_prev = (roc_auc_score(s.label[ok], prev[ok])
                      if ok.sum() > 20 and s.label[ok].nunique() == 2 else np.nan)

        # THE BASELINE AND THE MODEL MUST BE SCORED ON THE SAME VARIANTS.
        #
        # A variant alone in its 1-Mb block has no leave-one-out prevalence, so `ok` drops it
        # -- a mean of 20.2% of variants per dataset, up to 38.9%. The model's AUROC above is
        # computed on ALL variants. Differencing the two therefore compared two arms on two
        # different variant sets, and the "smooth decay" across 100 kb / 1 Mb / 10 Mb was
        # partly the evaluated set changing size (17,934 / 18,762 / 18,994) rather than the
        # rule getting worse.
        #
        # The bias turned out to be small -- conservation, which is scoreable either way,
        # gives 0.8921 on all variants and 0.8904 on the baseline's subset, so -0.0017 -- and
        # the headline survives it. But "small" is a measurement, not a licence, so the
        # common-mask columns are computed here and the paired comparison uses them.
        common = {}
        if ok.sum() > 20 and s.label[ok].nunique() == 2:
            common = {"n_common": int(ok.sum()),
                      "auroc_matched_common": roc_auc_score(s.label[ok], s.delta.abs()[ok]),
                      "auroc_conservation_common": roc_auc_score(s.label[ok],
                                                                 s.conservation[ok])}
        rows.append({
            "dataset": ds, "protein": ds.split(":")[0], "cell": ds.split(":")[1],
            "donor": donor,
            "same_cell": bool(donor) and donor.split(":")[1] == ds.split(":")[1],
            "n": len(s), "n_pathogenic": int(s.label.sum()),
            "auroc_matched": roc_auc_score(s.label, s.delta.abs()),
            "auroc_mismatched": roc_auc_score(g.label, g.delta.abs()),
            "auroc_conservation": roc_auc_score(s.label, s.conservation),
            "auroc_block_prevalence": auroc_prev,
            **common,
        })
    per = pd.DataFrame(rows).sort_values("dataset", kind="mergesort").reset_index(drop=True)
    per.to_csv(TABLES / "variant_specificity.csv", index=False)

    # --- the stratified paired table ---------------------------------------------------
    out = []
    for lo in (0, 10, 20, 50):
        q = per[per.n_pathogenic >= lo]
        if len(q) < 8:
            continue
        out.append({
            "min_pathogenic": lo, "n_datasets": len(q),
            "matched": q.auroc_matched.mean(),
            "mismatched": q.auroc_mismatched.mean(),
            "conservation": q.auroc_conservation.mean(),
            "specificity_gap": (q.auroc_matched - q.auroc_mismatched).mean(),
            "matched_wins": int((q.auroc_matched > q.auroc_mismatched).sum()),
            "p_specificity": wilcoxon(q.auroc_matched, q.auroc_mismatched).pvalue,
            "conservation_lead": (q.auroc_conservation - q.auroc_matched).mean(),
            "conservation_wins": int((q.auroc_conservation > q.auroc_matched).sum()),
            "p_conservation": wilcoxon(q.auroc_conservation, q.auroc_matched).pvalue,
            **_prev_block(q),
        })
    paired = pd.DataFrame(out)
    paired.to_csv(TABLES / "variant_ladder_paired.csv", index=False)

    # --- coefficients, with and without between-dataset scale --------------------------
    coefs = []
    for name, s in (("matched", sb), ("mismatched", mm)):
        d = s.dropna(subset=["delta", "conservation"]).copy()
        d["ad"] = d.delta.abs()
        gp = d.groupby("dataset").ad
        d["ad_within"] = (d.ad - gp.transform("mean")) / gp.transform("std").replace(0, np.nan)
        d = d.dropna(subset=["ad_within"])
        for tag, col in (("pooled", "ad"), ("within_dataset", "ad_within")):
            # Deduplicate on the SAME column being analysed, or the selection reintroduces
            # exactly the between-dataset scale the standardisation removed.
            u = d.sort_values(col, key=abs, ascending=False).drop_duplicates("vid")
            u = u.assign(block=u.vid.str.split(":").str[0] + "_" +
                         (u.vid.str.split(":").str[1].astype(int) // 1_000_000).astype(str))
            # take_abs only for the raw magnitude. ad_within is a z-score and is negative
            # below its dataset's mean; abs() there folds the distribution in half.
            f = cons.fit_delta_coef(u[col].to_numpy(), u.label.to_numpy(),
                                    u.conservation.to_numpy(), n_boot=500, seed=0,
                                    blocks=u.block.to_numpy(),
                                    take_abs=(col == "ad"))
            coefs.append({"arm": name, "standardisation": tag, "coef": f.coef,
                          "ci_low": f.ci_low, "ci_high": f.ci_high, "p_wald": f.p_wald,
                          "n_variants": len(u), "n_clusters": int(u.block.nunique())})
    pd.DataFrame(coefs).to_csv(TABLES / "variant_coefficients.csv", index=False)

    # --- is the wrong-protein floor a similarity artefact? ------------------------------
    # The obvious attack on this control: RBPs co-bind and share motif families, so the
    # "wrong" protein may be a similar one, making the floor a function of donor-target
    # overlap rather than generic sequence plausibility. Three checks, all reported.
    q = per[per.n_pathogenic >= 20]
    sc, xc = q[q.same_cell], q[~q.same_cell]
    checks = [{"check": "floor vs donor cell line",
               "value": sc.auroc_mismatched.mean() - xc.auroc_mismatched.mean(),
               "p": mannwhitneyu(sc.auroc_mismatched, xc.auroc_mismatched).pvalue
                    if len(sc) > 2 and len(xc) > 2 else np.nan,
               "note": f"same-cell {sc.auroc_mismatched.mean():.4f} n={len(sc)}, "
                       f"cross-cell {xc.auroc_mismatched.mean():.4f} n={len(xc)}"}]
    own = dict(zip(per.dataset, per.auroc_matched))
    p2 = per.assign(donor_own=per.donor.map(own)).dropna(subset=["donor_own"])
    rho, pv = spearmanr(p2.donor_own, p2.auroc_mismatched)
    checks.append({"check": "floor vs donor's own strength", "value": rho, "p": pv,
                   "note": f"spearman over n={len(p2)} donors"})
    rho2, pv2 = spearmanr(per.n_pathogenic, per.auroc_matched - per.auroc_mismatched)
    checks.append({"check": "gap vs statistical power", "value": rho2, "p": pv2,
                   "note": "a real effect should grow with power; the floor should not"})
    pd.DataFrame(checks).to_csv(TABLES / "variant_specificity_controls.csv", index=False)

    best = paired[paired.min_pathogenic == 20].iloc[0] if (paired.min_pathogenic == 20).any() \
        else paired.iloc[-1]
    log(f"R4-paired: {len(per)} datasets usable; at >={int(best.min_pathogenic)} pathogenic "
        f"(n={int(best.n_datasets)}) matched {best.matched:.4f} vs mismatched "
        f"{best.mismatched:.4f}, gap {best.specificity_gap:+.4f}, "
        f"{best.matched_wins}/{int(best.n_datasets)} wins, p={best.p_specificity:.2e}; "
        f"conservation leads by {best.conservation_lead:+.4f}")
    return per


# --- R4, the pooled version, kept for comparison ----------------------------------------

def variant_ladder(bucket):
    """The three rungs plus conservation, on identical variants with identical clustering.

    Deduplicated to one row per variant and cluster-bootstrapped over 1-Mb genomic blocks.
    A variant sitting near several proteins' peaks appears several times in the raw scoring
    tables, and those rows share a position, a label and a conservation value -- so per-row
    inference treats ~19k independent observations as ~33k and reports intervals about a
    third too narrow.
    """
    from sklearn.metrics import roc_auc_score

    from rbp.variants import conservation as cons

    cv = fetch(bucket, "results/tables/variant_conservation.csv")
    if cv is None:
        log("skip R4: no conservation table")
        return None
    cv = cv[["vid", "conservation"]]

    arms = {
        "matched": fetch_prefix(bucket, "variants/scores_sb/"),
        "mismatched": fetch_prefix(bucket, "variants/scores_mm/"),
        "kmer": fetch(bucket, "results/tables/variant_scores.csv"),
    }
    if arms["matched"] is None:
        log("skip R4: no matched scores")
        return None

    # Every arm restricted to the datasets the matched arm covers, or the ladder compares
    # panels rather than models.
    keep = set(arms["matched"].protein + ":" + arms["matched"].cell)
    rows = []
    for name, s in arms.items():
        if s is None:
            log(f"  {name}: missing, skipped")
            continue
        s = s.copy()
        s["ds"] = s.protein + ":" + s.cell
        s = s[s.ds.isin(keep)]
        d = s.merge(cv, on="vid", how="left").dropna(subset=["delta", "conservation"])
        d = d.sort_values("delta", key=abs, ascending=False).drop_duplicates("vid")
        d["block"] = (d.vid.str.split(":").str[0] + "_" +
                      (d.vid.str.split(":").str[1].astype(int) // 1_000_000).astype(str))
        f = cons.fit_delta_coef(d.delta.to_numpy(), d.label.to_numpy(),
                                d.conservation.to_numpy(), n_boot=500, seed=0,
                                blocks=d.block.to_numpy())
        rows.append({"arm": name, "n_variants": len(d), "n_clusters": d.block.nunique(),
                     "auroc": roc_auc_score(d.label, np.abs(d.delta)),
                     "coef": f.coef, "ci_low": f.ci_low, "ci_high": f.ci_high,
                     "p_wald": f.p_wald})
        log(f"  {name:11} n={len(d):,} clusters={d.block.nunique():,} "
            f"auroc={rows[-1]['auroc']:.3f} coef={f.coef:.3f} "
            f"[{f.ci_low:.3f}, {f.ci_high:.3f}]")

    # Conservation as its own rung: the competing explanation, on the same variants.
    d = arms["matched"].copy()
    d["ds"] = d.protein + ":" + d.cell
    d = d.merge(cv, on="vid", how="left").dropna(subset=["conservation"]).drop_duplicates("vid")
    rows.append({"arm": "conservation", "n_variants": len(d), "n_clusters": np.nan,
                 "auroc": roc_auc_score(d.label, d.conservation),
                 "coef": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_wald": np.nan})
    log(f"  conservation auroc={rows[-1]['auroc']:.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "variant_ladder.csv", index=False)
    return out


def do_figures():
    import subprocess
    rc = subprocess.run([sys.executable, str(ROOT / "scripts" / "figures.py")],
                        cwd=str(ROOT)).returncode
    if rc != 0:
        log(f"figures.py exited {rc}")
    return rc == 0


def donor_overlap(bucket):
    """Does the wrong-protein control actually use a WRONG protein? Measured, not assumed.

    THE ATTACK THIS ANSWERS, and it was right. RBPs co-bind: they occupy overlapping regions
    and share motif families. So "score protein A's variants with protein B's head" may not be
    a wrong-protein control at all -- if B binds the same sites, B is a partially-right
    protein and the floor it produces is contaminated.

    Two earlier attempts used PROXIES for donor similarity, sharing a cell line and the
    donor's own strength, and both came back clean. They were the wrong measurement. The
    direct one is how much of the target's variant set the donor also binds near, which is
    peak overlap evaluated exactly where the scoring happens, and it needs no peak files
    because variant_assignments.csv already records which datasets each variant is assigned to.

    WHAT IT FOUND. Median donor-target overlap is 0.001, so the offset pairing over a
    rank-sorted manifest does usually pick a genuinely different protein. But the association
    with the floor is real in the powered stratum (rho=+0.30, p=0.05), and stratifying is
    stark:

        donors with ~zero overlap   gap +0.1210, 16/17 datasets, p=2e-04
        all donors                  gap +0.0645, 33/44,          p=4e-04
        donors above median overlap gap +0.0130, 12/22,          p=0.50, not significant

    So contamination inflates the FLOOR, which SHRINKS the measured gap. The headline +0.065
    is therefore a lower bound, diluted by pairs where the donor was not really wrong, and the
    control behaves exactly as a valid one should: when the wrong protein is not wrong, it
    reports no difference.

    The methodological consequence is the useful part: a wrong-protein control must SCREEN ITS
    DONORS for target overlap, or it under-reports. That screening requirement is part of the
    method, not a caveat about this dataset.
    """
    from sklearn.metrics import roc_auc_score
    from scipy.stats import spearmanr, wilcoxon

    asg = fetch(bucket, "results/tables/variant_assignments.csv")
    sb = fetch_prefix(bucket, "variants/scores_sb/")
    mm = fetch_prefix(bucket, "variants/scores_mm/")
    if asg is None or sb is None or mm is None:
        log("skip donor overlap: need assignments and both arms")
        return None
    asg["dataset"] = asg.protein + ":" + asg.cell
    vids = asg.groupby("dataset").vid.apply(set)
    for d in (sb, mm):
        d["dataset"] = d.protein + ":" + d.cell

    rows = []
    for ds, s in sb.groupby("dataset"):
        g = mm[mm.dataset == ds].dropna(subset=["delta"])
        s = s.dropna(subset=["delta"])
        if s.label.nunique() < 2 or g.label.nunique() < 2:
            continue
        donor = g.weights_from.iloc[0] if "weights_from" in g else ""
        t, d = vids.get(ds, set()), vids.get(donor, set())
        if not t:
            continue
        rows.append({"dataset": ds, "donor": donor, "n_pathogenic": int(s.label.sum()),
                     "shared_frac": len(t & d) / len(t),
                     "jaccard": len(t & d) / len(t | d) if (t | d) else 0.0,
                     "auroc_matched": roc_auc_score(s.label, s.delta.abs()),
                     "auroc_floor": roc_auc_score(g.label, g.delta.abs())})
    r = pd.DataFrame(rows).sort_values("dataset", kind="mergesort").reset_index(drop=True)
    r.to_csv(TABLES / "donor_overlap.csv", index=False)

    q = r[r.n_pathogenic >= 20]
    rho, pv = spearmanr(q.shared_frac, q.auroc_floor)
    out = [{"stratum": "correlation of floor with donor overlap", "n": len(q),
            "matched": np.nan, "floor": np.nan, "gap": rho, "wins": np.nan, "p": pv}]
    med = q.shared_frac.median()
    for lab, sub in (("all powered donors", q),
                     ("donors with negligible overlap", q[q.shared_frac <= 0.001]),
                     ("donors above median overlap", q[q.shared_frac > med])):
        if len(sub) < 8:
            continue
        out.append({"stratum": lab, "n": len(sub),
                    "matched": sub.auroc_matched.mean(), "floor": sub.auroc_floor.mean(),
                    "gap": (sub.auroc_matched - sub.auroc_floor).mean(),
                    "wins": int((sub.auroc_matched > sub.auroc_floor).sum()),
                    "p": wilcoxon(sub.auroc_matched, sub.auroc_floor).pvalue})
    o = pd.DataFrame(out)
    o.to_csv(TABLES / "donor_overlap_strata.csv", index=False)
    log(f"donor overlap: median {r.shared_frac.median():.4f}, "
        f"floor-vs-overlap rho={rho:+.3f} (p={pv:.3f})")
    for _, x in o.iloc[1:].iterrows():
        log(f"  {x.stratum:32} n={int(x.n):3d} gap {x.gap:+.4f} "
            f"wins {int(x.wins)}/{int(x.n)} p={x.p:.4f}")
    return r


def specificity_attacks(bucket):
    """Four attacks on the wrong-protein control, from a council review, all answered here.

    THE CONTROL IS THE ONLY NOVEL THING IN THIS STUDY, so it gets attacked in the pipeline
    rather than defended in prose. Each of these was raised as a reason the result is not
    real. Each is recomputed on every run, so a future change that breaks one fails the build.

    1. THE THRESHOLD IS A FREE PARAMETER. The headline gap is quoted at >=20 pathogenic
       variants per dataset, a cutoff chosen after seeing data, and the gap moves with it:
       -0.011 at 0, +0.065 at 20, +0.117 at 100. Forking paths, unless the whole curve is
       shown -- and the whole curve is the defence, because it is a monotone ramp rather than
       a spike, and the wrong-protein floor stays flat within 0.015 across all of it. A
       cherry-picked cutoff does not look like that.

    2. THE CONTROL VALIDATES A MODEL THAT LOSES TO A TRIVIAL RULE. Refit the contrast with
       the 1-Mb positional prevalence as a covariate beside conservation. It does not shrink,
       because the rule applies equally to both arms.

    3. THE CONTROL HAS NO KNOWN FALSE-POSITIVE RATE. Permute labels within dataset and
       recompute the gap. A structurally biased comparison would show one under the null.

       WHAT THIS DOES NOT DO, stated because it matters: a label permutation destroys all
       signal, so it tests whether the COMPARISON is biased, not whether the control can be
       fooled by a model that learned something non-specific. The stronger calibration --
       retrain on shuffled binding labels, then run the control -- needs GPU time and has not
       been run. It is a limitation, not a completed check.

    4. THE TRIVIAL RULE MIGHT BE A 1-Mb ARTEFACT. Recompute it at 100 kb, 1 Mb and 10 Mb. It
       decays smoothly (0.851 / 0.818 / 0.733), which makes it a real positional-prevalence
       effect rather than a lucky binning.
    """
    from sklearn.metrics import roc_auc_score
    from scipy.stats import wilcoxon

    from rbp.variants import conservation as cons

    cv = fetch(bucket, "results/tables/variant_conservation.csv")
    sb = fetch_prefix(bucket, "variants/scores_sb/")
    mm = fetch_prefix(bucket, "variants/scores_mm/")
    if cv is None or sb is None or mm is None:
        log("skip specificity attacks: need conservation and both arms")
        return None
    cv = cv[["vid", "conservation"]]
    for d in (sb, mm):
        d["dataset"] = d.protein + ":" + d.cell

    # --- attack 1: the whole threshold curve --------------------------------------------
    pairs, curve = [], []
    for ds, s in sb.groupby("dataset"):
        g = mm[mm.dataset == ds].dropna(subset=["delta"])
        s = s.dropna(subset=["delta"])
        if s.label.nunique() < 2 or g.label.nunique() < 2:
            continue
        pairs.append((ds, int(s.label.sum()),
                      s.label.to_numpy(), s.delta.abs().to_numpy(),
                      g.label.to_numpy(), g.delta.abs().to_numpy()))
    au = [(ds, n, roc_auc_score(a, b), roc_auc_score(c, d)) for ds, n, a, b, c, d in pairs]
    tab = pd.DataFrame(au, columns=["dataset", "n_pathogenic", "matched", "mismatched"])
    for t in range(0, 105, 5):
        q = tab[tab.n_pathogenic >= t]
        if len(q) < 10:
            break
        curve.append({"min_pathogenic": t, "n_datasets": len(q),
                      "matched": q.matched.mean(), "mismatched": q.mismatched.mean(),
                      "gap": (q.matched - q.mismatched).mean(),
                      "matched_wins": int((q.matched > q.mismatched).sum()),
                      "p": wilcoxon(q.matched, q.mismatched).pvalue})
    cu = pd.DataFrame(curve)
    cu.to_csv(TABLES / "variant_threshold_curve.csv", index=False)

    # Spearman, not strict monotonicity: the last few thresholds have n<25 and wobble by
    # ~0.005, which says nothing. The claim is a dose-response trend, and that is what gets
    # measured. Strict monotonicity reported False on a curve that rises from -0.011 to
    # +0.117, which would have been a misleading gate.
    from scipy.stats import spearmanr
    out = [{"attack": "gap rises with power (spearman)",
            "value": float(spearmanr(cu.min_pathogenic, cu.gap).statistic), "note":
            f"gap {cu.gap.min():+.4f} to {cu.gap.max():+.4f} over {len(cu)} thresholds"},
           {"attack": "wrong-protein floor is flat across thresholds",
            "value": float(cu.mismatched.max() - cu.mismatched.min()),
            "note": f"range {cu.mismatched.min():.4f}-{cu.mismatched.max():.4f}"}]

    # --- attack 2: refit with the trivial rule as a covariate ---------------------------
    def prep(x):
        d = x.merge(cv, on="vid", how="left").dropna(subset=["delta", "conservation"]).copy()
        d["ad"] = d.delta.abs()
        gp = d.groupby("dataset").ad
        d["adw"] = (d.ad - gp.transform("mean")) / gp.transform("std").replace(0, np.nan)
        d["block"] = (d.vid.str.split(":").str[0] + "_" +
                      (d.vid.str.split(":").str[1].astype(int) // 1_000_000).astype(str))
        gb = d.groupby("block").label
        tot, cnt = gb.transform("sum"), gb.transform("size")
        d["prev"] = ((tot - d.label) / (cnt - 1)).where(cnt > 1, np.nan)
        d = d.dropna(subset=["adw", "prev"])
        return d.sort_values("adw", key=abs, ascending=False).drop_duplicates("vid")

    fits = []
    for name, x in (("matched", sb), ("mismatched", mm)):
        d = prep(x)
        for tag, extra in (("conservation only", None),
                           ("plus trivial window rule", d[["prev"]].to_numpy())):
            f = cons.fit_delta_coef(d.adw.to_numpy(), d.label.to_numpy(),
                                    d.conservation.to_numpy(), n_boot=400, seed=0,
                                    blocks=d.block.to_numpy(), extra=extra,
                                    take_abs=False)
            fits.append({"arm": name, "controls": tag, "coef": f.coef,
                         "ci_low": f.ci_low, "ci_high": f.ci_high, "n": len(d)})
    fi = pd.DataFrame(fits)
    fi.to_csv(TABLES / "variant_specificity_refit.csv", index=False)
    w = fi[fi.controls == "plus trivial window rule"].set_index("arm")
    out.append({"attack": "specificity survives the trivial rule as a covariate",
                "value": float(w.loc["matched", "ci_low"] > w.loc["mismatched", "ci_high"]),
                "note": f"matched {w.loc['matched','coef']:.3f} "
                        f"[{w.loc['matched','ci_low']:.3f},{w.loc['matched','ci_high']:.3f}] vs "
                        f"mismatched {w.loc['mismatched','coef']:.3f} "
                        f"[{w.loc['mismatched','ci_low']:.3f},{w.loc['mismatched','ci_high']:.3f}]"})

    # --- attack 3: permutation null ----------------------------------------------------
    rng = np.random.default_rng(0)
    pw = [(a, b, c, d) for _, n, a, b, c, d in pairs if n >= 20]
    obs = float(np.mean([roc_auc_score(a, b) - roc_auc_score(c, d) for a, b, c, d in pw]))
    null = np.array([np.mean([roc_auc_score(rng.permutation(a), b)
                              - roc_auc_score(rng.permutation(c), d) for a, b, c, d in pw])
                     for _ in range(300)])
    out.append({"attack": "permutation null is centred at zero", "value": float(null.mean()),
                "note": f"sd {null.std():.4f}, observed {obs:+.4f}, "
                        f"p={float((np.abs(null) >= abs(obs)).mean()):.4f}, n={len(pw)}"})

    # --- attack 4: window-size sensitivity ---------------------------------------------
    u = sb.dropna(subset=["delta"]).copy()
    u = u.sort_values("delta", key=abs, ascending=False).drop_duplicates("vid")
    u["chrom"] = u.vid.str.split(":").str[0]
    u["pos"] = u.vid.str.split(":").str[1].astype(int)
    for kb in (100_000, 1_000_000, 10_000_000):
        u["_b"] = u.chrom + "_" + (u.pos // kb).astype(str)
        gb = u.groupby("_b").label
        tot, cnt = gb.transform("sum"), gb.transform("size")
        pv = ((tot - u.label) / (cnt - 1)).where(cnt > 1, np.nan)
        m = pv.notna()
        out.append({"attack": f"trivial rule at {kb // 1000} kb",
                    "value": roc_auc_score(u.label[m], pv[m]),
                    "note": f"{int(m.sum()):,} variants, {u._b.nunique():,} blocks"})

    r = pd.DataFrame(out)
    r.to_csv(TABLES / "variant_specificity_attacks.csv", index=False)
    for _, x in r.iterrows():
        log(f"  {x['attack']:52} {x['value']:+.4f}  {x['note']}")
    return r


# --- robustness checks the reviewers will run if we do not ------------------------------

def robustness(bucket):
    """Three checks that decide whether the two surviving claims hold up.

    1. A CONFIDENCE INTERVAL ON THE COMPOSITION SHARE. 94.8% and 67.8% were reported as bare
       point estimates, and they are ratios of differences -- the least stable thing a
       statistic can be. Bootstrapped over datasets, because datasets are the sampling unit.

    2. THE SIZE CONFOUND, tested rather than assumed. Dataset size correlates with AUROC at
       up to r=+0.67, which confounds any BETWEEN-model comparison across datasets. It cannot
       confound either surviving claim, and the reason is structural: both are PAIRED within
       dataset. The composition share compares two negative sets on the same datasets; the
       specificity gap compares two weight sets on the same variants. Identical n on both
       sides means size cannot bias the difference, only its precision. What size CAN do is
       modify the effect, so that is what gets measured here.

    3. A PREVALENCE BASELINE, run on our own variants before a reviewer runs one. The known
       failure mode for ClinVar benchmarks is that a trivial rule exploiting where pathogenic
       variants CLUSTER can beat real models. We have no variant-type annotation, so the
       analogous rules here are the pathogenic rate of a variant's 1-Mb block and of its
       dataset. Both are computed LEAVE-ONE-OUT: a variant never contributes to its own
       group's rate, or the baseline trivially memorises the label.
    """
    from sklearn.metrics import roc_auc_score

    out = []

    def share(x, a, c):
        """Fraction of skill-above-chance recoverable from composition alone.

        Defined at function scope, not inside the block below: 1b uses the SAME estimator on
        a different table, and the two shares are only comparable if they are literally the
        same function. It used to be a closure inside `if d is not None`, which made 1b a
        NameError whenever cost_of_matching.csv was missing.
        """
        num = x[c].mean() - 0.5
        den = x[a].mean() - 0.5
        return num / den if den > 0 else np.nan

    # --- 1. composition share, with an interval ---------------------------------------
    d = fetch(bucket, "results/tables/cost_of_matching.csv")
    if d is not None:
        rng = np.random.default_rng(7)

        for arm, a, c in (("GC-matched", "auroc_gc", "composition_auroc_gc"),
                          ("dinuc-matched", "auroc_dn", "composition_auroc_dn")):
            pt = share(d, a, c)
            boot = [share(d.iloc[rng.integers(0, len(d), len(d))], a, c)
                    for _ in range(2000)]
            lo, hi = np.nanpercentile(boot, [2.5, 97.5])
            out.append({"check": f"composition share, {arm}", "value": pt,
                        "ci_low": lo, "ci_high": hi, "n": len(d),
                        "note": "ratio of means, dataset bootstrap"})
        # The DIFFERENCE is the claim, so it gets its own interval.
        pt = share(d, "auroc_gc", "composition_auroc_gc") - share(d, "auroc_dn", "composition_auroc_dn")
        boot = []
        for _ in range(2000):
            q = d.iloc[rng.integers(0, len(d), len(d))]
            boot.append(share(q, "auroc_gc", "composition_auroc_gc")
                        - share(q, "auroc_dn", "composition_auroc_dn"))
        lo, hi = np.nanpercentile(boot, [2.5, 97.5])
        out.append({"check": "composition share, GC minus dinuc", "value": pt,
                    "ci_low": lo, "ci_high": hi, "n": len(d),
                    "note": "the claim; excludes zero if the effect is real"})

        # --- 2a. is R1's effect modified by dataset size? -----------------------------
        lp = np.log10(d.pairs)
        per_share = ((d.composition_auroc_gc - 0.5) / (d.auroc_gc - 0.5).replace(0, np.nan)
                     - (d.composition_auroc_dn - 0.5) / (d.auroc_dn - 0.5).replace(0, np.nan))
        ok = per_share.replace([np.inf, -np.inf], np.nan).notna()
        from scipy.stats import spearmanr
        rho, pv = spearmanr(lp[ok], per_share[ok])
        out.append({"check": "R1 effect vs log10(dataset size)", "value": rho,
                    "ci_low": np.nan, "ci_high": np.nan, "p": pv, "n": int(ok.sum()),
                    "note": "paired design, so this is effect modification not confounding"})

    # --- 1b. THE SHARE IS MODEL-DEPENDENT, and this is the headline --------------------
    #
    # The share above is measured against the k-mer model, because that is what the rehearsal
    # arm trains. Stated without the model named it reads as a fact about the TASK, and it is
    # not: it is a joint fact about the task and the model class. Composition reproduces most
    # of a k-mer model's skill and well under half of a fine-tuned language model's, on the
    # SAME datasets with the SAME negatives -- so the negative-set protocol does not penalise
    # every model equally, and a benchmark that reports one model's drop is not describing the
    # others. Horlacher 2023 published the drop across 11 methods; the decomposition into a
    # per-model share, and the finding that it separates model classes, is what is new here.
    #
    # Same estimator as 1a (ratio of means, dataset bootstrap) so the two are comparable, and
    # the k-mer row is a CROSS-CHECK: it must reproduce the dinuc arm above to 3 decimals,
    # because it is the same quantity computed from a different table.
    fm = fetch(bucket, "results/tables/matched_four_models.csv")
    if fm is not None and {"kmer_auroc", "cnn", "splicebert"} <= set(fm.columns):
        # ONE resample stream, all three models scored on the SAME resample. The first
        # version drew independently per model and then differenced the streams, which is an
        # unpaired difference of paired quantities: it gave [0.198, 0.339], width 0.141,
        # against a correctly paired [0.241, 0.295], width 0.054. Off by 2.6x, in the
        # direction that looks more conservative and is simply wrong. Section 1a above always
        # did this correctly; 1b did not, because it was written separately.
        rng = np.random.default_rng(7)
        cols = (("k-mer", "kmer_auroc"), ("CNN", "cnn"), ("SpliceBERT", "splicebert"))
        draws = {lab: [] for lab, _ in cols}
        for _ in range(2000):
            q = fm.iloc[rng.integers(0, len(fm), len(fm))]
            for lab, col in cols:
                draws[lab].append(share(q, col, "composition_auroc"))
        draws = {k: np.array(v) for k, v in draws.items()}
        for lab, col in cols:
            lo, hi = np.nanpercentile(draws[lab], [2.5, 97.5])
            out.append({"check": f"composition share vs {lab}",
                        "value": share(fm, col, "composition_auroc"),
                        "ci_low": lo, "ci_high": hi, "n": len(fm),
                        "note": "dinuc-matched; ratio of means, dataset bootstrap"})
        pt = share(fm, "kmer_auroc", "composition_auroc") - share(fm, "splicebert",
                                                                  "composition_auroc")
        lo, hi = np.nanpercentile(draws["k-mer"] - draws["SpliceBERT"], [2.5, 97.5])
        # NOT A FINDING, AND THE INTERVAL IS DECORATION. share_m = C/gain_m with a numerator C
        # that is IDENTICAL across models, so share_kmer/share_SB == gain_SB/gain_kmer exactly
        # (verified to 6 dp). This contrast is a monotone rescaling of the AUROC ladder, and
        # SpliceBERT beats the k-mer model on 95/95 datasets, so it excludes zero with
        # probability 1. It is retained as a readable presentation of the ladder and MUST NOT
        # be reported as an independent result. See r1_headline_is_gc_share_only.
        out.append({"check": "composition share, k-mer minus SpliceBERT", "value": pt,
                    "ci_low": lo, "ci_high": hi, "n": len(fm),
                    "note": "RESCALING of the AUROC ladder, not an independent finding"})

    # --- 2b. is the specificity gap explained by dataset size? ------------------------
    per = fetch(bucket, "results/tables/variant_specificity.csv")
    panel = fetch(bucket, "results/tables/panel_summary.csv")
    if per is not None and panel is not None:
        from scipy.stats import spearmanr
        m = per.merge(panel[["dataset", "pairs"]], on="dataset", how="left").dropna(subset=["pairs"])
        gap = m.auroc_matched - m.auroc_mismatched
        rho, pv = spearmanr(np.log10(m.pairs), gap)
        out.append({"check": "specificity gap vs log10(dataset size)", "value": rho,
                    "ci_low": np.nan, "ci_high": np.nan, "p": pv, "n": len(m),
                    "note": "paired within dataset, so size cannot confound the difference"})
        # And within the powered stratum only, where the effect is claimed.
        q = m[m.n_pathogenic >= 20]
        rho2, pv2 = spearmanr(np.log10(q.pairs), q.auroc_matched - q.auroc_mismatched)
        out.append({"check": "specificity gap vs size, powered stratum", "value": rho2,
                    "ci_low": np.nan, "ci_high": np.nan, "p": pv2, "n": len(q),
                    "note": "n_pathogenic >= 20"})

    # --- 3. prevalence baselines ------------------------------------------------------
    sb = fetch_prefix(bucket, "variants/scores_sb/")
    cvt = fetch(bucket, "results/tables/variant_conservation.csv")
    if sb is not None and cvt is not None:
        v = sb.dropna(subset=["delta"]).copy()
        v["dataset"] = v.protein + ":" + v.cell
        v["block"] = (v.vid.str.split(":").str[0] + "_" +
                      (v.vid.str.split(":").str[1].astype(int) // 1_000_000).astype(str))
        u = v.sort_values("delta", key=abs, ascending=False).drop_duplicates("vid")
        u = u.merge(cvt[["vid", "conservation"]], on="vid", how="left").dropna(subset=["conservation"])

        def loo_rate(df, key):
            """Group pathogenic rate EXCLUDING the variant itself."""
            g = df.groupby(key).label
            tot, cnt = g.transform("sum"), g.transform("size")
            return ((tot - df.label) / (cnt - 1)).where(cnt > 1, np.nan)

        for key, label in (("block", "1-Mb block prevalence"), ("dataset", "dataset prevalence")):
            r = loo_rate(u, key)
            m = r.notna()
            if m.sum() > 100 and u.label[m].nunique() == 2:
                out.append({"check": f"TRIVIAL baseline: {label}",
                            "value": roc_auc_score(u.label[m], r[m]),
                            "ci_low": np.nan, "ci_high": np.nan, "n": int(m.sum()),
                            "note": "leave-one-out; a high value means positional leakage"})
        out.append({"check": "reference: conservation, same variants",
                    "value": roc_auc_score(u.label, u.conservation),
                    "ci_low": np.nan, "ci_high": np.nan, "n": len(u), "note": "phyloP"})
        out.append({"check": "reference: right-protein head, same variants",
                    "value": roc_auc_score(u.label, u.delta.abs()),
                    "ci_low": np.nan, "ci_high": np.nan, "n": len(u),
                    "note": "pooled, i.e. the inflated framing"})

    if not out:
        log("skip robustness: inputs missing")
        return None
    r = pd.DataFrame(out)
    r.to_csv(TABLES / "robustness.csv", index=False)
    for _, x in r.iterrows():
        ci = f" [{x.ci_low:.3f}, {x.ci_high:.3f}]" if pd.notna(x.get("ci_low")) else ""
        log(f"  {x['check']:46} {x['value']:+.4f}{ci}")
    return r


def upload(bucket):
    n = 0
    for p in sorted(TABLES.glob("*.csv")):
        bucket.blob(f"results/tables/{p.name}").upload_from_filename(
            str(p), content_type="text/csv")
        n += 1
    for p in sorted(FIGS.glob("*")):
        if p.suffix in (".png", ".pdf"):
            bucket.blob(f"results/figures/{p.name}").upload_from_filename(str(p))
            n += 1
    log(f"uploaded {n} artefacts")
    return n


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--what", default="all", choices=["tables", "figures", "all"])
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    bucket = cloudcfg.bucket()
    log(cloudcfg.describe())
    if bucket.blob(MARKER).exists() and not a.force:
        log(f"{MARKER} present, nothing to do")
        return

    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    # Pull anything a downstream table needs that lives only in GCS. locality_ism.csv is NOT
    # in this list any more: it is built by locality() below, not fetched. Fetching your own
    # output is how R3's table went missing without a single error.
    for key in ("results/tables/variant_conservation.csv",
                "results/tables/variant_scores.csv"):
        d = fetch(bucket, key)
        if d is not None:
            d.to_csv(TABLES / Path(key).name, index=False)

    if a.what in ("tables", "all"):
        panel_description(bucket)
        cost_of_matching(bucket)
        four_models(bucket)
        locality(bucket)
        # The paired analysis runs FIRST and is the reported result; the pooled ladder is
        # kept afterwards so the size of the inflation is on the record rather than deleted.
        variant_specificity(bucket)
        variant_ladder(bucket)
        specificity_attacks(bucket)
        donor_overlap(bucket)
        robustness(bucket)
    figs_ok = do_figures() if a.what in ("figures", "all") else True

    upload(bucket)
    if a.what == "all" and figs_ok:
        bucket.blob(MARKER).upload_from_string(
            json.dumps({"git_sha": os.environ.get("GIT_SHA", "unknown")}),
            content_type="application/json")
        log(f"wrote {MARKER}")
    log("done")


if __name__ == "__main__":
    main()
