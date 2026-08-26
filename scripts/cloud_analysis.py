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


# --- R4 -------------------------------------------------------------------------------

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
        variant_ladder(bucket)
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
