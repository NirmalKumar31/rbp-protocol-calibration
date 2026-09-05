"""Structural and configured constants the manuscript quotes, computed from their sources.

WHY THIS EXISTS. `audit_manuscript.py` now checks bare integers, and the first run left 18
orphans that were every one of them correct: the 101 nt window, the 256 4-mer columns, the
7089 CNN parameters, the 4000 bootstrap draws. None is an aggregate of a result column, so no
result table could ever source one, and they were the numbers "checked by hand" that the
Discussion used to confess to.

WHAT THIS DELIBERATELY IS NOT. It is not a transcription of the manuscript's integers into a
CSV so the audit passes. That would be a forgery of exactly the kind this repository has
already caught once: a table written from the paper cannot falsify the paper. Every row below
is DERIVED -- imported from the config, counted off the model object, or measured over the
committed evidence tree -- so a value that drifts from the source moves here first and the
manuscript then fails the audit against it. The `source` column names the derivation for each.
"""

import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval import baseline  # noqa: E402
from rbp.models.cnn import DeepBindCNN  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

OUT = ROOT / "results" / "tables" / "manuscript_constants.csv"
EVIDENCE = ROOT / "data" / "evidence"


def main():
    cfg = cfgmod.load()
    rows = []

    def add(name, value, source):
        rows.append({"name": name, "value": int(value), "source": source})

    size = int(cfg.windows["size"])
    add("window_size_nt", size, "config/params.yaml windows.size")
    add("cv_folds", cfg.cv["k"], "config/params.yaml cv.k")
    add("min_peak_distance_nt", cfg.negatives["min_peak_distance"],
        "config/params.yaml negatives.min_peak_distance")

    # the 4-mer block's shape follows from k and the window, so both are derived and neither
    # is a number anyone chose to write down
    k = baseline.fit_fold_models.__defaults__[0]
    add("kmer_k", k, "default of rbp.eval.baseline.fit_fold_models")
    add("kmer_columns", 4 ** k, f"4**k with k={k}")
    add("kmers_per_window", size - k + 1, f"window_size - k + 1 = {size} - {k} + 1")

    # LogisticRegression(max_iter=...) as the nested fit actually constructs it
    src = (ROOT / "src" / "rbp" / "stats.py").read_text()
    iters = sorted({int(t.split(")")[0]) for t in src.split("max_iter=")[1:]
                    if t.split(")")[0].isdigit()})
    for v in iters:
        add(f"logreg_max_iter_{v}", v, "literal in src/rbp/stats.py LogisticRegression(...)")

    # bootstrap draws, imported from the scripts that perform them rather than quoted
    boots = {}
    for name in ("cluster_intervals", "protocol_or_baseline", "recommendation_works",
                 "scale_sweep", "three_arm_models"):
        t = (ROOT / "scripts" / f"{name}.py").read_text()
        for line in t.splitlines():
            if line.startswith("N_BOOT"):
                boots[name] = int(line.split("=")[1].split("#")[0].strip())
    if len(set(boots.values())) == 1:
        add("bootstrap_draws", next(iter(boots.values())),
            f"N_BOOT, identical across {len(boots)} scripts: {', '.join(sorted(boots))}")
    else:
        for name, v in sorted(boots.items()):
            add(f"bootstrap_draws_{name}", v, f"N_BOOT in scripts/{name}.py")

    # parameter counts off the constructed objects, not off the paper
    cnn = DeepBindCNN()
    add("cnn_params", sum(p.numel() for p in cnn.parameters()),
        "sum of numel over DeepBindCNN().parameters()")

    # SpliceBERT's count is recorded by every run that trained it
    sb = sorted(EVIDENCE.glob("scores*/*/*/splicebert/fold*/metrics.json"))
    if sb:
        tot = {json.loads(p.read_text()).get("params_total") for p in sb[:200]}
        tot.discard(None)
        if len(tot) == 1:
            add("splicebert_params", tot.pop(),
                f"params_total, identical across {len(sb[:200])} committed metrics.json")

    # the committed evidence tree's own size
    runs = sorted(EVIDENCE.glob("scores*/*/*/*/fold*/scores.tsv.gz"))
    add("committed_fold_runs", len(runs),
        "count of scores.tsv.gz under data/evidence/scores*/")

    # The retrained datasets' size, which is why a count-stratified subset picks a different
    # 20. Read from the FROZEN list and not from fold_integrity.py, which by design reports
    # no flagged datasets once the retrain has landed.
    # THE STORE PATH IS A PARAMETER, and it was an absolute path under one author's home
    # directory. Guarded by exists(), so on any other machine it did not fail: it silently
    # skipped, and four manuscript constants disappeared from the table with no error and no
    # message. A missing input that says nothing is worse than one that stops.
    frozen = ROOT / "cloud" / "modal" / "retrain_dinuc_20.txt"
    store = Path(os.environ.get("RBP_STORE", ROOT.parent / "rbp-store"))
    man = store / "manifest" / "sweep_tasks_cnn_dinuc.tsv"
    if frozen.exists() and not man.exists():
        print(f"  note: {man} absent, so the retrain size constants are not derived. "
            f"Set RBP_STORE if the window store lives elsewhere.")
    if frozen.exists() and man.exists():
        leaky = {ln.strip() for ln in frozen.read_text().splitlines()
                 if ln.strip() and not ln.lstrip().startswith("#")}
        m = pd.read_csv(man, sep="\t")
        m["ds"] = m.protein + ":" + m.cell_line
        pr = m.groupby("ds").pairs.first()
        add("retrained_datasets", len(leaky), "rows of cloud/modal/retrain_dinuc_20.txt")
        add("retrained_median_pairs", pr[pr.index.isin(leaky)].median(),
            "median pairs over the retrained datasets, from the dinuc manifest")
        add("clean_median_pairs", pr[~pr.index.isin(leaky)].median(),
            "median pairs over the datasets that were not retrained")

    # The candidate pool the panel was selected from: one committed row per released ENCODE
    # eCLIP experiment in each cell line, which is what "139 and 105" counts.
    for f, cell in (("panel_full.tsv", "K562"), ("panel_full_HepG2.tsv", "HepG2")):
        q = ROOT / "config" / f
        if q.exists():
            d = pd.read_csv(q, sep="\t")
            add(f"encode_candidate_experiments_{cell}", d.experiment.nunique(),
                f"distinct experiments in config/{f}")

    # The eligible set the panel was halved out of: one committed row per dataset that cleared
    # the out-of-fold pair floor in the dinucleotide arm, across both cell lines.
    elig = sorted((ROOT / "config").glob("panel_final_*_dinuc.tsv"))
    if elig:
        n = sum(len(pd.read_csv(q, sep="\t")) for q in elig)
        add("eligible_datasets", n,
            "rows of " + " + ".join(f"config/{q.name}" for q in elig))

    # Per-arm pair totals and the gap between them. A column SUM is in the float haystack but
    # a difference of two sums is not, and the 1,264-pair gap between the dinucleotide arm and
    # the other two is a number the Methods has to account for.
    mq = ROOT / "results" / "tables" / "match_quality_per_dataset.csv"
    if mq.exists():
        t = pd.read_csv(mq).groupby("arm").pairs.sum()
        for a, v in t.items():
            add(f"pairs_{a}_arm", v, "sum of pairs over match_quality_per_dataset.csv")
        if {"dn", "gc"} <= set(t.index):
            add("pairs_dn_minus_gc", t["dn"] - t["gc"],
                "difference of those two sums, the arms' retained-positive gap")

    # What verify.py actually ran, recorded by verify.py rather than transcribed from its
    # console output, so a manuscript claiming a coverage it no longer has fails the audit.
    vs = ROOT / "results" / "tables" / "verify_summary.csv"
    if vs.exists():
        d = pd.read_csv(vs)
        for _, r in d.iterrows():
            add(f"verify_{r['name']}", r["value"], "written by scripts/verify.py")

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
