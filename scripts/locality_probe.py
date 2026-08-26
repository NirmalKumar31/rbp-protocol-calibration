"""Run the ISM locality probe on trained models, and compare model classes.

    python scripts/locality_probe.py --n 12                 # k-mer vs SpliceBERT
    python scripts/locality_probe.py --n 12 --models kmer    # k-mer only, no downloads

THE QUESTION. The composition control says a k-mer model gains only +0.048 AUROC over 19
composition features, while SpliceBERT gains +0.158. Two explanations fit that:

  (a) SpliceBERT reads LOCAL sequence features that a bag of k-mers cannot represent
  (b) SpliceBERT is simply a better global-composition estimator

Those imply opposite things about the paper. (a) means the confound is a property of weak
models and real motif signal exists above it; (b) means the confound goes all the way up.
AUROC alone cannot separate them -- both predict the same ranking.

The ISM probe can: it asks whether a model's sensitivity is CONCENTRATED at a few positions
(local) or spread evenly (global). Run on the same windows, for both models, the difference
in Gini answers the question directly.

WHY THIS IS THE JUSTIFICATION FOR TRAINING NEURAL MODELS AT ALL. A bag of k-mers is local by
construction and cannot be non-local, so it could never have answered this. See
src/rbp/eval/locality_ism.py for the two faults that made the previous probe unusable.

PAIRING. Both models are probed on the SAME held-out windows of the SAME fold, so the
comparison is within-dataset and within-window. Nothing about dataset difficulty, window
composition or fold assignment can drive a difference between them.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from rbp.data.splits import split_of_fold  # noqa: E402
from rbp.eval import locality_ism as loc  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402
from rbp.utils import cloud as cloudcfg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BUCKET = cloudcfg.derived_bucket()
CACHE = ROOT / ".cache" / "weights"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def sweep_datasets(n):
    """n datasets spanning the size range, taken from the SpliceBERT panel.

    Systematic by pair rank, same principle as the manifest subsetting: a size-selected
    sample would be confounded with the very thing being measured, since delta_auroc
    correlates with log(pairs) at r = +0.418.
    """
    from google.cloud import storage
    c = storage.Client(project=cloudcfg.project())
    rows = [ln.split("\t") for ln in c.bucket(BUCKET)
            .blob("manifest/sweep_tasks_splicebert.tsv").download_as_text()
            .strip().splitlines()[1:]]
    d = pd.DataFrame(rows, columns=["idx", "cell", "protein", "model", "fold", "pairs"])
    d = d.astype({"pairs": int}).drop_duplicates(subset=["cell", "protein"])
    d = d.sort_values("pairs").reset_index(drop=True)
    step = max(1, len(d) // n)
    return d.iloc[::step].head(n)


def fold0_split(df, k):
    roles = df.fold.map(lambda f: split_of_fold(int(f), 0, k))
    return df[roles == "train"], df[roles == "test"]


def kmer_probe(df, windows, fold=0):
    """The fold-`fold` k-mer model, built by the SAME function that produced the AUROC table.

    THIS USED TO FIT ITS OWN MODEL AND THAT WAS WRONG. It called LogisticRegression with
    C=0.01 on the train split only, under a comment claiming that matched the composition
    arm. Two things were off: baseline.C_DEFAULT is 1.0, not 0.01, and fit_fold_models
    trains on every fold except the held-out one (four folds) while the train split excludes
    the validation fold too (three). So the "k-mer model" in the locality figure was not the
    "k-mer model" in the AUROC table, and the whole claim is a comparison BETWEEN those two
    figures.

    Measured before fixing: the C difference moved mean k-mer Gini by +0.003 against a
    +0.058 SpliceBERT advantage, so the conclusion never depended on it. Fixed anyway,
    because "the same model" is a claim the paper makes in words and must be true in code.
    """
    from rbp.eval import baseline

    models, vec = baseline.fit_fold_models(df.seq_rna.tolist(), df.label.to_numpy(),
                                           df.fold.to_numpy())
    m = models.get(fold)
    if m is None:
        return None
    return loc.locality(loc.kmer_score_fn(m, vec), windows, max_windows=len(windows))


def splicebert_probe(cell, protein, windows, cfg, device):
    import torch
    from google.cloud import storage

    from rbp.models import registry
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{cell}_{protein}_splicebert_fold0.pt"
    if not dest.exists():
        c = storage.Client(project=cloudcfg.project())
        b = c.bucket(BUCKET).blob(
            f"runs/dinuc/{cell}/{protein}/splicebert/fold0/best.pt")
        if not b.exists():
            return None
        b.download_to_filename(str(dest))
    handle = registry.build("splicebert", cfg)
    handle.load(torch.load(dest, map_location="cpu", weights_only=False))
    handle.to(device)
    return loc.locality(loc.torch_score_fn(handle, device), windows,
                        max_windows=len(windows), batch=64)


def cloud_one(cfg, index, force):
    """One dataset, everything from GCS. Same probe, same seed, on a GPU.

    WHY THIS IS WORTH MOVING. ISM is 3*L mutants per window: 101 positions, 20 windows,
    two models, 95 datasets is about 570,000 forward passes. On this laptop's CPU that is
    ~90 minutes; on a T4 it is a few minutes, and unlike the variant job it is genuinely
    compute-bound rather than waiting on a download.

    Like the variant task it needs no FASTA -- the windows are already columns in
    dataset.tsv -- so the only inputs are that file and one checkpoint.
    """
    import io

    import torch
    from google.cloud import storage

    from rbp.models import registry

    bucket = storage.Client(project=cloudcfg.project()).bucket(BUCKET)
    man = pd.read_csv(io.StringIO(
        bucket.blob("manifest/locality_tasks.tsv").download_as_text()), sep="\t")
    if index >= len(man):
        sys.exit(f"index {index} beyond manifest of {len(man)}")
    r = man.iloc[index]
    # cell_line OR cell. The study panel writes cell_line; an earlier hand-built manifest
    # wrote cell. Accepting both costs one line and avoids a column rename propagating into
    # every consumer of the panel -- which is the sort of change that quietly invalidates a
    # file other stages already read.
    cell = getattr(r, "cell_line", None) or r.cell
    prot = r.protein
    out = f"runs/locality/{cell}_{prot}.json"
    if bucket.blob(out).exists() and not force:
        log(f"{prot}:{cell} already done")
        return

    df = pd.read_csv(io.StringIO(bucket.blob(
        f"processed/dinuc/{cell}/{prot}/dataset.tsv").download_as_text()), sep="\t")
    k = cfg.cv["k"]
    tr, te = fold0_split(df, k)
    windows = te[te.label == 1].seq_rna.tolist()[:20]
    if len(windows) < 8:
        log(f"{prot}:{cell} only {len(windows)} held-out positives, skipped")
        return

    row = {"dataset": f"{prot}:{cell}", "protein": prot, "cell": cell,
           "pairs": int(r.pairs), "n_windows": len(windows)}
    t0 = time.time()
    res = kmer_probe(df, windows)
    if res:
        # gini_sd is the window-to-window spread, and without it a per-dataset difference
        # cannot be told from measurement scatter. Three datasets came back with SpliceBERT
        # LESS local than the k-mer model; two of those three turned out to sit inside
        # 1 SE of zero, and establishing that needed a separate local rerun purely because
        # this line dropped the number. Carry it.
        row["kmer_gini"], row["kmer_top10"] = res["gini"], res["top10_frac"]
        row["kmer_gini_sd"] = res["gini_sd"]

    dest = Path("/tmp/fold0.pt")
    b = bucket.blob(f"runs/dinuc/{cell}/{prot}/splicebert/fold0/best.pt")
    if b.exists():
        b.download_to_filename(str(dest))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        handle = registry.build("splicebert", cfg)
        handle.load(torch.load(dest, map_location="cpu", weights_only=False))
        handle.to(device)
        res = loc.locality(loc.torch_score_fn(handle, device), windows,
                           max_windows=len(windows), batch=256)
        if res:
            row["sb_gini"], row["sb_top10"] = res["gini"], res["top10_frac"]
            row["sb_gini_sd"] = res["gini_sd"]
        row["accelerator"] = device.type
    row["seconds"] = round(time.time() - t0, 1)
    bucket.blob(out).upload_from_string(json.dumps(row), content_type="application/json")
    log(f"{prot}:{cell} kmer {row.get('kmer_gini', float('nan')):.3f} "
        f"sb {row.get('sb_gini', float('nan')):.3f} ({row['seconds']}s)")


def cloud_gather(out_path):
    """Pull every per-dataset JSON down into the same CSV the local path writes."""
    import json as _json

    from google.cloud import storage
    c = storage.Client(project=cloudcfg.project())
    rows = [_json.loads(b.download_as_text())
            for b in c.list_blobs(BUCKET, prefix="runs/locality/")]
    d = pd.DataFrame(rows)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(out_path, index=False)
    log(f"{len(d)} datasets -> {out_path}")
    return d


def summarise(d, out_path):
    print(f"\n{'':22} {'median':>8} {'mean':>8} {'min':>8} {'max':>8}")
    for c in ("kmer_gini", "sb_gini", "kmer_top10", "sb_top10"):
        if c in d and d[c].notna().any():
            v = d[c].dropna()
            print(f"{c:22} {v.median():8.3f} {v.mean():8.3f} {v.min():8.3f} {v.max():8.3f}")
    if {"kmer_gini", "sb_gini"} <= set(d.columns):
        both = d.dropna(subset=["kmer_gini", "sb_gini"])
        if len(both) >= 3:
            from scipy.stats import wilcoxon
            diff = both.sb_gini - both.kmer_gini
            print(f"\nSpliceBERT - k-mer gini: median {diff.median():+.3f}  "
                  f"more local in {int((diff > 0).sum())}/{len(both)} datasets")
            if len(both) >= 6:
                print(f"  paired Wilcoxon p = {wilcoxon(both.sb_gini, both.kmer_gini)[1]:.3g}")
    print(f"\nwrote {out_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=12, help="datasets to probe")
    p.add_argument("--windows", type=int, default=30, help="held-out windows per dataset")
    p.add_argument("--models", default="kmer,splicebert")
    p.add_argument("--out", default="results/tables/locality_ism.csv")
    p.add_argument("--cloud", action="store_true", help="score one dataset from GCS")
    p.add_argument("--gather", action="store_true", help="collect the cloud results")
    p.add_argument("--index", type=int, default=None)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    cfg = cfgmod.load()
    if a.cloud:
        idx = a.index if a.index is not None else int(os.environ.get("TASK_INDEX", 0))
        return cloud_one(cfg, idx, a.force)
    if a.gather:
        return summarise(cloud_gather(a.out), a.out)
    k, C = cfg.cv["k"], 0.01
    want = [m.strip() for m in a.models.split(",")]
    device = None
    if "splicebert" in want:
        import torch
        device = torch.device("cpu")     # deterministic; MPS op coverage is uneven here

    picked = sweep_datasets(a.n)
    log(f"{len(picked)} datasets, pairs {picked.pairs.min():,}-{picked.pairs.max():,}, "
        f"{a.windows} windows each, models={want}")

    out = []
    for r in picked.itertuples():
        f = panelmod.data_dir(r.cell, "dinuc") / r.protein / "dataset.tsv"
        if not f.exists():
            log(f"  skip {r.protein}:{r.cell} (no local dataset)")
            continue
        df = pd.read_csv(f, sep="\t")
        tr, te = fold0_split(df, k)
        # Held-out BOUND windows only. Probing negatives would measure sensitivity on
        # sequence the model was trained to score low, which is a different question.
        windows = te[te.label == 1].seq_rna.tolist()[:a.windows]
        if len(windows) < 8:
            log(f"  skip {r.protein}:{r.cell} (only {len(windows)} held-out positives)")
            continue

        row = {"dataset": f"{r.protein}:{r.cell}", "protein": r.protein, "cell": r.cell,
               "pairs": int(r.pairs), "n_windows": len(windows)}
        t0 = time.time()
        if "kmer" in want:
            # The fold-0 model from baseline.fit_fold_models, i.e. the exact model behind
            # the AUROC table, not a re-fit with different hyperparameters.
            res = kmer_probe(df, windows)
            if res:
                row["kmer_gini"] = res["gini"]
                row["kmer_top10"] = res["top10_frac"]
        if "splicebert" in want:
            res = splicebert_probe(r.cell, r.protein, windows, cfg, device)
            if res:
                row["sb_gini"] = res["gini"]
                row["sb_top10"] = res["top10_frac"]
        row["seconds"] = round(time.time() - t0, 1)
        out.append(row)
        log(f"  {row['dataset']:16} {row['pairs']:6,} pairs  "
            f"kmer gini {row.get('kmer_gini', float('nan')):.3f}  "
            f"sb gini {row.get('sb_gini', float('nan')):.3f}  ({row['seconds']:.0f}s)")

    d = pd.DataFrame(out)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out, index=False)
    print(f"\n{'':22} {'median':>8} {'mean':>8} {'min':>8} {'max':>8}")
    for c in ("kmer_gini", "sb_gini", "kmer_top10", "sb_top10"):
        if c in d and d[c].notna().any():
            v = d[c].dropna()
            print(f"{c:22} {v.median():8.3f} {v.mean():8.3f} {v.min():8.3f} {v.max():8.3f}")
    if {"kmer_gini", "sb_gini"} <= set(d.columns):
        both = d.dropna(subset=["kmer_gini", "sb_gini"])
        if len(both) >= 3:
            from scipy.stats import wilcoxon
            diff = both.sb_gini - both.kmer_gini
            print(f"\nSpliceBERT - k-mer gini: median {diff.median():+.3f}  "
                  f"more local in {int((diff > 0).sum())}/{len(both)} datasets")
            if len(both) >= 6:
                print(f"  paired Wilcoxon p = {wilcoxon(both.sb_gini, both.kmer_gini)[1]:.3g}")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
