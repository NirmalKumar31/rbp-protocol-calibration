"""The ClinVar arm, rerun with SpliceBERT instead of the k-mer rehearsal model.

WHY THIS EXISTS. scripts/rehearsal_variants.py says it in its own docstring: "Every number
here is a rehearsal, because the scores come from a k-mer model rather than a trained
network." That rehearsal reported delta AUROC 0.558 against conservation's 0.910, and the
obvious reading -- a binding model adds nothing over conservation -- is not supported by it,
because the scoring model was the WEAK one. On held-out binding the same k-mer model scores
0.688 where SpliceBERT scores 0.809. Asking whether binding models help on variants, using
the model that barely predicts binding, answers a different question.

So: same variants, same conservation, same test, same collapse rule. Only the score changes.

    --what score   SpliceBERT delta per variant   -> variant_scores_splicebert.csv
    --what test    conservation control + pooled  -> variant_results_splicebert.csv

COVERAGE. 94 of the 95 datasets with SpliceBERT weights carry ClinVar variants: 32,967
variant-dataset pairs, 8,716 pathogenic. Smaller than the rehearsal's 187 datasets because
SpliceBERT was swept on 95, and that difference is a caveat to state, not to hide -- the
k-mer numbers are recomputed on the same 94 by `test` so the comparison is like for like.

WEIGHTS ARE STREAMED. 94 datasets x 5 folds x 75 MB is 35 GB and the disk has 26 GB free.
Each dataset's folds are fetched, used and deleted, so peak usage is ~400 MB. Per-dataset
results are cached under data/interim/vsb/ and skipped on rerun, because a 90 minute job
that has to restart from zero is a job that never finishes.
"""

import argparse
import os
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.utils import config as cfgmod  # noqa: E402
from rbp.variants import assign  # noqa: E402
from rbp.utils import cloud as cloudcfg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
INTERIM = ROOT / "data" / "interim" / "vsb"
BUCKET = cloudcfg.derived_bucket()
PROJECT = cloudcfg.project()
# Under variants/, NOT manifest/.
#
# `--what tables` runs on the driver VM, which runs as rbp-analysis, and rbp-analysis may
# write results/, variants/ and driver/. manifest/ belongs to rbp-prep. Left at manifest/ this
# would have thrown 403 storage.objects.create AFTER cutting and uploading 164,835 variant
# windows -- the same failure mode as the variants marker, one stage later.
#
# Nothing else references this object, and rbp-modal has read access across the whole derived
# bucket, so moving it costs nothing.
MANIFEST = "variants/variant_tasks.tsv"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# The four-model table is written under one name and was read under another.
#
# cloud_analysis.four_models() writes matched_four_models.csv. This script read
# matched95_four_models.csv -- the earlier study's filename, which encoded the panel size in
# the name and which nothing in this pipeline produces. So `--what tables` would have died
# with FileNotFoundError before cutting a single window, and the failure was invisible to
# every test because the tests that touch that file SKIP when it is absent.
#
# verify.py already tried both names. Doing the same here rather than renaming the producer
# keeps the one script that gates the science as the authority on what the file is called.
DEEP_PANEL_NAMES = ("matched_four_models.csv", "matched95_four_models.csv")


def deep_panel():
    """The datasets that have all four models, i.e. the panel the variant arm scores."""
    for name in DEEP_PANEL_NAMES:
        p = TABLES / name
        if p.exists():
            return set(pd.read_csv(p).dataset)
    raise SystemExit(
        f"none of {DEEP_PANEL_NAMES} is present in {TABLES}. It is written by "
        f"cloud_analysis.py's four_models(), which needs results/rehearsal_binding_dinuc.csv "
        f"and results/sweep_dinuc.csv -- run `cloud_analysis.py --what tables` first.")


def fetch_folds(client, cell, protein, dest):
    """Pull every fold checkpoint for one dataset. Returns {fold: path}."""
    dest.mkdir(parents=True, exist_ok=True)
    got = {}
    for f in range(5):
        b = client.bucket(BUCKET).blob(
            f"runs/dinuc/{cell}/{protein}/splicebert/fold{f}/best.pt")
        if not b.exists():
            continue
        p = dest / f"fold{f}.pt"
        b.download_to_filename(str(p))
        got[f] = p
    return got


def score_seqs(handle, seqs, device, batch=64):
    import torch
    handle.model.eval()
    out = []
    with torch.no_grad():
        for s in range(0, len(seqs), batch):
            logits = handle.forward(handle.batch(seqs[s:s + batch], device))
            out.append(logits.float().cpu().numpy().ravel())
    return np.concatenate(out) if out else np.zeros(0)


def delta_for_dataset(handle, ckpts, table, device):
    """score(ref) - score(alt), each variant scored only by a fold that never saw it.

    Same contract as baseline.variant_delta: a variant whose fold has no checkpoint comes
    back NaN rather than being quietly scored by the wrong model.
    """
    import torch
    folds = table.fold.to_numpy(dtype=float)
    ref, alt = table.seq_ref.tolist(), table.seq_alt.tolist()
    out = np.full(len(table), np.nan)
    for f in sorted(ckpts):
        sel = folds == f
        if not sel.any():
            continue
        handle.load(torch.load(ckpts[f], map_location="cpu", weights_only=False))
        handle.to(device)
        idx = np.flatnonzero(sel)
        r = score_seqs(handle, [ref[i] for i in idx], device)
        a = score_seqs(handle, [alt[i] for i in idx], device)
        out[idx] = r - a
    return out


def stage_score(cfg, limit):
    """The local path: cut windows and score in one process, no GCS round trip."""
    import torch
    from google.cloud import storage
    from pyfaidx import Fasta

    device = torch.device("cpu")
    INTERIM.mkdir(parents=True, exist_ok=True)
    tmp = ROOT / ".cache" / "vsb_ckpt"

    a = pd.read_csv(TABLES / "variant_assignments.csv")
    have = deep_panel()
    a["ds"] = a.protein + ":" + a.cell
    a = a[a.ds.isin(have)]

    fasta = Fasta(str(ROOT / "data/raw/GRCh38.primary_assembly.genome.fa"))
    size, shifts = cfg.windows["size"], cfg["variants"]["shifts"]
    how = cfg["variants"]["delta"]

    client = storage.Client(project=PROJECT)
    from rbp.models import registry     # torch-importing; see module docstring
    handle = registry.build("splicebert", cfg)          # built once, reloaded per fold

    groups = list(a.groupby(["protein", "cell"]))
    if limit:
        groups = groups[:limit]
    log(f"{len(groups)} datasets, {len(a):,} variant rows, device={device}")

    t0 = time.time()
    for i, ((prot, cell), g) in enumerate(groups, 1):
        cache = INTERIM / f"{cell}_{prot}.csv"
        if cache.exists():
            continue
        try:
            ckpts = fetch_folds(client, cell, prot, tmp)
            if not ckpts:
                log(f"  [{i}/{len(groups)}] {prot}:{cell} no checkpoints, skipped")
                continue
            table, _ = assign.build_scoring_table(g.to_dict("records"), fasta, size, shifts)
            if not table:
                continue
            t = pd.DataFrame(table)
            d = delta_for_dataset(handle, ckpts, t, device)
            vids, deltas = assign.collapse_delta(t.vid.to_numpy(), d, how=how)
            lab, fld = dict(zip(t.vid, t.label)), dict(zip(t.vid, t.fold))
            pd.DataFrame([{"protein": prot, "cell": cell, "vid": v, "label": lab[v],
                           "fold": fld[v], "delta": dv}
                          for v, dv in zip(vids, deltas)]).to_csv(cache, index=False)
            el = time.time() - t0
            log(f"  [{i}/{len(groups)}] {prot}:{cell} {len(vids):,} variants, "
                f"{el:.0f}s elapsed, ~{el / i * (len(groups) - i) / 60:.0f}m left")
        except Exception as e:                          # unattended: one bad dataset
            log(f"  [{i}/{len(groups)}] {prot}:{cell} FAILED {type(e).__name__}: {e}")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    parts = sorted(INTERIM.glob("*.csv"))
    df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    df.to_csv(TABLES / "variant_scores_splicebert.csv", index=False)
    log(f"wrote variant_scores_splicebert.csv  {len(df):,} rows from {len(parts)} datasets")


def stage_tables(cfg, limit):
    """Cut every variant's ref/alt windows locally, upload them, so the cloud needs no FASTA.

    THE POINT. Scoring needs two things: the 3.1 GB genome, to cut windows, and the
    checkpoints, which already live in GCS. Uploading the genome to run this in the cloud
    would move 3.1 GB to save a step that costs seconds. Cutting the windows here instead
    turns the genome into ~30 MB of sequence, and the cloud task never learns what a FASTA
    is.
    """
    # Imported here, not at module scope. The cloud task never opens a FASTA -- that is
    # the whole point of --what tables -- so the Modal image does not ship pyfaidx, and a
    # module-level import would crash every cloud container on a dependency it never uses.
    import gzip

    from pyfaidx import Fasta

    from google.cloud import storage

    a = pd.read_csv(TABLES / "variant_assignments.csv")
    have = deep_panel()
    a["ds"] = a.protein + ":" + a.cell
    a = a[a.ds.isin(have)]

    fasta = Fasta(str(ROOT / "data/raw/GRCh38.primary_assembly.genome.fa"))
    size, shifts = cfg.windows["size"], cfg["variants"]["shifts"]
    bucket = storage.Client(project=PROJECT).bucket(BUCKET)

    groups = list(a.groupby(["protein", "cell"]))
    if limit:
        groups = groups[:limit]
    rows, t0 = [], time.time()
    for i, ((prot, cell), g) in enumerate(groups, 1):
        table, dropped = assign.build_scoring_table(g.to_dict("records"), fasta, size, shifts)
        if not table:
            log(f"  {prot}:{cell} no usable windows, skipped ({dropped})")
            continue
        t = pd.DataFrame(table)[["vid", "label", "fold", "seq_ref", "seq_alt"]]
        key = f"variants/tables/{cell}_{prot}.tsv.gz"
        bucket.blob(key).upload_from_string(
            gzip.compress(t.to_csv(sep="\t", index=False).encode("utf-8")),
            content_type="application/gzip")
        rows.append({"idx": len(rows), "cell": cell, "protein": prot,
                     "n_windows": len(t), "n_variants": int(t.vid.nunique())})
        if i % 20 == 0 or i == len(groups):
            log(f"  [{i}/{len(groups)}] {time.time() - t0:.0f}s")

    man = pd.DataFrame(rows)
    bucket.blob(MANIFEST).upload_from_string(man.to_csv(sep="\t", index=False),
                                             content_type="text/tab-separated-values")
    log(f"{len(man)} datasets, {man.n_windows.sum():,} windows, "
        f"{man.n_variants.sum():,} variants -> gs://{BUCKET}/{MANIFEST}")


def stage_cloud(cfg, index, force, mismatch=0):
    """One dataset, from the uploaded table, weights straight out of GCS. No FASTA, no disk.

    Same completion-marker discipline as everywhere else in this project: the result object
    IS the marker, so a task killed midway redoes its work rather than being skipped.
    """
    import gzip
    import io

    import torch
    from google.cloud import storage

    bucket = storage.Client(project=PROJECT).bucket(BUCKET)
    man = pd.read_csv(io.StringIO(bucket.blob(MANIFEST).download_as_text()), sep="\t")
    if index >= len(man):
        sys.exit(f"index {index} beyond manifest of {len(man)}")
    r = man.iloc[index]
    cell, prot = r.cell, r.protein

    # THE MISMATCHED-HEAD CONTROL. Score THIS dataset's variants with a DIFFERENT protein's
    # fine-tuned weights. If a mismatched head scores nearly as well as the matched one, the
    # signal is not binding-specific -- it is SpliceBERT noticing that a substitution makes
    # the sequence less plausible, which every protein's head would inherit from the shared
    # pretrained body. That is the single cheapest way to tell "this model knows where THIS
    # protein binds" from "this model knows what human RNA looks like", and it needs no new
    # training because all 95 fine-tuned checkpoints are already on disk.
    #
    # The offset walks the pair-rank-sorted manifest, so the donor is a genuinely different
    # protein rather than the same protein in the other cell line.
    wcell, wprot = cell, prot
    if mismatch:
        w = man.iloc[(index + mismatch) % len(man)]
        wcell, wprot = w.cell, w.protein
        if wprot == prot:
            w = man.iloc[(index + mismatch + 1) % len(man)]
            wcell, wprot = w.cell, w.protein
    sub = "scores_mm" if mismatch else "scores_sb"
    out = f"variants/{sub}/{cell}_{prot}.csv"
    if bucket.blob(out).exists() and not force:
        log(f"{prot}:{cell} already done")
        return

    raw = bucket.blob(f"variants/tables/{cell}_{prot}.tsv.gz").download_as_bytes()
    t = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    from rbp.models import registry     # torch-importing; see module docstring
    handle = registry.build("splicebert", cfg)
    tmp = Path("/tmp/ckpt")
    ckpts = fetch_folds(storage.Client(project=PROJECT), wcell, wprot, tmp)
    if not ckpts:
        sys.exit(f"no checkpoints for {wprot}:{wcell}")

    t0 = time.time()
    d = delta_for_dataset(handle, ckpts, t, device)
    vids, deltas = assign.collapse_delta(t.vid.to_numpy(), d, how=cfg["variants"]["delta"])
    lab, fld = dict(zip(t.vid, t.label)), dict(zip(t.vid, t.fold))
    df = pd.DataFrame([{"protein": prot, "cell": cell, "vid": v, "label": lab[v],
                        "fold": fld[v], "delta": dv} for v, dv in zip(vids, deltas)])
    df["platform"] = os.environ.get("PLATFORM", "local")
    df["accelerator"] = str(device.type)
    df["weights_from"] = f"{wprot}:{wcell}"
    bucket.blob(out).upload_from_string(df.to_csv(index=False), content_type="text/csv")
    log(f"{prot}:{cell} weights={wprot}:{wcell} {len(df):,} variants on {device.type} "
        f"in {time.time() - t0:.0f}s")


def stage_gather(cfg):
    """Pull the cloud results down and merge them with anything scored locally."""
    import io

    from google.cloud import storage

    client = storage.Client(project=PROJECT)
    parts = []
    for b in client.list_blobs(BUCKET, prefix="variants/scores_sb/"):
        parts.append(pd.read_csv(io.BytesIO(b.download_as_bytes())))
    cloud = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    log(f"{len(cloud):,} rows from {len(parts)} cloud datasets")

    local_parts = [pd.read_csv(p) for p in sorted(INTERIM.glob("*.csv"))]
    local = pd.concat(local_parts, ignore_index=True) if local_parts else pd.DataFrame()
    if len(local):
        local["platform"] = local.get("platform", "local")
        local["accelerator"] = local.get("accelerator", "cpu")
        log(f"{len(local):,} rows from {len(local_parts)} local datasets")

    # Where both exist, the CPU and GPU runs of the same weights on the same windows should
    # agree to numerical noise. Reported rather than assumed, because "the two platforms
    # agree" is a claim, and a free one to check.
    if len(cloud) and len(local):
        j = cloud.merge(local, on=["protein", "cell", "vid"], suffixes=("_gpu", "_cpu"))
        if len(j):
            diff = (j.delta_gpu - j.delta_cpu).abs()
            log(f"CPU/GPU agreement on {len(j):,} shared variants: "
                f"max|diff| {diff.max():.2e}, median {diff.median():.2e}, "
                f"corr {j.delta_gpu.corr(j.delta_cpu):.6f}")

    df = cloud if len(cloud) else local
    if len(cloud) and len(local):
        keep = local[~local.set_index(["protein", "cell"]).index.isin(
            cloud.set_index(["protein", "cell"]).index)]
        df = pd.concat([cloud, keep], ignore_index=True)
    df.to_csv(TABLES / "variant_scores_splicebert.csv", index=False)
    log(f"wrote variant_scores_splicebert.csv  {len(df):,} rows, "
        f"{df.groupby(['protein', 'cell']).ngroups} datasets")


def stage_test(cfg, n_boot):
    """The same conservation control the rehearsal ran, on the same 94 datasets.

    The k-mer arm is recomputed here rather than quoted from variant_results_pooled.csv,
    because that one pooled 187 datasets and 65,940 variants. Comparing a 94-dataset
    SpliceBERT number against a 187-dataset k-mer number would confound the model with the
    panel.
    """
    from rbp.variants import conservation as cons

    c = pd.read_csv(TABLES / "variant_conservation.csv")[["vid", "conservation"]]
    sb = pd.read_csv(TABLES / "variant_scores_splicebert.csv")
    km = pd.read_csv(TABLES / "variant_scores.csv")
    km["ds"] = km.protein + ":" + km.cell
    sb["ds"] = sb.protein + ":" + sb.cell
    km = km[km.ds.isin(set(sb.ds))]

    rows = {}
    for tag, s in (("splicebert", sb), ("kmer", km)):
        df = s.merge(c, on="vid", how="left").dropna(subset=["delta", "conservation"])
        df["dataset"] = df.protein + ":" + df.cell
        log(f"{tag}: {len(df):,} pairs, {int(df.label.sum()):,} pathogenic, "
            f"{df.dataset.nunique()} datasets")
        pooled = cons.run(df, ["delta"], group_col=None, n_boot=n_boot,
                          method=cfg.conservation["method"])
        pooled.insert(0, "score_model", tag)
        rows[tag] = pooled
        per = cons.run(df, ["delta"], group_col="dataset", n_boot=n_boot,
                       method=cfg.conservation["method"])
        per.insert(0, "score_model", tag)
        per.to_csv(TABLES / f"variant_results_{tag}_per.csv", index=False)
        ok = per[per.note == ""]
        log(f"  {len(ok)} testable, {int(ok.controlled_survives.sum())} survive control, "
            f"{int(ok.get('controlled_survives_fdr', pd.Series(dtype=bool)).sum())} after FDR")

    out = pd.concat(rows.values(), ignore_index=True)
    out.to_csv(TABLES / "variant_results_splicebert.csv", index=False)
    show = ["score_model", "n", "n_pathogenic", "conservation_auroc", "delta_auroc",
            "corr_delta_conservation", "alone_coef", "controlled_coef",
            "controlled_ci_low", "controlled_ci_high", "attenuation",
            "controlled_survives"]
    print("\n" + out[[x for x in show if x in out]].T.to_string())
    log("wrote variant_results_splicebert.csv")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--what", required=True,
                   choices=["score", "tables", "cloud", "gather", "test"])
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--index", type=int, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("--mismatch", type=int, default=0,
                   help="score with another protein's weights, offset N in the manifest")
    p.add_argument("--n-boot", type=int, default=500)
    a = p.parse_args()
    cfg = cfgmod.load()
    if a.what == "score":
        stage_score(cfg, a.limit)
    elif a.what == "tables":
        stage_tables(cfg, a.limit)
    elif a.what == "cloud":
        idx = a.index if a.index is not None else int(os.environ.get("TASK_INDEX", 0))
        stage_cloud(cfg, idx, a.force, a.mismatch)
    elif a.what == "gather":
        stage_gather(cfg)
    else:
        stage_test(cfg, a.n_boot)


if __name__ == "__main__":
    main()
