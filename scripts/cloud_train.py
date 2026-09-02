"""The GPU sweep: one (cell, protein, model, fold) run per Batch task.

    python scripts/cloud_train.py manifest              # freeze the run list (local)
    python scripts/cloud_train.py run                   # one run, by BATCH_TASK_INDEX
    python scripts/cloud_train.py aggregate             # pool the folds -> the result

WHY THIS EXISTS SEPARATELY FROM scripts/train.py. train.py is the training path and knows
nothing about buckets, which is what makes it testable on a laptop. This file is the part
that makes an unattended run on reclaimable machines safe, and it does exactly three
things train.py cannot:

  1. THE COMPLETION MARKER LIVES IN GCS. train.py skips a run whose metrics.json is on the
     local disk. After a preemption the task restarts on a machine that has never seen that
     disk, so finished work would be redone. Here the marker is an object, and objects
     outlive the VM.

  2. THE CHECKPOINT LIVES IN GCS. Same reason, but the consequence is worse: without it a
     preemption at epoch 9 of 12 throws away nine epochs of GPU time. The checkpoint is
     mirrored after every epoch and pulled back before the run starts.

  3. OUT-OF-FOLD SCORES ARE KEPT PER ROW. Summary metrics cannot be pooled across folds and
     cannot be compared against the composition arm. Every task uploads its held-out rows
     with their ids; `aggregate` concatenates the five folds into one score per row and
     runs DeLong against the composition arm's scores for the same rows.

Shape and machinery are the same as cloud_prep.py and cloud_rehearsal.py, which between
them have survived 677 tasks and one mid-flight job deletion.
"""

import argparse
import gzip
import io
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402
from rbp.utils import cloud as cloudcfg  # noqa: E402

WORK = Path(os.environ.get("WORK_DIR", "/tmp/rbp"))
# ONE MANIFEST PER JOB, keyed by a tag.
#
# WHY, AND IT COST A WRONG DESIGN TO LEARN. Cloud Batch does NOT hand out task indices in
# global order. Measured on a live 1,890-task job: completed indices were spread from 2 to
# 1,582 with no sequential pattern -- Batch partitions the index space across nodes. So
# ORDERING a manifest cannot control what finishes first, and an earlier "breadth ordering"
# design was built on a false premise.
#
# The only reliable control is what is IN the manifest. Two models therefore mean two
# manifests and two jobs, not one manifest ordered cleverly.
MANIFEST_TAG = os.environ.get("MANIFEST_TAG", "")
MANIFEST = f"manifest/sweep_tasks{MANIFEST_TAG}.tsv"
CELLS = ("K562", "HepG2")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_tsv(text):
    lines = text.strip().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


def bucket(a):
    """GCS, or a directory that answers to the same interface.

    STORE_DIR (or --store) selects the local store. That is not a testing convenience: the
    GCP project's billing account was closed, so the bucket this sweep was built against
    returns 403 on every object. See rbp.utils.localstore.
    """
    from rbp.utils import localstore
    store = getattr(a, "store", None) or os.environ.get("STORE_DIR")
    if store:
        return localstore.LocalBucket(
            store, getattr(a, "store_ro", None) or os.environ.get("STORE_RO"))
    from google.cloud import storage
    return storage.Client(project=a.project).bucket(a.derived)


def _where(a):
    """How to name the store in a log line, so it never claims gs:// for a directory."""
    store = getattr(a, "store", None) or os.environ.get("STORE_DIR")
    return store if store else f"gs://{a.derived}"


def _not_found_errors():
    """What "that object is not there" looks like, for whichever store is in use.

    google-cloud-storage is not importable at all in a container that no longer needs it,
    so this cannot be a module-level import.
    """
    from rbp.utils.localstore import NotFound
    try:
        from google.api_core import exceptions as gexc
        return (NotFound, gexc.NotFound)
    except ImportError:
        return (NotFound,)


def read_scores(blob):
    """Read a scores table, gzipped or not, deciding by CONTENT rather than by name.

    cloud_rehearsal.py writes the composition arm to `<name>.scores.tsv.gz` with
    content_type="application/gzip" but WITHOUT compressing the bytes. The name and the
    content type both claim gzip; the payload is plain text. That went unnoticed for a day
    because nothing read those files until now, and it surfaced as
    `BadGzipFile: Not a gzipped file (b'id')` -- the 'id' being the TSV header.
    The writer is fixed, but 189 objects already exist in the old form.

    Sniffing the two-byte gzip magic number is the honest fix: it is correct for both
    forms, needs no migration, and cannot be wrong about a file it is looking at.
    """
    raw = blob.download_as_bytes()
    comp = "gzip" if raw[:2] == b"\x1f\x8b" else None
    return pd.read_csv(io.BytesIO(raw), sep="\t", compression=comp)


def run_prefix(arm, cell, protein, model, fold):
    """Every axis that distinguishes a run appears here. Miss one and two runs collide."""
    return f"runs/{arm}/{cell}/{protein}/{model}/fold{fold}"


# --- mode: manifest ---------------------------------------------------------------------

def do_manifest(a):
    """Freeze the run list from the panels the cloud itself produced.

    Ordering is by estimated cost, largest first. A Batch job finishes when its slowest
    node finishes, so the long tasks have to start in the first wave; alphabetical ordering
    would leave a 30,000-pair SpliceBERT run to begin near the end and add its whole
    duration to the makespan. Cost is approximated as pairs x sqrt(params), which is not
    exact but only has to get the ordering roughly right.
    """
    cfg = cfgmod.load()
    b = bucket(a)
    models = {n: s for n, s in cfg["models"].items()}
    if a.models:
        want = [m.strip() for m in a.models.split(",")]
        unknown = [m for m in want if m not in models]
        if unknown:
            sys.exit(f"unknown model(s) {unknown}; config defines {sorted(models)}")
        models = {n: models[n] for n in want}
    weight = {n: max(float(s.get("params_m", 0.5)), 0.5) ** 0.5 for n, s in models.items()}

    rows, picked = [], []
    for cell in CELLS:
        blob = b.blob(f"panel/{a.arm}/panel_final_{cell}_{a.arm}.tsv")
        if not blob.exists():
            sys.exit(f"no panel for {cell} {a.arm}; run cloud_prep.py finalize first")
        for r in read_tsv(blob.download_as_text()):
            pairs = int(r["pairs"])
            if pairs < a.min_pairs:
                continue
            picked.append((cell, r["protein"], pairs))

    # THE STUDY PANEL DECIDES MEMBERSHIP, if one exists. This used to be the --every flag
    # below, and that was the root cause of the project carrying four different dataset
    # counts: the panel was an emergent property of a flag typed during one sweep, recorded
    # nowhere. scripts/select_panel.py now writes it once as an artefact and every stage
    # reads it, so the cheap and expensive arms cannot end up describing different
    # populations. --every is kept only for the case where no panel has been defined.
    study = cloudcfg.study_panel(b)
    if study is not None:
        before = len(picked)
        picked = [t for t in picked if cloudcfg.in_study_panel(study, t[0], t[1])]
        log(f"study panel: {len(picked)} of {before} datasets")
        if a.every and a.every > 1:
            log("ignoring --every: the study panel already defines membership")
            a.every = None

    # SYSTEMATIC SAMPLE BY PAIR RANK, when a budget will not cover the whole panel.
    #
    # Cost scales with pairs, so a fixed budget buys a fraction of the panel. WHICH fraction
    # is a scientific decision, not a convenience. Taking the cheapest datasets would be the
    # obvious move and is wrong: delta_auroc correlates with log(pairs) at r = +0.418 on this
    # panel, so a size-selected subset is confounded with the very thing being measured.
    #
    # Sorting by pairs and keeping every Nth is unbiased in size BY CONSTRUCTION -- the
    # subset spans the full range with the same shape as the whole, which is what makes it
    # reportable as "a systematic half of the panel" rather than "the small ones we could
    # afford".
    if a.every and a.every > 1:
        picked.sort(key=lambda t: t[2])
        picked = picked[::a.every]

    for cell, protein, pairs in picked:
        for m in models:
            for f in range(cfg.cv["k"]):
                rows.append((cell, protein, m, f, pairs, pairs * weight[m]))
    if a.order == "cost":
        # Most expensive first. A job ends when its slowest node ends, so the long tasks
        # must start in the first wave. Right when the whole run will complete.
        rows.sort(key=lambda r: (-r[5], r[0], r[1], r[2], r[3]))
    else:
        # BREADTH. Cheapest model first, then cheapest datasets first, so that a run cut
        # short at an arbitrary time has covered the largest possible NUMBER of datasets.
        #
        # These two orderings are not interchangeable and the difference is enormous when a
        # model is slow. Measured on this panel: in 9.5 hours RNABERT gets through ~8
        # datasets ordered by cost, or ~107 ordered by breadth. Same compute, same money,
        # thirteen times the coverage, because cost is dominated by a handful of huge
        # datasets while the science cares about how many proteins are covered.
        #
        # It carries a real caveat: stopping partway leaves a SIZE-BIASED subset, and
        # delta_auroc is known to correlate with log(pairs) at r = +0.418 on this panel. So
        # a partial breadth-ordered arm must be reported with its pair range stated, and
        # ideally calibrated against a model that did finish the full panel.
        order = list(models)
        rows.sort(key=lambda r: (order.index(r[2]), r[4], r[0], r[1], r[3]))
    body = "idx\tcell_line\tprotein\tmodel\tfold\tpairs\n" + "".join(
        f"{i}\t{c}\t{p}\t{m}\t{f}\t{n}\n"
        for i, (c, p, m, f, n, _) in enumerate(rows))
    b.blob(MANIFEST).upload_from_string(body)

    log(f"{len(rows)} runs -> {_where(a)}/{MANIFEST}")
    log(f"  arm={a.arm}  models={list(models)}  k={cfg.cv['k']}  min_pairs={a.min_pairs}")
    n_ds = len({(c, p) for c, p, _, _, _, _ in rows})
    log(f"  {n_ds} datasets x {len(models)} models x {cfg.cv['k']} folds")
    for cell in CELLS:
        log(f"  {cell:6} {sum(1 for r in rows if r[0] == cell)}")


# --- mode: run ----------------------------------------------------------------------------

def task_index(n_tasks):
    """Which manifest row this container is responsible for.

    SHARD/NSHARDS exist because the sweep is five separate Batch jobs, one per region,
    since a subnetwork is regional and V100 quota is one per region. Striding rather than
    slicing means each region gets an interleaved mix of expensive and cheap runs; a
    contiguous slice would hand the whole first block of 30,000-pair SpliceBERT runs to one
    region and leave it running hours after the others had drained.
    """
    shard = int(os.environ.get("SHARD", 0))
    nshards = int(os.environ.get("NSHARDS", 1))
    local = int(os.environ.get("BATCH_TASK_INDEX", 0))
    return local * nshards + shard


def restore(b, prefix, outdir):
    """Pull a previous attempt's checkpoint back onto this machine, if there is one."""
    from rbp.train import trainer
    got = []
    for name in (trainer.CHECKPOINT, trainer.BEST):
        blob = b.blob(f"ckpt/{prefix}/{name}")
        if blob.exists():
            blob.download_to_filename(str(outdir / name))
            got.append(name)
    return got


def syncer(b, prefix):
    from rbp.train import trainer

    def sync(epoch, outdir):
        # best.pt goes up FIRST, for the same reason it is written to disk first: a
        # checkpoint that names a best epoch whose weights were never uploaded is worse
        # than a best.pt that is ahead of its checkpoint. The latter costs one repeated
        # epoch; the former restores the wrong model.
        for name in (trainer.BEST, trainer.CHECKPOINT):
            p = outdir / name
            if p.exists():
                b.blob(f"ckpt/{prefix}/{name}").upload_from_filename(str(p))
    return sync


def do_run(a):
    import torch

    from rbp.models import registry
    from rbp.train import data as tdata
    from rbp.train import trainer

    cfg = cfgmod.load()
    b = bucket(a)
    tasks = read_tsv(b.blob(MANIFEST).download_as_text())
    idx = a.index if a.index is not None else task_index(len(tasks))
    if idx >= len(tasks):
        log(f"index {idx} beyond manifest of {len(tasks)}: nothing to do")
        return
    t = tasks[idx]
    cell, protein, model, fold = t["cell_line"], t["protein"], t["model"], int(t["fold"])
    prefix = run_prefix(a.arm, cell, protein, model, fold)
    marker = f"{prefix}/metrics.json"
    log(f"task {idx}: {protein} {cell} {model} fold{fold}  ({int(t['pairs']):,} pairs)")

    if b.blob(marker).exists() and not a.force:
        log("already complete, nothing to do")
        return
    if a.dry_run:
        log("dry run: manifest resolved, stopping before any work")
        return

    src = f"processed/{a.arm}/{cell}/{protein}/dataset.tsv"
    local = WORK / src
    local.parent.mkdir(parents=True, exist_ok=True)
    if not local.exists():
        b.blob(src).download_to_filename(str(local))

    outdir = WORK / prefix
    outdir.mkdir(parents=True, exist_ok=True)
    if not a.no_resume:
        got = restore(b, prefix, outdir)
        if got:
            log(f"  restored {', '.join(got)} from a previous attempt")

    tcfg = cfg["train"]
    device = registry.device_of(a.device)
    # REFUSE TO RUN ON CPU ON A GPU NODE. registry.device_of falls back to CPU when cuda is
    # unavailable, which is right on a laptop and disastrous here: the node bills at V100
    # rates either way, the run is perhaps a hundred times slower, and nothing in the logs
    # says so. If the driver install failed or the container cannot see the device, that is
    # a job to fix, not a job to run slowly. --device cpu is the deliberate override.
    if a.device is None and device.type != "cuda":
        sys.exit(f"no CUDA device visible (torch sees {device}); refusing to run a GPU "
                 f"task on CPU. Pass --device cpu if that is genuinely what you want.")
    log(f"  device {device}  {torch.cuda.get_device_name(0) if device.type == 'cuda' else ''}")
    # SEED BEFORE BUILD, and this was wrong for the entire study.
    #
    # trainer.train() calls set_seed, but it does so AFTER this function has already
    # constructed the network -- so every weight was drawn from an unseeded RNG and
    # torch.manual_seed(7) only ever governed dropout and batch order. Three fresh processes
    # produced three different first-layer sums. All 945 deep-model fold-runs in this project
    # were therefore trained from an uncontrolled, unrecorded initialisation, and the paper's
    # "identical seed" was false.
    #
    # Measured consequence, from a referee's replication: per-dataset SD of the nested
    # contribution across independent training runs is 0.006 (CNN) to 0.010 (SpliceBERT).
    # That is immaterial for a panel mean over 94 datasets, where it induces ~0.001 against a
    # reported CI half-width of ~0.008, but it is fatal for exact reproducibility and for any
    # per-dataset count.
    trainer.set_seed(cfg.seed)
    handle = registry.build(model, cfg)
    dl = tdata.loaders(local, tcfg["batch_size"], seed=cfg.seed, fold=fold, k=cfg.cv["k"])
    ids = tdata.test_ids(local, fold=fold, k=cfg.cv["k"])

    metrics = trainer.train(
        handle, dl, outdir, epochs=a.epochs or tcfg["epochs"],
        patience=tcfg["patience"], seed=cfg.seed, device=device,
        weight_decay=float(tcfg["weight_decay"]), resume=not a.no_resume,
        on_epoch=syncer(b, prefix), test_ids=ids, log=log)

    metrics.update({"dataset": f"{protein}:{cell}", "protein": protein, "cell": cell,
                    "arm": a.arm, "fold": fold, "pairs": int(t["pairs"]),
                    "git_sha": os.environ.get("GIT_SHA", "unknown"),
                    "image_digest": os.environ.get("IMAGE_DIGEST", "unknown"),
                    # WHERE this ran. The three large models cannot run on GCP -- the
                    # project's GPU quota is zero and unraisable -- so they execute on Modal
                    # while cnn and rnabert ran on Cloud Batch CPU. Different hardware means
                    # a reader is entitled to know which rows came from where, and a results
                    # table that cannot answer that is hiding a real methodological detail.
                    "platform": os.environ.get("PLATFORM", "gcp-batch"),
                    # The real device, not a two-way guess. This said "cpu" for anything
                    # that was not CUDA, so a run on Apple MPS -- which is how the CNN arm
                    # is trained now that there is no GCP -- recorded hardware it never
                    # touched, in the one field a reader would use to check exactly that.
                    "accelerator": (torch.cuda.get_device_name(0)
                                    if device.type == "cuda" else device.type)})

    # BENCH MODE UPLOADS NOTHING, and this exists because it already went wrong once.
    #
    # A GPU timing probe was run with --epochs 2 --force to compare three accelerators on
    # one dataset. Each overwrote the last, and the survivor landed in the results table as
    # a legitimate-looking row: BCLAF1 fold1, 2 epochs, test AUROC 0.9152, indistinguishable
    # from a real result except by its epoch count. A deliberately under-trained model
    # sitting in the results is exactly the sort of thing that survives every downstream
    # check and reaches a paper.
    #
    # Timing needs the clock, not the artefact. So --bench trains, prints, and returns.
    if a.bench:
        log(f"  BENCH: {metrics['epochs_run']} epochs in {metrics['seconds']:.0f}s "
            f"({metrics['seconds']/max(metrics['epochs_run'],1):.1f} s/epoch) -- nothing uploaded")
        return

    npz = np.load(outdir / "test_predictions.npz", allow_pickle=True)
    sc = pd.DataFrame({"id": npz["id"], "label": npz["label"].astype(int),
                       "fold": fold, "score": npz["prob"]})
    b.blob(f"{prefix}/scores.tsv.gz").upload_from_string(
        gzip.compress(sc.to_csv(sep="\t", index=False).encode()),
        content_type="application/gzip")

    # PERSIST THE TRAINED WEIGHTS. Not optional, and it was missing.
    #
    # trainer.train writes best.pt to the VM's disk and cloud_train deletes the mirrored
    # copy under ckpt/ once the run finishes, so the first 945 CNN runs produced scores and
    # metrics and NO MODEL. That is fine for the AUROC comparison and fatal for the locality
    # probe, which has to feed perturbed sequences through the trained network -- the whole
    # scientific justification for training a neural model here rather than stopping at the
    # k-mer baseline is that a CNN *could* be non-local and a bag of k-mers structurally
    # cannot. You cannot probe a model you threw away.
    #
    # Cost of keeping them: the CNN is 7,089 parameters, about 28 KB, so 945 folds is 26 MB.
    # RNABERT is 496,697 trainable, roughly 2 MB, so 1.9 GB. Both are cents per month. A
    # 19.7M-parameter full fine-tune would be ~75 GB and would deserve a second thought.
    best_local = outdir / trainer.BEST
    if best_local.exists():
        b.blob(f"{prefix}/{trainer.BEST}").upload_from_filename(str(best_local))

    # Written LAST, so a task killed between the two uploads redoes the work rather than
    # being skipped with no scores. Same rule as cloud_prep and cloud_rehearsal.
    b.blob(marker).upload_from_string(json.dumps(metrics),
                                      content_type="application/json")
    # The run is finished, so the checkpoint is dead weight. Left behind, 4,725 of them
    # would sit in the bucket at up to 240 MB each. Missing is fine and expected: a run
    # that never resumed and finished inside one attempt still uploaded them, but a
    # forced re-run may not have.
    missing = _not_found_errors()
    for name in (trainer.CHECKPOINT, trainer.BEST):
        try:
            b.blob(f"ckpt/{prefix}/{name}").delete()
        except missing:
            pass
    log(f"  val {metrics['val_auroc']:.4f}  test {metrics['test_auroc']:.4f}  "
        f"in {metrics['seconds']:.0f}s over {metrics['epochs_run']} epochs")


# --- mode: aggregate ------------------------------------------------------------------

def do_aggregate(a):
    from rbp.eval import delong

    b = bucket(a)
    client = b.client
    rows = [json.loads(x.download_as_text())
            for x in client.list_blobs(a.derived, prefix=f"runs/{a.arm}/")
            if x.name.endswith("metrics.json")]
    if not rows:
        sys.exit(f"no results under runs/{a.arm}/")
    per_fold = pd.DataFrame(rows)

    # Pool the folds. This is the number params.yaml calls primary: one score per row,
    # every row scored by a model that never trained on it, one AUROC over all of them.
    # It is NOT the mean of the five per-fold AUROCs, which weights a small fold as
    # heavily as a large one.
    out = []
    for (dataset, model), g in per_fold.groupby(["dataset", "model"]):
        protein, cell = dataset.split(":")
        folds = sorted(g.fold)
        parts = []
        for f in folds:
            blob = b.blob(f"{run_prefix(a.arm, cell, protein, model, f)}/scores.tsv.gz")
            if not blob.exists():
                continue
            parts.append(read_scores(blob))
        if not parts:
            continue
        s = pd.concat(parts, ignore_index=True)
        if s.id.duplicated().any():
            raise ValueError(f"{dataset} {model}: a row is scored twice; the folds do "
                             "not partition the data")
        auroc, lo, hi = delong.auc_ci(s.score.to_numpy(), s.label.to_numpy())
        rec = {"dataset": dataset, "protein": protein, "cell": cell, "model": model,
               "arm": a.arm, "n_folds": len(folds), "n": len(s),
               "pooled_auroc": auroc, "ci_low": lo, "ci_high": hi,
               "mean_fold_auroc": float(g.test_auroc.mean()),
               "seconds": float(g.seconds.sum()),
               "git_sha": g.git_sha.iloc[0]}

        # Against the composition arm, on exactly the same rows. This is the comparison
        # the whole study is built on, and it is only possible because both sides carry
        # row ids.
        comp = b.blob(f"rehearsal/{a.arm}/{cell}/{protein}.scores.tsv.gz")
        if comp.exists():
            c = read_scores(comp)
            j = s.merge(c[["id", "score"]], on="id", suffixes=("", "_kmer"))
            if len(j) == len(s):
                d = delong.delong_test(j.score.to_numpy(), j.score_kmer.to_numpy(),
                                       j.label.to_numpy())
                rec.update({"kmer_auroc": d["auc_b"], "delta_vs_kmer": d["diff"],
                            "delong_z": d["z"], "delong_p": d["p"]})
        out.append(rec)

    res = pd.DataFrame(out).sort_values(["dataset", "model"]).reset_index(drop=True)
    body = res.to_csv(index=False)
    b.blob(f"results/sweep_{a.arm}.csv").upload_from_string(body)
    Path("results/tables").mkdir(parents=True, exist_ok=True)
    Path(f"results/tables/sweep_{a.arm}.csv").write_text(body)

    print(f"\n{len(per_fold)} runs -> {len(res)} (dataset, model) results   [arm={a.arm}]")
    print(f"{'model':12} {'n':>4} {'median':>8} {'mean':>8} {'vs kmer':>9} {'gpu-h':>7}")
    for model, g in res.groupby("model"):
        beats = (g.delong_p < 0.05).sum() if "delong_p" in g else 0
        print(f"{model:12} {len(g):4d} {g.pooled_auroc.median():8.4f} "
              f"{g.pooled_auroc.mean():8.4f} {beats:9d} {g.seconds.sum()/3600:7.2f}")
    incomplete = res[res.n_folds < 5]
    if len(incomplete):
        print(f"\n{len(incomplete)} incomplete (fewer than 5 folds):")
        print(incomplete[["dataset", "model", "n_folds"]].to_string(index=False))
    print(f"\nwrote {_where(a)}/results/sweep_{a.arm}.csv and results/tables/")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["manifest", "run", "aggregate"])
    p.add_argument("--derived", default=os.environ.get("DERIVED_BUCKET"))
    p.add_argument("--project", default=os.environ.get("GOOGLE_CLOUD_PROJECT"))
    p.add_argument("--store", default=os.environ.get("STORE_DIR"),
                   help="use this directory instead of GCS. See rbp.utils.localstore.")
    p.add_argument("--store-ro", default=os.environ.get("STORE_RO"),
                   help="read-only fallback root for --store, e.g. a mounted Volume.")
    p.add_argument("--arm", default=os.environ.get("ARM"),
                   choices=sorted(panelmod.ARMS))
    p.add_argument("--min-pairs", type=int, default=None)
    p.add_argument("--order", default="cost", choices=["cost", "breadth"],
                   help="cost: expensive first, minimises total wall time. breadth: cheap "
                        "first. NOTE: Batch does not run tasks in manifest order, so this "
                        "affects makespan only, never which tasks complete first. To "
                        "control that, restrict --models and use --tag.")
    p.add_argument("--every", type=int, default=None,
                   help="keep every Nth dataset after sorting by pairs. A systematic, "
                        "size-unbiased sample for when a budget cannot cover the panel.")
    p.add_argument("--tag", default=None,
                   help="suffix for the manifest object, so one job per model set can have "
                        "its own frozen list. Workers read it via MANIFEST_TAG.")
    p.add_argument("--models", default=os.environ.get("MODELS"),
                   help="comma-separated subset for `manifest`. The manifest IS the run "
                        "list, so restricting it here is how a partial sweep is expressed. "
                        "Completion markers are keyed by run, not by manifest, so a later "
                        "fuller manifest skips whatever is already done.")
    p.add_argument("--epochs", type=int, default=None)
    p.add_argument("--device", default=None)
    p.add_argument("--index", type=int, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("--bench", action="store_true",
                   help="train and time it, upload nothing. For accelerator comparisons.")
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if not a.derived and not a.store:
        sys.exit("--derived or DERIVED_BUCKET required (or --store for a local store)")
    cfg = cfgmod.load()
    a.arm = panelmod.arm_of(cfg, a.arm)
    if a.min_pairs is None:
        a.min_pairs = cfg.cv["min_pairs"]
    if a.tag is not None:
        global MANIFEST
        MANIFEST = f"manifest/sweep_tasks{a.tag}.tsv"
    {"manifest": do_manifest, "run": do_run, "aggregate": do_aggregate}[a.mode](a)


if __name__ == "__main__":
    main()
