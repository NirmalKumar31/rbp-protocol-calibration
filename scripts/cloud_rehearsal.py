"""Stage 5 in the cloud: the composition control, one dataset per Batch task.

WHAT THIS MEASURES, AND WHY IT IS THE POINT OF THE STUDY. A published RBP binding model
reports an AUROC. Some of that number is the protein's actual sequence preference and some
of it is nothing but nucleotide composition -- binding sites are, on average, compositionally
unlike the rest of the transcriptome, and a model can score well by noticing only that. This
step fits both arms out-of-fold on identical folds and reports the difference:

    composition_auroc   what 19 composition features alone achieve
    auroc               what a k-mer sequence model achieves
    delta_auroc         what the sequence model adds OVER composition

The median delta across the panel is the headline. It needs no GPU, which is why it runs
here rather than waiting on the sweep.

WHY ONE DATASET PER TASK. Every dataset is independent, and the per-dataset cost is
dominated by a cluster bootstrap that is itself serial. Same shape as cloud_prep.py, and it
reuses the same manifest, staging, and completion-marker machinery -- which has now survived
488 tasks and one mid-flight job deletion.

    python scripts/cloud_rehearsal.py manifest    # freeze the dataset list (run locally)
    python scripts/cloud_rehearsal.py run         # one dataset, by BATCH_TASK_INDEX
    python scripts/cloud_rehearsal.py aggregate   # collect rows -> the panel result
"""

import argparse
import gzip
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.eval import baseline, nested  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402
from rbp.utils import cloud as cloudcfg  # noqa: E402

WORK = Path(os.environ.get("WORK_DIR", "/tmp/rbp"))
MANIFEST = "manifest/rehearsal_tasks.tsv"
CELLS = ("K562", "HepG2")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_tsv(text):
    lines = text.strip().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


def buckets(a):
    from google.cloud import storage
    c = storage.Client(project=a.project)
    return c.bucket(a.derived)


# --- mode: manifest ---------------------------------------------------------------------

def do_manifest(a):
    """Freeze the dataset list from the panel the cloud itself produced.

    Read from GCS rather than from config/ so the list is whatever `finalize` decided,
    not whatever happens to be on the laptop.
    """
    derived = buckets(a)
    # THE STUDY PANEL, if one has been defined, decides membership. Without this filter the
    # rehearsal would silently run on every dataset the arm produced while the expensive
    # arms ran on the sampled subset, and the four-model comparison would then be built
    # from two different populations. See docs/PANELS.md and scripts/select_panel.py.
    study = cloudcfg.study_panel(derived)
    if study is not None:
        log(f"filtering to the study panel: {len(study)} datasets")
    # EVERY ARM IN ONE MANIFEST, with the arm carried on the row.
    #
    # This used to be one manifest per arm, which forced one Batch job per arm and therefore
    # a LOCAL process to wait for the first and submit the second. That makes the laptop part
    # of the pipeline: close the lid between arms and the second never starts. cloud_prep.py
    # already solved this by putting `arm` in the row and letting each task read its own, so
    # both arms are one job. Same pattern here.
    #
    # The point is not tidiness. A pipeline whose sequencing lives in a shell loop on someone's
    # machine is not a cloud pipeline, it is a cloud-assisted manual process.
    # PINNED to the two composition-matched arms rather than sorted(ARMS). This stage is the
    # k-mer/composition rehearsal and its whole design is the paired two-arm contrast; when
    # neg2 joined ARMS for the neural sweep, iterating ARMS here would have silently added a
    # third arm to a job that has no paired interpretation for it. Pass --arm neg2 --arm-only
    # to run it deliberately.
    arms = [a.arm] if a.arm_only else sorted(panelmod.COMPOSITION_MATCHED_ARMS)
    rows = []
    skipped = 0
    for arm in arms:
        for cell in CELLS:
            blob = derived.blob(f"panel/{arm}/panel_final_{cell}_{arm}.tsv")
            if not blob.exists():
                sys.exit(f"no panel for {cell} {arm}; run cloud_prep.py finalize first")
            for r in read_tsv(blob.download_as_text()):
                pairs = int(r["pairs"])
                if pairs < a.min_pairs:
                    continue
                if not cloudcfg.in_study_panel(study, cell, r["protein"]):
                    skipped += 1
                    continue
                rows.append((cell, r["protein"], arm, pairs))
    if skipped:
        log(f"  {skipped} dataset-arms outside the study panel, skipped")
    # Biggest first: the bootstrap is linear in rows, so the longest tasks start first or the
    # job ends when the unluckiest node finishes. Same reasoning as prep.
    rows.sort(key=lambda r: (-r[3], r[0], r[1], r[2]))
    body = "idx\tcell_line\tprotein\tarm\tpairs\n" + "".join(
        f"{i}\t{c}\t{p}\t{m}\t{n}\n" for i, (c, p, m, n) in enumerate(rows))
    derived.blob(f"{MANIFEST}").upload_from_string(body)
    log(f"{len(rows)} dataset-arms -> gs://{a.derived}/{MANIFEST}")
    log(f"  arms={','.join(arms)}  min_pairs={a.min_pairs}  "
        f"pairs {min(r[3] for r in rows):,} to {max(r[3] for r in rows):,}")
    for arm in arms:
        log(f"  {arm:6} {sum(1 for r in rows if r[2] == arm)}")


# --- mode: run --------------------------------------------------------------------------

def do_run(a):
    derived = buckets(a)
    tasks = read_tsv(derived.blob(MANIFEST).download_as_text())
    idx = a.index if a.index is not None else int(os.environ.get("BATCH_TASK_INDEX", 0))
    if os.environ.get("TASK_LIST"):
        idx = [int(x) for x in os.environ["TASK_LIST"].split(",")][idx]
    if idx >= len(tasks):
        sys.exit(f"index {idx} beyond manifest of {len(tasks)}")
    t = tasks[idx]
    cell, name = t["cell_line"], t["protein"]
    # The arm comes from the ROW. A task must not depend on how the job was invoked, or the
    # same manifest index means different work depending on an environment variable.
    arm = t.get("arm", a.arm)
    out = f"rehearsal/{arm}/{cell}/{name}.json"
    log(f"task {idx}: {name} {cell} {arm}  ({int(t['pairs']):,} pairs)")

    if derived.blob(out).exists() and not a.force:
        log("already present, nothing to do")
        return
    if a.dry_run:
        log("dry run: manifest resolved, stopping before any work")
        return

    src = f"processed/{arm}/{cell}/{name}/dataset.tsv"
    local = WORK / src
    local.parent.mkdir(parents=True, exist_ok=True)
    derived.blob(src).download_to_filename(str(local))
    df = pd.read_csv(local, sep="\t")

    t0 = time.time()
    y, folds, seqs = df.label.to_numpy(), df.fold.to_numpy(), df.seq_rna.tolist()
    res = baseline.evaluate(df, k=a.k)
    g = nested.gain_over_composition(seqs, res["scores"], y, folds)
    nt = nested.test_score(seqs, res["scores"], y, n_boot=a.n_boot, seed=7)

    row = {
        "dataset": f"{name}:{cell}", "protein": name, "cell": cell, "arm": arm,
        "pairs": int(len(df) // 2), "n": int(res["n"]),
        "auroc": res["auroc"], "ci_low": res["ci_low"], "ci_high": res["ci_high"],
        "composition_auroc": g.auroc_composition,
        "with_score_auroc": g.auroc_with_score,
        "delta_auroc": g.delta, "delta_ci_low": g.delta_ci_low,
        "delta_ci_high": g.delta_ci_high, "delta_p": g.delta_p, "helps": bool(g.helps),
        "coef": nt.coef, "coef_ci_low": nt.ci_low, "coef_ci_high": nt.ci_high,
        "lr_p": nt.lr_p,
        "k": a.k, "n_boot": a.n_boot, "seconds": round(time.time() - t0, 1),
        "git_sha": os.environ.get("GIT_SHA", "unknown"),
    }
    # Out-of-fold scores are written alongside the summary. Without them the DeLong
    # comparisons and any re-analysis would need the whole run repeating.
    sc = pd.DataFrame({"id": df.id, "label": y, "fold": folds, "score": res["scores"]})
    # ACTUALLY COMPRESS IT. This said content_type="application/gzip" and uploaded plain
    # bytes, so 189 objects are named .gz, declare themselves gzip, and are not. Nothing
    # read them for a day, then cloud_train.py aggregate died on
    # `BadGzipFile: Not a gzipped file (b'id')`. A filename is not a format; if you claim
    # a content type, produce it. The reader sniffs the magic number so the old objects
    # still work -- see cloud_train.read_scores.
    derived.blob(f"rehearsal/{arm}/{cell}/{name}.scores.tsv.gz").upload_from_string(
        gzip.compress(sc.to_csv(sep="\t", index=False).encode("utf-8")),
        content_type="application/gzip")

    # The summary is written LAST and is the completion marker, for the same reason as in
    # cloud_prep.py: a task preempted between the two uploads must redo the work, not skip.
    derived.blob(out).upload_from_string(json.dumps(row), content_type="application/json")
    log(f"  auroc {row['auroc']:.3f}  composition {row['composition_auroc']:.3f}  "
        f"gain {row['delta_auroc']:+.4f}  in {row['seconds']}s")


# --- mode: aggregate ----------------------------------------------------------------------

def do_aggregate(a):
    derived = buckets(a)
    rows = [json.loads(b.download_as_text())
            for b in derived.client.list_blobs(a.derived, prefix=f"rehearsal/{a.arm}/")
            if b.name.endswith(".json")]
    if not rows:
        sys.exit(f"no results under rehearsal/{a.arm}/")
    res = pd.DataFrame(rows).sort_values("dataset").reset_index(drop=True)

    body = res.to_csv(index=False)
    derived.blob(f"results/rehearsal_binding_{a.arm}.csv").upload_from_string(body)
    outdir = Path("results/tables")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"rehearsal_binding_{a.arm}.csv").write_text(body)

    print(f"\n{'':=<74}")
    print(f"{len(res)} datasets, {res.protein.nunique()} distinct proteins, "
          f"{res.pairs.sum():,} pairs   [arm={a.arm}]")
    print(f"{'':=<74}\n")
    print(f"{'':24} {'median':>8} {'mean':>8} {'min':>8} {'max':>8}")
    for col, lab in (("auroc", "baseline AUROC"),
                     ("composition_auroc", "composition alone"),
                     ("delta_auroc", "gain over composition")):
        c = res[col]
        print(f"{lab:24} {c.median():8.4f} {c.mean():8.4f} {c.min():8.4f} {c.max():8.4f}")

    print("\ngain over composition, distribution:")
    for lo, hi, lab in ((-1, 0.005, "<0.005 (composition explains it)"),
                        (0.005, 0.02, "0.005-0.02"), (0.02, 0.05, "0.02-0.05"),
                        (0.05, 1, ">0.05 (substantial)")):
        n = int(((res.delta_auroc >= lo) & (res.delta_auroc < hi)).sum())
        print(f"  {lab:36} {n:4d} ({100 * n / len(res):4.1f}%)")
    print(f"\n  interval excludes zero: {int(res.helps.sum())}/{len(res)}")

    # Per cell line, as MEANS, because that is how doc 24 tabulates them and the point of
    # rerunning is to be able to compare directly.
    #
    # NOTE ON A NUMBER THIS DOES NOT COMPUTE. The published "cost of proper matching",
    # -0.0975 in K562 and -0.0998 in HepG2, is the model's AUROC on the GC-matched arm
    # minus its AUROC on the dinucleotide-matched arm. It needs BOTH arms and therefore
    # cannot come out of a single-arm run. An earlier version of this function printed
    # composition_auroc - auroc under that heading, which is a different quantity entirely.
    print(f"\nper cell line (means, comparable with doc 24)   [arm={a.arm}]")
    print(f"  {'':6} {'model':>9} {'composition':>12} {'gain':>9}   n")
    for cell, sub in res.groupby("cell"):
        print(f"  {cell:6} {sub.auroc.mean():9.4f} {sub.composition_auroc.mean():12.4f} "
              f"{sub.delta_auroc.mean():+9.4f}   {len(sub)}")

    lp = np.log(res.pairs)
    print("\nconfound check, correlation with log(pairs):")
    for c in ("auroc", "composition_auroc", "delta_auroc"):
        print(f"  {c:20} r = {np.corrcoef(lp, res[c])[0, 1]:+.3f}")
    print(f"\nwrote gs://{a.derived}/results/rehearsal_binding_{a.arm}.csv "
          f"and results/tables/")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["manifest", "run", "aggregate"])
    p.add_argument("--derived", default=os.environ.get("DERIVED_BUCKET"))
    p.add_argument("--project", default=os.environ.get("GOOGLE_CLOUD_PROJECT"))
    p.add_argument("--arm", default=os.environ.get("ARM", "dinuc"),
                   choices=sorted(panelmod.ARMS))
    p.add_argument("--k", type=int, default=int(os.environ.get("KMER_K", 4)))
    p.add_argument("--n-boot", type=int, default=int(os.environ.get("N_BOOT", 2000)))
    p.add_argument("--min-pairs", type=int, default=None)
    p.add_argument("--index", type=int, default=None)
    p.add_argument("--arm-only", action="store_true",
                   help="manifest for --arm alone. Default builds EVERY arm in one job so "
                        "no local process has to sequence them.")
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if not a.derived:
        sys.exit("--derived or DERIVED_BUCKET required")
    if a.min_pairs is None:
        a.min_pairs = cfgmod.load().cv["min_pairs"]
    {"manifest": do_manifest, "run": do_run, "aggregate": do_aggregate}[a.mode](a)


if __name__ == "__main__":
    main()
