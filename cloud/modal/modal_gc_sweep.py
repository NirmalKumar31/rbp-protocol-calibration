"""The GC arm of the sweep, on Modal, with no GCS anywhere.

    modal run cloud/modal/modal_gc_sweep.py::upload                  # push inputs, once
    modal run cloud/modal/modal_gc_sweep.py::sweep --datasets 12     # the pilot
    modal run cloud/modal/modal_gc_sweep.py::sweep --yes             # the full 94
    modal run cloud/modal/modal_gc_sweep.py::status

WHY A SECOND SWEEP FILE. `modal_sweep.py` ran the dinucleotide arm and hardcodes
`"ARM": "dinuc"`, but that is not the reason it cannot be reused. It reads its datasets from
`gs://rbp-repro-2026-derived` and writes every score back there, and that project's billing
account has been closed: every object returns 403. The inputs all still exist on local disk,
so this file keeps the training path byte-for-byte and replaces only the storage.

WHAT CHANGED, AND WHAT DELIBERATELY DID NOT.

  * Inputs live on a Modal Volume, mounted READ-ONLY. Nothing writes to shared storage
    during a sweep, so there is no commit to conflict over and no way for one task to
    corrupt another's input.
  * Outputs come back as the function's RETURN VALUE. A run produces about 4 KB of scores
    and metrics, so 470 runs is ~2 MB. The local store is the only source of truth.
  * Weights are NOT returned. `cloud_train.py` still writes best.pt into the container's
    own store and it is discarded with the container. The dinucleotide arm kept them for
    the locality probe; that result is cut, and 470 full fine-tunes would be ~37 GB.
  * Resume is driver-side: an index whose metrics.json is already in the local store is
    never dispatched. Same rule as the GCS completion marker, same guarantee.
  * A10G, 10 at a time, retries=2 -- all unchanged, because they were measured. See the
    note in modal_sweep.py: A10G is 1.98x a T4 for 1.42x the price, and an A100 measured
    only 2.89x because a 20M-parameter model does not saturate one.

THE COST GUARD. `--yes` is required before anything is dispatched, and the estimate printed
first is built from the dinucleotide arm's actual bill rather than from published prices --
an estimate assembled by summing three list prices once came out 44% high, and reached $6
out of pocket before anyone compared it to the invoice.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

# WHERE "the repo" IS depends on which side of the wire this module is imported on. Locally
# it is three levels up from cloud/modal/. On Modal the file is mounted at
# /root/modal_gc_sweep.py, which has no third parent at all, and the image puts config, src
# and scripts under /app. modal_sweep.py hid this behind three nested dirname() calls, which
# quietly return "/" instead of raising; that is not better, it just fails later.
_here = Path(__file__).resolve()
_local = _here.parents[2] if len(_here.parents) > 2 else None
REPO = _local if (_local and (_local / "scripts").is_dir()) else Path("/app")
if (REPO / "src").is_dir():
    sys.path.insert(0, str(REPO / "src"))
import modal  # noqa: E402

APP = "rbp-gc-sweep"
VOLUME = "rbp-gc-store"
ARM = "gc"
STORE = REPO.parent / "rbp-store"

# The dinucleotide SpliceBERT sweep: 475 runs over 458,336 pairs, billed ~$31.
# Cost scales with pairs, so this is the only honest unit to estimate a new arm in.
COST_PER_PAIR = 31.0 / 458336
MAX_CONTAINERS = 10

vol = modal.Volume.from_name(VOLUME, create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("torch==2.13.0", index_url="https://download.pytorch.org/whl/cu126")
    .pip_install(
        "numpy==2.5.2", "pandas==3.0.5", "scipy==1.18.0", "scikit-learn==1.9.0",
        "PyYAML==6.0.3", "requests==2.34.2",
        "transformers==5.15.1", "tokenizers==0.22.2", "safetensors==0.8.0",
        "huggingface_hub==1.28.0", "accelerate==1.14.0",
        "multimolecule==0.2.1", "peft==0.20.0",
    )
    .env({"HF_HOME": "/opt/hf", "PYTHONUNBUFFERED": "1"})
    .add_local_dir(f"{REPO}/config", "/app/config", copy=True)
    .add_local_dir(f"{REPO}/src", "/app/src", copy=True)
    .add_local_file(f"{REPO}/docker/bake_weights.py", "/app/docker/bake_weights.py",
                    copy=True)
    .workdir("/app")
    .env({"PYTHONPATH": "/app/src"})
    .run_commands("python docker/bake_weights.py fetch")
    # Offline from here on, so a missing weight fails loudly instead of reaching out.
    .env({"HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1"})
    .add_local_dir(f"{REPO}/scripts", "/app/scripts")
)

app = modal.App(APP, image=image)


# --- the local side: the store is the source of truth -----------------------------------

def manifest_rows(model):
    f = STORE / "manifest" / f"sweep_tasks_{model}.tsv"
    if not f.exists():
        sys.exit(f"no manifest at {f}. Run:\n"
                 f"  python scripts/cloud_train.py manifest --arm {ARM} "
                 f"--models {model} --tag _{model} --store {STORE}")
    lines = f.read_text().strip().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


def done(rows):
    """Indices already finished, by the same rule the GCS marker enforced."""
    out = set()
    for r in rows:
        p = (STORE / "runs" / ARM / r["cell_line"] / r["protein"] / r["model"]
             / f"fold{r['fold']}" / "metrics.json")
        if p.exists():
            out.add(int(r["idx"]))
    return out


def stratified(rows, n_datasets):
    """A size-unbiased subset of whole datasets, never a subset of windows.

    Whole datasets, because the dinucleotide arm was trained on full datasets and capping
    the GC arm would confound protocol with training-set size. Size-STRATIFIED rather than
    cheapest-first, because the contrast is known to grow with dataset size (rho +0.307,
    p=0.0026), so the cheap end is the one subset guaranteed to be biased.
    """
    by_ds = {}
    for r in rows:
        by_ds.setdefault(f"{r['protein']}:{r['cell_line']}", int(r["pairs"]))
    order = sorted(by_ds, key=lambda d: (by_ds[d], d))
    if n_datasets >= len(order):
        return set(order)
    step = len(order) / n_datasets
    return {order[int(i * step)] for i in range(n_datasets)}


# --- one task, on a GPU -----------------------------------------------------------------

@app.function(gpu="A10G", cpu=2.0, memory=4096, timeout=60 * 60 * 2,
              max_containers=MAX_CONTAINERS,
              retries=modal.Retries(max_retries=2),
              volumes={"/vol": vol})
def task(idx: int, model: str, epochs: int = 0):
    """Train one fold and hand the two small artefacts back to the driver."""
    env = dict(os.environ)
    env.update({
        "WORK_DIR": "/tmp/rbp",
        "STORE_DIR": "/tmp/out",       # writable, container-local, thrown away
        "STORE_RO": "/vol",            # the inputs, read-only
        "ARM": ARM,
        "MANIFEST_TAG": f"_{model}",
        "PLATFORM": "modal",
        # One BLAS thread, as on every other platform this study has run on: thread count
        # changes summation order inside BLAS, so this is reproducibility, not speed.
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    })
    cmd = [sys.executable, "scripts/cloud_train.py", "run", "--index", str(idx),
           "--arm", ARM, "--store", "/tmp/out", "--store-ro", "/vol"]
    if epochs:
        cmd += ["--epochs", str(epochs), "--force", "--bench"]
    rc = subprocess.run(cmd, env=env, cwd="/app").returncode

    out = []
    root = Path("/tmp/out/runs")
    if root.exists():
        for p in sorted(root.rglob("*")):
            if p.is_file() and p.name in ("metrics.json", "scores.tsv.gz"):
                out.append((str(p.relative_to("/tmp/out")), p.read_bytes()))
    return idx, rc, out


# --- entrypoints ------------------------------------------------------------------------

@app.local_entrypoint()
def upload(model: str = "splicebert"):
    """Push the inputs the sweep reads: panels, the study panel, manifests, datasets."""
    rows = manifest_rows(model)
    want = {(r["cell_line"], r["protein"]) for r in rows}
    files = []
    for key in ("manifest/study_panel.tsv", f"manifest/sweep_tasks_{model}.tsv"):
        files.append((STORE / key, key))
    for cell in ("K562", "HepG2"):
        key = f"panel/{ARM}/panel_final_{cell}_{ARM}.tsv"
        files.append((STORE / key, key))
    for cell, protein in sorted(want):
        key = f"processed/{ARM}/{cell}/{protein}/dataset.tsv"
        files.append((STORE / key, key))

    total = sum(p.resolve().stat().st_size for p, _ in files)
    print(f"uploading {len(files)} files, {total/1e6:.0f} MB -> volume {VOLUME}")
    with vol.batch_upload(force=True) as up:
        for local, key in files:
            up.put_file(local.resolve(), f"/{key}")
    print(f"done. {len(want)} datasets on the volume.")


@app.local_entrypoint()
def sweep(model: str = "splicebert", datasets: int = 0, yes: bool = False):
    rows = manifest_rows(model)
    if datasets:
        keep = stratified(rows, datasets)
        rows = [r for r in rows if f"{r['protein']}:{r['cell_line']}" in keep]
    finished = done(rows)
    todo = [r for r in rows if int(r["idx"]) not in finished]

    n_ds = len({(r["cell_line"], r["protein"]) for r in rows})
    pairs = sum(int(r["pairs"]) for r in todo) / 5      # each fold trains on one dataset
    est = pairs * COST_PER_PAIR
    print(f"\n{model} / arm={ARM}: {n_ds} datasets, {len(rows)} runs, "
          f"{len(finished)} already done, {len(todo)} to run")
    print(f"  estimated ${est:.2f} at {MAX_CONTAINERS} x A10G "
          f"(${MAX_CONTAINERS*1.10:.2f}/h, so roughly {est/(MAX_CONTAINERS*1.10):.1f} h)")
    print("  basis: the dinucleotide arm's actual bill, $31 for 458,336 pairs")
    if not todo:
        return
    if not yes:
        print("\nnothing dispatched. Re-run with --yes to spend this.")
        return

    ok = failed = 0
    for idx, rc, files in task.map([int(r["idx"]) for r in todo],
                                   [model] * len(todo)):
        if rc != 0 or not files:
            failed += 1
            print(f"  task {idx} FAILED rc={rc}")
            continue
        for rel, blob in files:
            dst = STORE / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(blob)
        ok += 1
        if ok % 25 == 0:
            print(f"  {ok}/{len(todo)} done")
    print(f"\ndone. {ok} ok, {failed} failed. Store: {STORE}/runs/{ARM}")


@app.local_entrypoint()
def status(model: str = "splicebert"):
    rows = manifest_rows(model)
    finished = done(rows)
    print(f"{model} / arm={ARM}: {len(finished)} / {len(rows)} runs complete")
    if finished:
        secs = 0.0
        for r in rows:
            p = (STORE / "runs" / ARM / r["cell_line"] / r["protein"] / r["model"]
                 / f"fold{r['fold']}" / "metrics.json")
            if p.exists():
                secs += json.loads(p.read_text()).get("seconds", 0.0)
        print(f"  {secs/3600:.2f} GPU-h of training so far")
