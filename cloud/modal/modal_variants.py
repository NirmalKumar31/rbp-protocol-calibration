"""Score the ClinVar variants with SpliceBERT on Modal, because GCP still has no GPU.

    modal run cloud/modal/modal_variants.py::probe     # ONE dataset, measure first
    modal run cloud/modal/modal_variants.py::sweep     # all 94
    modal run cloud/modal/modal_variants.py::status    # how far along

WHY MODAL AND NOT CLOUD BATCH. Same wall as the sweep: `GPUS_ALL_REGIONS` is 0 and
`CPUS_ALL_REGIONS` is 12, so GCP gives three e2-standard-4 nodes, and transformer inference
on those CPUs measured 4.9x slower than the laptop. Three nodes at 1/4.9 speed is slower
than one Mac, so "run it in the cloud" on GCP would have been slower than not bothering.
Modal caps neither, and the work is done before the GCP image would have finished building.

WHY THIS NEEDED NO FASTA. Scoring wants two things: the genome, to cut each variant's
ref/alt windows, and the checkpoints, which are already in GCS. Uploading a 3.1 GB genome
to save a step that costs 27 seconds locally would be absurd, so
`variant_splicebert.py --what tables` cuts the windows on the laptop and uploads 164,835 of
them as ~30 MB of sequence. The cloud task never learns what a FASTA is.

THE IMAGE IS REUSED FROM modal_sweep.py ON PURPOSE. Same pinned versions, same baked
weights, same digest. If this file declared its own image, the variant scores could come
from a different build of SpliceBERT than the binding scores, and the two would no longer
be comparable -- which is the entire point of running them.
"""

import os
import subprocess
import sys

import modal

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from modal_sweep import DERIVED, PROJECT, SECRET  # noqa: E402
from modal_sweep import image as _base  # noqa: E402

# THE sys.path LINE ABOVE ONLY FIXES THE CLIENT. Modal ships the entrypoint file to the
# container as /root/modal_variants.py and nothing else, so `import modal_sweep` resolved
# here and then failed inside every container with ModuleNotFoundError. The containers
# crash-looped at import, five seconds each, doing no work and producing no error the
# client surfaced until the whole run was inspected.
#
# The fix ships the dependency too, rather than duplicating the image spec in this file:
# a second spec would drift, and variant scores produced by a different build of SpliceBERT
# than the binding scores would not be comparable, which is the entire point of running
# them. Same image object, one extra file in it.
image = _base.add_local_file(f"{HERE}/modal_sweep.py", "/root/modal_sweep.py")

APP = "rbp-variants"

# TASK COUNTS COME FROM THE MANIFEST, NEVER FROM A TYPED NUMBER.
#
# N_TASKS was hardcoded to 94, the earlier study's variant panel. This pipeline's manifest has
# 95 rows, so range(N_TASKS) dispatched indices 0..93 and silently skipped index 94, K562
# ZNF800 -- 288 variants dropped from R4 with nothing in any log to say so. The driver's own
# gate was `>= 94`, so it would have been satisfied by 94 of 95 and moved on.
#
# cloud/submit.sh already derives every Batch task count from manifest_rows() for exactly this
# reason. Modal was the one place that did not, because it is not driven by submit.sh.
MANIFEST_OBJ = "variants/variant_tasks.tsv"
N_LOCALITY = 95


def n_tasks():
    """How many datasets the variant manifest actually lists."""
    import io

    import pandas as pd
    from google.cloud import storage

    b = storage.Client(project=PROJECT).bucket(DERIVED).blob(MANIFEST_OBJ)
    if not b.exists():
        raise SystemExit(
            f"gs://{DERIVED}/{MANIFEST_OBJ} is absent. Run "
            f"`cloud_variants.py --what windows` before scoring.")
    return len(pd.read_csv(io.StringIO(b.download_as_text()), sep="\t"))

# T4, not the A10G the sweep used. The sweep was compute-bound -- fine-tuning 20M parameters
# for twelve epochs -- so a card at 1.98x for 1.42x the price was the cheaper unit of work.
# This is the opposite shape: ~3,500 forward passes per dataset against 375 MB of checkpoint
# download, so the task is network-bound and the GPU idles either way. Paying A10G rates to
# wait on a download is paying for the wrong thing. T4 is Modal's cheapest.
#
# max_containers is the budget control, exactly as in modal_sweep.py, and it matters more
# here because the $30 credit is gone and every second is out of pocket. 10 containers of
# T4 burns roughly $6/hour; this job should take about ten minutes.
MAX_CONTAINERS = 10

app = modal.App(APP, image=image)


def _env():
    return {
        "GOOGLE_CLOUD_PROJECT": PROJECT,
        "DERIVED_BUCKET": DERIVED,
        "PLATFORM": "modal",
        "GOOGLE_APPLICATION_CREDENTIALS": "/tmp/sa.json",
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    }


def _run_one(idx: int, force: bool = False, mismatch: int = 0) -> int:
    import json
    import pathlib

    raw = os.environ["SERVICE_ACCOUNT_JSON"]
    pathlib.Path("/tmp/sa.json").write_text(raw)
    json.loads(raw)                      # fail now, loudly, if it is not valid JSON

    cmd = [sys.executable, "scripts/variant_splicebert.py", "--what", "cloud",
           "--index", str(idx)]
    if force:
        cmd.append("--force")
    if mismatch:
        cmd += ["--mismatch", str(mismatch)]
    return subprocess.run(cmd, env={**os.environ, **_env()}, cwd="/app").returncode


@app.function(gpu="T4", cpu=2.0, memory=4096, timeout=60 * 30,
              max_containers=MAX_CONTAINERS,
              secrets=[modal.Secret.from_name(SECRET)],
              retries=modal.Retries(max_retries=2))
def task(idx: int, force: bool = False, mismatch: int = 0) -> int:
    return _run_one(idx, force, mismatch)


@app.local_entrypoint()
def probe(index: int = 0):
    """One dataset, so the cost of the other 93 is measured rather than guessed."""
    import time
    t0 = time.time()
    rc = task.remote(index, True)
    el = time.time() - t0
    print(f"rc={rc} in {el:.0f}s")
    n = n_tasks()
    print(f"projected {n} tasks / {MAX_CONTAINERS} containers: "
          f"{el * n / MAX_CONTAINERS / 60:.1f} min wall, "
          f"${el * n / 3600 * 0.59:.2f} at T4 rates")


@app.local_entrypoint()
def sweep(force: bool = False):
    n = n_tasks()
    rcs = list(task.map(range(n), kwargs={"force": force}))
    bad = [i for i, rc in enumerate(rcs) if rc != 0]
    print(f"{len(rcs) - len(bad)}/{len(rcs)} ok" + (f", failed {bad}" if bad else ""))


@app.local_entrypoint()
def mismatch_sweep(offset: int = 47):
    """The negative control: every dataset scored with ANOTHER protein's fine-tuned head.

    Offset 47 is about half the manifest, which is sorted by pair rank, so a dataset is
    paired with one of a very different size as well as a different protein -- the donor
    cannot be a near-duplicate assay of the same factor.
    """
    n = n_tasks()
    rcs = list(task.map(range(n), kwargs={"force": True, "mismatch": offset}))
    bad = [i for i, rc in enumerate(rcs) if rc != 0]
    print(f"{len(rcs) - len(bad)}/{len(rcs)} ok" + (f", failed {bad}" if bad else ""))


@app.local_entrypoint()
def status():
    from google.cloud import storage
    c = storage.Client(project=PROJECT)
    n = sum(1 for _ in c.list_blobs(DERIVED, prefix="variants/scores_sb/"))
    print(f"{n}/{n_tasks()} datasets scored")


# --- the ISM locality probe, same image, same discipline ---------------------------------
#
# Unlike the variant task this one is genuinely COMPUTE-bound: 3*L mutants per window, 101
# positions, 20 windows, two models, 95 datasets is ~570,000 forward passes against a single
# 75 MB checkpoint download. That is the shape a GPU actually helps with -- ~90 minutes on
# the laptop's CPU, a few minutes here -- which is why it is worth moving and the variant
# job was borderline.


def _run_locality(idx: int, force: bool = False) -> int:
    import json
    import pathlib as _p

    raw = os.environ["SERVICE_ACCOUNT_JSON"]
    _p.Path("/tmp/sa.json").write_text(raw)
    json.loads(raw)

    cmd = [sys.executable, "scripts/locality_probe.py", "--cloud", "--index", str(idx)]
    if force:
        cmd.append("--force")
    return subprocess.run(cmd, env={**os.environ, **_env()}, cwd="/app").returncode


@app.function(gpu="T4", cpu=2.0, memory=4096, timeout=60 * 30,
              max_containers=MAX_CONTAINERS,
              secrets=[modal.Secret.from_name(SECRET)],
              retries=modal.Retries(max_retries=2))
def locality_task(idx: int, force: bool = False) -> int:
    return _run_locality(idx, force)


@app.local_entrypoint()
def locality_probe(index: int = 0):
    """One dataset first, same rule as everywhere else: measure before committing."""
    import time
    t0 = time.time()
    rc = locality_task.remote(index, True)
    el = time.time() - t0
    print(f"rc={rc} in {el:.0f}s")
    print(f"projected {N_LOCALITY} / {MAX_CONTAINERS}: "
          f"{el * N_LOCALITY / MAX_CONTAINERS / 60:.1f} min wall, "
          f"${el * N_LOCALITY / 3600 * 0.59:.2f} at T4 rates")


@app.local_entrypoint()
def locality_sweep(force: bool = False):
    rcs = list(locality_task.map(range(N_LOCALITY), kwargs={"force": force}))
    bad = [i for i, rc in enumerate(rcs) if rc != 0]
    print(f"{len(rcs) - len(bad)}/{len(rcs)} ok" + (f", failed {bad}" if bad else ""))
