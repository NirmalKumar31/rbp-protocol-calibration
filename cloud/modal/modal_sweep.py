"""Run the sweep on Modal, because this GCP project cannot create a GPU.

    modal run cloud/modal/modal_sweep.py::probe        # ONE task, measure before committing
    modal run cloud/modal/modal_sweep.py::sweep        # fan out
    modal run cloud/modal/modal_sweep.py::status       # how far along

WHY THIS EXISTS. `GPUS_ALL_REGIONS` is 0 on this project and the increase request is
auto-denied with NOT_ENOUGH_USAGE_HISTORY, for 8, for 4, and for 1. AWS returns 0 on all
four GPU families too (adjustable, but a request away). Azure forbids GPU quota on a free
trial outright. Modal gates nothing.

WHAT MAKES THIS CHEAP TO BUILD. A training task's only cloud dependency is
google-cloud-storage: it reads the manifest and its dataset from GCS, writes scores,
metrics and weights back, and takes its index from an environment variable with a default.
No Batch API calls, no metadata-server assumptions. So `scripts/cloud_train.py` runs here
UNCHANGED, and `cloud_train.py aggregate` cannot tell the difference afterwards.

THE PART THAT ACTUALLY MATTERS, AND IT IS NOT THE FREE CREDIT. GCP caps this project at
CPUS_ALL_REGIONS = 12, which is why the sweep runs 12 tasks at a time and why SpliceBERT
would take 155 hours there. Modal has no equivalent cap, so the same work fans out:

    SpliceBERT, GCP CPU, 12 vCPU cap    ~155 h   (6.5 days)
    SpliceBERT, Modal, 10 concurrent      ~2 h

The constraint that shaped every decision for two days simply is not here.

CREDENTIALS. Off GCP there is no metadata server, so this needs a real key. It uses
`rbp-modal`, an identity created solely for this, whose write access is bounded by an IAM
condition to `runs/` and `ckpt/`. Fully compromised it cannot alter a dataset, touch the
raw bucket, or delete anything outside those two prefixes. See cloud/terraform/iam.tf.
"""

import os
import subprocess
import sys

# src/ on the path before the rbp import: `modal run` executes this file directly, so the
# package is not importable unless we say where it lives.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "src"))
import modal

from rbp.utils import cloud as cloudcfg  # noqa: E402

APP = "rbp-sweep"
PROJECT = cloudcfg.project()
DERIVED = cloudcfg.derived_bucket()
SECRET = "rbp-gcp"          # holds SERVICE_ACCOUNT_JSON for rbp-modal

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# --- the image ------------------------------------------------------------------------
#
# Built by Modal rather than pulled from our Artifact Registry, for two concrete reasons:
# the GCP image is ~7 GB and pulling it out of GCS costs egress on every cold machine, and
# pulling a private GCP registry would need a SECOND credential here purely to fetch an
# image. Modal builds this once and caches it.
#
# Versions are pinned to exactly what requirements-gpu.txt pins, so the numbers are
# comparable with the GCP runs. torch comes from the CUDA wheel index.
#
# LAYER ORDER IS THE SAME DISCIPLINE AS THE DOCKERFILE: dependencies, then the 1.2 GB of
# model weights, then our source last, so editing a script does not re-download the weights.
image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("torch==2.13.0", index_url="https://download.pytorch.org/whl/cu126")
    .pip_install(
        "numpy==2.5.2", "pandas==3.0.5", "scipy==1.18.0", "scikit-learn==1.9.0",
        "PyYAML==6.0.3", "requests==2.34.2",
        "transformers==5.15.1", "tokenizers==0.22.2", "safetensors==0.8.0",
        "huggingface_hub==1.28.0", "accelerate==1.14.0",
        "multimolecule==0.2.1", "peft==0.20.0",
        "google-cloud-storage==2.19.0",
    )
    .env({"HF_HOME": "/opt/hf", "PYTHONUNBUFFERED": "1"})
    # Bake the weights, same reasoning as docker/Dockerfile.gpu: a run must not depend on
    # huggingface.co being up, and weights fetched at run time are not pinned to anything.
    .add_local_dir(f"{REPO}/config", "/app/config", copy=True)
    .add_local_dir(f"{REPO}/src", "/app/src", copy=True)
    .add_local_file(f"{REPO}/docker/bake_weights.py", "/app/docker/bake_weights.py", copy=True)
    .workdir("/app")
    .env({"PYTHONPATH": "/app/src"})
    .run_commands("python docker/bake_weights.py fetch")
    # Offline from here on, so a missing weight fails loudly instead of reaching out.
    .env({"HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1"})
    # Source last: it changes constantly.
    .add_local_dir(f"{REPO}/scripts", "/app/scripts")
)

app = modal.App(APP, image=image)


def _env(model: str):
    """Everything scripts/cloud_train.py reads from the environment.

    MANIFEST_TAG is how scope is expressed: Cloud Batch does not run tasks in manifest
    order, so which model runs is decided by WHICH MANIFEST, never by ordering. The same
    rule holds here even though Modal's map() is explicit, because both platforms share
    one code path and one convention.
    """
    return {
        "GOOGLE_CLOUD_PROJECT": PROJECT,
        "DERIVED_BUCKET": DERIVED,
        "WORK_DIR": "/tmp/rbp",
        "ARM": "dinuc",
        "MANIFEST_TAG": f"_{model}",
        "PLATFORM": "modal",
        "GOOGLE_APPLICATION_CREDENTIALS": "/tmp/sa.json",
        # One BLAS thread per task, same as the GCP image. Thread count changes the
        # summation order inside BLAS, so this is about reproducibility as well as speed.
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    }


def _run_one(idx: int, model: str, gpu: bool = True, epochs: int = 0) -> int:
    """Materialise the credential, then hand off to the unmodified pipeline."""
    import json
    import pathlib

    # The key arrives as a Modal secret in the environment; the google client wants a file.
    # /tmp only, never baked into the image, gone when the container dies.
    raw = os.environ["SERVICE_ACCOUNT_JSON"]
    pathlib.Path("/tmp/sa.json").write_text(raw)
    json.loads(raw)                      # fail now, loudly, if it is not valid JSON

    env = {**os.environ, **_env(model)}
    cmd = [sys.executable, "scripts/cloud_train.py", "run", "--index", str(idx)]
    if not gpu:
        cmd += ["--device", "cpu"]
    if epochs:
        # Benchmarking. Caps epochs so a timing probe on the largest dataset costs cents
        # instead of forty minutes; --force to re-run a dataset that already has a marker;
        # --bench so NOTHING is uploaded. That last flag was added after a 2-epoch timing
        # run landed in the results table looking like a real one.
        cmd += ["--epochs", str(epochs), "--force", "--bench"]
    return subprocess.run(cmd, env=env, cwd="/app").returncode


# --- one task, on a GPU ---------------------------------------------------------------
#
# T4 is the cheapest GPU Modal offers and is plenty for a 20M-parameter model over
# 200-token sequences.
#
# MAX_CONTAINERS IS A BUDGET CONTROL, NOT A PERFORMANCE KNOB, AND THIS IS THE ONE NUMBER
# TO UNDERSTAND BEFORE RUNNING ANYTHING HERE.
#
# On GCP the spend rate was capped by quota whether we liked it or not: CPUS_ALL_REGIONS=12
# meant three e2-standard-4 nodes and $0.24/hour, so even total abandonment took six days
# to reach the $40 killswitch. Modal removes that cap, which is exactly why it is useful --
# and it removes the accidental cost ceiling along with it. Fan-out multiplies burn rate:
#
#     concurrency x A10G price    = burn rate      $30 credit lasts
#      1 x $1.10                  = $1.10/h        ~27 h
#     10 x $1.10                  = $11.00/h       ~2.7 h
#     50 x $1.10                  = $55.00/h       ~33 min
#
# Note the price is the GPU alone. cpu= and memory= on a GPU function set the allocation but
# are NOT billed on top -- verified against Modal's dashboard after an estimate built by
# summing three published prices came out 44% high.
#
# So 10 is deliberate: fast enough to finish SpliceBERT in a couple of hours, slow enough
# that a mistake costs a few dollars rather than the whole credit before anyone notices.
# There is no Modal equivalent of the billing killswitch built in killswitch.tf, so this
# cap IS the guardrail.
MAX_CONTAINERS = 10

# timeout is generous because the largest dataset is 32,387 pairs. retries=2 because
# preemption and transient GCS errors both happen, and a retry resumes from the GCS
# checkpoint rather than restarting -- same machinery as Cloud Batch.
# cpu=2 is not decoration. Modal's default CPU allocation is a fraction of a core, and the
# GPU is fed by a CPU-side dataloader that tokenises 200-nucleotide strings one batch at a
# time. Starve that and the T4 idles waiting for input -- you pay GPU rates for CPU work,
# which is the same failure the device guard in cloud_train.py exists to prevent, arriving
# through a different door. Two cores per task, one BLAS thread, matching the GCP shape.
# A10G, measured. Not T4 and not A100: at 1.98x T4 for 1.42x the price it is the
# cheapest per unit of work, because Modal bills CPU and memory per hour alongside
# the GPU so halving the hours more than pays for a dearer card. A100 measured only
# 2.89x T4 -- a 20M-parameter model does not saturate one -- so it costs more.
@app.function(gpu="A10G", cpu=2.0, memory=4096, timeout=60 * 60 * 2,
              max_containers=MAX_CONTAINERS,
              secrets=[modal.Secret.from_name(SECRET)],
              retries=modal.Retries(max_retries=2))
def task(idx: int, model: str, epochs: int = 0) -> int:
    return _run_one(idx, model, gpu=True, epochs=epochs)


@app.function(timeout=60 * 20, secrets=[modal.Secret.from_name(SECRET)])
def manifest_len(model: str) -> int:
    import pathlib

    from google.cloud import storage
    pathlib.Path("/tmp/sa.json").write_text(os.environ["SERVICE_ACCOUNT_JSON"])
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/tmp/sa.json"
    b = storage.Client(project=PROJECT).bucket(DERIVED)
    txt = b.blob(f"manifest/sweep_tasks_{model}.tsv").download_as_text()
    return len(txt.strip().splitlines()) - 1


@app.function(timeout=60 * 20, secrets=[modal.Secret.from_name(SECRET)])
def done_count(model: str) -> int:
    import pathlib

    from google.cloud import storage
    pathlib.Path("/tmp/sa.json").write_text(os.environ["SERVICE_ACCOUNT_JSON"])
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/tmp/sa.json"
    c = storage.Client(project=PROJECT)
    return sum(1 for x in c.list_blobs(DERIVED, prefix="runs/dinuc/")
               if x.name.endswith("metrics.json") and f"/{model}/" in x.name)


# --- entrypoints ----------------------------------------------------------------------

@app.local_entrypoint()
def probe(model: str = "splicebert", index: int = 0):
    """ONE task, timed. Run this before any fan-out.

    Every cost estimate in this project that came from extrapolation has been wrong -- the
    RNABERT slowdown was 4.9x where the CNN's was 1.65x, and a 36-hour projection was
    really 26. The GPU speedup here is a guess between 100x and 200x, so it gets measured
    on one real task, for a few cents of the free credit, before anything is committed.
    """
    import time
    t0 = time.time()
    rc = task.remote(index, model)
    print(f"\nindex {index} of {model}: rc={rc} in {time.time()-t0:.0f}s")
    if rc != 0:
        raise SystemExit(f"probe FAILED (rc={rc}); do not fan out")


@app.local_entrypoint()
def sweep(model: str = "splicebert", limit: int = 0):
    """Fan out. Already-complete runs cost one GCS existence check each and return.

    `limit` caps how many indices are dispatched, for a cautious first batch. Concurrency
    is fixed at MAX_CONTAINERS because it is a budget control -- see the note there.
    """
    n = manifest_len.remote(model)
    idx = list(range(n if limit <= 0 else min(limit, n)))
    print(f"{model}: manifest {n} runs, dispatching {len(idx)}, "
          f"{MAX_CONTAINERS} at a time (~${MAX_CONTAINERS*0.59:.2f}/h), "
          f"already-done are skipped")
    bad = sum(1 for rc in task.map(idx, [model] * len(idx)) if rc != 0)
    print(f"\ndone. {len(idx)-bad} ok, {bad} failed")
    print(f"completed {model} runs now in GCS: {done_count.remote(model)}")


@app.local_entrypoint()
def status(model: str = "splicebert"):
    print(f"{model}: {done_count.remote(model)} / {manifest_len.remote(model)} complete")


@app.local_entrypoint()
def bench(model: str = "splicebert", index: int = 1, epochs: int = 2,
          gpus: str = "T4,A10G,A100"):
    """Time the SAME dataset on several GPUs, then price the sweep from measurement.

    WHY THIS EXISTS AS ITS OWN ENTRYPOINT. Modal charges per-hour for CPU and memory
    alongside the GPU, and every run pays a fixed cost to load a model and fetch a dataset.
    Both amortise over wall time, so a faster GPU can be CHEAPER overall even at a higher
    hourly rate -- on paper an A100 came out at $28 for the full SpliceBERT sweep against
    $73 on a T4.

    That conclusion rested on guessed relative speeds (A100 = 6x T4). SpliceBERT is only
    20M parameters and may not saturate an A100 at all, in which case the guess is wildly
    wrong in the expensive direction. Measuring three GPUs on one dataset for two epochs
    costs a few cents and replaces the guess with a number.

    Every extrapolation in this project has been wrong: the CNN cloud penalty was 1.65x and
    RNABERT's was 4.9x; the GPU speedup was guessed at 100-200x and measured at 29.6x.
    """
    import time
    out = []
    for g in [x.strip() for x in gpus.split(",")]:
        f = task.with_options(gpu=g)
        t0 = time.time()
        rc = f.remote(index, model, epochs)
        el = time.time() - t0
        out.append((g, rc, el))
        print(f"  {g:6} rc={rc}  {el:6.0f}s wall for {epochs} epochs "
              f"(includes cold start + model load + dataset fetch)")
    print("\n  wall time includes fixed overhead, so read the RATIOS not the absolutes;")
    print("  per-epoch cost is in the GCS checkpoint; rates are in docs/COST.md.")
    return out
