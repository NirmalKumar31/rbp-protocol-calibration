#!/usr/bin/env bash
# Submit the GPU sweep: one Batch job per region, because that is the only way to get
# more than one GPU.
#
#   bash cloud/submit_sweep.sh smoke        # 5 runs, one per region, then stop
#   bash cloud/submit_sweep.sh half         # the first half of the manifest
#   bash cloud/submit_sweep.sh all          # everything not already done
#   bash cloud/submit_sweep.sh cancel       # delete every sweep job in every region
#
# WHY FIVE JOBS AND NOT ONE. PREEMPTIBLE_NVIDIA_V100_GPUS is 1 per region and V100 hardware
# exists in exactly five regions. A Batch job is submitted to one location, and a subnet is
# regional, so five GPUs means five jobs, five subnets, five submissions. The manifest is
# shared and each job strides through it with SHARD/NSHARDS, so the work is interleaved
# rather than sliced -- a contiguous slice would give one region every large SpliceBERT run.
#
# `set -e` is deliberately absent, same as submit_prep.sh: a failed submission in region 3
# must not stop regions 4 and 5 from being reported.
# Project and buckets come from the environment, never from a literal. A hardcoded id is
# how a pipeline ends up only running on its author's account. Override with:
#   export GOOGLE_CLOUD_PROJECT=your-project
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(.venv/bin/python -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())')}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"

set -uo pipefail

PROJECT=${PROJECT_ID}
RAW=${RAW}
DERIVED=${DERIVED}
# rbp-train, NOT rbp-prep. Separate identities per stage is the whole point of iam.tf:
# this one can read the datasets, write only under runs/ and ckpt/, and append to
# BigQuery. It cannot touch raw/ or overwrite a dataset.
SA=rbp-train@${PROJECT}.iam.gserviceaccount.com
OUT=cloud/jobs/rendered
ARM=${ARM:-dinuc}

# Must match var.gpu_regions in cloud/terraform/variables.tf, and each must have a
# rbp-gpu-<region> subnet. Checked below rather than assumed.
REGIONS=(us-central1 us-east1 us-west1 europe-west4 asia-east1)
NSHARDS=${#REGIONS[@]}

# n1 is the only family V100 attaches to. 2 vCPU because CPUS_ALL_REGIONS is 12 globally
# and 5 x 2 = 10 leaves headroom; the GPU is the bottleneck, not the host.
MACHINE=n1-standard-2
GPU_TYPE=nvidia-tesla-v100
BOOT_GB=60          # the GPU image is ~7 GB and the driver install needs room

MODE=${1:?usage: submit_sweep.sh smoke|half|all|cancel}

if [ "$MODE" = cancel ]; then
  for r in "${REGIONS[@]}"; do
    gcloud batch jobs list --project="$PROJECT" --location="$r" \
      --filter="name:sweep-" --format='value(name)' 2>/dev/null |
      while read -r j; do
        echo "deleting $j"
        gcloud batch jobs delete "$j" --location="$r" --quiet
      done
  done
  echo; echo "confirm nothing survives:"
  gcloud compute instances list --project="$PROJECT"
  exit 0
fi

DIGEST=$(gcloud storage cat gs://${PROJECT}-artifacts/images/gpu_digest.txt 2>/dev/null)
if [ -z "$DIGEST" ]; then
  echo "no GPU image digest published; build it first:" >&2
  echo "  gcloud builds submit --config docker/cloudbuild.gpu.yaml \\" >&2
  echo "    --substitutions=_GIT_SHA=\$(git rev-parse --short HEAD) ." >&2
  exit 1
fi
IMAGE=us-central1-docker.pkg.dev/${PROJECT}/rbp/gpu@${DIGEST}

TOTAL=$(gcloud storage cat "gs://${DERIVED}/manifest/sweep_tasks.tsv" 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
if [ -z "$TOTAL" ] || [ "$TOTAL" -eq 0 ]; then
  echo "no manifest; run: python scripts/cloud_train.py manifest" >&2
  exit 1
fi

case "$MODE" in
  smoke) LIMIT=$NSHARDS ;;
  half)  LIMIT=$((TOTAL / 2)) ;;
  all)   LIMIT=$TOTAL ;;
  *) echo "unknown mode: $MODE" >&2; exit 2 ;;
esac

mkdir -p "$OUT"
STAMP=$(date +%m%d-%H%M%S)
echo "manifest $TOTAL runs, submitting $LIMIT across $NSHARDS regions, arm=$ARM"
echo "image $IMAGE"
echo

# Per-region task counts, computed before the prompt so what is being approved is visible.
for i in "${!REGIONS[@]}"; do
  r=${REGIONS[$i]}
  n=$(( (LIMIT - i + NSHARDS - 1) / NSHARDS ))
  [ "$n" -lt 0 ] && n=0
  echo "  $r  shard $i  $n tasks"
done
echo
read -r -p "submit ${NSHARDS} jobs as sweep-${MODE}-${STAMP}? [y/N] " ok
[ "$ok" = y ] || { echo "not submitted"; exit 0; }

for i in "${!REGIONS[@]}"; do
  r=${REGIONS[$i]}
  n=$(( (LIMIT - i + NSHARDS - 1) / NSHARDS ))
  [ "$n" -le 0 ] && continue
  JOB=sweep-${MODE}-${STAMP}-s${i}
  SPEC=${OUT}/${JOB}.json

  cat > "$SPEC" <<JSON
{
  "taskGroups": [{
    "taskCount": ${n},
    "parallelism": 1,
    "taskCountPerNode": 1,
    "taskSpec": {
      "computeResource": {"cpuMilli": 1800, "memoryMib": 6000},
      "maxRunDuration": "10800s",
      "maxRetryCount": 3,
      "runnables": [{
        "script": {"text": "mkdir -p /var/lib/rbp && chmod 1777 /var/lib/rbp"}
      }, {
        "container": {
          "imageUri": "${IMAGE}",
          "entrypoint": "python",
          "commands": ["scripts/cloud_train.py", "run"],
          "options": "--gpus all",
          "volumes": ["/var/lib/rbp:/mnt/disks/work:rw"]
        },
        "environment": {"variables": {
          "RAW_BUCKET": "${RAW}",
          "DERIVED_BUCKET": "${DERIVED}",
          "GOOGLE_CLOUD_PROJECT": "${PROJECT}",
          "WORK_DIR": "/mnt/disks/work",
          "ARM": "${ARM}",
          "SHARD": "${i}",
          "NSHARDS": "${NSHARDS}",
          "IMAGE_DIGEST": "${DIGEST}",
          "OMP_NUM_THREADS": "1",
          "OPENBLAS_NUM_THREADS": "1",
          "MKL_NUM_THREADS": "1",
          "NUMEXPR_NUM_THREADS": "1"
        }}
      }]
    }
  }],
  "allocationPolicy": {
    "instances": [{
      "installGpuDrivers": true,
      "policy": {
        "machineType": "${MACHINE}",
        "provisioningModel": "SPOT",
        "bootDisk": {"sizeGb": ${BOOT_GB}},
        "accelerators": [{"type": "${GPU_TYPE}", "count": 1}]
      }
    }],
    "network": {
      "networkInterfaces": [{
        "network": "projects/${PROJECT}/global/networks/rbp-net",
        "subnetwork": "projects/${PROJECT}/regions/${r}/subnetworks/rbp-gpu-${r}",
        "noExternalIpAddress": true
      }]
    },
    "serviceAccount": {"email": "${SA}"}
  },
  "logsPolicy": {"destination": "CLOUD_LOGGING"}
}
JSON

  # Keep the error. A submission that fails in region 3 must not stop 4 and 5, but
  # swallowing the reason turns a one-line diagnosis into a rerun.
  if err=$(gcloud batch jobs submit "$JOB" --project="$PROJECT" --location="$r" \
             --config="$SPEC" 2>&1 >/dev/null); then
    echo "  submitted $JOB  ($r, $n tasks)"
  else
    echo "  FAILED     $JOB  ($r)" >&2
    echo "$err" | sed 's/^/      /' >&2
  fi
done

echo
echo "watch:   bash cloud/watch_sweep.sh"
echo "stop:    bash cloud/submit_sweep.sh cancel"
