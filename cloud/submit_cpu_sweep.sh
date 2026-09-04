#!/usr/bin/env bash
# The sweep, on CPU, because this project has no GPU quota.
#
#   bash cloud/submit_cpu_sweep.sh smoke 2    # two tasks, prove the path
#   bash cloud/submit_cpu_sweep.sh all        # everything in the manifest
#   bash cloud/submit_cpu_sweep.sh cancel
#
# WHY A SEPARATE FILE FROM submit_sweep.sh. That one is GPU-shaped in every dimension:
# n1 machines because V100 requires them, one task per node because there is one GPU per
# node, five regions because GPU quota is one per region, and five separate jobs striding a
# shared manifest as a result. None of that applies here. Bending it with a flag would leave
# both paths harder to read than two honest files.
#
# THE SHAPE HERE. GPUS_ALL_REGIONS is 0 and cannot be raised, so training runs on CPU.
# CPUS_ALL_REGIONS is 12 globally, so 3 x e2-standard-4 is the entire budget: 12 vCPU, one
# thread per task, four tasks per node. One region is enough because nothing is
# region-scarce; us-central1 has the data and the registry.
#
# It uses the GPU IMAGE. That image is the one with torch and the baked model weights, and
# a CUDA build of torch runs perfectly on a machine with no GPU -- torch.cuda.is_available()
# simply returns False. Building a second near-identical image to avoid pulling 7 GB three
# times would cost more in build minutes than the pull costs in time.
# Project and buckets come from the environment, never from a literal. A hardcoded id is
# how a pipeline ends up only running on its author's account. Override with:
#   export GOOGLE_CLOUD_PROJECT=your-project
# THE INTERPRETER IS A PARAMETER, and it has to be. This said `.venv/bin/python`, which does
# not exist in this repository: the working environment is a venv in a sibling project, so the
# substitution produced an EMPTY project id and the next line then reported "no image digest
# published", which points at the image rather than at the interpreter. Fall back to python3 and
# fail loudly if the id cannot be resolved.
PY="${PY:-python3}"
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$("$PY" -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())' 2>/dev/null)}"
if [ -z "${PROJECT_ID}" ]; then
  echo "cannot resolve the GCP project id. Set GOOGLE_CLOUD_PROJECT, or point PY at an" >&2
  echo "interpreter that can import src/rbp (PY=/path/to/python $0 ...)." >&2
  exit 1
fi
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"

set -uo pipefail

PROJECT=${PROJECT_ID}
DERIVED=${DERIVED}
SA=rbp-train@${PROJECT}.iam.gserviceaccount.com
OUT=cloud/jobs/rendered
REGION=us-central1
ARM=${ARM:-dinuc}

MACHINE=e2-standard-4
NODES=3               # x 4 vCPU = 12 = CPUS_ALL_REGIONS, the hard global cap
PER_NODE=4            # one task per vCPU, one BLAS thread per task
PAR=$((NODES * PER_NODE))
BOOT_GB=50            # ~7 GB image + transient checkpoints and datasets

MODE=${1:?usage: submit_cpu_sweep.sh smoke [n]|all|cancel}
# Which frozen manifest this job reads. One job per model set, because Batch does
# not run tasks in manifest order so scope must be expressed by CONTENTS.
TAG=${TAG:-}

if [ "$MODE" = cancel ]; then
  gcloud batch jobs list --project="$PROJECT" --location="$REGION" \
    --filter="name:cpusweep-" --format='value(name)' 2>/dev/null |
    while read -r j; do echo "deleting $j"; gcloud batch jobs delete "$j" --location="$REGION" --quiet; done
  echo; gcloud compute instances list --project="$PROJECT"
  exit 0
fi

DIGEST=$(gcloud storage cat gs://${PROJECT}-artifacts/images/gpu_digest.txt 2>/dev/null)
[ -z "$DIGEST" ] && { echo "no image digest published; build docker/cloudbuild.gpu.yaml first" >&2; exit 1; }
IMAGE=us-central1-docker.pkg.dev/${PROJECT}/rbp/gpu@${DIGEST}

TOTAL=$(gcloud storage cat "gs://${DERIVED}/manifest/sweep_tasks${TAG}.tsv" 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
[ -z "$TOTAL" ] || [ "$TOTAL" -eq 0 ] && { echo "no manifest sweep_tasks${TAG}.tsv; run scripts/cloud_train.py manifest --tag '${TAG}'" >&2; exit 1; }

case "$MODE" in
  smoke) COUNT=${2:-2}; PAR=$COUNT; PER_NODE=$COUNT ;;
  all)   COUNT=$TOTAL ;;
  *) echo "unknown mode: $MODE" >&2; exit 2 ;;
esac

mkdir -p "$OUT"
JOB=cpusweep${TAG//_/-}-${MODE}-$(date +%m%d-%H%M%S)
SPEC=${OUT}/${JOB}.json
NODES_USED=$(( (PAR + PER_NODE - 1) / PER_NODE ))

cat > "$SPEC" <<JSON
{
  "taskGroups": [{
    "taskCount": ${COUNT},
    "parallelism": ${PAR},
    "taskCountPerNode": ${PER_NODE},
    "taskSpec": {
      "computeResource": {"cpuMilli": 900, "memoryMib": 3500},
      "maxRunDuration": "14400s",
      "maxRetryCount": 3,
      "runnables": [{
        "script": {"text": "mkdir -p /var/lib/rbp && chmod 1777 /var/lib/rbp"}
      }, {
        "container": {
          "imageUri": "${IMAGE}",
          "entrypoint": "python",
          "commands": ["scripts/cloud_train.py", "run", "--device", "cpu"],
          "volumes": ["/var/lib/rbp:/mnt/disks/work:rw"]
        },
        "environment": {"variables": {
          "DERIVED_BUCKET": "${DERIVED}",
          "GOOGLE_CLOUD_PROJECT": "${PROJECT}",
          "WORK_DIR": "/mnt/disks/work",
          "ARM": "${ARM}",
          "MANIFEST_TAG": "${TAG}",
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
      "policy": {
        "machineType": "${MACHINE}",
        "provisioningModel": "SPOT",
        "bootDisk": {"sizeGb": ${BOOT_GB}}
      }
    }],
    "network": {
      "networkInterfaces": [{
        "network": "projects/${PROJECT}/global/networks/rbp-net",
        "subnetwork": "projects/${PROJECT}/regions/${REGION}/subnetworks/rbp-workers",
        "noExternalIpAddress": true
      }]
    },
    "serviceAccount": {"email": "${SA}"}
  },
  "logsPolicy": {"destination": "CLOUD_LOGGING"}
}
JSON

echo "manifest $TOTAL runs | submitting $COUNT | $NODES_USED x $MACHINE spot | $PAR concurrent | arm=$ARM"
echo "image  $IMAGE"
echo "spec   $SPEC"
echo
if [ "${CONFIRM:-ask}" = ask ]; then
  read -r -p "submit as ${JOB}? [y/N] " ok
  [ "$ok" = y ] || { echo "not submitted"; exit 0; }
fi

if err=$(gcloud batch jobs submit "$JOB" --project="$PROJECT" --location="$REGION" \
           --config="$SPEC" 2>&1 >/dev/null); then
  echo "submitted $JOB"
  echo
  echo "watch:  bash cloud/watch_sweep.sh"
  echo "stop:   bash cloud/submit_cpu_sweep.sh cancel"
else
  echo "FAILED to submit $JOB" >&2
  echo "$err" | sed 's/^/    /' >&2
  exit 1
fi
