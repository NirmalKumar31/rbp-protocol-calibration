#!/usr/bin/env bash
# One Batch submitter for every stage. Replaces seven near-identical scripts.
#
#   ./cloud/submit.sh ingest        one VM, EXTERNAL IP (needs public internet)
#   ./cloud/submit.sh panel         one VM, external IP (ENCODE API)
#   ./cloud/submit.sh prep          fan out over every candidate, sealed
#   ./cloud/submit.sh rehearsal     fan out, sealed
#   ./cloud/submit.sh variants      one VM, EXTERNAL IP (UCSC phyloP over HTTP)
#   ./cloud/submit.sh analysis      one VM, sealed
#
# THE BUG THIS DESIGN EXISTS TO PREVENT. submit_prep.sh hardcoded COUNT=189. The gc arm has
# 187 datasets, so the job dispatched two tasks past the end of its manifest, and Batch
# reported the whole job FAILED even though every real task had succeeded. Task count is
# ALWAYS read from the manifest in GCS here. A count you typed is a count that will be wrong
# the first time the panel changes size.
#
# NETWORK POSTURE, and it is not uniform. Workers get Private Google Access and no external
# IP, so they can reach *.googleapis.com and nothing else. Three stages genuinely need the
# public internet -- ingest (ENCODE, GENCODE, NCBI), panel (ENCODE API) and variants (UCSC
# phyloP) -- and those get an external IP on a single short-lived VM. There is deliberately
# no Cloud NAT: it is billed per VM-hour plus per GB and would hand internet access to every
# other worker for no reason.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

JOB_TYPE="${1:?usage: submit.sh <ingest|panel|prep|rehearsal|sweep|variants|analysis>}"

PY="${PY:-python3}"
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
PROJECT="${GOOGLE_CLOUD_PROJECT:?set GOOGLE_CLOUD_PROJECT}"
DERIVED="${DERIVED_BUCKET:-${PROJECT}-derived}"
RAW="${RAW_BUCKET:-${PROJECT}-raw}"
REGION="${REGION:-us-central1}"
ARM="${ARM:-dinuc}"

# --- how many tasks? From the manifest, never from a literal ---------------------------
manifest_rows() {
  local key=$1
  local n
  n=$(gcloud storage cat "gs://${DERIVED}/${key}" 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
  [ -n "$n" ] && [ "$n" -gt 0 ] || { echo "EMPTY_MANIFEST"; return 1; }
  echo "$n"
}

# --- per-stage shape -------------------------------------------------------------------
# EXTERNAL=1 means the VM gets a public IP because the stage must reach the internet.
case "$JOB_TYPE" in
  ingest)
    SA=rbp-ingest
    SCRIPT="scripts/cloud_ingest.py"; ARGS="--bucket ${RAW}"
    COUNT=1; PAR=1; PER_NODE=1; MACHINE=e2-standard-4; CPU=4000; MEM=16384
    EXTERNAL=1; DISK=200; TIMEOUT=10800 ;;
  panel)
    SA=rbp-ingest
    SCRIPT="scripts/build_panel.py"; ARGS="--all"
    COUNT=1; PAR=1; PER_NODE=1; MACHINE=e2-standard-2; CPU=2000; MEM=8192
    EXTERNAL=1; DISK=50; TIMEOUT=3600 ;;
  prep)
    # prep_tasks.tsv, NOT study_panel.tsv. Prep runs before the panel exists -- `pairs` is
    # a result of preprocessing, so the size-ranked sample cannot be taken until prep is
    # done. Pointing this at the study panel would have made stage 5 unrunnable on a fresh
    # project, which is exactly the class of error a fresh project is meant to catch.
    # No --arm: each task reads its own arm from its manifest row (cloud_prep.py:173), so
    # this is ONE job covering both negative arms, not one job per arm.
    SA=rbp-prep
    SCRIPT="scripts/cloud_prep.py"; ARGS="prep"
    COUNT=$(manifest_rows "manifest/prep_tasks.tsv") || exit 1
    PAR=12; PER_NODE=4; MACHINE=e2-standard-4; CPU=900; MEM=3500
    EXTERNAL=0; DISK=100; TIMEOUT=7200 ;;
  rehearsal)
    SA=rbp-train
    SCRIPT="scripts/cloud_rehearsal.py"; ARGS="run --arm ${ARM}"
    COUNT=$(manifest_rows "manifest/rehearsal_tasks.tsv") || exit 1
    PAR=12; PER_NODE=4; MACHINE=e2-standard-4; CPU=900; MEM=3500
    EXTERNAL=0; DISK=50; TIMEOUT=7200 ;;
  variants)
    SA=rbp-ingest
    SCRIPT="scripts/cloud_variants.py"; ARGS="--what all"
    COUNT=1; PAR=1; PER_NODE=1; MACHINE=e2-standard-4; CPU=4000; MEM=16384
    EXTERNAL=1; DISK=200; TIMEOUT=14400 ;;
  sweep)
    # The CNN arm. Task count is models x datasets x folds and comes from the manifest that
    # cloud_train.py wrote, never from arithmetic done here.
    SA=rbp-train
    SCRIPT="scripts/cloud_train.py"; ARGS="run --arm ${ARM}"
    COUNT=$(manifest_rows "manifest/sweep_tasks${MANIFEST_TAG:-}.tsv") || exit 1
    PAR=12; PER_NODE=4; MACHINE=e2-standard-4; CPU=900; MEM=3500
    EXTERNAL=0; DISK=100; TIMEOUT=14400 ;;
  analysis)
    SA=rbp-analysis
    SCRIPT="scripts/cloud_analysis.py"; ARGS="--what all"
    COUNT=1; PAR=1; PER_NODE=1; MACHINE=e2-standard-4; CPU=4000; MEM=16384
    EXTERNAL=0; DISK=100; TIMEOUT=7200 ;;
  *) echo "unknown job type: $JOB_TYPE" >&2; exit 1 ;;
esac

# --- the image, pinned BY DIGEST -------------------------------------------------------
# A tag is mutable. Pinning by digest is what makes "which image produced this result?"
# answerable, and it is the same reason the model weights are baked into the image.
DIGEST=$(gcloud storage cat "gs://${PROJECT}-artifacts/images/cpu_digest.txt" 2>/dev/null | tr -d '[:space:]')
if [ -z "$DIGEST" ]; then
  echo "no image digest at gs://${PROJECT}-artifacts/images/cpu_digest.txt -- run stage 2" >&2
  exit 1
fi
IMAGE="${REGION}-docker.pkg.dev/${PROJECT}/rbp/cpu@${DIGEST}"

# NESTING AND NAMES BOTH MATTERED HERE. Batch wants
# allocationPolicy.network.networkInterfaces, not allocationPolicy.networkInterfaces -- the
# API rejects the latter outright with `Unknown name "networkInterfaces"`. And the network is
# rbp-net, the one Terraform creates, not "default": the whole point of network.tf is that
# workers sit on a subnet with Private Google Access and no route to the internet.
#
# The external case omits the block entirely, which puts the VM on the default network with
# an external IP. That is deliberate for the three stages that must reach ENCODE, GENCODE,
# NCBI and UCSC, and it is why there is no Cloud NAT.
if [ "$EXTERNAL" = "1" ]; then
  NETWORK=""
else
  NETWORK='"network": {"networkInterfaces": [{"network": "projects/'"${PROJECT}"'/global/networks/rbp-net", "subnetwork": "projects/'"${PROJECT}"'/regions/'"${REGION}"'/subnetworks/rbp-workers", "noExternalIpAddress": true}]},'
fi

JOB="${JOB_TYPE}-$(date +%m%d-%H%M%S)"
mkdir -p cloud/jobs/rendered
SPEC="cloud/jobs/rendered/${JOB}.json"

cat > "$SPEC" <<JSON
{
  "taskGroups": [{
    "taskCount": ${COUNT},
    "parallelism": ${PAR},
    "taskCountPerNode": ${PER_NODE},
    "taskSpec": {
      "maxRetryCount": 2,
      "maxRunDuration": "${TIMEOUT}s",
      "computeResource": {"cpuMilli": ${CPU}, "memoryMib": ${MEM}},
      "runnables": [{
        "container": {
          "imageUri": "${IMAGE}",
          "entrypoint": "python",
          "commands": ["${SCRIPT}"$(for a in ${ARGS}; do printf ', "%s"' "$a"; done)]
        },
        "environment": {"variables": {
          "GOOGLE_CLOUD_PROJECT": "${PROJECT}",
          "DERIVED_BUCKET": "${DERIVED}",
          "RAW_BUCKET": "${RAW}",
          "ARM": "${ARM}",
          "MANIFEST_TAG": "${MANIFEST_TAG:-}",
          "MODELS": "${MODELS:-}",
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
        "bootDisk": {"sizeGb": ${DISK}, "type": "pd-balanced"}
      }
    }],
    ${NETWORK}
    "serviceAccount": {"email": "${SA}@${PROJECT}.iam.gserviceaccount.com"}
  },
  "logsPolicy": {"destination": "CLOUD_LOGGING"}
}
JSON

echo "job=${JOB} script=${SCRIPT} ${ARGS}"
echo "tasks=${COUNT} parallelism=${PAR} per_node=${PER_NODE} nodes=$(( (PAR + PER_NODE - 1) / PER_NODE ))"
echo "machine=${MACHINE} spot external_ip=${EXTERNAL} image=@${DIGEST:0:19}..."

if [ "${RBP_YES:-0}" != "1" ]; then
  read -r -p "submit? [y/N] " ok
  [ "$ok" = "y" ] || { echo "not submitted"; exit 1; }
fi

gcloud batch jobs submit "$JOB" --project="$PROJECT" --location="$REGION" --config="$SPEC" \
  || { echo "submit failed" >&2; exit 1; }
echo "submitted ${JOB}"
echo "watch: gcloud batch jobs describe ${JOB} --project=${PROJECT} --location=${REGION} --format='value(status.state,status.taskGroups)'"
