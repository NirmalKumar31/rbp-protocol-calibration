#!/usr/bin/env bash
# Submit a CPU Batch job. One code path, several presets.
#
#   bash cloud/submit_prep.sh index          # build regions.pkl from the GTF (1 task)
#   bash cloud/submit_prep.sh smoke 12,301   # two named manifest rows, the md5 gate
#   bash cloud/submit_prep.sh all            # the full 488 preprocessing tasks
#   bash cloud/submit_prep.sh rerun 7,44,91  # whatever spot preemption ate
#   bash cloud/submit_prep.sh rehearse       # the composition control, 189 datasets
#   bash cloud/submit_prep.sh rehearse-smoke 0,188
#
# The rendered JSON is written to disk and printed before anything is submitted, so what
# runs is inspectable rather than buried in a shell variable.
#
# `set -e` is deliberately absent: a failed gcloud call should still let the script print
# where the job spec went, and every command that matters is checked explicitly.
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
REGION=us-central1
RAW=${RAW}
DERIVED=${DERIVED}
SA=rbp-prep@${PROJECT}.iam.gserviceaccount.com
OUT=cloud/jobs/rendered

MODE=${1:?usage: submit_prep.sh index|smoke|all|rerun [task_list]}
TASK_LIST=${2:-}

# WHY HIGHMEM RATHER THAN THE CHEAPER e2-standard-4. The GC negative matcher pulls sequence
# for candidate after candidate until one lands inside the GC tolerance, so it is thousands
# of small random reads into a 3.1 GB FASTA. Measured on the laptop: 7s per dataset with the
# genome in page cache, and the 971s recorded for the same protein on the same inputs when
# it was cold. 32 GB holds the genome in cache alongside 4 resident tasks; 16 GB does not.
# Catalog price is $0.1085/node-hour against $0.0804, so the insurance costs $0.028/hour.
MEM=3500
CPU=900
MACHINE=e2-highmem-4
SCRIPT=scripts/cloud_prep.py
CMD=prep
# Rehearsal tasks never touch the genome, so they do not need the page cache that made
# highmem worth paying for; e2-standard-4 is $0.0804/node-hour against $0.1085.
REHEARSE_MACHINE=e2-standard-4
case "$MODE" in
  index) COUNT=1;   PAR=1;  PER_NODE=1; MEM=8192; CPU=2000; MACHINE=e2-standard-4 ;;
  smoke) COUNT=$(tr -cd ',' <<<"$TASK_LIST" | wc -c); COUNT=$((COUNT + 1)); PAR=$COUNT; PER_NODE=$COUNT ;;
  rerun) COUNT=$(tr -cd ',' <<<"$TASK_LIST" | wc -c); COUNT=$((COUNT + 1)); PAR=$((COUNT < 28 ? COUNT : 28)); PER_NODE=4 ;;
  # 7 nodes x 4 tasks. Three separate quotas bind here and only the third is obvious:
  #   INSTANCES      8      -> 7 leaves one spare
  #   CPUS          32      -> 7 x 4 vCPU = 28
  #   SSD_TOTAL_GB 250      -> 7 x 30 GB boot = 210
  #   IN_USE_ADDRESSES 4    -> no longer binds: workers take no external IP at all,
  #                            see cloud/terraform/network.tf. Before that they did, and
  #                            this job silently ran on 3 nodes instead of 7.
  #   CPUS_ALL_REGIONS 12   -> THE REAL CEILING, and global rather than regional. 3 nodes
  #                            of 4 vCPU is the most that can ever run at once, so PAR=28
  #                            is aspirational; it degrades to 3 nodes rather than failing.
  #                            Cannot be raised: NOT_ENOUGH_USAGE_HISTORY on a new account.
  # SSD_TOTAL_GB is why there is no second disk. pd-balanced counts against SSD_TOTAL_GB,
  # NOT against DISKS_TOTAL_GB (2048), so an earlier 50 GB boot + 100 GB work disk was
  # 150 GB per node and only ONE node fit. Batch does not fail in that situation: it logged
  # 52 instance-creation errors as OPERATIONAL_INFO, kept the job RUNNING, and quietly
  # delivered a seventh of the throughput. Scratch now lives on the boot disk, shared
  # between the four tasks on a node by a bind mount.
  all)   COUNT=488; PAR=28; PER_NODE=4 ;;
  rehearse)
    # READ the manifest length, do not hardcode it. The gc arm has 187 datasets and the
    # dinuc arm has 189 -- two proteins clear min_pairs on one arm and not the other. With
    # COUNT hardcoded at 189, the gc run dispatches two tasks past the end of the manifest,
    # cloud_rehearsal exits non-zero on them, Batch retries each three times and marks the
    # whole job FAILED even though all 187 real results landed.
    COUNT=$(gcloud storage cat "gs://${DERIVED}/manifest/rehearsal_tasks.tsv" 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
    [ -z "$COUNT" ] || [ "$COUNT" -eq 0 ] && { echo "no rehearsal manifest; run cloud_rehearsal.py manifest" >&2; exit 1; }
    PAR=12; PER_NODE=4; MACHINE=$REHEARSE_MACHINE
    SCRIPT=scripts/cloud_rehearsal.py; CMD=run ;;
  rehearse-smoke)
    COUNT=$(tr -cd ',' <<<"$TASK_LIST" | wc -c); COUNT=$((COUNT + 1))
    PAR=$COUNT; PER_NODE=$COUNT; MACHINE=$REHEARSE_MACHINE
    SCRIPT=scripts/cloud_rehearsal.py; CMD=run ;;
  *) echo "unknown mode: $MODE" >&2; exit 2 ;;
esac

DIGEST=$(gcloud storage cat gs://${PROJECT}-artifacts/images/cpu_digest.txt 2>/dev/null)
if [ -z "$DIGEST" ]; then echo "no image digest published; build the CPU image first" >&2; exit 1; fi
IMAGE=${REGION}-docker.pkg.dev/${PROJECT}/rbp/cpu@${DIGEST}

[ "$MODE" = index ] && CMD=index

mkdir -p "$OUT"
JOB=prep-${MODE}-$(date +%m%d-%H%M%S)
SPEC=${OUT}/${JOB}.json

cat > "$SPEC" <<JSON
{
  "taskGroups": [{
    "taskCount": ${COUNT},
    "parallelism": ${PAR},
    "taskCountPerNode": ${PER_NODE},
    "taskSpec": {
      "computeResource": {"cpuMilli": ${CPU}, "memoryMib": ${MEM}},
      "maxRunDuration": "14400s",
      "maxRetryCount": 3,
      "runnables": [{
        "script": {"text": "mkdir -p /var/lib/rbp && chmod 1777 /var/lib/rbp"}
      }, {
        "container": {
          "imageUri": "${IMAGE}",
          "entrypoint": "python",
          "commands": ["${SCRIPT}", "${CMD}"],
          "volumes": ["/var/lib/rbp:/mnt/disks/work:rw"]
        },
        "environment": {"variables": {
          "RAW_BUCKET": "${RAW}",
          "DERIVED_BUCKET": "${DERIVED}",
          "WORK_DIR": "/mnt/disks/work",
          "TASK_LIST": "${TASK_LIST}",
          "ARM": "${ARM:-dinuc}",
          "N_BOOT": "${N_BOOT:-2000}",
          "KMER_K": "${KMER_K:-4}",
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
        "bootDisk": {"sizeGb": 30}
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

echo "spec: $SPEC"
cat "$SPEC"
echo
echo "mode=$MODE script=$SCRIPT machine=$MACHINE tasks=$COUNT parallelism=$PAR per_node=$PER_NODE nodes=$(( (PAR + PER_NODE - 1) / PER_NODE ))"
read -r -p "submit as ${JOB}? [y/N] " ok
[ "$ok" = y ] || { echo "not submitted"; exit 0; }

gcloud batch jobs submit "$JOB" --project="$PROJECT" --location="$REGION" --config="$SPEC" \
  && echo && echo "watch:  gcloud batch jobs describe $JOB --location=$REGION --format='value(status.state)'"
