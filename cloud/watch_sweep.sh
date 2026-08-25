#!/usr/bin/env bash
# Progress across all five sweep regions, from the only source that cannot lie.
#
# Job state is the scheduler's opinion and task state is the container's exit code. Neither
# proves work exists. Completion markers in the bucket do, so that is what gets counted
# first and everything else is context.
# Project and buckets come from the environment, never from a literal. A hardcoded id is
# how a pipeline ends up only running on its author's account. Override with:
#   export GOOGLE_CLOUD_PROJECT=your-project
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(.venv/bin/python -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())')}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"

set -uo pipefail

PROJECT=${PROJECT_ID}
DERIVED=${DERIVED}
ARM=${ARM:-dinuc}
REGIONS=(us-central1 us-east1 us-west1 europe-west4 asia-east1)

total=$(gcloud storage cat "gs://${DERIVED}/manifest/sweep_tasks.tsv" 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
done_n=$(gcloud storage ls "gs://${DERIVED}/runs/${ARM}/**/metrics.json" 2>/dev/null | wc -l | tr -d ' ')

echo "=========================================================="
echo " SWEEP  $(date '+%H:%M:%S')  arm=$ARM"
echo "=========================================================="
echo "  complete: ${done_n:-0} / ${total:-?}"

echo
echo "RUNNING GPUs (this is what bills)"
vms=$(gcloud compute instances list --project "$PROJECT" \
        --format="value(name,zone.basename(),status)" 2>/dev/null)
[ -z "$vms" ] && echo "  none" || echo "$vms" | sed 's/^/  /'

echo
echo "JOBS"
for r in "${REGIONS[@]}"; do
  s=$(gcloud batch jobs list --project="$PROJECT" --location="$r" \
        --filter="name:sweep-" --format='value(name.basename(),status.state)' 2>/dev/null)
  [ -z "$s" ] && s="  (none)"
  echo "  $r"
  echo "$s" | sed 's/^/    /'
done

# Checkpoints still in flight. A number that stays constant while `complete` also stays
# constant means work has stalled, which no job state will tell you.
# NOTE THE PATH. run_prefix() already begins with "runs/", so checkpoints land at
# ckpt/runs/<arm>/... not ckpt/<arm>/... . This looked in the wrong place and therefore
# always reported zero -- silently disabling the one diagnostic that distinguishes "work is
# progressing" from "work keeps starting and never finishing".
ck=$(gcloud storage ls "gs://${DERIVED}/ckpt/runs/${ARM}/**/checkpoint.pt" 2>/dev/null | wc -l | tr -d ' ')
echo
echo "  in-flight checkpoints: ${ck:-0}"
echo "  stop everything:  bash cloud/submit_sweep.sh cancel"
