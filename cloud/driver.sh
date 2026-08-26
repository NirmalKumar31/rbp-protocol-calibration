#!/usr/bin/env bash
# The remaining pipeline, sequenced FROM A GCP VM instead of from a laptop.
#
#   ./cloud/driver.sh launch    start it on a VM and return; the laptop is then free
#   ./cloud/driver.sh watch     tail its progress from anywhere
#   ./cloud/driver.sh stop      delete the VM
#
# WHY THIS EXISTS. Every individual step already survives a sleeping laptop: Cloud Build runs
# server-side, Batch jobs run server-side, and Modal apps survive with --detach. What did NOT
# survive was the TRANSITIONS -- the shell loop that waits for one thing and starts the next.
# That loop lived on the laptop, so closing the lid stopped the pipeline between stages even
# though nothing was actually running there.
#
# This is the same defect, twice already fixed at a smaller scale and twice reappearing:
#   * stage 7 sequenced two Batch jobs from a shell loop  -> fixed by one job for both arms
#   * modal run without --detach tied 475 GPU tasks to a  -> fixed by --detach
#     local client
# Both times the work was in the cloud and the CONTROL was not. This moves the control too.
#
# WHY A PLAIN VM AND NOT WORKFLOWS. Cloud Workflows is the "right" GCP answer and would mean
# learning a YAML dialect, granting it permissions, and expressing a shell pipeline as a state
# machine. A g1-small running the script we already have costs about a cent an hour and is
# debuggable with ssh. For a pipeline that runs a handful of times, the boring option wins.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

PROJECT="${GOOGLE_CLOUD_PROJECT:?set GOOGLE_CLOUD_PROJECT}"
DERIVED="${DERIVED_BUCKET:-${PROJECT}-derived}"
ZONE="${ZONE:-us-central1-c}"
VM="rbp-driver"

launch() {
  # The driver needs to submit Batch jobs and read/write GCS, so it runs as rbp-analysis --
  # the identity that already has exactly those rights. It does NOT need to be able to create
  # VMs or change IAM.
  gcloud compute instances create "$VM" \
    --project="$PROJECT" --zone="$ZONE" \
    --machine-type=e2-small \
    --service-account="rbp-analysis@${PROJECT}.iam.gserviceaccount.com" \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --image-family=debian-12 --image-project=debian-cloud \
    --boot-disk-size=20GB \
    --metadata-from-file=startup-script=cloud/driver_startup.sh \
    --metadata="rbp-project=${PROJECT},rbp-derived=${DERIVED}" \
    || { echo "could not create the driver VM" >&2; exit 1; }
  echo "driver VM ${VM} created in ${ZONE}"
  echo "it writes progress to gs://${DERIVED}/driver/driver.log"
  echo "watch with: ./cloud/driver.sh watch"
}

watch() {
  while :; do
    clear 2>/dev/null || true
    echo "=== driver log $(date '+%H:%M:%S') ==="
    gsutil cat "gs://${DERIVED}/driver/driver.log" 2>/dev/null | tail -25 || echo "(no log yet)"
    if gsutil -q stat "gs://${DERIVED}/driver/DONE" 2>/dev/null; then
      echo; echo "DRIVER FINISHED"; break
    fi
    sleep 30
  done
}

stop() {
  gcloud compute instances delete "$VM" --project="$PROJECT" --zone="$ZONE" --quiet \
    && echo "driver VM deleted"
}

case "${1:-}" in
  launch) launch ;;
  watch)  watch ;;
  stop)   stop ;;
  *) sed -n '2,12p' "$0" ;;
esac
