#!/usr/bin/env bash
# Is anything billing right now? Run this any time, trust it over anyone's word.
#
# Checks every category of resource in this project that can charge money, not just the
# ones we happen to use. Billing reports lag up to 24h, so they cannot answer "is money
# being spent right now". A resource list can.
#
# The rule this encodes: with no VM, no Cloud Run service and no reserved IP, the only
# charge is storage at cents per month. Everything expensive requires a running instance.
# Project and buckets come from the environment, never from a literal. A hardcoded id is
# how a pipeline ends up only running on its author's account. Override with:
#   export GOOGLE_CLOUD_PROJECT=your-project
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(.venv/bin/python -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())')}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"

set -uo pipefail

P=${PROJECT:-${PROJECT_ID}}
REGIONS=(us-central1 us-east1 us-west1 europe-west4 asia-east1)
bad=0

say() { printf '  %-34s %s\n' "$1" "$2"; }
check() {                       # name, command
  local out; out=$(eval "$2" 2>/dev/null)
  if [ -z "$out" ]; then say "$1" "none"; else
    say "$1" "*** ACTIVE ***"; echo "$out" | sed 's/^/      /'; bad=1
  fi
}

echo "=============================================================="
echo " BILLING AUDIT   $(date '+%Y-%m-%d %H:%M %Z')   project=$P"
echo "=============================================================="

check "VMs (all zones)" \
  "gcloud compute instances list --project=$P --format='value(name,zone.basename(),status)'"

for r in "${REGIONS[@]}"; do
  check "Batch active, $r" \
    "gcloud batch jobs list --project=$P --location=$r \
       --filter='status.state:(QUEUED OR SCHEDULED OR RUNNING)' \
       --format='value(name.basename(),status.state)'"
done

check "Vertex AI jobs" \
  "gcloud ai custom-jobs list --project=$P --region=us-central1 \
     --filter='state:(JOB_STATE_PENDING OR JOB_STATE_RUNNING)' --format='value(displayName)'"

# Cloud Run, EXCLUDING the killswitch.
#
# A gen2 Cloud Function IS a Cloud Run service, so billing-killswitch shows up here
# permanently. Left in, this script would report "SOMETHING IS ACTIVE" every single time and
# stop being a usable signal -- the classic way a monitor gets ignored. It scales to zero
# (maxScale=1, no minimum) and is invoked about three times an hour against a 2M/month free
# tier, so its cost is zero. It is reported separately as the guardrail it is.
check "Cloud Run services" \
  "gcloud run services list --project=$P --format='value(name)' | grep -v '^billing-killswitch$'"
check "Cloud Build running" \
  "gcloud builds list --project=$P --ongoing --format='value(id)'"

# UNATTACHED disks only. A boot disk of a running VM is already counted by the VM check
# above, and Batch creates them with autoDelete=true so they vanish with it. An ORPHANED
# disk is the real risk: it bills at ~$0.10/GB/month forever with nothing to show it.
check "Unattached disks" \
  "gcloud compute disks list --project=$P --filter='-users:*' --format='value(name,sizeGb)'"
check "Reserved IPs" \
  "gcloud compute addresses list --project=$P --format='value(name,status)'"

# The guardrail, reported but never alarmed on.
ks=$(gcloud functions describe billing-killswitch --project="$P" --region=us-central1 --gen2 \
       --format='value(state,serviceConfig.environmentVariables.DRY_RUN)' 2>/dev/null)
if [ -n "$ks" ]; then
  case "$ks" in
    *false*) say "killswitch" "ACTIVE and ARMED (kills at \$40)" ;;
    *)       say "killswitch" "present but NOT ARMED -- DRY_RUN is true" ;;
  esac
else
  say "killswitch" "*** NOT DEPLOYED ***"
fi

echo
if [ "$bad" -eq 0 ]; then
  echo "  CLEAR. Nothing can bill except stored data (~\$0.15/month)."
else
  echo "  SOMETHING IS ACTIVE. To stop everything immediately:"
  echo "    bash cloud/audit.sh --kill"
fi

# --kill: stop every possible source of charge, in one command. Safe to run when idle.
if [ "${1:-}" = "--kill" ]; then
  echo
  echo "KILLING EVERYTHING"
  for r in "${REGIONS[@]}"; do
    gcloud batch jobs list --project="$P" --location="$r" --format='value(name)' 2>/dev/null |
      while read -r j; do echo "  delete batch $j"; gcloud batch jobs delete "$j" --location="$r" --quiet; done
  done
  gcloud compute instances list --project="$P" --format='value(name,zone.basename())' 2>/dev/null |
    while read -r n z; do echo "  delete vm $n"; gcloud compute instances delete "$n" --zone="$z" --quiet; done
  gcloud builds list --project="$P" --ongoing --format='value(id)' 2>/dev/null |
    while read -r b; do echo "  cancel build $b"; gcloud builds cancel "$b" --quiet; done
  echo
  echo "re-checking:"
  gcloud compute instances list --project="$P" 2>&1 | sed 's/^/  /'
fi
