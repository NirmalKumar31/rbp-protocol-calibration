#!/usr/bin/env bash
# What have we actually spent, and what is running right now?
#
# Run this before and after every pipeline stage. Budget alerts are necessary but they are
# LAGGING -- GCP billing data can be several hours behind, so an alert at 25% tells you
# about money already gone. This checks the two things that are immediate:
#
#   1. what compute is running RIGHT NOW (the only thing that can still be stopped)
#   2. what is stored (the only recurring charge when nothing is running)
#
# Usage: bash cloud/cost.sh
# Project and buckets come from the environment, never from a literal. A hardcoded id is
# how a pipeline ends up only running on its author's account. Override with:
#   export GOOGLE_CLOUD_PROJECT=your-project
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(.venv/bin/python -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())')}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"

set -uo pipefail

PROJECT="${PROJECT:-${PROJECT_ID}}"
# NO DEFAULT. A billing account ID is a credential-adjacent identifier and must not live in a
# public repository; this used to carry the real one as a default so the script "just worked".
BILLING="${BILLING:?set BILLING to your billing account ID, e.g. export BILLING=XXXXXX-XXXXXX-XXXXXX}"

echo "=========================================================="
echo " COST CHECK  $(date '+%Y-%m-%d %H:%M')  project=$PROJECT"
echo "=========================================================="

# --- 1. Anything running that costs money by the second -------------------------------
echo
echo "RUNNING COMPUTE (this is what can still be stopped)"

vms=$(gcloud compute instances list --project "$PROJECT" \
        --format="value(name,machineType.basename(),status)" 2>/dev/null)
if [ -z "$vms" ]; then echo "  VMs:              none"; else echo "$vms" | sed 's/^/  VM: /'; fi

batch=$(gcloud batch jobs list --project "$PROJECT" --location us-central1 \
          --filter="status.state:(QUEUED OR SCHEDULED OR RUNNING)" \
          --format="value(name.basename(),status.state)" 2>/dev/null)
if [ -z "$batch" ]; then echo "  Batch jobs:       none active"; else echo "$batch" | sed 's/^/  Batch: /'; fi

vertex=$(gcloud ai custom-jobs list --project "$PROJECT" --region us-central1 \
           --filter="state:(JOB_STATE_PENDING OR JOB_STATE_RUNNING)" \
           --format="value(displayName,state)" 2>/dev/null)
if [ -z "$vertex" ]; then echo "  Vertex jobs:      none active"; else echo "$vertex" | sed 's/^/  Vertex: /'; fi

# --- 2. Stored data, the recurring charge when idle -----------------------------------
echo
echo "STORED DATA (recurring; ~\$0.02/GB/month for standard storage)"
total=0
for b in raw derived artifacts tfstate; do
  bytes=$(gcloud storage du -s "gs://${PROJECT}-${b}" 2>/dev/null | awk '{print $1}')
  bytes=${bytes:-0}
  total=$((total + bytes))
  printf "  %-12s %8.2f GB\n" "$b" "$(echo "$bytes / 1073741824" | bc -l)"
done
printf "  %-12s %8.2f GB  -> ~\$%.2f/month\n" "TOTAL" \
  "$(echo "$total / 1073741824" | bc -l)" \
  "$(echo "$total / 1073741824 * 0.02" | bc -l)"

imgs=$(gcloud artifacts docker images list \
         "us-central1-docker.pkg.dev/${PROJECT}/rbp" --format="value(IMAGE)" 2>/dev/null | wc -l | tr -d ' ')
echo "  container images: ${imgs:-0}  (~\$0.10/GB/month)"

# --- 3. Where to see actual billed spend ----------------------------------------------
echo
echo "BILLED SPEND (lags by up to ~24h -- the console is authoritative)"
echo "  https://console.cloud.google.com/billing/${BILLING}/reports?project=${PROJECT}"
echo "  Credit remaining:"
echo "  https://console.cloud.google.com/billing/${BILLING}/credits"
echo
echo "Budget: \$100 cap, alerts at \$25 / \$50 / \$80 / \$100 + forecast."
echo "Estimated total for the whole pipeline: ~\$34 (docs/26)."
echo
echo "PANIC BUTTON -- stop everything that can bill:"
echo "  gcloud batch jobs list --project $PROJECT --location us-central1 \\"
echo "    --format='value(name)' | xargs -I{} gcloud batch jobs delete {} --quiet"
echo "  gcloud ai custom-jobs list --project $PROJECT --region us-central1 \\"
echo "    --format='value(name)' | xargs -I{} gcloud ai custom-jobs cancel {} --quiet"
