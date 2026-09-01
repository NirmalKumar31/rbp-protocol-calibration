"""Disable billing on the project when spend crosses the budget. The actual hard stop.

WHY THIS EXISTS. A GCP budget alert sends an email. It does not stop anything. If a job
runs away at 3am the alert arrives, nobody reads it, and the spend continues until someone
wakes up. The only mechanism that genuinely halts charges is detaching the project from its
billing account, and nothing in GCP does that for you.

HOW IT IS WIRED. The budget publishes a JSON message to a Pub/Sub topic every time it
re-evaluates -- roughly every 20-30 minutes, not only when a threshold trips. This function
is triggered by that topic, compares actual spend against the limit, and detaches billing
if it is over.

WHAT HAPPENS WHEN IT FIRES. Every VM stops, Batch jobs fail, Cloud Run stops serving.
Buckets and their contents SURVIVE -- storage is not deleted, it just becomes inaccessible
until billing is re-attached. So the worst case is an interrupted run and a manual
re-enable, never lost data.

RE-ENABLING IS DELIBERATELY MANUAL:

    gcloud billing projects link $PROJECT_ID \
      --billing-account=$BILLING_ACCOUNT

If this could re-enable itself the guardrail would be pointless.

TESTED BY DRY RUN, NOT BY FIRING IT. Setting KILL_THRESHOLD_USD low enough to trip on
current spend proves the whole path -- Pub/Sub delivery, parsing, permissions, the decision
-- while DRY_RUN=true stops it at the last step. Firing it for real would take the project
down to prove it can take the project down.
"""

import base64
import json
import os

import functions_framework
from googleapiclient import discovery

PROJECT = os.environ["TARGET_PROJECT"]
LIMIT = float(os.environ["KILL_THRESHOLD_USD"])
DRY_RUN = os.environ.get("DRY_RUN", "true").lower() == "true"
PROJECT_NAME = f"projects/{PROJECT}"


def billing():
    return discovery.build("cloudbilling", "v1", cache_discovery=False)


def is_enabled(api):
    info = api.projects().getBillingInfo(name=PROJECT_NAME).execute()
    return bool(info.get("billingAccountName", ""))


def disable(api):
    """Detach the billing account. This is the stop.

    An empty billingAccountName is how the API expresses "no billing account", and it is
    the same operation as clicking Disable Billing in the console.
    """
    return api.projects().updateBillingInfo(
        name=PROJECT_NAME, body={"billingAccountName": ""}).execute()


@functions_framework.cloud_event
def handle(event):
    payload = json.loads(base64.b64decode(event.data["message"]["data"]).decode())
    # costAmount is absent from the very first message of a budget period, when nothing
    # has been spent yet. Treat missing as zero rather than crashing, or the function
    # errors on every new month.
    cost = float(payload.get("costAmount", 0) or 0)
    budget = float(payload.get("budgetAmount", 0) or 0)
    name = payload.get("budgetDisplayName", "?")

    print(f"budget={name} cost={cost:.2f} budget={budget:.2f} "
          f"kill_at={LIMIT:.2f} dry_run={DRY_RUN}", flush=True)

    if cost < LIMIT:
        return
    api = billing()
    if not is_enabled(api):
        print("billing already disabled, nothing to do", flush=True)
        return
    if DRY_RUN:
        print(f"DRY RUN: would disable billing on {PROJECT} "
              f"(cost {cost:.2f} >= {LIMIT:.2f})", flush=True)
        return
    print(f"DISABLING BILLING on {PROJECT}: cost {cost:.2f} >= {LIMIT:.2f}", flush=True)
    print(disable(api), flush=True)
