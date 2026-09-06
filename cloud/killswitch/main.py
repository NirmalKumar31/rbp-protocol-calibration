"""Disable billing on the project when spend crosses the budget. The actual hard stop.

WHY THIS EXISTS. A GCP budget alert sends an email. It does not stop anything. If a job
runs away at 3am the alert arrives, nobody reads it, and the spend continues until someone
wakes up. The only mechanism that genuinely halts charges is detaching the project from its
billing account, and nothing in GCP does that for you.

HOW IT IS WIRED. The budget publishes a JSON message to a Pub/Sub topic every time it
re-evaluates -- roughly every 20-30 minutes, not only when a threshold trips. This function
is triggered by that topic, compares actual spend against the limit, and detaches billing
if it is over.

WHAT HAPPENS WHEN IT FIRES, STATED THE WAY GOOGLE STATES IT. Every VM stops, Batch jobs
fail, Cloud Run stops serving. Google's own documentation warns that disabling billing may
DELETE some resources and that the deletion can be non-recoverable; it does not promise that
buckets survive, and this docstring used to, ending on the words "never lost data". That was a
guarantee we were not in a position to give.

What is true: charges already incurred remain payable, billing reports lag by up to a day so
the figure that triggered this is not the final one, and anything whose loss would matter must
be backed up OUTSIDE the project this can fire on. Treat a firing as data loss until you have
checked otherwise.

RE-ENABLING IS DELIBERATELY MANUAL:

    gcloud billing projects link $PROJECT_ID \
      --billing-account=$BILLING_ACCOUNT

If this could re-enable itself the guardrail would be pointless.

WHAT THE DRY RUN PROVES, AND WHAT IT DOES NOT. Setting KILL_THRESHOLD_USD low enough to trip
on current spend exercises Pub/Sub delivery, message parsing, the threshold decision, and the
getBillingInfo READ. It stops before updateBillingInfo, so it does NOT prove the service
account may perform the write that is the actual stop.

That distinction was missing here for the whole project, and the docstring claimed the dry run
proved "permissions" without qualification. It cannot: billing-account access and the
project-level updateBillingInfo permission are granted separately in GCP, and this function
only ever reads. The two are checked apart now -- `permission_check()` below asks the IAM API
whether the caller holds the write, which is a read-only question with a real answer.

Firing it for real would take the project down to prove it can take the project down, so the
remaining honest options are the permission check below or a rehearsal in a disposable
project. Neither is a substitute for the other: the check can pass while a deny policy or an
org constraint blocks the call.
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


def _verdict(x):
    return {True: "HELD", False: "NOT HELD", None: "UNKNOWN"}[x]


def permission_check(api=None):
    """Does this identity actually hold the permission that constitutes the stop?

    testIamPermissions is a read: it answers "would this call be allowed" without making it.
    Called at cold start so the answer appears in the logs of every deployment, rather than
    being discovered at 3am by the one invocation that needed to work.

    A False here means the killswitch is decorative. It does not raise, because a killswitch
    that refuses to start is worse than one that starts and warns: the threshold logging and
    the alerting path still have value, and an exception at import would take those too.

    AND A TRUE HERE IS STILL NOT A REHEARSAL. testIamPermissions answers about IAM. An
    organisation policy, a deny policy or a billing-account link lock can block the write with
    every permission held. The only thing that proves the stop works is firing it in a
    disposable project.
    """
    # THE RIGHT PERMISSION, ON THE RIGHT RESOURCE, AND THERE ARE TWO ROUTES. The first version
    # of this asked for resourcemanager.projects.updateBillingInfo, which is the name of the API
    # METHOD and not of a permission Google grants. Unlinking a project is authorised by either
    # billing.resourceAssociations.delete on the billing account, or
    # resourcemanager.projects.deleteBillingAssignment on the project. Testing the wrong string
    # returns "not held" against an identity that can in fact make the call, which is the worse
    # of the two possible errors: it teaches you to ignore the check.
    want = "resourcemanager.projects.deleteBillingAssignment"
    held = None
    try:
        crm = discovery.build("cloudresourcemanager", "v1", cache_discovery=False)
        got = crm.projects().testIamPermissions(
            resource=PROJECT, body={"permissions": [want]}).execute()
        held = want in (got.get("permissions") or [])
    except Exception as e:                                    # noqa: BLE001
        print(f"project-route permission check failed to run: {e}", flush=True)
    print(f"permission {want}: {_verdict(held)}", flush=True)
    if held:
        return True

    # The billing-account route. Checked second because it needs the account id, which the
    # function only has when the project is still linked.
    try:
        api = billing()
        acct = api.projects().getBillingInfo(
            name=PROJECT_NAME).execute().get("billingAccountName", "")
        if not acct:
            print("no billing account linked, so the account route cannot be tested",
                  flush=True)
            return held
        want2 = "billing.resourceAssociations.delete"
        got = api.billingAccounts().testIamPermissions(
            resource=acct, body={"permissions": [want2]}).execute()
        held2 = want2 in (got.get("permissions") or [])
        print(f"permission {want2} on {acct}: {_verdict(held2)}", flush=True)
        return bool(held or held2)
    except Exception as e:                                    # noqa: BLE001
        print(f"billing-account route check failed to run: {e}", flush=True)
        return held


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
        held = permission_check()
        print(f"DRY RUN: would disable billing on {PROJECT} "
              f"(cost {cost:.2f} >= {LIMIT:.2f}); write permission "
              f"{ {True: 'held', False: 'NOT HELD', None: 'unknown'}[held] }", flush=True)
        return
    print(f"DISABLING BILLING on {PROJECT}: cost {cost:.2f} >= {LIMIT:.2f}", flush=True)
    print(disable(api), flush=True)
