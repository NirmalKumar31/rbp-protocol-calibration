"""Cost guard for Modal. The thing killswitch.tf is for GCP, built by hand because Modal
has no budget API.

    python cloud/modal/guard.py                 # report once
    python cloud/modal/guard.py --watch         # loop, and STOP the app if over budget

WHY THIS IS NEEDED HERE SPECIFICALLY. On GCP the spend rate was capped by quota whether we
liked it or not -- CPUS_ALL_REGIONS=12 meant three nodes and $0.24/hour, so even total
abandonment took six days to reach the $40 killswitch. Modal's whole value is that it
removes that cap, and it removes the accidental cost ceiling along with it. Ten A10G
containers burn $11.00/hour, which is roughly two days of the GCP sweep every sixty minutes.
(That figure read $15.80 for months. It was MAX_CONTAINERS times the pre-correction rate of
$1.58, left behind by the very commit that establishes forty lines below why $1.58 was wrong
by 44%. A prose number beside a constant is a second copy of that constant.)

Modal exposes no spend figure over the CLI and no budget limit, so the guard has to be
built from what is observable: how long the app has been up, and how much work has landed
in GCS.

TWO ESTIMATES, AND THE CONSERVATIVE ONE IS THE ONE THAT ACTS.

  upper bound   elapsed x MAX_CONTAINERS x rate. Assumes every container busy every second.
                Always >= reality. This is what triggers the stop, because a guard that
                under-estimates is not a guard.

  lower bound   GPU-seconds actually recorded by finished runs, priced at the same rate.
                Cannot include work still in flight, so it always UNDER-states.

The truth sits between them, and the gap is itself a diagnostic: a wide gap means containers
are idle, which on Modal you still pay for.

AN EARLIER VERSION ADDED A PER-TASK 43-SECOND OVERHEAD to the work figure and produced a
number LARGER than the supposed upper bound -- $14.75 against $13.64. That was the giveaway
that the model was wrong, not the measurement: 43s covers cold start plus model load, and
Modal REUSES warm containers, so it is paid once per container and not once per task. With
87 tasks and ~10 containers the overhead was over-counted roughly eightfold. Removed.

An estimate that exceeds its own upper bound is telling you the model is broken. Believe it.
"""

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import time

# src/ on the path before the rbp import: `modal run` executes this file directly, so the
# package is not importable unless we say where it lives. This block sat ABOVE the imports
# until 2026-08-26, so the file raised NameError on line 3 and had never once run.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "src"))
from rbp.utils import cloud as cloudcfg  # noqa: E402

# A10G per hour. JUST the GPU price -- CPU and memory are NOT additive.
#
# THIS WAS WRONG AND THE ERROR WAS 44%. The first version summed three published prices:
# A10G $1.10 + 2 CPU x $0.192 + 4 GiB x $0.024 = $1.580/h. That produced a SpliceBERT
# estimate of $36.89 against an actual of about $31.50 on Modal's own dashboard.
#
# Back-solving from the real figure across ~29 billed container-hours gives ~$1.09/h, which
# lands on the A10G list price alone. A GPU container bundles a baseline CPU and RAM
# allocation; the per-core and per-GiB prices apply to CPU-only functions, not on top of a
# GPU. Reading three numbers off a pricing page and assuming they sum was the mistake --
# the same shape as every other bad estimate in this project: asserting an unverified number.
#
# The lesson that generalises: a cost model has TWO parts, the work estimate and the rate.
# Ours drifted upward all run and I kept blaming the work estimate. It was the rate. When an
# estimate is biased in one direction rather than noisy, suspect the constant, not the model.
RATE = 1.10
MAX_CONTAINERS = 10
PROJECT = cloudcfg.project()
DERIVED = cloudcfg.derived_bucket()


def _epoch(ts):
    """Modal's created_at is an ISO string with offset; we need an epoch."""
    try:
        return dt.datetime.fromisoformat(ts).timestamp()
    except Exception:
        return None


# EVERY SWEEP APP, NOT ONE. This matched the literal description "rbp-sweep" while the apps
# that actually ran the GC, bias-aware and region-matched arms are named rbp-gc-sweep,
# rbp-neg2-sweep and rbp-neg2-rm-sweep -- modal_gc_sweep.py derives the name from RBP_ARM. So
# the guard reported "no running rbp-sweep app" and returned success for three of the four
# sweeps it exists to guard, on the most expensive arms, and looked healthy doing it.
APP_NAME = re.compile(r"^rbp-[a-z0-9-]*-?sweep$")
LIVE = ("ephemeral", "running", "deployed")


class Unobservable(RuntimeError):
    """The Modal CLI could not be read. NOT the same as nothing running."""


def app_state():
    """The running sweep app, if any: (app_id, started_epoch).

    RAISES rather than returning (None, None) when Modal cannot be read. Those two states used
    to be identical: the returncode was never checked, so an expired token, a network failure
    or a CLI upgrade that changed the output shape all produced an empty parse, which the
    caller printed as "no running app" and exited zero. A cost guard that cannot see must fail
    closed, because the failure it guards against is invisible by definition.
    """
    r = subprocess.run(["modal", "app", "list", "--json"], capture_output=True, text=True)
    if r.returncode != 0:
        raise Unobservable(f"modal app list exited {r.returncode}: "
                           f"{(r.stderr or r.stdout).strip()[:200]}")
    try:
        apps = json.loads(r.stdout)
    except json.JSONDecodeError as e:
        raise Unobservable(f"modal app list returned unparseable output: {e}") from e
    for a in apps:
        # Modal's own field names have varied; be liberal about which key holds what.
        state = str(a.get("State") or a.get("state") or "").lower()
        desc = a.get("Description") or a.get("description") or ""
        if APP_NAME.match(desc) and state in LIVE:
            return (a.get("app_id") or a.get("App ID"),
                    _epoch(a.get("created_at") or a.get("Created at") or ""))
    return None, None


# ARMS TO COUNT. The prefix was the literal "runs/dinuc/", so work landing under runs/gc/,
# runs/neg2/ or runs/neg2_rm/ contributed zero to the lower bound -- the guard's only
# measurement of what has actually been paid for. Combined with the app-name bug above, the
# guard was blind at both ends for every arm except the first one written.
ARMS = ("dinuc", "gc", "neg2", "neg2_rm")


def work_done(model, arms=ARMS):
    """(runs, gpu_seconds) from GCS, over every arm. The receipt, not the promise."""
    from google.cloud import storage
    c = storage.Client(project=PROJECT)
    n, secs = 0, 0.0
    for arm in arms:
        for b in c.list_blobs(DERIVED, prefix=f"runs/{arm}/"):
            if b.name.endswith("metrics.json") and f"/{model}/" in b.name:
                m = json.loads(b.download_as_text())
                if m.get("platform") == "modal":
                    n += 1
                    secs += float(m.get("seconds", 0))
    return n, secs


def report(model, budget, started_epoch):
    n, secs = work_done(model)
    elapsed_h = (time.time() - started_epoch) / 3600 if started_epoch else 0.0
    upper = elapsed_h * MAX_CONTAINERS * RATE
    lower = secs / 3600 * RATE
    print(f"[{dt.datetime.now():%H:%M:%S}] {model}: {n} runs done on modal")
    print(f"    upper bound  ${upper:6.2f}   (elapsed {elapsed_h:.2f} h x {MAX_CONTAINERS}"
          f" x ${RATE:.3f})")
    print(f"    lower bound  ${lower:6.2f}   ({secs/3600:.2f} GPU-h recorded)")
    print(f"    budget       ${budget:6.2f}")
    return upper, lower


def stop(app_id):
    print(f"    STOPPING APP {app_id}", flush=True)
    r = subprocess.run(["modal", "app", "stop", app_id], capture_output=True, text=True)
    print("   ", (r.stdout or r.stderr).strip()[:300])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", default="splicebert")
    # 40 = the $30 free credit plus an authorised $10 buffer, and the buffer is the point.
    #
    # It started at 32 and that was the wrong instinct. At 46% of the SpliceBERT sweep the
    # cost model projected $28.52; at 80% it projected $32.29, because the model runs ~16%
    # optimistic. A ceiling of 32 would therefore have stopped the run at roughly 95%
    # complete -- and stopping near the end is the worst possible outcome, not a near miss:
    # the manifest keeps a dataset's five folds adjacent, so an interruption leaves two or
    # three datasets holding 3 or 4 folds, and a partial fold set cannot be pooled at all.
    # Those datasets are simply lost.
    #
    # So the guard exists to stop a RUNAWAY, not to shave a budget. Set it far enough above
    # the projection that normal drift cannot trip it -- here 1.24x - and let the bounded
    # task count do the actual cost limiting. A guard that fires on ordinary estimation
    # error destroys work for no benefit.
    p.add_argument("--budget", type=float, default=40.0)
    p.add_argument("--watch", action="store_true")
    p.add_argument("--interval", type=int, default=120)
    a = p.parse_args()

    try:
        app_id, started = app_state()
    except Unobservable as e:
        print(f"CANNOT OBSERVE MODAL: {e}")
        print("Refusing to report 'nothing running'. Fix the CLI or the credentials, or stop "
              "the app by hand at modal.com/apps.")
        sys.exit(2)
    if not app_id:
        print("no running sweep app")
        report(a.model, a.budget, None)
        return
    print(f"guarding app {app_id}, budget ${a.budget:.2f}, "
          f"burn ${MAX_CONTAINERS*RATE:.2f}/h at full fan-out\n")
    if not started:               # fall back to now if the timestamp did not parse
        started = time.time()

    while True:
        upper, lower = report(a.model, a.budget, started)
        if upper >= a.budget:
            print(f"    OVER BUDGET on the upper bound (${upper:.2f} >= ${a.budget:.2f})")
            stop(app_id)
            sys.exit(1)
        if not a.watch:
            return
        try:
            app_id2, _s = app_state()
        except Unobservable as e:
            # Mid-watch, this is the dangerous moment: the app may still be burning. Keep
            # watching rather than exiting, and say so, so a transient network blip does not
            # silently end the only supervision the run has.
            print(f"    cannot observe modal ({e}); still watching, app NOT assumed finished")
            time.sleep(a.interval)
            continue
        if not app_id2:
            print("    app finished on its own")
            return
        time.sleep(a.interval)


if __name__ == "__main__":
    main()
