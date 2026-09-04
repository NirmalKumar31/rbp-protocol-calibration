"""Stage 0. Refuse to start until the environment can actually finish.

WHY THIS EXISTS, AND WHY IT IS THE MOST VALUABLE FILE HERE. Every failure of the original
three-day build was one of these, and every one was discovered mid-run, after money had been
spent:

  * an API not enabled            -> a job submits and dies minutes later
  * GPU quota 0                   -> two days of work planned around hardware that cannot
                                     be allocated, and the quota request auto-denied
  * a budget reporting $0         -> credit_types_treatment INCLUDE_ALL_CREDITS, so every
                                     alert threshold is unreachable while credit remains
  * an unauthenticated client     -> gcloud returns an EMPTY LIST, not an error, so a
                                     missing scope reads exactly like a missing result
  * a bucket name already taken   -> buckets are globally unique; terraform fails at apply
  * Modal credit exhausted        -> the one stage that costs real money silently bills

None of these need a single cent to detect. All of them cost hours to discover late.

    python scripts/preflight.py                # check, report, exit non-zero on any FAIL
    python scripts/preflight.py --skip-modal   # GCP-only stages

Exit code 0 means every gate passed and the pipeline may start. Anything else means stop.
WARN does not fail the run but is printed loudly, because some checks are advisory (a
budget that exists but is generous) and some are fatal (no billing account).
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

REQUIRED_APIS = [
    "batch.googleapis.com",
    "compute.googleapis.com",
    "storage.googleapis.com",
    "artifactregistry.googleapis.com",
    "cloudbuild.googleapis.com",
    "billingbudgets.googleapis.com",
    "cloudresourcemanager.googleapis.com",
]

# The sweep needs 12 vCPU to run 3 e2-standard-4 concurrently. Below that it still works,
# just slower; at 0 something is wrong with the project.
MIN_CPUS_ALL_REGIONS = 8

results = []


def check(name, ok, detail="", fatal=True):
    status = "PASS" if ok else ("FAIL" if fatal else "WARN")
    results.append((status, name, detail))
    print(f"  [{status}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    return ok


def sh(cmd, timeout=90):
    """Run a command, return (rc, stdout). Never raises; a broken tool is a finding."""
    try:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except Exception as e:
        return 1, str(e)


# --- section 1: is anything even configured? -------------------------------------------

def section_config():
    print("\nconfiguration")
    from rbp.utils import cloud as cloudcfg
    try:
        proj, derived, raw = cloudcfg.project(), cloudcfg.derived_bucket(), cloudcfg.raw_bucket()
    except RuntimeError as e:
        check("project resolves", False, str(e))
        return None
    check("project resolves", True, f"project={proj}")
    check("bucket names derived", True, f"{derived} / {raw}")
    return proj, derived, raw


# --- section 2: authentication and the trap that looks like an empty result -------------

def section_auth(proj):
    print("\nauthentication")
    rc, out = sh("gcloud auth list --filter=status:ACTIVE --format='value(account)'")
    acct = out.strip().splitlines()[0] if rc == 0 and out.strip() else ""
    check("gcloud has an active account", bool(acct), acct)

    # THE TRAP. `gcloud <verb> list` with no project returns an empty list rather than an
    # error, so "not configured" is indistinguishable from "nothing exists". Assert the
    # project is reachable by NAME, which does fail loudly.
    rc, out = sh(f"gcloud projects describe {proj} --format='value(projectId)'")
    check("project exists and is visible", rc == 0 and proj in out,
          out.strip().splitlines()[-1] if out.strip() else "")

    rc, out = sh("gcloud config get-value project 2>/dev/null")
    cfgproj = out.strip()
    check("gcloud default project set", cfgproj == proj,
          f"config={cfgproj or '(unset)'} expected={proj}; "
          "an unset default makes list commands return empty, not error", fatal=False)
    return bool(acct)


# --- section 3: billing, and the budget that reports zero ------------------------------

def section_billing(proj):
    print("\nbilling")
    rc, out = sh(f"gcloud billing projects describe {proj} --format=json")
    linked = False
    if rc == 0:
        try:
            linked = json.loads(out).get("billingEnabled", False)
        except Exception:
            pass
    check("billing account linked", linked,
          "without this every job submission fails at creation")

    rc, out = sh("gcloud billing budgets list --format=json 2>/dev/null || true")
    budgets = []
    try:
        budgets = json.loads(out) if out.strip().startswith("[") else []
    except Exception:
        pass
    check("a budget exists", bool(budgets),
          f"{len(budgets)} found; without one there is no killswitch trigger", fatal=False)

    # The subtle one. A budget whose creditTypesTreatment includes credits reports $0 spend
    # for as long as free credit lasts, so no threshold ever fires and the killswitch is
    # decorative.
    bad = [b.get("displayName", "?") for b in budgets
           if (b.get("budgetFilter", {}) or {}).get("creditTypesTreatment")
           == "INCLUDE_ALL_CREDITS"]
    check("budgets exclude credits", not bad,
          f"INCLUDE_ALL_CREDITS on {bad}: these report $0 while credit remains, so no "
          f"alert can fire" if bad else "", fatal=not not bad)


def section_apis(proj):
    print("\nAPIs")
    rc, out = sh(f"gcloud services list --enabled --project={proj} "
                 "--format='value(config.name)' --limit=500")
    enabled = set(out.split())
    if rc != 0 or not enabled:
        check("can list enabled services", False, "cannot verify APIs")
        return
    for api in REQUIRED_APIS:
        check(f"api {api}", api in enabled,
              "" if api in enabled else f"enable with: gcloud services enable {api}")


def section_quota(proj):
    print("\nquota")
    rc, out = sh(f"gcloud compute project-info describe --project={proj} --format=json")
    if rc != 0:
        check("can read quota", False, out.strip()[:120])
        return
    try:
        q = {x["metric"]: x for x in json.loads(out).get("quotas", [])}
    except Exception:
        check("can read quota", False, "unparseable")
        return
    cpus = q.get("CPUS_ALL_REGIONS", {}).get("limit", 0)
    check("CPUS_ALL_REGIONS sufficient", cpus >= MIN_CPUS_ALL_REGIONS,
          f"limit={cpus:g}, need >= {MIN_CPUS_ALL_REGIONS} to run 3 workers concurrently")

    # Informational ONLY. GPU quota is expected to be zero on a new project and is not a
    # blocker: the GPU stage runs on Modal precisely because of this. Reported so nobody
    # spends a day trying to fix it, and NOT filtered to limit>0 -- filtering that way is
    # how the original survey hid GPUS_ALL_REGIONS=0 from its own output, twice.
    gpus = q.get("GPUS_ALL_REGIONS", {}).get("limit", 0)
    check("GPUS_ALL_REGIONS (informational)", True,
          f"limit={gpus:g}" + (" -- expected on a new project; stage 8 uses Modal"
                               if gpus == 0 else ""), fatal=False)


def section_buckets(derived, raw):
    print("\nbuckets")
    for b in (derived, raw):
        rc, out = sh(f"gcloud storage buckets describe gs://{b} --format='value(name)' 2>&1")
        if rc == 0:
            check(f"gs://{b}", True, "exists and is readable by us")
        elif "does not have storage.buckets.get" in out or "403" in out:
            # Globally unique namespace: a 403 means SOMEBODY ELSE owns this name.
            check(f"gs://{b}", False,
                  "name is taken by another project. Bucket names are global; pick a "
                  "different project id or set DERIVED_BUCKET/RAW_BUCKET explicitly")
        else:
            check(f"gs://{b}", True, "does not exist yet; terraform will create it",
                  fatal=False)


def section_tools():
    print("\nlocal tools")
    for tool, why in [("gcloud", "GCP control plane"), ("gsutil", "bucket listing"),
                      ("terraform", "stage 1"), ("docker", "not needed; Cloud Build is used")]:
        rc, out = sh(f"command -v {tool}")
        check(f"{tool} on PATH", rc == 0, why, fatal=tool in {"gcloud", "terraform"})


def section_modal():
    print("\nModal")
    rc, out = sh("command -v modal")
    if not check("modal CLI on PATH", rc == 0, "pip install modal"):
        return
    rc, out = sh("modal profile current 2>&1")
    check("modal authenticated", rc == 0 and "not" not in out.lower(), out.strip()[:80])
    rc, out = sh("modal secret list 2>&1")
    check("secret 'rbp-gcp' exists", "rbp-gcp" in out,
          "create it from the rbp-modal service account key: "
          "modal secret create rbp-gcp SERVICE_ACCOUNT_JSON=\"$(cat key.json)\"")
    # Modal does not expose a balance via CLI. The stage that costs real money is stage 8
    # (~$31), so this is a deliberate human gate rather than a silent assumption.
    check("Modal credit confirmed by a human", False,
          "CANNOT be checked automatically. Open modal.com and confirm >= $35 available "
          "BEFORE stage 8. Re-run with --modal-credit-ok once verified.", fatal=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--skip-modal", action="store_true")
    p.add_argument("--modal-credit-ok", action="store_true",
                   help="you have checked the Modal dashboard yourself")
    a = p.parse_args()

    print("=" * 72)
    print("PREFLIGHT -- nothing is spent until every gate passes")
    print("=" * 72)

    cfg = section_config()
    if cfg is None:
        print("\nSTOP: nothing else can be checked without a project.")
        return 2
    proj, derived, raw = cfg

    section_tools()
    if section_auth(proj):
        section_billing(proj)
        section_apis(proj)
        section_quota(proj)
        section_buckets(derived, raw)
    if not a.skip_modal:
        section_modal()
        if a.modal_credit_ok:
            results[:] = [(("PASS" if "credit confirmed" in n else s), n, d)
                          for s, n, d in results]

    fails = [(n, d) for s, n, d in results if s == "FAIL"]
    warns = [(n, d) for s, n, d in results if s == "WARN"]
    print("\n" + "=" * 72)
    print(f"{sum(1 for s, _, _ in results if s == 'PASS')} passed, "
          f"{len(warns)} warnings, {len(fails)} failures")
    for n, d in fails:
        print(f"  FAIL  {n}: {d}")
    print("=" * 72)
    if fails:
        print("STOP. Fix the failures above. No stage should run.")
        return 1
    print("Preflight clear. Safe to run stage 1.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
