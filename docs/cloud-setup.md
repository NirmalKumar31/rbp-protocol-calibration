# 53. Cloud from zero, taught three ways

Every concept here is explained **three times**, deliberately:

- **CLI / Terraform** — what to type, and what each flag does
- **UI** — the same thing clicking through the console, step by step
- **BACKEND** — what Google actually does when you do it, because that is what lets you debug

Assume no prior GCP knowledge. Nothing is skipped. If a term appears, it gets defined.

---

# Chapter 1. The mental model

Before any commands, four ideas. Everything else is detail.

## 1.1 A project is a billing and permission boundary

A GCP **project** is not a folder. It is the unit that:

- has its own set of enabled APIs (about 200 services, **all off by default**)
- has its own quotas
- is linked to exactly one billing account
- owns its own IAM policy

Two projects on the same billing account share the **credit balance** but nothing else. This
matters: `rbp-repro-2026` draws from the same ~$291 of free trial credit as
`rbp-composition-2026`, but has entirely separate quotas, buckets and permissions.

**Why "all APIs off by default" matters.** A brand-new project cannot create a VM, cannot
store an object, cannot run a container. Every one of those is an API you must enable first,
and the error when you forget is *not* "enable the API" — it is a job that submits fine and
dies two minutes later. This is the single most common wasted hour for a beginner.

## 1.2 There are three planes, and they fail differently

| plane | what it is | how it fails |
|---|---|---|
| **control plane** | the API you talk to: "create a job", "make a bucket" | fast, loud, returns a specific error code |
| **data plane** | the machines actually doing work | slow, quiet, fails hours in |
| **billing plane** | what you are charged, and budgets | **lags by hours**, which is why budgets are not a safety net |

Nearly every confusing GCP experience is a control-plane success followed by a data-plane
failure. `gcloud batch jobs submit` returning `successfully submitted` means the *control
plane* accepted your JSON. It says nothing about whether the work will run.

## 1.3 Global versus regional versus zonal

Three scopes, and mixing them up produces errors that read like nonsense.

- **global**: buckets (names are unique across all of Google Cloud, not per project),
  networks, IAM
- **regional**: subnets, some quotas, Batch job locations (`us-central1`)
- **zonal**: actual VMs (`us-central1-c`), GPU availability

**A concrete trap.** Quota exists at two scopes simultaneously:

```
CPUS_ALL_REGIONS   limit 12    <- global, and this is the binding one
CPUS (us-central1) limit 32    <- regional, irrelevant while the global cap is lower
```

You can have plenty of regional quota and still be unable to start a VM. The binding
constraint is whichever is lower, and beginners read the regional number because it is the
one the console shows first.

## 1.4 Quota is not availability

Quota is *permission* to allocate. Availability is whether the hardware exists in that zone
right now. They are independent:

- quota 8 GPUs, zero V100s in the region → you get nothing, and the error is about capacity
- quota 0 → you get nothing, and the error is about quota

Checking one tells you nothing about the other. `gcloud compute accelerator-types list`
answers availability; `gcloud compute project-info describe` answers quota.

---

# Chapter 2. Creating the project

## 2.1 CLI

```bash
gcloud projects create rbp-repro-2026 --name="RBP reproducible rebuild"
```

- `rbp-repro-2026` is the **project ID**: globally unique, immutable, 6–30 chars, lowercase
  letters, digits and hyphens. This is what every API call uses.
- `--name` is the display name. Mutable, cosmetic.

Then link billing, because a project without billing can do almost nothing:

```bash
gcloud billing accounts list
# ACCOUNT_ID            NAME         OPEN
# XXXXXX-XXXXXX-XXXXXX  RBP project  True

gcloud billing projects link rbp-repro-2026 --billing-account=XXXXXX-XXXXXX-XXXXXX
```

Verify, and verify by reading the field rather than trusting the absence of an error:

```bash
gcloud billing projects describe rbp-repro-2026 --format="value(billingEnabled)"
# True
```

## 2.2 UI, step by step

1. Go to `console.cloud.google.com`.
2. Click the **project picker** in the top blue bar (it shows the current project name).
3. Click **NEW PROJECT**, top right of the dialog.
4. **Project name**: `RBP reproducible rebuild`. Note the grey text underneath showing the
   generated **Project ID** — click **EDIT** and set it to `rbp-repro-2026`. Do this now; the
   ID is immutable.
5. Leave **Location** as *No organisation* for a personal account.
6. **CREATE**. Watch the bell icon top-right for completion, ~30 seconds.
7. Now link billing: hamburger menu (☰) → **Billing**. If it says *This project has no
   billing account*, click **LINK A BILLING ACCOUNT** → choose your account → **SET ACCOUNT**.

**A UI advantage worth knowing.** The console's project picker *always* shows which project
you are in. The CLI keeps that in a config file you never look at, which produces the failure
in §3.3. This is one of the few places the UI is genuinely safer.

## 2.3 Backend

`gcloud projects create` calls `cloudresourcemanager.projects.create`. That returns a
long-running **Operation** immediately; the project does not exist yet. `gcloud` then polls
the operation until done, which is why you see `Waiting for [operations/create_project...]`.

Google allocates a **project number** (`PROJECT_NUMBER` here) alongside your chosen ID. The
number is what internal systems actually use — you will see it in service account emails,
log resource names and error messages, and it is *not* interchangeable with the ID in those
places.

Linking billing writes a `ProjectBillingInfo` record. Until that exists, most APIs reject
calls with `FAILED_PRECONDITION: Billing must be enabled`.

---

# Chapter 3. Authentication, and the trap that looks like an empty result

## 3.1 CLI

```bash
gcloud auth login                 # your human identity, opens a browser
gcloud auth application-default login   # what SDKs use, a separate credential
gcloud config set project rbp-repro-2026
```

Those are **three different things** and confusing them is common:

| | what it is | who uses it |
|---|---|---|
| `auth login` | your user credential for the `gcloud` CLI | `gcloud` commands |
| `auth application-default login` | Application Default Credentials (ADC) on disk | Python `google.cloud.*` libraries |
| `config set project` | ambient default project for the CLI | `gcloud` commands that omit `--project` |

## 3.2 UI

Authentication is implicit — you are logged into the browser. To get the CLI-equivalent, use
**Cloud Shell** (the `>_` icon top right): a container in your browser, already authenticated,
with `gcloud`, `gsutil` and `terraform` installed. For learning, Cloud Shell removes an entire
class of setup problem.

## 3.3 The trap: emptiness is not absence

This cost real time on this project. Watch:

```bash
$ gcloud batch jobs list --location=us-central1
$                       # <- nothing. No error. Exit code 0.
```

There were fourteen jobs. The command returned nothing because **no project was configured**,
and `gcloud <noun> list` with no project resolves to an empty collection rather than an error.
An empty collection is a *valid response to a well-formed query*, so nothing is wrong from the
API's point of view.

**A missing scope is indistinguishable from a missing result.**

The defence, used in `scripts/preflight.py`, is to assert the project **by name**, which does
fail loudly:

```python
rc, out = sh(f"gcloud projects describe {proj} --format='value(projectId)'")
check("project exists and is visible", rc == 0 and proj in out, ...)
```

**BACKEND:** `list` calls take a `parent` field, `projects/{p}/locations/{l}`. With no project
the client library sends a parent it cannot fill and the server returns an empty page. The
emptiness is generated **at the client**, before the request is even meaningful.

---

# Chapter 4. Enabling APIs

## 4.1 CLI

```bash
for api in \
  cloudresourcemanager.googleapis.com \
  compute.googleapis.com \
  storage.googleapis.com \
  batch.googleapis.com \
  artifactregistry.googleapis.com \
  cloudbuild.googleapis.com \
  billingbudgets.googleapis.com
do
  gcloud services enable "$api" --project=rbp-repro-2026
done
```

What each one is for, since a list of domain names teaches nothing:

| API | why this project needs it |
|---|---|
| `cloudresourcemanager` | read and modify the project itself; Terraform needs it constantly |
| `compute` | VMs, networks, subnets, **and the quota system** |
| `storage` | GCS buckets and objects |
| `batch` | the managed batch-job service that runs our fan-out |
| `artifactregistry` | private container image registry |
| `cloudbuild` | builds those images, remotely |
| `billingbudgets` | budgets, which the killswitch listens to |

## 4.2 UI

1. ☰ → **APIs & Services** → **Library**.
2. Search the API name, click the card, click **ENABLE**.
3. Repeat. There is no bulk enable in the Library UI, which is why the CLI loop is worth
   learning.

To audit: **APIs & Services** → **Enabled APIs & services**.

## 4.3 Backend

Enabling an API flips a per-project flag *and* often creates a **service agent** — a
Google-managed service account like
`service-PROJECT_NUMBER@gcp-sa-batch.iam.gserviceaccount.com` — that the service uses to act on
your behalf. This is why enabling can take a minute and why it is a long-running operation.

**The failure mode when you forget.** The control plane accepts your request (the API surface
exists) and the data plane cannot act. You get a job that submits successfully and fails
minutes later with a permission error mentioning a service account you have never heard of.
That is the signature of a missing API or a missing service agent.

---

# Chapter 5. Storage: buckets and objects

## 5.1 The model

GCS has **no directories**. `processed/dinuc/K562/QKI/dataset.tsv` is one flat object name
containing slashes. Tools display it as a tree; the store has no such concept.

Three consequences that matter daily:

1. **A wrong path is indistinguishable from missing data.** There is no schema to violate.
   `gsutil ls gs://bucket/runz/` returns nothing, exactly as `runs/` would if it were empty.
   *Always list the parent prefix before concluding data is absent.* This cost time on this
   project: I looked for `runs/gc/` when the rehearsal writes to `rehearsal/gc/`, and read
   the emptiness as a failed job.
2. **Listing is a prefix scan.** `ls gs://b/a/b/c/` is a filtered scan, not a directory
   lookup, and costs proportionally to matches.
3. **Names are globally unique.** Not per project — global. `rbp-composition-2026-derived`
   can exist exactly once on Earth, which is why reproducers cannot reuse it and why the
   convention is `{project_id}-derived`.

## 5.2 CLI

```bash
gcloud storage buckets create gs://rbp-repro-2026-derived \
  --project=rbp-repro-2026 --location=us-central1 --uniform-bucket-level-access

gsutil ls gs://rbp-repro-2026-derived/                 # top-level prefixes
gsutil ls "gs://rbp-repro-2026-derived/processed/**"   # recursive
gsutil du -sh gs://rbp-repro-2026-raw                  # total size
gsutil cat gs://.../manifest/study_panel.tsv | head    # read without downloading
gsutil stat gs://.../results/tables/x.csv              # metadata only
```

`--uniform-bucket-level-access` turns off per-object ACLs so IAM is the only permission
system. **Required** for IAM conditions (Chapter 7), and simply better: two permission
systems that can disagree is a bug generator.

## 5.3 UI

1. ☰ → **Cloud Storage** → **Buckets** → **CREATE**.
2. **Name**: globally unique. It will tell you if taken.
3. **Location type**: *Region* → `us-central1`. Region matters: co-locating buckets with
   compute makes transfer free and fast. Cross-region egress is billed.
4. **Storage class**: Standard.
5. **Access control**: **Uniform**. Not Fine-grained.
6. **CREATE**.

Browsing: click the bucket, and the console *renders prefixes as folders*. Remember that is a
UI fiction.

## 5.4 Backend

An object write is a single atomic operation — an object either fully exists or does not.
There are no partial objects and **no cross-object transactions**.

That last point drives a pattern used everywhere in this pipeline. A task writes several
objects and cannot make them atomic together, so **write the cheap summary LAST**:

```python
# 1. the payload (large, slow)
bucket.blob(f"{prefix}/scores.tsv.gz").upload_from_string(...)
# 2. the marker (tiny, fast) -- written LAST, on purpose
bucket.blob(f"{prefix}/metrics.json").upload_from_string(...)
```

Now consider a crash between the two. The marker is absent, so the next run redoes the task.
Costly, and correct. Reverse the order and a crash leaves a marker with no payload — the task
is skipped forever and the gap is permanent and silent. **Order the writes so the survivable
failure is the one that happens.**

---

# Chapter 6. Identity: service accounts and IAM

## 6.1 The model

A **service account** is a non-human identity with an email address. Code runs *as* one.

IAM is `(who, what, where)`:

- **who** = a member: `user:you@gmail.com`, `serviceAccount:x@proj.iam.gserviceaccount.com`
- **what** = a role: a named bundle of permissions, e.g. `roles/storage.objectViewer`
- **where** = the resource the binding is attached to: project, bucket, single object

Roles are bundles, not single permissions. `roles/storage.objectAdmin` includes create,
delete, get and list. There are **predefined** roles (Google's), **basic** roles (Owner,
Editor, Viewer — far too broad, avoid) and **custom** roles.

## 6.2 This project's five identities, and why five

| service account | job | what it can do |
|---|---|---|
| `rbp-ingest` | download from ENCODE, GENCODE, NCBI, UCSC | write `raw/` |
| `rbp-prep` | preprocessing fan-out | read `raw/`, write `processed/` |
| `rbp-train` | rehearsal and training sweeps | read `processed/`, write `runs/` |
| `rbp-analysis` | aggregation and figures | read most, write `results/` |
| `rbp-modal` | the **only credential that leaves Google's network** | write `runs/`, `ckpt/`, `variants/` only |

**Why not one account?** Because the point of separation is that a preprocessing task
*cannot* write a model, and an ingest task *cannot* touch results. A single identity discards
that for no benefit. This project shipped a bug where `submit.sh` ran every job as
`rbp-train`, which silently threw the whole scheme away — see the chronicle.

## 6.3 CLI

```bash
gcloud iam service-accounts create rbp-train \
  --project=rbp-repro-2026 --display-name="Training jobs"

# grant on a BUCKET, not the project: the narrowest scope that works
gcloud storage buckets add-iam-policy-binding gs://rbp-repro-2026-derived \
  --member="serviceAccount:rbp-train@rbp-repro-2026.iam.gserviceaccount.com" \
  --role="roles/storage.objectViewer"

# audit
gcloud storage buckets get-iam-policy gs://rbp-repro-2026-derived --format=json
```

## 6.4 UI

**Creating one:** ☰ → **IAM & Admin** → **Service Accounts** → **CREATE SERVICE ACCOUNT** →
name → **CREATE AND CONTINUE** → optionally grant project roles → **DONE**.

**Granting on a bucket (the better habit):** **Cloud Storage** → click bucket →
**PERMISSIONS** tab → **GRANT ACCESS** → paste the service account email → pick a role →
**SAVE**.

**Note the UI nudge.** The service-account creation flow offers *project-level* roles, which
is the broad option. Granting on the bucket is narrower and better, and the UI does not lead
you there. Knowing that is the difference between following a wizard and understanding the
model.

## 6.5 Backend

Every API call carries a token. The server resolves the caller, gathers every IAM policy on
the resource **and its ancestors** (object → bucket → project → folder → org), and checks
whether any binding grants the required permission. Deny is the default; there is no implicit
allow.

**Policies are eventually consistent.** A fresh binding can take seconds to tens of seconds
to take effect. A 403 immediately after granting is often just propagation, and retrying is
the right move — which is genuinely confusing the first time.

---

# Chapter 7. IAM conditions: the sharpest tool here

## 7.1 The problem

`rbp-modal`'s key **leaves Google's network** — it sits in a Modal secret on a third-party
platform. If that key leaks, what can the holder do? With plain `roles/storage.objectAdmin`
on the bucket: delete every result, overwrite every dataset, destroy the study.

## 7.2 CLI / Terraform

A **condition** is a CEL expression attached to a binding. The binding applies only when it
evaluates true.

```hcl
resource "google_storage_bucket_iam_member" "modal_derived_write" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.modal.email}"

  condition {
    title      = "runs-checkpoints-and-variants-only"
    expression = <<-EOT
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/runs/") ||
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/ckpt/") ||
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/variants/")
    EOT
  }
}
```

Fully compromised, that key cannot alter a dataset under `processed/`, cannot touch the raw
bucket, and cannot delete anything outside those three prefixes.

**This worked in practice, and I have the receipt.** The first Modal ClinVar probe scored its
dataset correctly and then died on upload:

```
rbp-modal@... does not have storage.objects.create access to ...
  /objects/variants/scores_sb/K562_AATF.csv
```

That is the condition doing exactly its job: a **new write path is denied until somebody
widens it deliberately.** I widened it in Terraform, not with a console click, because a
hand-made binding is invisible to Terraform and gets reverted by the next apply.

## 7.3 UI

**Cloud Storage** → bucket → **PERMISSIONS** → **GRANT ACCESS** → member and role → then
**ADD IAM CONDITION**. In the dialog: give it a title, switch to the **CONDITION EDITOR** tab
and paste the CEL. The visual builder cannot express prefix-OR chains, so the editor is the
only workable path.

**Requires uniform bucket-level access.** If the bucket uses fine-grained ACLs, conditions
are unavailable and the option is greyed out with an unhelpful tooltip.

## 7.4 Backend

The condition is evaluated **per request**, at authorisation time, against a context object
containing `resource.name`, `resource.type`, `request.time` and more. A condition that
references something not in the context for a given call type evaluates false, and the
binding silently does not apply. That is why conditions are best kept to `resource.name`
prefix tests: they are the reliably-present field.

---

# Chapter 8. Networking, and why it is deliberately not uniform

## 8.1 The design

```
rbp-net  (custom VPC, no auto subnets)
├── rbp-workers          us-central1   Private Google Access ON
│     Batch workers: prep, rehearsal, sweep, analysis
│     NO external IP -> can reach *.googleapis.com and NOTHING else
└── rbp-gpu-us-central1  us-central1   (present, unused: GPU quota is 0)

default VPC
      ingest, panel, variants -> EXTERNAL IP
      because these must reach ENCODE, GENCODE, NCBI and UCSC
```

## 8.2 Private Google Access, precisely

A VM with **no external IP** normally cannot send a packet off the machine. PGA is a subnet
flag that routes traffic destined for Google's public API ranges over Google's internal
network anyway.

**What it does not do:** reach anything that is not a Google API. `huggingface.co` is not a
Google API. `hgdownload.soe.ucsc.edu` is not a Google API.

This is exactly why the model weights are **baked into the container image**. A worker cannot
`from_pretrained("multimolecule/splicebert")` — the DNS lookup fails, or worse, the library
retries until `maxRunDuration`. Baking makes the weights part of the artefact the digest
identifies, so "which weights produced this result?" and "which image produced this result?"
have the same answer.

## 8.3 Why no Cloud NAT

NAT would give every worker internet access. Three reasons not to:

1. It is billed per gateway-hour **plus** per GB processed. For an occasional download that
   costs more than the download.
2. It hands general internet access to 488 workers that have no reason to have it.
3. It makes every run depend on third-party sites being up.

Instead the three stages that genuinely need the internet run on a **single short-lived VM
with a public IP**. `cloud/submit.sh` encodes this:

```bash
if [ "$EXTERNAL" = "1" ]; then
  NETWORK=""          # default network, gets an external IP
else
  NETWORK='"network": {"networkInterfaces": [{"network": "projects/'"${PROJECT}"'/global/networks/rbp-net", "subnetwork": "...rbp-workers", "noExternalIpAddress": true}]},'
fi
```

## 8.4 CLI

```bash
gcloud compute networks create rbp-net --subnet-mode=custom --project=rbp-repro-2026

gcloud compute networks subnets create rbp-workers \
  --network=rbp-net --region=us-central1 --range=10.0.0.0/20 \
  --enable-private-ip-google-access --project=rbp-repro-2026
```

`--subnet-mode=custom` matters: `auto` creates a subnet in **every** region, which is a lot
of surface you did not ask for.

## 8.5 UI

1. ☰ → **VPC network** → **VPC networks** → **CREATE VPC NETWORK**.
2. Name `rbp-net`. **Subnet creation mode: Custom**.
3. Add a subnet: name `rbp-workers`, region `us-central1`, IPv4 range `10.0.0.0/20`, and
   **Private Google Access: On**. That toggle is the whole point of the page.
4. **CREATE**.

## 8.6 Backend

PGA works by installing routes for Google's API IP ranges pointing at an internal next hop.
Traffic never traverses the public internet. DNS for `*.googleapis.com` resolves to those
ranges (`199.36.153.x` and similar) from inside the VPC.

**The debugging signature:** a worker that hangs on an outbound connection to a non-Google
host, rather than failing fast. There is no route, so the SYN goes nowhere and you wait for a
TCP timeout. If a task times out doing something that should be quick, suspect the network
before suspecting the code.

---

# Chapter 9. Terraform: infrastructure as code

## 9.1 Why, in one sentence

Because "what does the infrastructure look like?" should be answerable by reading a file, not
by clicking through twelve console pages and hoping you remember what you changed.

## 9.2 The core loop

```bash
terraform init      # download providers, configure the state backend
terraform plan      # compute the diff between config and reality. READ THIS.
terraform apply     # execute the diff
```

**State** is the crux. Terraform records what it believes it created in a state file. `plan`
diffs three things: your config, the state, and (with a refresh) reality.

## 9.3 The mistake that nearly destroyed the original study

`main.tf` had:

```hcl
backend "gcs" {
  bucket = "rbp-composition-2026-tfstate"     # HARDCODED
  prefix = "terraform/state"
}
```

**Terraform backends cannot use variables.** So every checkout of this repo shared one state
file. I pointed `config/params.yaml` at the new project, ran `terraform plan`, and got:

```
Plan: 63 to add, 1 to change, 63 to destroy.
```

The 63 destroys were **the original study's buckets, containing every result**. Terraform was
behaving perfectly: the state said those resources existed and the config now described a
different project, so the correct diff is "destroy those, create these". And `run.sh` ran
`terraform apply -auto-approve`.

**Three fixes, and each is worth internalising.**

**(a) Partial backend configuration.** Remove the bucket from the file; supply it at init:

```hcl
backend "gcs" {
  prefix = "terraform/state"      # no bucket
}
```

```bash
terraform init -reconfigure -backend-config="bucket=${PROJECT_ID}-tfstate"
```

`-reconfigure` is load-bearing: without it, init reuses whatever backend a previous checkout
cached in `.terraform/`, which is precisely how another project's state leaks in.

**(b) A destroy guard.** A first apply on an empty project is additive **by definition**. Any
destroy means the state does not describe this project:

```bash
terraform plan -out=tfplan.new
DESTROYS=$(terraform show -no-color tfplan.new | grep -cE "^  # .* will be destroyed" || true)
if [ "${DESTROYS:-0}" -gt 0 ]; then
  die "plan contains ${DESTROYS} DESTROY actions. On a fresh project it must contain none."
fi
terraform apply tfplan.new
```

Note `apply tfplan.new` applies the **saved plan**, so what you inspected is exactly what
runs. `apply -auto-approve` without a saved plan re-plans and can do something different.

**(c) The state bucket is created outside Terraform**, by `run.sh`, before init. Storing
state *about* the bucket that stores the state is a bootstrap paradox.

## 9.4 The `for_each` trap

```hcl
# BROKEN
resource "google_artifact_registry_repository_iam_member" "pullers" {
  for_each = toset([
    google_service_account.prep.email,      # "known only after apply"
    ...
  ])
}
```

`for_each` **keys must be known at plan time**. A resource attribute of something not yet
created is not. Terraform refuses to plan *at all*, which also blocks `terraform import` of
anything else in the configuration — so you cannot even adopt existing resources to get
unstuck.

The fix is to key on something known statically and build the derived value:

```hcl
for_each = toset(["rbp-prep", "rbp-ingest", "rbp-train", "rbp-analysis"])
depends_on = [google_service_account.prep, ...]
member = "serviceAccount:${each.value}@${var.project_id}.iam.gserviceaccount.com"
```

`var.project_id` is known at plan time. `depends_on` preserves ordering. Identical grants, and
the configuration is now plannable from empty. **`for_each` over computed values is an
anti-pattern, not a style preference.**

## 9.5 Importing what already exists

I created the project and the state bucket by hand, so Terraform's first apply hit
`Error 409: Requested entity already exists`. Adopt them:

```bash
terraform import google_project.rbp rbp-repro-2026
terraform import google_storage_bucket.tfstate rbp-repro-2026-tfstate
```

Import writes the resource into state without changing infrastructure. Note the second import
took state from 2 entries to 78 — Terraform refreshed and discovered the resources a partial
apply had already created.

## 9.6 UI equivalent

There is none, and that is the point. The console has no "apply this file" and no state. To
build this by hand you would create, in order: project, billing link, 7 APIs, 4 buckets, 5
service accounts, ~20 IAM bindings (3 with CEL conditions), 1 VPC, 2 subnets, 1 Artifact
Registry repo, 1 budget, 1 Pub/Sub topic, 1 Cloud Function. Roughly 45 console pages, no
record of what you did, and no way to tear it down reliably.

Do it once by hand to learn the objects. Never do it twice.

---

# Chapter 10. Containers: Artifact Registry and Cloud Build

## 10.1 Why two images

| | size | contents | used by |
|---|---|---|---|
| `cpu` | ~1.2 GB | numpy, pandas, scikit-learn, `google-cloud-storage`. **No torch.** | ingest, panel, prep, rehearsal, variants, analysis |
| `gpu` | ~6 GB | the above **plus** torch, transformers, multimolecule, and baked model weights | the CNN sweep |

Preprocessing fans out to 488 tasks. Using one 6 GB image would pull 4.8 GB of unused CUDA
into every worker — more time spent pulling than preprocessing.

**And a trap that follows directly.** The CNN is a torch model. `submit.sh` originally pinned
*every* job to the cpu image, so all 475 sweep tasks would have died on `import torch`. The
fix:

```bash
IMAGE_KIND=cpu
[ "$JOB_TYPE" = "sweep" ] && IMAGE_KIND=gpu
```

The gpu image runs fine on a CPU machine — it just carries CUDA it will not touch, which is
the cheaper mistake.

## 10.2 Layer order is the cost story

Docker caches per instruction, and **everything below a changed layer is rebuilt**:

```dockerfile
COPY docker/requirements-cpu.txt /tmp/          # changes rarely
RUN pip install -r /tmp/requirements-cpu.txt    # slow
COPY src/ /app/src/                             # changes constantly
COPY scripts/ /app/scripts/                     # changes constantly
```

Reversed, every one-character edit to a script would reinstall scipy.

## 10.3 CLI

```bash
gcloud artifacts repositories create rbp \
  --repository-format=docker --location=us-central1 --project=rbp-repro-2026

gcloud builds submit --project=rbp-repro-2026 \
  --config=docker/cloudbuild.cpu.yaml \
  --substitutions="_IMAGE=us-central1-docker.pkg.dev/rbp-repro-2026/rbp/cpu,_ARTIFACTS=rbp-repro-2026-artifacts,_GIT_SHA=$(git rev-parse --short HEAD)" .
```

**Why pass substitutions explicitly?** Because Cloud Build does **not** recursively expand
`$PROJECT_ID` inside a user-defined substitution's *default value*. Written as

```yaml
substitutions:
  _IMAGE: us-central1-docker.pkg.dev/$PROJECT_ID/rbp/cpu
```

the literal string reaches docker, which rejects it:

```
invalid reference format: repository name must be lowercase
```

— because `$PROJECT_ID` contains capitals. Passing them at submit time makes the values
concrete.

## 10.4 The `.gitignore` trap, which is subtle and cost a build

`gcloud builds submit .` uploads the source directory, **honouring `.gitignore`**. My
`.gitignore` had:

```
data/
```

That is **unanchored**, so it matches a directory of that name at *any depth* — including
`src/rbp/data/`. Git stopped tracking the package, the upload excluded it, and the image was
built without it. The image's own test step caught it:

```
ModuleNotFoundError: No module named 'rbp.data'
```

Anchor with a leading slash: `/data/` matches only at the repo root. Verify with
`git check-ignore -v <path>`, which names the offending rule and line.

## 10.5 Test gates: a floor, not an equality

The build runs the test suite inside the image and asserts a count, because a suite that
silently collects fewer tests still exits 0.

It originally asserted **equality**, and failed:

```
collected 548 tests (expected 460)
FAIL: expected 460 tests, collected 548
```

It failed because tests were **added**. An exact count punishes the right behaviour and
trains people to edit the expected number until it goes green. Under-collection is the
failure worth catching, so:

```bash
if [ "$$n" -lt "${_EXPECTED_TESTS}" ]; then ... fi
```

Note `$$n` — in a Cloud Build `bash` step, `$$` escapes a literal `$` so the shell expands it
rather than Cloud Build's substitution engine.

## 10.6 Push before recording the digest

Cloud Build's `images:` block pushes only after **every** step finishes. An earlier version
recorded the digest in a step that ran before the push, and printed an empty digest. The build
now pushes explicitly, then records:

```bash
gcloud storage cp <(echo "$DIGEST") gs://${_ARTIFACTS}/images/cpu_digest.txt
```

## 10.7 Pin by digest, never by tag

```bash
DIGEST=$(gcloud storage cat "gs://${PROJECT}-artifacts/images/${IMAGE_KIND}_digest.txt")
IMAGE="${REGION}-docker.pkg.dev/${PROJECT}/rbp/${IMAGE_KIND}@${DIGEST}"
```

A tag is mutable — `:latest` today is not `:latest` tomorrow. A digest is content-addressed.
"Which image produced this result?" is only answerable by digest.

## 10.8 UI

**Registry:** ☰ → **Artifact Registry** → **CREATE REPOSITORY** → format Docker, region
`us-central1`.

**Builds:** ☰ → **Cloud Build** → **History**. Click a build to see every step, its duration
and its log. This is genuinely better than the CLI for diagnosis: failed steps are
colour-coded and you can jump straight to the failing step's output instead of scrolling.

**Images:** Artifact Registry → repository → click an image to see tags, digest and size.

---

# Chapter 11. Cloud Batch

## 11.1 The model

You submit a **job** describing N identical **tasks**. Batch provisions VMs in a Managed
Instance Group, runs an agent on each, and hands out tasks. Each task gets
`BATCH_TASK_INDEX`, which is how it knows which slice of work is its own.

```
job (taskCount 488, parallelism 8, taskCountPerNode 4)
 └── MIG: 2 x e2-standard-4          <- 8 vCPU total
       ├── agent, 4 tasks: BATCH_TASK_INDEX = 149, 167, 32, 78
       └── agent, 4 tasks: BATCH_TASK_INDEX = 22, 30, 394, 418
```

## 11.2 The index-space rule that shapes the whole pipeline

**Batch partitions the index space across nodes and does NOT honour manifest order.** Those
indices above are not contiguous.

Consequence: **you cannot express "run the CNN first" by ordering the manifest.** Scope is
expressed by *which manifest*:

```python
MANIFEST_TAG = os.environ.get("MANIFEST_TAG", "")
MANIFEST = f"manifest/sweep_tasks{MANIFEST_TAG}.tsv"
```

One manifest per model. `MANIFEST_TAG=_cnn` selects it. Ordering within a manifest is only
useful for *scheduling* — longest-first, so the job does not end when the unluckiest node
finishes.

## 11.3 The job spec, annotated

```json
{
  "taskGroups": [{
    "taskCount": 488,            // from the manifest. NEVER a literal.
    "parallelism": 8,            // how many tasks may run at once
    "taskCountPerNode": 4,       // tasks packed per VM
    "taskSpec": {
      "maxRetryCount": 2,        // transient GCS errors and preemption both happen
      "maxRunDuration": "7200s", // a hung task must die, not bill forever
      "computeResource": {"cpuMilli": 900, "memoryMib": 3500},
      "runnables": [{
        "container": {
          "imageUri": "us-central1-docker.pkg.dev/PROJ/rbp/cpu@sha256:...",
          "entrypoint": "python",
          "commands": ["scripts/cloud_prep.py", "prep"]
        },
        "environment": {"variables": {"OMP_NUM_THREADS": "1", ...}}
      }]
    }
  }],
  "allocationPolicy": {
    "instances": [{"policy": {"machineType": "e2-standard-4",
                              "provisioningModel": "STANDARD",
                              "bootDisk": {"sizeGb": 100, "type": "pd-balanced"}}}],
    "network": {"networkInterfaces": [{"network": "projects/PROJ/global/networks/rbp-net",
                                       "subnetwork": "...rbp-workers",
                                       "noExternalIpAddress": true}]},
    "serviceAccount": {"email": "rbp-prep@PROJ.iam.gserviceaccount.com"}
  },
  "logsPolicy": {"destination": "CLOUD_LOGGING"}
}
```

**Two schema traps I hit:**

1. `networkInterfaces` goes under `allocationPolicy.network.networkInterfaces`, **not**
   `allocationPolicy.networkInterfaces`. The latter is rejected with
   `Unknown name "networkInterfaces" at 'job.allocation_policy'`.
2. `cpuMilli: 900` with a 4000-milli machine gives 4 tasks per node. It also **caps each
   task at 0.9 CPU**, which is why `OMP_NUM_THREADS=1` matters — see §11.5.

## 11.4 Never type a task count

```bash
manifest_rows() {
  n=$(gcloud storage cat "gs://${DERIVED}/${key}" | tail -n +2 | wc -l | tr -d ' ')
  [ -n "$n" ] && [ "$n" -gt 0 ] || { echo "EMPTY_MANIFEST"; return 1; }
  echo "$n"
}
```

The old script hardcoded `COUNT=189`. The gc arm has 187 datasets, so the job dispatched two
tasks past the end of its manifest and Batch reported the whole job **FAILED** — with every
real task having succeeded. That is an exceptionally confusing failure: complete output, red
status.

## 11.5 The thread-count trap, measured

numpy and scipy size their thread pools from the cores they can **see**, which is the host's
count, not the container's cgroup limit. Four tasks share a 4-core node at 0.9 CPU each, so
the default is four processes each spawning four threads: sixteen threads fighting over four
cores.

Measured on this project: NCBP2 is a 406-pair dataset that computes in **7 seconds** on a
laptop. Unpinned in a container it ran **over 20 minutes without finishing**. Pinned to one
thread: **26.7 seconds**.

```dockerfile
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
```

A 189-task sweep would have looked like a hang rather than an error. It also makes numerics
reproducible: thread count changes summation order inside BLAS, which moves iterative solvers
in the last few decimals.

## 11.6 Quota: the ceiling you cannot touch

```
parallelism 8 ÷ taskCountPerNode 4 = 2 nodes × 4 vCPU = 8 vCPU
CPUS_ALL_REGIONS = 12
```

Why 8 and not 12? Because **VM creation fails AT the limit, not approaching it.** Asking for
exactly 12 gets:

```
OPERATIONAL_INFO: CODE_GCE_QUOTA_EXCEEDED
```

Batch then retries the third node forever, the job runs 8-wide anyway, and the events log
fills with an alarming error that changes nothing.

**I misdiagnosed this.** I saw vCPU usage drop 8 → 4, called it spot preemption, switched to
on-demand, and measured no improvement. The job's events contained the quota error and **zero
preemption events**. The dip was a node cycling out after finishing its 4 tasks with its
replacement refused. Lesson: read the evidence before choosing the cause, especially when you
already have a story that fits.

## 11.7 Submit does not mean run

```bash
gcloud batch jobs submit "$JOB" --config="$SPEC"
# Job ... was successfully submitted.       <- CONTROL PLANE ONLY
```

That is acceptance of your JSON. It says nothing about the work. Stage 5 submitted 488 tasks
and its very next line was `cloud_prep.py finalize`, which would have run seconds later
against an empty bucket and written **a panel of zero datasets that every downstream stage
would have trusted.**

So `submit.sh` polls:

```bash
while :; do
  STATE=$(gcloud batch jobs describe "$JOB" --format="value(status.state)")
  case "$STATE" in
    SUCCEEDED) exit 0 ;;
    FAILED|CANCELLED) exit 2 ;;      # not 1: a failed job is not always a failed run
    *) echo "$STATE $COUNTS" ;;
  esac
  sleep 60
done
```

Exit **2** rather than 1 for FAILED, because a failed job is not necessarily a failed run:
preemption leaves a scatter of resumable failures, and a task count past the end of a manifest
once failed a job whose every task succeeded. The caller decides.

## 11.8 UI

☰ → **Batch** → **Job list**. Click a job for:

- **Details**: the spec as submitted
- **Task groups**: per-state counts, refreshing live
- **Events**: `STATUS_CHANGED` and `OPERATIONAL_INFO`. **Read these first on any failure** —
  this is where `CODE_GCE_QUOTA_EXCEEDED` appeared
- **Logs**: straight into Cloud Logging, filtered to the job

**The UI is better than the CLI here.** Per-state counts refresh without re-running a command,
and the Events tab surfaces `OPERATIONAL_INFO` that is easy to miss in a `--format` string.

## 11.9 Backend

`batch.jobs.create` validates and persists the spec, then a Batch-owned controller creates a
MIG. Each VM boots a COS image with a Batch agent, which registers, polls for task
assignments, runs the container via `docker run`, and reports status. `logsPolicy:
CLOUD_LOGGING` makes the agent forward container stdout/stderr to Cloud Logging labelled with
`job_uid`, which is how you query one job's output:

```bash
gcloud logging read 'labels.job_uid="<uid>"' --limit=20 --format="value(textPayload)"
```

**Warning:** that returns agent state dumps too, which are enormous protobufs. Filter to your
own log lines, or you will drown.

---

# Chapter 12. Budgets, and why a budget is not a brake

## 12.1 The lag

Budgets are evaluated against **reported** spend, and reporting lags by hours. A runaway job
can spend a day's budget before a single alert fires. **A budget is a smoke alarm, not a
sprinkler.**

The real brakes are the ones that act at submission time: `parallelism`, `max_containers` on
Modal, `maxRunDuration` per task, and quota itself.

## 12.2 The default that makes budgets useless

```hcl
credit_types_treatment = "EXCLUDE_ALL_CREDITS"
```

The default is `INCLUDE_ALL_CREDITS`, which **subtracts free credit from reported spend**. On
a project with $300 of trial credit, reported spend is **$0** until the credit is gone. Every
threshold is unreachable. The budget exists, the alerts are configured, and nothing can ever
fire.

`scripts/preflight.py` checks this explicitly, because it is invisible otherwise:

```python
bad = [b.get("displayName") for b in budgets
       if (b.get("budgetFilter", {}) or {}).get("creditTypesTreatment") == "INCLUDE_ALL_CREDITS"]
check("budgets exclude credits", not bad, ...)
```

## 12.3 The killswitch

Budget → Pub/Sub topic → Cloud Function → **unlink billing from the project**. Unlinking is
the only true stop: it halts everything at once, and it is reversible with one command.

```
gcloud billing projects link $PROJECT_ID --billing-account=$BILLING_ACCOUNT
```

Verified against a real budget message on the original project: `cost=5.49 ... dry_run=False`.

**Careful:** unlinking billing stops *everything*, including things you wanted. It is a
last-resort brake, which is why it sits at a threshold you do not expect to reach.

## 12.4 UI

☰ → **Billing** → **Budgets & alerts** → **CREATE BUDGET**.

- **Scope**: this project.
- **Amount**: your cap.
- **Actions**: thresholds at 25/50/80/100%.
- **Manage notifications** → **Connect a Pub/Sub topic to this budget**. This is the part
  people miss; email alerts alone cannot trigger automation.
- **Credits**: under Scope there is a checkbox for discounts/credits. **Uncheck it** — that is
  the UI equivalent of `EXCLUDE_ALL_CREDITS`.

---

# Chapter 13. Modal, and why the pipeline needs a second cloud

## 13.1 The wall

```
GCP:   GPUS_ALL_REGIONS 0, increase auto-denied (NOT_ENOUGH_USAGE_HISTORY) for 8, 4, and 1
AWS:   0 on all four GPU families
Azure: GPU quota forbidden on free trial
Modal: no quota gate at all
```

## 13.2 What made the migration cheap

A training task's only cloud dependency is `google-cloud-storage`. It reads its manifest and
dataset from GCS, writes scores, metrics and weights back, and takes its index from an
environment variable. **No Batch API calls, no metadata-server assumptions.** So
`scripts/cloud_train.py` runs on Modal *unchanged*, and `aggregate` cannot tell afterwards
which platform produced a row.

This is a design property worth copying: keep the task's interface to the world as narrow as
possible and it becomes portable for free.

## 13.3 max_containers is a budget control, not a performance knob

On GCP, quota capped the burn rate whether we liked it or not. Modal removes the cap — which
is why it is useful, and which also removes the accidental cost ceiling.

```
concurrency × A10G price = burn rate      $30 credit lasts
 1 × $1.10                = $1.10/h        ~27 h
10 × $1.10                = $11.00/h       ~2.7 h
50 × $1.10                = $55.00/h       ~33 min
```

`MAX_CONTAINERS = 10` is deliberate: fast enough to finish in hours, slow enough that a
mistake costs a few dollars rather than the whole credit before anyone notices. **There is no
Modal equivalent of the billing killswitch, so this cap IS the guardrail.**

## 13.4 Choosing the GPU by measuring, not by price list

- **A10G for training**: measured 1.98× a T4 for 1.42× the price, so it is the cheaper unit of
  work. A100 measured only 2.89× a T4 — a 20M-parameter model does not saturate one — so it
  costs more per unit.
- **T4 for the variant scoring**: that job is **network-bound**, 375 MB of checkpoint per
  dataset against ~3,500 forward passes. Paying A10G rates to wait on a download is paying for
  the wrong thing.
- **T4 for the locality probe**: that job *is* compute-bound — ~570,000 forward passes against
  one 75 MB download — and it went from ~90 minutes of laptop CPU to ~8 minutes.

Also worth recording: `cpu=` and `memory=` on a GPU function set the allocation but are **not
billed on top** of the GPU. An earlier estimate built by summing three published prices came
out **44% high**.

## 13.5 The three Modal failures, in order

**(a) `ModuleNotFoundError: No module named 'modal_sweep'`.** Modal ships **only the
entrypoint file** to the container. `import modal_sweep` resolved on my machine because I had
added the directory to `sys.path`; inside the container the file did not exist. Containers
crash-looped at import, ~5 seconds each, doing no work. Fix: ship the dependency in the image
rather than duplicating the image spec, because two specs drift and variant scores from a
different SpliceBERT build than the binding scores would defeat the purpose.

```python
image = _base.add_local_file(f"{HERE}/modal_sweep.py", "/root/modal_sweep.py")
```

**(b) `ModuleNotFoundError: No module named 'pyfaidx'`.** The cloud task never opens a FASTA —
that is the entire point of precomputing windows — but the import sat at module scope, so it
crashed every container on a dependency it never used. Fix: import inside the two functions
that actually cut windows.

**(c) 403 on upload.** The IAM condition permitted `runs/` and `ckpt/` only. The task scored
its dataset correctly and then could not write `variants/scores_sb/...`. That is the guardrail
working. Widened in Terraform.

**All three were caught by a single-task probe costing about a cent**, before the 94-task
sweep. That is what a probe is for:

```python
@app.local_entrypoint()
def probe(index: int = 0):
    t0 = time.time(); rc = task.remote(index, True); el = time.time() - t0
    print(f"projected {N_TASKS} / {MAX_CONTAINERS}: {el*N_TASKS/MAX_CONTAINERS/60:.1f} min, "
          f"${el*N_TASKS/3600*0.59:.2f} at T4 rates")
```

## 13.6 Setting up a Modal account for this project

```bash
modal token set --token-id ak-... --token-secret as-... --profile=NAME
modal profile activate NAME

# AFTER terraform has created the service account:
gcloud iam service-accounts keys create /tmp/k.json \
  --iam-account=rbp-modal@$GOOGLE_CLOUD_PROJECT.iam.gserviceaccount.com
modal secret delete rbp-gcp 2>/dev/null || true      # Modal will NOT overwrite
modal secret create rbp-gcp SERVICE_ACCOUNT_JSON="$(cat /tmp/k.json)"
rm /tmp/k.json
```

**Two things that bite.** `modal` must be on `PATH`, not merely installed in a virtualenv —
preflight failed on exactly this. And if a `rbp-gcp` secret already exists from a previous
project, **delete it first**: Modal will not overwrite, and the stage would silently write
into the old project's bucket.

**Terraform deliberately contains no `google_service_account_key` resource.** That would put a
live private key in Terraform state in plaintext. Keys are minted by hand, used, and deleted.

---

# Chapter 14. Doing all of this in the UI, end to end

If you had no CLI and no Terraform, here is the whole build. Worth doing once to learn the
objects; never worth doing twice.

| # | console path | what to set |
|---|---|---|
| 1 | project picker → **NEW PROJECT** | name, and **EDIT the Project ID** |
| 2 | ☰ **Billing** → LINK A BILLING ACCOUNT | your account |
| 3 | ☰ **APIs & Services** → Library | enable the 7 APIs from Ch.4, one at a time |
| 4 | ☰ **IAM & Admin** → Quotas | filter `CPUS_ALL_REGIONS`; note the limit. Request an increase here if needed |
| 5 | ☰ **Cloud Storage** → CREATE (×4) | `-raw`, `-derived`, `-artifacts`, `-tfstate`; region `us-central1`; **Uniform** access |
| 6 | ☰ **IAM & Admin** → Service Accounts → CREATE (×5) | `rbp-ingest`, `rbp-prep`, `rbp-train`, `rbp-analysis`, `rbp-modal` |
| 7 | Storage → bucket → **PERMISSIONS** → GRANT ACCESS | per-bucket roles. For `rbp-modal`: objectAdmin **+ ADD IAM CONDITION**, CEL from Ch.7 |
| 8 | ☰ **VPC network** → CREATE VPC NETWORK | `rbp-net`, **Custom** subnets; subnet `rbp-workers` `10.0.0.0/20` with **Private Google Access ON** |
| 9 | ☰ **Artifact Registry** → CREATE REPOSITORY | `rbp`, format Docker, `us-central1` |
| 10 | ☰ **Cloud Build** → History | watch builds; click a failed step for its log |
| 11 | ☰ **Batch** → CREATE JOB | script or container, machine type, parallelism, service account, network |
| 12 | ☰ **Billing** → Budgets & alerts → CREATE | amount, thresholds, **uncheck credits**, **connect a Pub/Sub topic** |
| 13 | ☰ **Cloud Functions** → CREATE FUNCTION | Pub/Sub trigger on that topic; the killswitch code |

**Where the UI genuinely wins:** Cloud Build step logs, Batch's Events tab, and the IAM
policy troubleshooter (paste a principal and a resource, get an explanation of why access was
denied). Use it for diagnosis even when you build with Terraform.

**Where the UI genuinely loses:** anything you must do twice, anything you must review before
applying, and anything you must be able to tear down.

---

# Chapter 15. The debugging ladder

In order. Do not skip a rung; each is cheaper than the next.

1. **Did the control plane accept it?** `gcloud batch jobs describe <job>`. Rejected specs
   fail here with a precise field path.
2. **What do the Events say?** `status.statusEvents`. `OPERATIONAL_INFO` is where quota,
   capacity and preemption appear. This is the rung most often skipped, and it is where
   `CODE_GCE_QUOTA_EXCEEDED` was sitting while I theorised about preemption.
3. **Are there VMs?** `gcloud compute instances list`. No VMs plus a RUNNING job means
   allocation is failing.
4. **What is quota doing?** `gcloud compute project-info describe`. Compare **usage** against
   **limit**, and never filter to `limit > 0` — that is how a survey hid
   `GPUS_ALL_REGIONS = 0` from its own output. Twice.
5. **Did the container start?** Cloud Logging, filtered by `job_uid`. Look for
   `Runnable command line:` — if present, the image pulled and docker ran.
6. **What did the code say?** Your own log lines. Filter out agent state dumps.
7. **Did it write anything?** `gsutil ls` the output prefix. Compare against the completion
   marker: payload without marker means it died between writes, which is the *designed*
   failure.
8. **Is the number right?** `scripts/verify.py`. Everything above can pass while the science
   is wrong.

## 15.1 The failure signatures worth memorising

| symptom | almost always |
|---|---|
| a `list` returns nothing, exit 0 | no project configured. Emptiness ≠ absence |
| job submits, dies in ~2 min | API not enabled, or service agent missing |
| job RUNNING, no VMs, quota error in Events | asking for exactly the quota limit |
| container crash-loops in seconds, 0 inputs | import error. Check what actually got shipped |
| task hangs then hits `maxRunDuration` | outbound connection to a non-Google host with no route |
| `FAILED` job with complete output | task count past the end of the manifest |
| budget shows $0 spent | `INCLUDE_ALL_CREDITS` |
| 403 immediately after granting IAM | policy propagation; retry |
| `repository name must be lowercase` | an unexpanded `$VAR` reached docker |
| suite exits 0 but collected too few | a collection error dropped whole files |
