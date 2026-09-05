# Reproduce the whole study, from raw inputs to verified results

Cloud only. The laptop submits jobs and reads results; it never computes. Roughly
**$5 of GCP credit** and **~$32 on Modal** against a $50 budget.

The study has been run; this is the procedure for running it again from raw inputs. To check
every published number without running any of it, see
[`../SUBMISSION.md`](../SUBMISSION.md):
`python scripts/verify.py --local results/tables` needs only a clone.

---

## What you need before you start

| | |
|---|---|
| A **new, empty GCP project** with billing linked | |
| A **Modal account** with at least **$35** available | stage 9 is 95% of the spend |
| `gcloud`, `terraform`, `modal` on PATH | `python3 -m pip install modal` |
| ~10 GB free disk **in the container**, not locally | Batch provisions it |

Bucket names are **globally unique across all of Google Cloud**. They are derived as
`{project_id}-derived` and `{project_id}-raw`, so a project id nobody has used gives you
bucket names nobody has used. Preflight checks this and tells you if a name is taken.

---

## Step 0. Two things that will bite immediately

**`modal` must be on PATH, not just installed in a venv.** Preflight failed on exactly this
during setup: `modal` lived in a virtualenv while preflight ran under system python. Either
activate the venv or export its bin directory.

**The Modal secret must be recreated for THIS project.** `rbp-gcp` holds a service-account
key, and a key from another project authenticates against the wrong buckets. Terraform
creates the `rbp-modal` service account but deliberately contains no key resource -- a
`google_service_account_key` would put a live private key in Terraform state in plaintext.
So after stage 1:

```bash
gcloud iam service-accounts keys create /tmp/rbp-modal.json \
  --iam-account=rbp-modal@$GOOGLE_CLOUD_PROJECT.iam.gserviceaccount.com
modal secret create rbp-gcp SERVICE_ACCOUNT_JSON="$(cat /tmp/rbp-modal.json)"
rm /tmp/rbp-modal.json
```

If a `rbp-gcp` secret already exists from a previous project, **delete it first** -- Modal
will not overwrite, and stage 9 would silently write into the old project's bucket.

---

## Step 1. Point the pipeline at your project

Two ways. Either edit `config/params.yaml`:

```yaml
cloud:
  project_id: your-new-project-id
```

or just export it, which wins over the file:

```bash
export GOOGLE_CLOUD_PROJECT=your-new-project-id
```

That single value redirects everything: buckets, service accounts, job specs, Modal
secrets, the killswitch. There is **no hardcoded project id anywhere in the source** and
`tests/unit/test_no_hardcoded_project.py` fails the build if one reappears.

Then fill in Terraform's inputs:

```bash
cp cloud/terraform/terraform.tfvars.example cloud/terraform/terraform.tfvars
# set project_id and billing_account
```

---

## Step 2. Preflight. This spends nothing and it is the most important step

```bash
./run.sh preflight
```

It refuses to let anything run until every gate is green. It checks the exact things that
went wrong the first time, all of which were discovered *after* money had been spent:

- required APIs enabled (a missing one lets a job submit and die minutes later)
- billing actually linked
- **a budget that excludes credits.** The default `INCLUDE_ALL_CREDITS` reports $0 spend
  while free credit remains, so no alert threshold can ever fire and the killswitch is
  decorative
- `CPUS_ALL_REGIONS` high enough to run the workers
- **that you are authenticated.** `gcloud <verb> list` with no project returns an *empty
  list*, not an error, so a missing scope reads exactly like a missing result. Preflight
  asserts the project by name, which does fail loudly
- bucket names not already taken by someone else
- Modal authenticated and the `rbp-gcp` secret present

GPU quota is reported but **is not a gate**. It will be 0 on a new project, it cannot be
raised (`NOT_ENOUGH_USAGE_HISTORY`), and that is precisely why the GPU stages run on Modal.

One gate cannot be automated: Modal does not expose a balance via CLI. Check the dashboard
yourself, then:

```bash
./run.sh preflight   # after: PREFLIGHT_ARGS=--modal-credit-ok
```

---

## Step 3. Run the stages

```bash
./run.sh status        # what exists so far
./run.sh stage 1       # one stage
./run.sh from 6        # stage 6 onward
./run.sh all           # everything, still pausing at each paid gate
```

| # | stage | where | cost | produces |
|---|---|---|---|---|
| 0 | preflight | local | $0 | permission to continue |
| 1 | terraform | GCP | $0 | buckets, service accounts, IAM, budget, killswitch |
| 2 | images | Cloud Build | ~$0.50 | CPU + GPU images, weights baked in |
| 3 | ingest **(public internet)** | Batch | ~$0.20 | genome, GENCODE, ClinVar, ENCODE peaks |
| 4 | panel **(public internet)** | Batch | ~$0.10 | candidate datasets, all arms |
| 5 | preprocess **all candidates** + finalize | Batch | ~$2 | matched datasets, all arms, and the pair counts |
| 6 | **select panel** | local | $0 | `manifest/study_panel.tsv` — *the* panel |
| 7 | rehearsal | Batch | ~$0.60 | **R1** |
| 8 | CNN | Batch | ~$3 | **R2** |
| 9 | **SpliceBERT** | **Modal** | **~$31** | **R2** |
| 10 | locality | Modal | ~$0.30 | **R3** |
| 11 | variants **(public internet)** | Batch | ~$0.30 | assignments + phyloP |
| 12 | ClinVar + mismatch control | Modal | ~$0.60 | **R4** |
| 13 | aggregate + figures | Batch | ~$0.10 | tables, figures |
| 14 | **verify** | local | $0 | pass/fail against golden numbers |

Paid stages ask before spending. `RBP_YES=1` skips the prompt once you have decided.

---

## The three things about this design worth understanding

**1. The panel is an artefact, not a flag.** Stage 6 writes `manifest/study_panel.tsv` once
and every later stage reads it. In the original build the panel was an emergent property of
a `--every 2` typed during one sweep, which is why four different dataset counts ended up in
circulation. Stage 6 refuses to redefine an existing panel without `--force`, because
redefining it mid-study invalidates everything downstream.

**Selection cannot precede preprocessing, and that order is forced.** The panel is a
size-ranked sample, and the size (`pairs`) counts the positives that could actually be
matched to a negative -- so it is a RESULT of preprocessing, not an input. Prep therefore
runs on every candidate first. Nothing is saved by trying to invert this: full prep is ~$2,
and the savings from a smaller panel come from the GPU stages, not from prep.

Sampling is **systematic by pair rank**, not a size threshold. AUROC correlates with dataset
size at r = +0.53 to +0.67, so "keep the biggest N" would confound the panel with the
quantity being measured. Stage 6 asserts the kept set spans the 5th and 95th percentile of
the full distribution and refuses to write a size-biased panel.

**2. Task counts are read from manifests, never typed.** `submit_prep.sh` once hardcoded
`COUNT=189`; the gc arm has 187 datasets, so the job dispatched two tasks past the end of its
manifest and Batch reported the whole job FAILED while every real task had succeeded.
`cloud/submit.sh` derives the count from the manifest in GCS every time.

**3. The network posture is not uniform, on purpose.** Workers get Private Google Access and
no external IP, so they reach `*.googleapis.com` and nothing else. Three stages genuinely
need the public internet — ingest (ENCODE, GENCODE, NCBI), panel (ENCODE API) and variants
(UCSC phyloP by HTTP range request) — and those run on a single short-lived VM with a public
IP. There is deliberately **no Cloud NAT**: it bills per VM-hour plus per GB and would hand
internet access to every other worker for no reason.

---

## Step 4. Verify

```bash
./run.sh stage 14
```

Every claim in `config/golden.yaml` is asserted with an explicit tolerance. Tolerances
absorb what legitimately varies — panel size, BLAS thread order, bootstrap seed, GPU-vs-CPU
inference (measured at max 1.1e-4 per variant) — and nothing more.

Where a claim is about unanimity or ordering, the **count** is checked rather than the mean,
because that is what the paper asserts:

- R1: every dataset must fall, and the gain over composition must **at least double** under
  the harder control. If that ratio inverts, the paper's thesis is false regardless of how
  close the other numbers land.
- R2: the four-model ordering must hold exactly.
- R3: SpliceBERT more concentrated on ≥85%, and **significantly reversed on zero**.
- R4: the ladder must be monotone, and matched must exceed mismatched by ≥0.60, or the
  signal is not binding-specific and R4 must be withdrawn.

**Exit 0 means the science reproduced.** Exit 1 names the claim that broke. A pipeline that
completes and quietly produces different science is worse than one that crashes, because
nobody diffs a plausible table.

---

## If something fails

| symptom | cause |
|---|---|
| a `list` command returns nothing | no project configured. Emptiness is not absence |
| budget shows $0 spent | `INCLUDE_ALL_CREDITS`. Preflight catches this |
| GPU quota request denied | expected. `NOT_ENOUGH_USAGE_HISTORY`. Stages 9/10/12 use Modal |
| terraform fails on a bucket | the name is taken globally. Change the project id |
| a Batch job FAILED with complete output | a task count past the end of a manifest |
| Modal containers crash-loop in seconds | an import error. Modal ships only the entrypoint file |
| first GCS call after a laptop wake fails DNS | processes resume before the network |

Every stage is resumable and writes its completion marker **last**, after the payload, so a
stage killed midway redoes its work rather than being skipped by a marker that arrived early.
Rerunning a finished stage costs seconds.
