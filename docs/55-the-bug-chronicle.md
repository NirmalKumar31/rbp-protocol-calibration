# 55. The bug chronicle: sixteen bugs, in the order they happened

Every bug found while rebuilding the study as a cloud pipeline and running it on a genuinely
new GCP project. For each one: the **symptom**, the **mechanism**, **what it would have cost**,
the **fix**, and — the part worth actually learning — **what class of thinking finds it**.

**Every one of these was invisible in the working project.** That is the argument for the
exercise in a single sentence.

---

## How to read this

Bugs fall into five families, and recognising the family is more useful than remembering the
instance:

| family | the shape | how many here |
|---|---|---|
| **Portability** | works because of ambient state on one machine or account | 4 |
| **Silent wrong answer** | completes successfully, produces different science | 4 |
| **Schema / contract** | one component's idea of a name or shape differs from another's | 4 |
| **Diagnosis** | I chose a cause that fit my story rather than the evidence | 2 |
| **Guard** | a safety mechanism that was absent, mis-scoped, or hiding another failure | 2 |

---

# Family: portability

## Bug 1 — Eighteen files hardcoded one GCP project

**Symptom.** None, until the project changes. Then either 403s that look transient, or writes
into the previous study.

**Mechanism.** `storage.Client(project="rbp-composition-2026")` in eighteen files.

**Cost.** The rebuild would not have started. Worse plausible outcome: with cached credentials
it would have *appended to the original study* and nothing would have said so.

**Fix.** `src/rbp/utils/cloud.py` resolves in one place — explicit argument, then environment,
then `config/params.yaml` — with **no default**, because a default is how you write into
someone else's bucket. Buckets derive as `{project_id}-derived`.

**What finds it.** Asking "what does this depend on that is not in the repo?" A grep for your
own project id is thirty seconds and should be a reflex.

---

## Bug 3 — `.gitignore` rule `data/` matched `src/rbp/data/`

**Symptom.** Cloud Build failed the image's own test step:

```
ModuleNotFoundError: No module named 'rbp.data'
```

**Mechanism.** A gitignore pattern without a leading slash matches a directory of that name at
**any depth**. `data/` was meant for the top-level data directory; it also matched
`src/rbp/data/`. Git stopped tracking the package. `gcloud builds submit` honours `.gitignore`.
So the upload excluded a source package and the image was built without it.

**Cost.** Every container missing a module used by preprocessing, variant assignment and
analysis. Caught here; had the image had no test step it would have surfaced as a task failure
hours later.

**Fix.** Anchor every rule: `/data/`, `/results/`, `/logs/`. Verified with
`git check-ignore -v src/rbp/data/annotation.py`, which names the offending rule *and line*.
Confirmed nothing else was hidden.

**What finds it.** Two habits. First, `git check-ignore -v` rather than reasoning about
patterns. Second, **the image test step**: running the suite inside the artefact is what turned
a silent packaging bug into a red build.

---

## Bug 4 — `cloudbuild.*.yaml` still named the old project

**Symptom.** Build step 0 (`pull-cache`):

```
denied: Permission 'artifactregistry.repositories.downloadArtifacts' denied on resource
'//artifactregistry.googleapis.com/projects/rbp-composition-2026/...'
```

**Mechanism.** `_IMAGE` and `_ARTIFACTS` substitution defaults still contained the old project
id, so the build went looking for a cache layer in a registry it had no access to.

**Cost.** Small in itself. Significant because it survived a test written to prevent exactly
this — see bug 6.

**Fix.** Use `$PROJECT_ID`, and pass concrete values from `run.sh` (see bug 5).

---

## Bug 16 — Modal secret from a different project

**Symptom.** Would have been none. Stage 9 would have written into the previous project's
bucket.

**Mechanism.** `rbp-gcp` holds a service-account key. A key from another project authenticates
against that project's buckets. **Modal will not overwrite an existing secret**, so creating it
again is a no-op that looks like success.

**Cost.** ~$31 of GPU work written into the wrong project, and a "reproduction" that was
actually reading the original's outputs.

**Fix.** `modal secret delete rbp-gcp` before creating, and the runbook now says so
explicitly. Verified by reading `project_id` out of the key JSON before uploading it.

**What finds it.** Asking "which project does this credential actually name?" rather than "does
a credential exist?"

---

# Family: silent wrong answer

## Bug 2 — The dataset panel was a command-line flag

**Symptom.** Four dataset counts in circulation (189, 187, 95, 94) with no document reconciling
them.

**Mechanism.** The 95-dataset panel came from `--every 2` typed during one sweep. Nothing
recorded it, so it could not be reproduced, audited, or explained.

**Cost.** Already paid, repeatedly, in confusion. Forward cost: a rerun with a different flag
silently produces a different study.

**Fix.** `scripts/select_panel.py` writes `manifest/study_panel.tsv` once; every later stage
reads it; it **refuses to redefine an existing panel** without `--force`; and it asserts the
sample spans the 5th–95th percentile so it cannot be size-biased. `docs/PANELS.md` plus
`tests/unit/test_panel_counts.py` pin the four counts and their nesting.

**What finds it.** Asking of every number in the paper: "where is this written down?" If the
answer is shell history, it is not a parameter.

---

## Bug 12 — `submit.sh` did not wait, so aggregation ran against an empty bucket

**Symptom.** None observed, because I checked before trusting it.

**Mechanism.** `gcloud batch jobs submit` returns when the **control plane accepts the JSON**.
Stage 5's next line was `cloud_prep.py finalize`.

**Cost.** This is the worst bug in the list. Finalize would have run seconds after submission,
seen almost no processed datasets, and written **a panel of zero (or a handful of) datasets**.
Stage 6 would have sampled it, stages 7–12 would have run on it, and stage 14 would have failed
with numbers that looked like a science problem rather than a plumbing one.

**Fix.** `submit.sh` polls to completion by default, printing task counts, with `NO_WAIT=1`
for the parallel track. FAILED exits **2**, not 1, because a failed job is not necessarily a
failed run — preemption leaves resumable failures, and bug 8 failed a job whose every real
task succeeded.

**What finds it.** Reading "successfully submitted" as the control-plane statement it is.
Whenever a command returns instantly for work that cannot be instant, ask what it actually
promised.

---

## Bug 9 — Two figures read table names nothing wrote

**Symptom.** Silent. Figures skipped.

**Mechanism.** `figures.py` needed `matched95_four_models.csv` and
`variant_results_splicebert.csv`; `cloud_analysis.py` writes `matched_four_models.csv` and
`variant_ladder.csv`.

**Cost.** A paper missing two of its five figures, with the pipeline reporting success.

**Why it was silent, and this is the interesting part.** `need()` is *designed* to skip missing
tables quietly, so `figures.py` can run mid-pipeline while other stages are still producing
inputs. That is a correct and useful design. It also hid this bug completely. **A tolerant gate
needs a strict counterpart somewhere** — here, `verify.py`.

**Fix.** Names aligned; `f4` rewritten against `variant_ladder.csv`, which is better science
anyway (three rungs where the claim is the gaps).

**What finds it.** Listing which artefacts each script **writes** against which it **reads**,
as sets, and diffing. Not reading the code — diffing the sets.

---

## Bug 10 — `verify.py` crashed instead of falling back

**Symptom.** `ValueError: The truth value of a DataFrame is ambiguous.`

**Mechanism.** `d = T.get(a) or T.get(b)`. `or` evaluates the left operand's truthiness, and
pandas refuses to define that for a DataFrame.

**Cost.** The verifier — the single stage whose whole job is to fail loudly and *correctly* —
would have died with an unrelated error, which reads as "the verifier is broken" rather than
"the science did not reproduce".

**Fix.** Explicit `if d is None:`.

**What finds it.** Executing it. This is pure syntax-valid, type-invalid Python; no amount of
reading catches it, and it is why `verify.py` was then run against the original run's real
tables (32/32) rather than merely imported.

---

# Family: schema and contract

## Bug 5 — `$PROJECT_ID` is not expanded inside a substitution default

**Symptom.**

```
invalid argument "us-central1-docker.pkg.dev/$PROJECT_ID/rbp/cpu" for "-t, --tag" flag:
  invalid reference format: repository name must be lowercase
```

**Mechanism.** Cloud Build substitutes `$PROJECT_ID` in **step arguments**, but does **not**
recursively expand it inside a user-defined substitution's *default value*. The literal string
reached docker, which rejected it because `$PROJECT_ID` contains capital letters.

**Fix.** `run.sh` passes `--substitutions=_IMAGE=...,_ARTIFACTS=...,_GIT_SHA=...` with concrete
values.

**What finds it.** Reading the error literally. It says the repository name must be lowercase,
and there are capitals visible in the string — the message was exactly right.

---

## Bug 8 — Task counts typed rather than read

**Symptom.** A Batch job reported **FAILED** with complete, correct output.

**Mechanism.** `COUNT=189` hardcoded. The gc arm has 187 datasets, so two tasks were dispatched
with indices past the end of the manifest, exited non-zero, and failed the job.

**Cost.** Already paid once on the original project: a completed run that looked failed, and
the hour spent working out which.

**Fix.** `manifest_rows()` reads the count from the manifest in GCS, every time. **A count you
typed is a count that will be wrong the first time the panel changes size.**

---

## Bug 11 — `networkInterfaces` at the wrong nesting level, and the wrong network

**Symptom.**

```
INVALID_ARGUMENT: Unknown name "networkInterfaces" at 'job.allocation_policy'
```

**Mechanism.** Batch expects `allocationPolicy.network.networkInterfaces`, one level deeper.
Separately, I had named the `default` network rather than `rbp-net`.

**Cost.** The nesting error is loud and free. **The network name was the dangerous half**: it
would have worked, and quietly placed every worker on the default network with a route to the
internet — discarding the entire point of `network.tf`.

**Fix.** Correct nesting; `rbp-net` and the `rbp-workers` subnet with `noExternalIpAddress`;
external-IP stages omit the block entirely.

**What finds it.** For the nesting, the API. For the network name, comparing against the
working submitter rather than writing from memory. **The bug that would have worked is the
dangerous one.**

---

## Bug 14 — `test_train_folds.py` needs torch transitively

**Symptom.** GPU-less CPU image failed its own test step with
`ModuleNotFoundError: No module named 'torch'`.

**Mechanism.** The test imports `rbp.train.data` for one helper. That module imports torch at
line 10. The CPU image has no torch **by design** — 1.2 GB against the GPU image's 6 GB,
because 488 preprocessing tasks should not pull CUDA they never use. Only `test_models.py` was
excluded, because it is the only test that imports torch *directly*.

**This is a latent bug in the original repo.** The CPU image predates the test, so that gate
has been failing since the day the test was added, and nobody rebuilt the image to find out.

**Fix.** Exclude `test_train_folds.py` too, with a comment explaining that transitive imports
are the reason exclusion lists rot.

**What finds it.** Building the artefact. A test-exclusion list maintained by filename drifts
the moment an import chain changes, and only a rebuild reveals it.

---

# Family: diagnosis

## Bug 15 — I called a quota problem "spot preemption"

**Symptom.** vCPU usage dropped 8 → 4. Throughput 2.2 datasets/min against an expected ~10.

**What I did.** Concluded spot preemption, switched every job to on-demand, wrote a commit
message about it, and told the user it was fixed.

**What was actually true.** The job's events said:

```
OPERATIONAL_INFO: CODE_GCE_QUOTA_EXCEEDED
```

and contained **zero** preemption events. `parallelism 12 ÷ 4 tasks-per-node = 3 nodes ×
4 vCPU = 12 vCPU`, which is exactly `CPUS_ALL_REGIONS`. **VM creation fails at the limit, not
approaching it.** Batch retried the third node forever; the job ran 8-wide regardless; the dip
to 4 was a node cycling out with its replacement refused.

**How I fooled myself.** My evidence was a single grep hit from
`grep -icE "preempt|QUOTA_EXCEEDED|FAILED"` — a pattern that matches the quota error too. I
counted one hit and read it as the thing I already believed.

**Proof the fix was not a fix.** I measured after switching to on-demand: still 2.2/min.

**Real fix.** `parallelism 8` on every fan-out stage: identical throughput, no failed
creations, no alarming error that changes nothing. Real speedup needs a quota increase.

**What finds it.** Two things. Grep patterns with alternation cannot tell you *which* branch
matched — print the lines. And **measure after fixing**: had I not, a wrong diagnosis would be
in the documentation as fact.

---

## Bug 13 — Estimating from a single sample

**Symptom.** I told the user prep would take 45–50 minutes. It took ~3.5 hours.

**Mechanism.** I extrapolated from one log line — a task that finished in 47.5s — without
noticing it was a small dataset, and that the manifest is deliberately sorted **biggest-first**
so the early tasks are the *slowest*. Measured over a 5-minute window: 2.2/min.

**Cost.** Trust. A wrong ETA is worse than no ETA, because plans get made against it.

**Fix.** Sample twice and divide:

```bash
A=$(count); sleep 300; B=$(count)
rate=$(( (B-A) / 5 ))   # per minute, measured
```

**What finds it.** Knowing your own scheduling policy. Longest-first is a deliberate choice
documented in the manifest builder, and it makes early throughput unrepresentative *by
design*.

---

# Family: guards

## Bug 7 — The Terraform backend was hardcoded to the original project

**This is the one that would have destroyed the study.**

**Symptom.**

```
Plan: 63 to add, 1 to change, 63 to destroy.
```

**Mechanism.** `main.tf` had `backend "gcs" { bucket = "rbp-composition-2026-tfstate" }`.
**Terraform backends cannot use variables**, so every checkout shares one state file. Pointing
the config at the new project made `plan` load the **original** project's state. Terraform then
computed the correct diff for what it was told: destroy those 63 resources, create these 63.
The 63 destroys were the original study's buckets, holding every result — including 485
SpliceBERT checkpoints representing ~$31 of GPU time.

And `run.sh` stage 1 ran `terraform apply -auto-approve`. **One command, no prompt.**

**What saved it.** Reading the plan instead of trusting it. Only `init` and `plan` ever ran,
both read-only; the original state file's timestamp confirms nothing wrote to it.

**Fix, three parts:**

1. **Partial backend config.** Remove the bucket from the file; pass
   `-backend-config="bucket=${PROJECT_ID}-tfstate"` at init, with `-reconfigure` so a cached
   backend from another checkout cannot leak in.
2. **A destroy guard.** Save the plan, count destroy actions, refuse to apply if any exist. A
   first apply on a fresh project is additive by definition, so a destroy proves the state
   describes a different project.
3. Apply the **saved plan file**, so what was inspected is exactly what runs.
   `apply -auto-approve` without a saved plan re-plans and may differ.

**What finds it.** Reading `terraform plan` output, every time, and specifically counting
destroys. `-auto-approve` on a plan nobody read is the actual defect; the hardcoded bucket was
merely the trigger.

---

## Bug 6 — My own hardcoded-project test never searched `docker/`

**Symptom.** Bug 4 passed a test written to prevent bug 4.

**Mechanism.**

```python
SEARCH = ["scripts", "src", "cloud"]     # docker/ absent
```

**Cost.** False confidence, which is worse than no test — I had cited it as evidence the
codebase was portable.

**Fix.** Add `docker`. Then, when a comment edit broke a YAML file, add a test that **parses**
both cloudbuild files rather than only grepping them.

**What finds it.** Asking what a passing test actually certifies. **A test that does not look
everywhere only certifies the places it looked**, and the gap is invisible from the green tick.

---

# The remaining two, for completeness

## Bug 17 — I broke a YAML file editing a comment into it

Slicing the leading `# ` off the first line of a comment block turned it into a bare YAML key.
`gcloud builds submit` failed before uploading anything. Fixed, and a test now parses both
cloudbuild files: **a YAML file that looks fine in a diff is not a YAML file that parses.**

## Bug 18 — The GPU test-count gate was still an equality

Fixed the CPU file's exact-count gate to a floor and forgot the GPU file. It then failed with
`collected 548 tests (expected 460)` — it failed *because tests were added*. An exact count
punishes the right behaviour and trains people to edit the number until it goes green.
Under-collection is the real failure mode; a floor catches it.

---

# What the sixteen have in common

**Three quarters of them cannot happen on a machine that already works.** They need a fresh
project, a fresh account, a fresh clone. That is not an argument for being careful; it is an
argument that **carefulness is not a substitute for a clean-room run**.

The four that could have produced wrong science without failing — the panel-as-a-flag, the
non-waiting submit, the mismatched table names, and the wrong-project Modal secret — are all
the same shape: **a component that succeeds while doing the wrong thing**. Every guard in this
pipeline exists because of that shape:

- completion markers written last, so a crash redoes work rather than skipping it
- task counts read, never typed
- the panel written once and refusing redefinition
- images pinned by digest
- a verify stage that asserts numbers rather than exit codes

And the two diagnosis bugs are a different lesson entirely, one no guard fixes: **I chose a
cause that fit the story I already had, twice.** The defence is mechanical, not intellectual —
print the matching line rather than counting matches, and measure again after the fix.

---

# The 2026-08-26 batch: seven more, and the worst four were mine

The bugs above are infrastructure. These are worse, because six of the seven produced
**green checks on wrong answers** rather than failures.

## Bug 25 — A stage did ninety minutes of correct work and 403'd on the last line

Stage 11 ran assign, score and phyloP over 66,010 variants, uploaded three correct tables,
then died writing `variants-complete.json` at the bucket root. `rbp-analysis` is scoped by an
IAM condition to `results/`, `variants/` and `driver/`; the bucket root is outside every
allowed prefix.

Four consecutive jobs recomputed identical results to reach the same failing line, because
`run_existing()` raised before `stage_out()` could save anything. The guard was right and the
path was wrong — the third time in this project an IAM condition caught a real mistake.

**Fix:** marker moved under `results/`; finished work now uploads even when a later sub-stage
fails. **And the general fix:** `tests/unit/test_write_paths_are_permitted.py` resolves every
static blob key in the codebase, including function locals, and checks it against the
Terraform conditions — treating `variant_splicebert.py` as dual-identity because Modal shells
out to it, so its writes must satisfy the *intersection* of two identities.

That test immediately found **three more** of the same bug waiting: `analysis-complete.json`
at the root, the phyloP cache to `interim/` (which I had just introduced), and the variant
task manifest to `manifest/`. All three would have failed at the end of a completed stage.

## Bug 26 — `tar --exclude='data'` deleted a source package

Packaging the repo for the driver VM, I wrote `--exclude='data'`. That matches a path
component **at any depth**, so it silently dropped `src/rbp/data/` — a Python package, not a
data directory. Modal uploaded that tree to 188 containers and every one died with
`ModuleNotFoundError: No module named 'rbp.data'`.

This is **the same bug as the `.gitignore data/` entry** already recorded above, in a
different tool. Knowing the story did not stop me typing the unanchored version.

**Fix:** `cloud/package_repo.sh` replaces the inline tar and *proves* the archive before
upload — every package under `src/rbp/` present, entrypoints present, and all nine modules
actually importable from an unpacked copy. It refuses to upload otherwise.

## Bug 27 — A task count typed instead of read, in the one place not driven by `submit.sh`

`modal_variants.py` had `N_TASKS = 94`, the earlier study's variant panel. This pipeline's
manifest lists 95, so `range(N_TASKS)` skipped index 94 — K562 ZNF800, 288 variants — with
nothing in any log. The driver's own gate was `-ge 94`, so it would have been *satisfied* by
94 of 95 and proceeded to the analysis.

`submit.sh` derives every Batch count from `manifest_rows()` for exactly this reason. Modal
was the one path that did not.

## Bug 28 — The analysis fetched its own output

`cloud_analysis.main()` fetched `results/tables/locality_ism.csv` from GCS — an object nothing
ever created. Stage 10 wrote 95 per-dataset JSONs and nothing aggregated them. So R3 quietly
had no table, `do_figures()` skipped the R3 figure without complaint, and **stage 14 was the
first thing to notice**, after the whole pipeline had run green.

Every other result had a producer. R3's was a fetch of itself.

## Bug 29 — The verifier certified a number nobody claimed

Stage 14 computed R1's "gain over composition" as `auroc - composition_auroc`: the difference
between two *separately fitted* models, with no confidence interval and no p-value. The claim
is the **nested** gain — composition alone versus composition plus the sequence score — which
the rehearsal already computes as `delta_auroc` with a CI and a per-dataset `helps` flag.

Naive: 3.94x. Nested: 2.50x. The gate was blessing 3.94x while the write-up quoted the nested
figure. **A gate that certifies a number nobody claims is worse than no gate, because it reads
as confirmation.**

## Bug 30 — The headline statistic was the wrong statistic, and 33/33 passed

R4 pooled ~19k variants across 95 datasets into one AUROC per arm. Mean |delta| per dataset
correlates with that dataset's pathogenic rate at **+0.73** and spans **10.4x**, so the pooled
number partly measured *which dataset a variant came from*. Matched 0.829 pooled against 0.755
paired; the specificity gap +0.149 against +0.065.

Conservation was the only arm immune, being on a fixed external scale — **so the artefact was
invisible precisely because the uninflatable arm was winning anyway.**

This was found only because a council reviewer challenged the wrong-protein control as
possibly contaminated. The contamination test came back clean; the pooling problem fell out
of the same analysis. **The bug was found by attacking a different claim.**

## Bug 31 — A trivial baseline beat the model, and nothing was checking

"What fraction of the OTHER variants in this 1-Mb window are pathogenic", leave-one-out, no
sequence and no model, reaches 0.8139 where the model reaches 0.7553. The model does not beat
it in any stratum. Nothing in the pipeline had ever asked.

`golden.yaml` now asserts `model_minus_prevalence_max: 0.02`, so a future run **fails** if
anyone claims the model beat the trivial rule, and asserts the all-datasets stratum that shows
nothing (gap −0.011, p=0.87) so it cannot be quietly dropped.

---

# What this batch has in common, and it is not what the first batch had

The first twenty-four were mostly **infrastructure failing loudly**. These seven are
**verification succeeding on the wrong thing**:

- three write paths that would 403 only after the work was done (25)
- a package that imported fine locally and not in the container (26)
- a count that silently dropped one dataset in 95 (27)
- a result table nothing produced, hidden by a tolerant figure guard (28)
- a gate certifying a quantity the paper never claimed (29)
- a gate certifying an inflated statistic, 33/33 green (30)
- a trivial baseline nobody had run (31)

**The pattern: every check I wrote confirmed the thing it was derived from.** The write-path
test, the golden values, the figure guards — each was built from the code it was checking, so
each agreed with it. The three that actually caught something came from *outside*: an IAM
condition written from a security model rather than from the code, a council reviewer with no
stake in the result, and a baseline borrowed from a paper that disagreed with us.

The defence is not more tests. It is **at least one check per claim that was not derived from
the code producing it.**
