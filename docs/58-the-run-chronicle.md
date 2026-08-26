# 58. The run chronicle: building and running it, in order, with everything that broke

This is the narrative record of 2026-08-25: the day the study was rebuilt as a cloud pipeline
and run end to end on a brand-new GCP project and a brand-new Modal account.

**Twenty-four bugs.** Every one is here with its symptom, its mechanism, what it would have
cost, and — the part worth actually learning — what *kind* of thinking finds it.

Read [55-the-bug-chronicle.md](55-the-bug-chronicle.md) for the same bugs grouped by family.
This file is chronological, because the order matters: several bugs were only reachable once
an earlier one was fixed, and two of my own fixes created new bugs.

---

## How to read this

Each episode follows the same shape:

> **What I was trying to do** → **what happened** → **why** → **what it would have cost** →
> **the fix** → **the lesson**

Skip to the lesson if you are in a hurry. Read the mechanism if you want to be able to find
that class of bug yourself.

---

# Act I. Building the thing (morning)

## Episode 1. The audit that started it

The study was scientifically finished. Four results, verified, two adversarial attacks
tested and survived. So why rebuild anything?

Because I asked one question: **"what does this depend on that is not in the repo?"**

```bash
$ grep -rl "rbp-composition-2026" --include="*.py" --include="*.sh" . | wc -l
18
```

Eighteen files hardcoded one GCP project. That is invisible for as long as you only ever run
in that project, and a total failure the first time anyone runs it anywhere else.

**Why it is worse than it looks.** A hardcoded project does not fail loudly. It fails by
*writing into the wrong account*, or by *reading a stale bucket and finding nothing*. Neither
says "this is not portable."

Three more defects followed from the same question:

- the 95-dataset panel existed only as `--every 2`, typed once during a sweep
- four analysis stages ran only on a laptop
- nothing verified that a rerun produced the same science

**Lesson.** "It works" and "it is reproducible" are different claims, and only the first was
true. The audit question is not "is the code correct" but **"what is not written down."**

## Episode 2. `rbp.utils.cloud`, and the decision not to have a default

```python
def _resolve(explicit, env_key, cfg_key, what):
    if explicit: return explicit
    if os.environ.get(env_key): return os.environ[env_key]
    v = _from_config().get(cfg_key)
    if v: return v
    raise RuntimeError(f"{what} is not configured. ...")
```

Four decisions in nine lines:

1. **Most-specific first**: argument beats environment beats config. That is the order of
   increasing generality, so a caller can always override.
2. **No default.** A default is how you write into somebody else's bucket. If the environment
   is not configured, that is a bug in the environment and the program should stop.
3. **Buckets derive from the project id.** Bucket names are globally unique *across all of
   Google Cloud*, so a project id nobody has used gives you bucket names nobody has used, for
   free, and one variable configures the system.
4. **`@lru_cache` on the config read**, because this is called in loops.

And a test, parametrised **per file** so a failure names the offender rather than dumping a
list.

## Episode 3. Three cross-file bugs, found by diffing sets rather than reading code

Before running anything, I listed which tables each script **writes** against which it
**reads**, as sets, and diffed them.

**Bug 1: `verify.py` crashed instead of falling back.**

```python
d = T.get("matched_four_models.csv") or T.get("matched95_four_models.csv")
```

`or` evaluates the left operand's truthiness, and pandas raises `ValueError` for a DataFrame.
The one stage whose entire job is to fail *loudly and correctly* would have died with an
unrelated error — which reads as "the verifier is broken" rather than "the science did not
reproduce."

**Bug 2: two figures read table names nothing writes.** `figures.py` wanted
`matched95_four_models.csv`; `cloud_analysis.py` writes `matched_four_models.csv`.

**Why this one was silent, and it is the interesting part.** `need()` is *designed* to skip
missing tables quietly so `figures.py` can run mid-pipeline while other stages are still
producing inputs. That design is correct and useful. It also hid this bug completely.

> **A tolerant gate needs a strict counterpart somewhere.** Here it is `verify.py`.

**Bug 3: `f4` read a table that no longer exists in that shape.** Rewritten against
`variant_ladder.csv` — which is better science anyway: four bars where the claim is the
*gaps* between them.

**Lesson.** Reading code finds bugs in logic. **Diffing the set of artefacts produced against
the set consumed** finds bugs in contracts, and contracts are where multi-stage systems break.

---

# Act II. The near-catastrophe (midday)

## Episode 4. `Plan: 63 to add, 1 to change, 63 to destroy`

I pointed the config at the new project and ran `terraform plan`. It proposed **destroying 63
resources** — the original study's buckets, holding every result, including 485 SpliceBERT
checkpoints representing ~$31 of GPU time.

**The mechanism:**

```hcl
backend "gcs" {
  bucket = "rbp-composition-2026-tfstate"    # HARDCODED
}
```

**Terraform backends cannot use variables.** So every checkout of the repo shares one state
file. Pointing the config at a new project made `plan` load the *old* project's state, and
Terraform then computed the correct diff for what it had been told: destroy those, create
these.

**And `run.sh` stage 1 ran `terraform apply -auto-approve`.** One command, no prompt.

**What saved it:** reading the plan instead of trusting it. Only `init` and `plan` ever ran,
both read-only. The original state file's timestamp confirms nothing wrote to it.

**Three fixes, each worth internalising:**

1. **Partial backend config.** Remove the bucket from the file; pass it at init with
   `-backend-config="bucket=${PROJECT_ID}-tfstate"` **and `-reconfigure`** — without the
   latter, init reuses whatever backend a previous checkout cached in `.terraform/`, which is
   exactly how another project's state leaks in.
2. **A destroy guard.** Save the plan, count destroys, refuse to apply if any exist. A first
   apply on a fresh project is additive *by definition*.
3. **Apply the saved plan file**, so what you inspected is exactly what runs.
   `apply -auto-approve` without a saved plan re-plans and can differ.

**Lesson.** The hardcoded bucket was the trigger. **The actual defect was `-auto-approve` on
a plan nobody read.**

## Episode 5. `for_each` over values that do not exist yet

```hcl
for_each = toset([google_service_account.prep.email, ...])   # "known only after apply"
```

`for_each` keys must be known at **plan** time. Terraform refused to plan *at all*, which also
blocked `terraform import` of anything else — so I could not even adopt existing resources to
get unstuck.

Fixed by keying on static account ids and building the email from `var.project_id`, which *is*
known at plan time. Identical grants; configuration now plannable from empty.

**Lesson.** `for_each` over computed values is an anti-pattern, not a style preference.

---

# Act III. The image builds (afternoon)

Four consecutive failures, each a different lesson.

## Episode 6. A `.gitignore` rule deleted a source package

```
ModuleNotFoundError: No module named 'rbp.data'
```

`.gitignore` had `data/`. **Unanchored, so it matches a directory of that name at any depth** —
including `src/rbp/data/`. Git stopped tracking the package. `gcloud builds submit` honours
`.gitignore`. The image shipped without it.

**The image's own test step caught it**, which is that gate doing exactly its job.

Fix: anchor every rule (`/data/`). Verify with `git check-ignore -v <path>`, which names the
offending rule *and line*.

## Episode 7. The cloudbuild files still named the old project — and my test did not look there

Step 0 took a permission denial fetching a cache layer from `rbp-composition-2026`'s registry.

Worse: **it survived a test written to prevent exactly this.**

```python
SEARCH = ["scripts", "src", "cloud"]     # docker/ absent
```

**Lesson.** A test that does not look everywhere only certifies the places it looked, and the
gap is invisible from the green tick.

## Episode 8. `$PROJECT_ID` is not expanded inside a substitution default

```
invalid reference format: repository name must be lowercase
```

Cloud Build substitutes `$PROJECT_ID` in step *arguments* but does **not** recursively expand
it inside a user-defined substitution's *default value*. The literal string reached docker,
which rejected it for containing capitals.

**Lesson.** Read the error literally. It said the name must be lowercase, and there were
capitals visible in the string. The message was exactly right.

## Episode 9. I broke the YAML editing a comment into it

Slicing the leading `# ` off a comment's first line turned it into a bare YAML key.
`gcloud builds submit` failed before uploading anything.

Fix: a test that **parses** both cloudbuild files. **A YAML file that looks fine in a diff is
not a YAML file that parses.**

## Episode 10. A transitive torch import, and a gate that punished good behaviour

The CPU image has no torch **by design** — 1.2 GB against the GPU image's 6 GB, because 488
preprocessing tasks should not pull CUDA they never use. Only `test_models.py` was excluded,
because it is the only test that imports torch *directly*.

`test_train_folds.py` imports `rbp.train.data` for one helper, and that module imports torch at
line 10.

**This is a latent bug in the original repo**: that gate has been failing since the day the
test was added, and nobody rebuilt the image to find out.

Then the GPU build failed too:

```
collected 548 tests (expected 460)
FAIL: expected 460 tests, collected 548
```

**It failed because tests were ADDED.** An exact count punishes the right behaviour and trains
people to edit the number until it goes green. Changed to a floor — under-collection is the
real failure mode.

---

# Act IV. Running it (afternoon into evening)

## Episode 11. Stage 3 rejected outright, and the dangerous half

```
INVALID_ARGUMENT: Unknown name "networkInterfaces" at 'job.allocation_policy'
```

Batch wants `allocationPolicy.network.networkInterfaces`, one level deeper.

**And separately I had named the `default` network instead of `rbp-net`.** That one would have
*worked* — quietly placing every worker on a network with a route to the internet, discarding
the whole point of `network.tf`.

**Lesson.** The bug that would have worked is the dangerous one. The loud one is free.

Also fixed here: every job was running as `rbp-train`. Four service accounts exist so a
preprocessing task cannot write a model. One identity discards that for no gain.

## Episode 12. Ingest and panel reproduce exactly

First real evidence the pipeline reproduces rather than merely runs:

| check | original | reproduced |
|---|---|---|
| raw bucket total | 3.97 GiB | **3.97 GiB** |
| candidate panel | 139 K562 + 105 HepG2 | **139 + 105** |

## Episode 13. `submit.sh` did not wait — the worst bug in the list

`gcloud batch jobs submit` returns when the **control plane accepts the JSON**. Stage 5's very
next line was `cloud_prep.py finalize`.

**What that would have done:** finalize runs seconds after submission, sees almost no processed
datasets, and writes **a panel of zero datasets**. Stage 6 samples it. Stages 7–12 run on it.
Stage 14 fails with numbers that look like a science problem rather than a plumbing one.

Fix: poll to completion, with `NO_WAIT=1` for the fire-and-forget case. **FAILED exits 2, not
1**, because a failed job is not necessarily a failed run.

**Lesson.** Whenever a command returns instantly for work that cannot be instant, ask what it
actually promised.

## Episode 14. I diagnosed a quota problem as spot preemption

vCPU usage dropped 8 → 4. Throughput 2.2 datasets/min against an expected ~10.

I concluded spot preemption, switched every job to on-demand, wrote a commit message about it,
and told the user it was fixed.

**Then I measured: still 2.2/min.**

The job's events had said, the whole time:

```
OPERATIONAL_INFO: CODE_GCE_QUOTA_EXCEEDED
```

with **zero** preemption events. `parallelism 12 ÷ 4 tasks-per-node = 3 nodes × 4 vCPU = 12
vCPU`, which is exactly `CPUS_ALL_REGIONS`. **VM creation fails AT the limit, not approaching
it.** Batch retried the third node forever; the job ran 8-wide regardless; the dip to 4 was a
node cycling out with its replacement refused.

**How I fooled myself:** my evidence was one hit from
`grep -icE "preempt|QUOTA_EXCEEDED|FAILED"` — a pattern that matches the quota error too. I
counted a match and read it as the thing I already believed.

**Two mechanical defences, not intellectual ones:**

- grep patterns with alternation cannot tell you which branch matched. **Print the lines.**
- **measure again after fixing.** Had I not, a wrong diagnosis would be in the documentation
  as fact.

## Episode 15. I gave a 3× wrong ETA

I told the user prep would take 45–50 minutes, extrapolated from one log line showing a 47.5s
task. It took ~3.5 hours.

The manifest is sorted **biggest-first by design**, so early tasks are the *slowest* and early
throughput is unrepresentative. Measured over a five-minute window: 2.28/min.

**Lesson.** One sample is not a rate. A wrong ETA is worse than no ETA, because plans get made
against it.

## Episode 16. The panel was not deterministic

Prep finished 488/488 and the counts reproduced exactly: dinuc 189 (88 + 101), gc 187 (88 +
99), matching the original including the per-cell split.

But the 95-dataset panel had **93 of 95 shared, 2 swapped.**

**Cause:** `sort_values("pairs")` defaults to **quicksort, which is not stable**, and three
pairs of datasets share a pair count exactly — 539, 3640, 7988. Which member of a tied pair
landed on an even index was arbitrary, so `[::2]` kept a different one. Two of the three
flipped.

**Is it an issue?** Not for the science — a tie means identical size, so the panel's size
distribution and its size-unbiasedness do not move. **It is fatal for the claim**: "the 95" is
supposed to name one specific set; that is the entire reason for writing the panel down.

Fixed with `sort_values(["pairs","dataset"], kind="mergesort")`. Verified: membership now
**identical**, 95 datasets, 79 proteins, 94 in both arms.

**Lesson.** This is a bug no code review finds. It needs a second run against real data that
happens to contain ties.

## Episode 17. Two least-privilege gaps, exposed by a 403

Stage 7 failed: every rehearsal task fitted its models correctly and then took a 403 uploading
`rehearsal/{arm}/{cell}/{name}.json`. `rbp-train`'s IAM condition allowed `runs/` and `ckpt/`
only.

**The condition was working; my service-account assignment was the mistake.**

And diagnosing it exposed why the earlier build never hit this: **`rbp-prep` held
unconditional `objectAdmin` on the entire derived bucket** — quietly the broadest identity in
the project, able to delete every trained model while only ever needing to write preprocessing
output. The old build ran the rehearsal under that account and was simply allowed to write
anywhere.

Result: no unconditional writers remain.

```
rbp-modal    runs | ckpt | variants
rbp-train    runs | ckpt | rehearsal | results
rbp-prep     processed | panel | manifest | interim
```

**And the apply exposed a gap in my own destroy guard.** It reported `0 destroys`, then the
apply printed `2 added, 0 changed, 2 destroyed` — because **Terraform describes a replacement
as "must be replaced", never "will be destroyed"**. Harmless here (IAM bindings re-created
after a title change), but a replacement of a data-bearing resource is a destroy with a
friendlier name, and a bucket replacement is data loss. The guard now counts both.

## Episode 18. The laptop was still part of the pipeline

The rehearsal wrote one manifest per arm, which forced one Batch job per arm, which forced a
**local process** to wait for the first and submit the second. Close the lid between arms and
the second never starts.

`cloud_prep.py` had already solved this by putting `arm` on the manifest row. Same pattern
applied: **189 dataset-arms (95 dinuc + 94 gc) in one job.**

The task now takes its arm from the **row**, not from `--arm`. That matters beyond
convenience: with the arm coming from an environment variable, the same manifest index means
different work depending on how the job was invoked, and `BATCH_TASK_INDEX` is only meaningful
against a fixed interpretation.

**Lesson.** Work that *runs* in the cloud but is *controlled* from a machine that can go away
is not cloud work.

## Episode 19. A SUCCEEDED job that did nothing

The combined job reported **189/189 SUCCEEDED**. `rehearsal/gc/` was empty.

```
image built at:          22e54a0
arm-from-row landed in:  79ca391
```

**The container was running code from before the change.** All 189 tasks took `ARM=dinuc` from
the environment, the 94 gc tasks found a dinuc marker already present, logged "already present,
nothing to do", and exited 0. Green job, zero gc work.

The completion-marker design worked perfectly — it just protected the wrong thing, because the
code deciding *which* marker to check was stale.

**Lesson.** After any code change, confirm the digest changed:
`gsutil cat gs://${PROJECT}-artifacts/images/cpu_digest.txt`.

## Episode 20. The device guard was right; the job spec was wrong

The CNN sweep had 12 tasks exit 1 immediately.

```python
if a.device is None and device.type != "cuda":
    sys.exit("no CUDA device visible; refusing to run a GPU task on CPU.")
```

That guard exists because a GPU node silently falling back to CPU bills at GPU rates, runs
~100× slower, and says nothing in the logs. On this project GPU quota is 0, so the CNN
genuinely runs on CPU — **and the spec has to say so rather than let a safety check fire.**

The image stays the GPU one, because that is where torch lives. Carrying unused CUDA is the
cheap mistake; having no torch is fatal.

## Episode 21. Modal died because the app was ephemeral

The SpliceBERT sweep stopped at 207 of 475 fold-runs.

```
modal.exception.ConflictError: App state is APP_STATE_STOPPED
RemoteError: Function call was cancelled by user or a failure.
```

No cost limit was close: $3 of a $30 free credit, $0 out of pocket against a $20 cap.

**`modal run` without `--detach` ties the app's lifetime to the local client.** The shell
exiting, a network blip, or the lid closing stops the app and every running container. Four
hundred and seventy-five GPU tasks were hanging off one `nohup`'d process on a laptop.

**This is Episode 18 again, on the other cloud**, three hours later. I fixed the pattern on
GCP and did not think to look for it on Modal.

Resumed rather than restarted — completion markers skipped the 207 already done.

**Does the interruption affect validity?** No, and the design anticipated it. Each fold-run is
an independent task; `metrics.json` is written **last**; killed tasks had no marker so they were
redone; finished tasks were skipped. Every epoch also checkpoints model, **optimiser state**,
epoch, best-so-far and history, so a resumed task continues the same trajectory. A resumed run
is scientifically equivalent, not bit-identical — the same magnitude of variation as BLAS
thread ordering, which `golden.yaml` already budgets for.

## Episode 22. I miscounted with a glob and reported it as fact

I told the user SpliceBERT was at 474/475 and "essentially done". It was at **195/475**.

I had counted `runs/dinuc/**/metrics.json`, which matches **both** models. The number was CNN
plus SpliceBERT. The CNN was the one nearly finished.

**Lesson.** A glob that matches more than you think reports a number that looks right. Count
per model. I had written this exact warning into doc 56 and then broke it.

## Episode 23. Stage 11 failed three times for three unrelated reasons

**(a) Wrong identity.** Ran as `rbp-ingest`, which has no access to the derived bucket at all.
I had conflated two unrelated things: **the external IP is a *network* property of the job**,
set by `EXTERNAL=1`, and says nothing about which identity it runs as. I picked the account
that can reach the internet instead of the one that can write the output.

The right fix was not a broader account but the correct one: variant assignment writes result
tables, so it is analysis work. `rbp-analysis` could *read* the derived bucket and not write to
it, which made every result-producing stage unassignable.

**(b) Could not read raw.** It then took a 403 on the raw bucket, needing the genome and
ClinVar VCF. Granted read-only there.

> **The trap in the middle of this:** fixing an access error by granting exactly the next thing
> that failed is how an identity quietly accumulates everything. The discipline is to ask what
> the stage *legitimately* needs. Stage 11 reads two immutable inputs and writes result tables.
> Nothing else.

**(c) A log line took the stage down.**

```python
b = raw_bucket.blob(name)
if not b.exists():          # True — but fetches NO metadata
    ...
log(f"downloading {name} ({b.size / 1e9:.2f} GB)")   # TypeError: NoneType / float
```

`bucket.blob()` builds a **local reference** and fetches nothing; `exists()` issues a HEAD and
returns a bool **without hydrating the object**. Reproduced against the live API: `exists()
True, size None`.

`get_blob()` does one GET and returns a hydrated blob or `None`, which also collapses the
existence check into the same call. The size is guarded with `(b.size or 0)` regardless,
because **a progress message must never be able to fail the thing it is reporting on.**

## Episode 24. The stage-in list was a contract I wrote from memory, twice

**First:** all 95 locality tasks failed with `'Series' object has no attribute 'cell'`. The
study panel writes `cell_line`; my task read `cell`.

**Second:** variants got much further — staged 95 datasets — then `no peak file for AATF`. My
`RAW_OBJECTS` list had the genome, its index and the ClinVar VCF: the three files the *phyloP*
half needs. But the *assign* half reads the ENCODE peaks first.

**Third:** I fixed that by staging peaks for the 95 study-panel datasets. It then failed on
`ADAT1`, which is **not in the study panel** — because `--what assign` walks the *full*
candidate panel, since a variant's nearest binding site can belong to a protein the study never
scores.

The peaks are **19.7 MiB** against a 3.1 GB genome already being downloaded. Being selective
saved nothing and cost two runs.

**Lesson, twice in one evening.** A stage-in list is a contract with **every function the stage
calls**, not just the one you had in mind while writing it. Drive it from what the code
*reads*, not from what seems relevant.

---

# Act V. The results

Three of four reproduced and verified as of this writing.

## R1 — the protocol effect

```
n = 94 datasets present in both arms
AUROC   GC-matched      0.7981
        dinuc-matched   0.6886
cost of proper matching  -0.1095
datasets that fall       94/94
paired Wilcoxon p        3.81e-17

THE FINDING: nested gain over composition
        GC-matched      +0.0265
        dinuc-matched   +0.0662     (2.50x)

COMPOSITION SHARE of skill above chance
        GC-matched      94.8%   [92.1, 97.4]
        dinuc-matched   67.8%   [62.1, 73.8]
        drop            27.0 points [23.1, 30.9]
```

**The finding is not the drop.** Under the harder control the model's gain over a
composition-only baseline **two-and-a-half times**. Rigour normally deflates deep learning;
here it reveals that the easy control was hiding most of what the model learned.

**CORRECTED 2026-08-26.** This section first reported +0.0154 -> +0.0607, a 3.94x ratio,
and the verifier certified it. Those were the difference between two SEPARATELY fitted
models, a quantity with no confidence interval and no p-value that nothing in the write-up
claims. The claim is the NESTED gain -- composition alone against composition plus the
sequence score -- which the rehearsal already computes as `delta_auroc` with a bootstrap CI
and a per-dataset `helps` flag. It is 2.50x. Smaller, and the only version that is defensible.

The clearer statement of the same result, added later, is the composition SHARE: under GC
matching a 19-feature composition model recovers 94.8% of the k-mer model's skill above
chance, so that benchmark is very nearly a composition test. Under dinucleotide matching it
recovers 67.8% -- better, and still most of it.

## R2 — four models on identical splits

| model | reproduced | reference |
|---|---|---|
| composition | 0.6279 | 0.628 |
| k-mer | 0.6875 | 0.688 |
| CNN | 0.7063 | 0.708 |
| SpliceBERT | 0.8091 | 0.809 |

Ordering intact. 95 datasets, 5 folds each, **0 NaN, 0 at chance, no duplicates**.

## R3 — positional concentration

```
k-mer Gini median       0.295
SpliceBERT Gini median  0.351
difference median      +0.064
more local in           91/95
significantly reversed  0/95
paired Wilcoxon p       3.94e-17
```

## Against the reference tolerances: 8 of 8 pass

Including the strong form of R3 — **reversed on none** — which is the claim that makes the
result more than a ranking.

## R4 — complete, and it took three corrections to become honest

The ClinVar arm finished on 2026-08-26: 95/95 matched and 95/95 mismatched on Modal, both
arms in sixteen minutes. What it first reported was wrong in a way that passed every check.

**As first computed (the pooled ladder):**

```
k-mer            0.5519
wrong protein    0.6797
right protein    0.8294
phyloP           0.9078
conservation-controlled coefficients   0.201 / 0.700 / 1.605
```

**Correction 1: pooling inflated it.** Concatenating ~19k variants across 95 datasets into
one AUROC per arm partly measures *which dataset a variant came from*. Mean |delta| per
dataset correlates with that dataset's pathogenic rate at Spearman **+0.73** and spans
**10.4x**. Paired within dataset, the matched arm is 0.755, not 0.829, and the specificity
gap is +0.065, not +0.149.

Conservation was the only arm immune, because phyloP is on a fixed external scale with no
between-dataset scale to contribute. **That is why it stayed invisible: the arm that could
not be inflated was winning anyway, so nothing looked out of place.**

**Correction 2: a trivial baseline beats the model.** "What fraction of the OTHER variants in
this 1-Mb window are pathogenic", leave-one-out, no sequence and no model:

| stratum | model | trivial rule | model wins | p |
|---|---|---|---|---|
| ≥20 pathogenic (n=44) | 0.7553 | **0.8139** | 15/44 | 0.007 *against* the model |
| ≥50 pathogenic (n=31) | 0.7997 | 0.8020 | 14/31 | 0.53 |

The model never beats it. **Absolute AUROC on peak-proximal ClinVar variants is therefore
uninformative about model utility**, and no version of this work may report 0.755 as evidence
of usefulness.

**What survives, and why it is better.** The specificity contrast is unaffected, because the
positional baseline applies equally to both arms: the right protein's head beats a wrong
protein's head by **+0.0645**, 33/44 datasets, **p=3.9e-04**, rising to +0.1031 (27/31) at
≥50 pathogenic. Within-dataset conservation-controlled coefficients **0.716** [0.640, 0.801]
against **0.445** [0.365, 0.523], non-overlapping.

So the result stops being "our model scores 0.83" and becomes something with a demonstrated
need: **AUROC on this task cannot tell you whether a model learned anything protein-specific,
and a wrong-protein control can.**

**Correction 3: the control was attacked and held.** Three checks, all negative for
contamination — the floor does not depend on the donor sharing a cell line (p=0.83), barely
tracks the donor's own strength (rho=+0.23), and stays flat at ~0.69 across power strata
while the matched arm climbs 0.66 → 0.80.

## The gate passed while all of this was wrong

Stage 14 reported **33/33** on the inflated numbers, and separately certified an R1 gain ratio
of 3.94x that the write-up never claimed. Both failures are the same: **it checked the value
the code produced, not whether that value was the right value.** It is now 56/56 and asserts
the unflattering results too — the all-datasets stratum that shows nothing, and a ceiling on
`model_minus_prevalence` that fails the build if anyone ever claims the model beat the
trivial rule.

---

# What the twenty-four bugs have in common

**Three quarters could not happen on a machine that already works.** They need a fresh
project, a fresh account, a fresh clone. That is not an argument for being careful; it is an
argument that **carefulness is not a substitute for a clean-room run**.

**Four could have produced wrong science without failing:** the panel-as-a-flag, the
non-waiting submit, the mismatched table names, and the stale image. All the same shape — **a
component that succeeds while doing the wrong thing.** Every guard in this pipeline exists
because of that shape.

**Five were safety mechanisms firing correctly on wrong instructions from me:** an IAM
condition (twice), a device guard, a completion marker, and an app-lifetime default. Those are
not failures of the system. They are the system working, and each one cost minutes instead of
a corrupted result.

**Three were my own reasoning, and no guard fixes those:**

- I chose a cause that fit the story I already had (Episode 14)
- I extrapolated a rate from one sample (Episode 15)
- I counted with a glob that matched more than I thought, and reported it as fact (Episode 22)

The defences there are mechanical, not intellectual: print the matching line rather than
counting matches; sample twice and divide; scope every glob as narrowly as the claim.
